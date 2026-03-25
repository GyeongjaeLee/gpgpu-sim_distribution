#include "booksim.hpp"
#include <vector>
#include <sstream>
#include <cassert>
#include <iostream>
#include <queue>
#include <map>
#include <limits>

#include "hbmnet.hpp"
#include "misc_utils.hpp"
#include "random_utils.hpp"

// ============================================================
//  Global tables for routing
// ============================================================
vector<vector<int>> gHBMNetDistBaseline;
vector<vector<int>> gHBMNetDistFabric;
vector<vector<pair<int,int>>> gHBMNetAdj;
vector<int> gHBMNetRouterConc;
vector<int> gHBMNetRouterFirstNode;

// Reverse mapping: node -> router, node -> port
static vector<int> gHBMNetNodeToRouter;
static vector<int> gHBMNetNodeToPort;

// Topology parameters cached for routing functions
static int gHBMNetNumSMs;
static int gHBMNetNumL2;
static int gHBMNetK;          // num_hbm_stacks
static int gHBMNetL;          // l2_per_hbm
static int gHBMNetInterleave;
static int gHBMNetUGALThreshold;
static int gHBMNetUGALIntmSelect;  // 0 = random, 1 = least-cost
static double gHBMNetBaselineRatio; // 0.0 = pure adaptive, 1.0 = pure baseline
static int gHBMNetHybridRouting;    // 0=baseline,1=min_adaptive,2=ugal,3=valiant,4=min_oblivious

// UGAL routing decision counters (reset each simulation run)
long long gUGALMinDecisions    = 0;  // source chose minimal path
long long gUGALNonMinDecisions = 0;  // source chose non-minimal (via intermediate)

// Escape VC usage: count flits that arrive at their destination on vc==0.
long long gEscapeVCEjects = 0;
long long gTotalEjects    = 0;

void hbmnet_reset_stats() {
  gUGALMinDecisions    = 0;
  gUGALNonMinDecisions = 0;
  gEscapeVCEjects      = 0;
  gTotalEjects         = 0;
}

int hbmnet_node_to_router(int node) { return gHBMNetNodeToRouter[node]; }
int hbmnet_node_to_port(int node)   { return gHBMNetNodeToPort[node]; }

// ============================================================
//  Router ID helpers
//    Router 0: Crossbar Partition 0
//    Router 1: Crossbar Partition 1
//    Router 2 .. K+1: Xbar HBM 0 .. K-1
// ============================================================
static inline int crossbar_router_id(int partition) { return partition; }
static inline int hbm_router_id(int hbm_idx) { return hbm_idx + 2; }
static inline bool is_crossbar_router(int router) { return router < 2; }
static inline int hbm_idx_from_router(int router) { return router - 2; }

// 2D layout helpers (need K)
static inline int hbm_col(int hbm_idx, int K) { return hbm_idx / (K / 2); }
static inline int hbm_row(int hbm_idx, int K) { return hbm_idx % (K / 2); }
static inline int hbm_partition(int hbm_idx, int K) {
  return (hbm_row(hbm_idx, K) < K / 4) ? 0 : 1;
}

// ============================================================
//  Constructor
// ============================================================
HBMNet::HBMNet( const Configuration &config, const string & name )
  : Network( config, name )
{
  _ComputeSize( config );
  _Alloc();
  _BuildNet( config );
  _BuildRoutingTable();
}

// ============================================================
//  _ComputeSize
// ============================================================
void HBMNet::_ComputeSize( const Configuration &config )
{
  _num_sms        = config.GetInt("num_sms");
  _num_l2_slices  = config.GetInt("num_l2_slices");
  _num_hbm_stacks = config.GetInt("num_hbm_stacks");
  _l2_per_hbm     = _num_l2_slices / _num_hbm_stacks;
  _l2_interleave  = config.GetInt("l2_interleave");
  _is_fabric      = config.GetInt("is_fabric");

  assert(_num_sms % 2 == 0);           // must split evenly across 2 partitions
  assert(_num_l2_slices % _num_hbm_stacks == 0);  // M divisible by K
  assert(_num_hbm_stacks % 4 == 0);    // K divisible by 4 for 2-col 2-partition layout

  _xbar_xbar_latency   = config.GetInt("xbar_xbar_latency");
  _xbar_hbm_latency    = config.GetInt("xbar_hbm_latency");
  _hbm_hbm_latency     = config.GetInt("hbm_hbm_latency");
  _xbar_xbar_bandwidth = config.GetInt("xbar_xbar_bandwidth");
  _xbar_hbm_bandwidth  = config.GetInt("xbar_hbm_bandwidth");
  _hbm_hbm_bandwidth   = config.GetInt("hbm_hbm_bandwidth");

  // Cache for routing functions
  gHBMNetNumSMs = _num_sms;
  gHBMNetNumL2  = _num_l2_slices;
  gHBMNetK      = _num_hbm_stacks;
  gHBMNetL      = _l2_per_hbm;
  gHBMNetInterleave = _l2_interleave;
  gHBMNetUGALThreshold  = config.GetInt("ugal_threshold");
  gHBMNetUGALIntmSelect = config.GetInt("ugal_intm_select");
  gHBMNetBaselineRatio  = config.GetFloat("baseline_ratio");

  string hr = config.GetStr("hybrid_routing");
  if      (hr == "min_adaptive") gHBMNetHybridRouting = 1;
  else if (hr == "ugal")         gHBMNetHybridRouting = 2;
  else if (hr == "valiant")      gHBMNetHybridRouting = 3;
  else if (hr == "min_oblivious") gHBMNetHybridRouting = 4;
  else                           gHBMNetHybridRouting = 0;

  // Total routers: 2 crossbars + K HBM routers
  _size = 2 + _num_hbm_stacks;

  // Build per-router concentration and node mapping
  gHBMNetRouterConc.resize(_size);
  gHBMNetRouterFirstNode.resize(_size);

  // Crossbar 0: SM nodes 0 .. N/2-1
  gHBMNetRouterConc[0] = _num_sms / 2;
  gHBMNetRouterFirstNode[0] = 0;

  // Crossbar 1: SM nodes N/2 .. N-1
  gHBMNetRouterConc[1] = _num_sms / 2;
  gHBMNetRouterFirstNode[1] = _num_sms / 2;

  // HBM routers: L2 slice nodes
  // We need to figure out which L2 slice nodes go to which HBM router.
  // Node IDs for L2 slices: N .. N+M-1
  // L2 slice index j (0-based): node ID = N + j
  //
  // Sequential: HBM h gets L2 slices j = h*L .. (h+1)*L-1
  // Interleaved: HBM h gets L2 slices j = h, h+K, h+2K, ..., h+(L-1)*K
  //
  // For concentration, each HBM router always has L nodes.
  // But we need contiguous node IDs per router for _Alloc to work.
  //
  // Since accel-sim expects node 0..N-1 = SM, N..N+M-1 = L2 slice,
  // we assign L2 slice nodes to HBM routers. The node ordering is
  // fixed (node N+j = L2 slice j), but the router assignment varies.

  // First, assign all HBM routers concentration = L
  for (int h = 0; h < _num_hbm_stacks; h++) {
    gHBMNetRouterConc[hbm_router_id(h)] = _l2_per_hbm;
  }

  // For _Alloc, we need contiguous node blocks per router.
  // The L2 nodes N..N+M-1 must be assigned to HBM routers.
  // We'll reorder: HBM router h gets nodes at positions determined by mapping.
  //
  // To keep node IDs contiguous per router, we assign:
  //   HBM router h -> nodes [N + h*L, N + (h+1)*L - 1]
  // The L2 interleave mapping is handled by the accel-sim side
  // when translating between "L2 slice index" and "booksim node ID".
  //
  // That is, for accel-sim integration:
  //   Sequential:   L2 slice j -> booksim node (N + j)
  //                 -> HBM router h = j / L
  //   Interleaved:  L2 slice j -> booksim node (N + h*L + offset)
  //                 where h = j % K, offset = j / K
  //                 This reorders L2 slice node IDs so that
  //                 contiguous blocks map to each HBM router.

  for (int h = 0; h < _num_hbm_stacks; h++) {
    gHBMNetRouterFirstNode[hbm_router_id(h)] = _num_sms + h * _l2_per_hbm;
  }

  _nodes = _num_sms + _num_l2_slices;

  // Build reverse mapping: node -> (router, port)
  gHBMNetNodeToRouter.resize(_nodes);
  gHBMNetNodeToPort.resize(_nodes);
  for (int rtr = 0; rtr < _size; rtr++) {
    int first = gHBMNetRouterFirstNode[rtr];
    int conc  = gHBMNetRouterConc[rtr];
    for (int p = 0; p < conc; p++) {
      gHBMNetNodeToRouter[first + p] = rtr;
      gHBMNetNodeToPort[first + p]   = p;
    }
  }

  // Count channels
  // Crossbar-Crossbar: 2 directions * bandwidth
  _xbar_xbar_channels = _xbar_xbar_bandwidth * 2;

  // Crossbar-HBM: each HBM connects to its partition's crossbar
  // K HBM routers, each with _xbar_hbm_bandwidth ports * 2 directions
  _xbar_hbm_channels = _num_hbm_stacks * _xbar_hbm_bandwidth * 2;

  // HBM-HBM (fabric)
  // Physical layout: HBMs are on left/right shorelines of the GPU die.
  // Column 0 = left shoreline, Column 1 = right shoreline.
  // Only vertical (same-column) adjacent HBMs can be connected.
  // No horizontal links (left↔right) — GPU die is in between.
  if (_is_fabric) {
    // Vertical pairs: (K/2 - 1) per column * 2 columns
    int vert_pairs = 2 * (_num_hbm_stacks / 2 - 1);
    _hbm_hbm_channels = vert_pairs * _hbm_hbm_bandwidth * 2;
  } else {
    _hbm_hbm_channels = 0;
  }

  _channels = _xbar_xbar_channels + _xbar_hbm_channels + _hbm_hbm_channels;

  cout << "HBMNet Config:"
       << " num_sms=" << _num_sms
       << " num_l2=" << _num_l2_slices
       << " K=" << _num_hbm_stacks
       << " L=" << _l2_per_hbm
       << " interleave=" << _l2_interleave
       << " is_fabric=" << _is_fabric << endl;
  cout << "  xbar_xbar_bw=" << _xbar_xbar_bandwidth
       << " xbar_hbm_bw=" << _xbar_hbm_bandwidth
       << " hbm_hbm_bw=" << _hbm_hbm_bandwidth << endl;
  cout << "  routers=" << _size
       << " nodes=" << _nodes
       << " channels=" << _channels << endl;
}

// ============================================================
//  _BuildNet
// ============================================================
void HBMNet::_BuildNet( const Configuration& config )
{
  ostringstream name;
  int chan_idx = 0;

  // Compute router degrees
  vector<int> degree(_size, 0);

  // Crossbar 0 & 1: concentration + HBM links + inter-partition links
  for (int p = 0; p < 2; p++) {
    degree[p] = gHBMNetRouterConc[p];  // SM eject ports
    // Links to HBM routers in this partition
    for (int h = 0; h < _num_hbm_stacks; h++) {
      if (hbm_partition(h, _num_hbm_stacks) == p) {
        degree[p] += _xbar_hbm_bandwidth;
      }
    }
    // Inter-partition link
    degree[p] += _xbar_xbar_bandwidth;
  }

  // HBM routers: concentration + crossbar link + fabric links
  for (int h = 0; h < _num_hbm_stacks; h++) {
    int rid = hbm_router_id(h);
    degree[rid] = gHBMNetRouterConc[rid];  // L2 eject ports
    degree[rid] += _xbar_hbm_bandwidth;   // link to crossbar

    if (_is_fabric) {
      int row = hbm_row(h, _num_hbm_stacks);
      int rows_per_col = _num_hbm_stacks / 2;

      // Vertical neighbors only (same shoreline column)
      // No horizontal links — GPU die separates left and right shorelines
      if (row > 0) degree[rid] += _hbm_hbm_bandwidth;
      if (row < rows_per_col - 1) degree[rid] += _hbm_hbm_bandwidth;
    }
  }

  // Create routers
  for (int rtr = 0; rtr < _size; rtr++) {
    if (is_crossbar_router(rtr)) {
      name << "router_xbar" << rtr;
    } else {
      int h = hbm_idx_from_router(rtr);
      name << "router_hbm" << h
           << "_c" << hbm_col(h, _num_hbm_stacks)
           << "_r" << hbm_row(h, _num_hbm_stacks);
    }
    _routers[rtr] = Router::NewRouter(config, this, name.str(), rtr,
                                       degree[rtr], degree[rtr]);
    _timed_modules.push_back(_routers[rtr]);
    name.str("");
  }

  // Inject/Eject
  int inject_lat = config.GetInt("inject_latency");
  int eject_lat  = config.GetInt("eject_latency");
  for (int n = 0; n < _nodes; n++) {
    int router = gHBMNetNodeToRouter[n];
    _routers[router]->AddInputChannel(_inject[n], _inject_cred[n]);
    _routers[router]->AddOutputChannel(_eject[n], _eject_cred[n]);
    _inject[n]->SetLatency(inject_lat);
    _eject[n]->SetLatency(eject_lat);
  }

  // Build adjacency
  gHBMNetAdj.clear();
  gHBMNetAdj.resize(_size);
  vector<int> opc(_size);
  for (int rtr = 0; rtr < _size; rtr++) {
    opc[rtr] = gHBMNetRouterConc[rtr];  // network ports start after eject ports
  }

  auto add_link = [&](int u, int v, int latency) {
    assert(chan_idx < _channels);
    _chan[chan_idx]->SetLatency(latency);
    _chan_cred[chan_idx]->SetLatency(latency);
    _routers[u]->AddOutputChannel(_chan[chan_idx], _chan_cred[chan_idx]);
    _routers[v]->AddInputChannel(_chan[chan_idx], _chan_cred[chan_idx]);
    gHBMNetAdj[u].push_back(make_pair(v, opc[u]++));
    chan_idx++;
  };

  // 1) Crossbar-HBM links
  for (int h = 0; h < _num_hbm_stacks; h++) {
    int part = hbm_partition(h, _num_hbm_stacks);
    int xbar = crossbar_router_id(part);
    int hbm  = hbm_router_id(h);
    for (int b = 0; b < _xbar_hbm_bandwidth; b++) {
      add_link(xbar, hbm, _xbar_hbm_latency);
      add_link(hbm, xbar, _xbar_hbm_latency);
    }
  }

  // 2) Crossbar-Crossbar inter-partition links
  for (int b = 0; b < _xbar_xbar_bandwidth; b++) {
    add_link(0, 1, _xbar_xbar_latency);
    add_link(1, 0, _xbar_xbar_latency);
  }

  // 3) HBM-HBM fabric links
  if (_is_fabric) {
    int rows_per_col = _num_hbm_stacks / 2;

    // Vertical links (same column, adjacent rows)
    for (int col = 0; col < 2; col++) {
      for (int row = 0; row < rows_per_col - 1; row++) {
        int h1 = col * rows_per_col + row;
        int h2 = col * rows_per_col + row + 1;
        int r1 = hbm_router_id(h1);
        int r2 = hbm_router_id(h2);
        for (int b = 0; b < _hbm_hbm_bandwidth; b++) {
          add_link(r1, r2, _hbm_hbm_latency);
          add_link(r2, r1, _hbm_hbm_latency);
        }
      }
    }

    // No horizontal links — left and right shorelines are separated by GPU die
  }

  assert(chan_idx == _channels);
}

// ============================================================
//  BFS helper — compute shortest distances on a given adjacency
// ============================================================
static void bfs_distances(int size,
                          const vector<vector<pair<int,int>>> &adj,
                          vector<vector<int>> &dist)
{
  dist.assign(size, vector<int>(size, -1));
  for (int src = 0; src < size; src++) {
    dist[src][src] = 0;
    queue<int> q;
    q.push(src);
    while (!q.empty()) {
      int u = q.front(); q.pop();
      for (size_t i = 0; i < adj[u].size(); i++) {
        int v = adj[u][i].first;
        if (dist[src][v] == -1) {
          dist[src][v] = dist[src][u] + 1;
          q.push(v);
        }
      }
    }
  }
}

// ============================================================
//  _BuildRoutingTable
//
//  Builds two distance tables:
//    gHBMNetDistBaseline — BFS on baseline-only links (Xbar-HBM, Xbar-Xbar)
//    gHBMNetDistFabric   — BFS on all links including HBM-HBM fabric
// ============================================================
void HBMNet::_BuildRoutingTable()
{
  // --- Baseline adjacency (no HBM-HBM links) ---
  // We extract only the Xbar-HBM and Xbar-Xbar edges from gHBMNetAdj.
  vector<vector<pair<int,int>>> baselineAdj(_size);
  for (int r = 0; r < _size; r++) {
    for (size_t i = 0; i < gHBMNetAdj[r].size(); i++) {
      int nb = gHBMNetAdj[r][i].first;
      int port = gHBMNetAdj[r][i].second;
      // Include edge only if at least one endpoint is a crossbar
      if (is_crossbar_router(r) || is_crossbar_router(nb)) {
        baselineAdj[r].push_back(make_pair(nb, port));
      }
    }
  }

  bfs_distances(_size, baselineAdj, gHBMNetDistBaseline);
  bfs_distances(_size, gHBMNetAdj, gHBMNetDistFabric);

  // Debug output
  cout << "=== Baseline Distance Table ===" << endl;
  for (int s = 0; s < _size; s++) {
    for (int d = 0; d < _size; d++) {
      cout << gHBMNetDistBaseline[s][d] << " ";
    }
    cout << endl;
  }

  cout << "=== Fabric Distance Table ===" << endl;
  for (int s = 0; s < _size; s++) {
    for (int d = 0; d < _size; d++) {
      cout << gHBMNetDistFabric[s][d] << " ";
    }
    cout << endl;
  }

  cout << "=== Adjacency Table ===" << endl;
  for (int r = 0; r < _size; r++) {
    cout << "  Router " << r << " (";
    if (is_crossbar_router(r))
      cout << "Crossbar" << r;
    else {
      int h = hbm_idx_from_router(r);
      cout << "HBM" << h << " c" << hbm_col(h, gHBMNetK) << "r" << hbm_row(h, gHBMNetK);
    }
    cout << ", conc=" << gHBMNetRouterConc[r] << "): ";
    map<int, int> nb_count;
    for (size_t i = 0; i < gHBMNetAdj[r].size(); i++) {
      nb_count[gHBMNetAdj[r][i].first]++;
    }
    for (auto &p : nb_count) {
      int nb = p.first;
      cout << "R" << nb << "(";
      if (is_crossbar_router(nb))
        cout << "Xbar" << nb;
      else {
        int h = hbm_idx_from_router(nb);
        cout << "HBM" << h;
      }
      cout << ")x" << p.second << " ";
    }
    cout << endl;
  }

  // Print L2 mapping info
  cout << "=== L2 Slice Mapping ===" << endl;
  cout << "  Mode: " << (gHBMNetInterleave ? "interleaved" : "sequential") << endl;
  for (int h = 0; h < gHBMNetK; h++) {
    int rid = hbm_router_id(h);
    cout << "  HBM" << h << " (Router " << rid
         << ", partition " << hbm_partition(h, gHBMNetK) << "): nodes ["
         << gHBMNetRouterFirstNode[rid] << ", "
         << gHBMNetRouterFirstNode[rid] + gHBMNetRouterConc[rid] - 1 << "]";
    if (gHBMNetInterleave) {
      cout << " -> L2 slices: ";
      for (int p = 0; p < gHBMNetL; p++) {
        cout << (h + p * gHBMNetK);
        if (p < gHBMNetL - 1) cout << ",";
      }
    } else {
      cout << " -> L2 slices: " << h * gHBMNetL << "-" << (h+1) * gHBMNetL - 1;
    }
    cout << endl;
  }
}

// ============================================================
//  RegisterRoutingFunctions
// ============================================================
void HBMNet::RegisterRoutingFunctions()
{
  gRoutingFunctionMap["baseline_hbmnet"]      = &hbmnet_baseline;
  gRoutingFunctionMap["min_adaptive_hbmnet"]  = &hbmnet_min_adaptive;
  gRoutingFunctionMap["ugal_hbmnet"]          = &hbmnet_ugal;
  gRoutingFunctionMap["valiant_hbmnet"]       = &hbmnet_valiant;
  gRoutingFunctionMap["min_oblivious_hbmnet"] = &hbmnet_min_oblivious;
  gRoutingFunctionMap["hybrid_hbmnet"]        = &hbmnet_hybrid;
}

// ============================================================
//  Routing helpers
// ============================================================
static inline int eject_port(int dest_node) {
  return hbmnet_node_to_port(dest_node);
}

// ============================================================
//  Determine baseline next-hop router (deterministic)
//
//  Routing tree:
//    SM -> own Crossbar -> [inter-partition if needed] -> dest Crossbar -> dest HBM -> L2
//    L2 -> HBM -> own Crossbar -> [inter-partition if needed] -> dest Crossbar -> SM
//    L2 -> HBM -> own Crossbar -> [inter-partition if needed] -> dest Crossbar -> dest HBM -> L2
//
//  In all cases: src goes up to Crossbar, possibly crosses partitions,
//  then goes down to destination.
// ============================================================
static int get_baseline_next_hop(int cur_router, int dest_router) {
  int K = gHBMNetK;

  if (is_crossbar_router(cur_router)) {
    // Current is a crossbar router
    if (is_crossbar_router(dest_router)) {
      // Dest is the other crossbar (inter-partition)
      return dest_router;
    }
    // Dest is an HBM router
    int dest_h = hbm_idx_from_router(dest_router);
    int dest_part = hbm_partition(dest_h, K);
    if (dest_part == cur_router) {
      // Dest HBM is in our partition -> go directly to it
      return dest_router;
    } else {
      // Dest HBM is in other partition -> go to other crossbar first
      return crossbar_router_id(1 - cur_router);
    }
  } else {
    // Current is an HBM router -> always go up to own crossbar
    int cur_h = hbm_idx_from_router(cur_router);
    int cur_part = hbm_partition(cur_h, K);
    return crossbar_router_id(cur_part);
  }
}

// ============================================================
//  Common routing preamble: inject, eject, and escape channel
//
//  All routing functions share this for fair VC comparison and
//  deadlock freedom via Duato's protocol:
//    VC 0:       escape channel (baseline tree routing, pri=0)
//    VCs 1..V-1: data channels (routing-specific)
//
//  Returns true if routing is complete (inject or eject handled).
//  Otherwise, adds the escape channel and sets cur_router/dest_router
//  for the caller to add data channel outputs on VCs 1..V-1.
// ============================================================
static bool routing_preamble(const Router *r, const Flit *f,
                             OutputSet *outputs, bool inject,
                             int &cur_router, int &dest_router)
{
  outputs->Clear();

  if (inject) {
    outputs->AddRange(-1, 1, gNumVCs - 1);
    return true;
  }

  cur_router  = r->GetID();
  dest_router = hbmnet_node_to_router(f->dest);

  if (cur_router == dest_router) {
    ++gTotalEjects;
    if (f->vc == 0) ++gEscapeVCEjects;
    outputs->AddRange(eject_port(f->dest), 0, gNumVCs - 1);
    return true;
  }

  assert(gNumVCs >= 2);

  // Escape channel (VC 0): deterministic baseline tree routing
  // Pick one random port among parallel links to baseline next-hop.
  int baseline_next = get_baseline_next_hop(cur_router, dest_router);
  vector<int> esc_ports;
  for (size_t i = 0; i < gHBMNetAdj[cur_router].size(); i++) {
    if (gHBMNetAdj[cur_router][i].first == baseline_next)
      esc_ports.push_back(gHBMNetAdj[cur_router][i].second);
  }
  assert(!esc_ports.empty());
  outputs->AddRange(esc_ports[RandomInt(esc_ports.size() - 1)], 0, 0, 0);

  return false;
}

// ============================================================
//  BASELINE ROUTING — deterministic tree
//  L2 -> HBM -> Crossbar -> [Crossbar] -> Crossbar -> HBM -> L2
//
//  Data channels (VCs 1..V-1) also use baseline tree routing.
//  Escape VC is redundant here (tree is already deadlock-free)
//  but included for fair VC-count comparison across algorithms.
// ============================================================
void hbmnet_baseline( const Router *r, const Flit *f, int in_channel,
                        OutputSet *outputs, bool inject )
{
  int cur_router, dest_router;
  if (routing_preamble(r, f, outputs, inject, cur_router, dest_router))
    return;


  // Data channels (VCs 1..V-1): deterministic baseline tree routing
  int next_router = get_baseline_next_hop(cur_router, dest_router);
  vector<int> base_ports;
  for (size_t i = 0; i < gHBMNetAdj[cur_router].size(); i++) {
    if (gHBMNetAdj[cur_router][i].first == next_router)
      base_ports.push_back(gHBMNetAdj[cur_router][i].second);
  }
  assert(!base_ports.empty());
  outputs->AddRange(base_ports[RandomInt(base_ports.size() - 1)], 1, gNumVCs - 1, 1);
}

// ============================================================
//  MINIMAL OBLIVIOUS ROUTING — random minimal path selection
//
//  Collects all minimal next-hop ports and picks one uniformly
//  at random, ignoring congestion. Provides natural load balancing
//  across parallel ports without herding.
// ============================================================
void hbmnet_min_oblivious( const Router *r, const Flit *f, int in_channel,
                                    OutputSet *outputs, bool inject )
{
  int cur_router, dest_router;
  if (routing_preamble(r, f, outputs, inject, cur_router, dest_router))
    return;


  // Collect all minimal ports
  int cur_dist = gHBMNetDistFabric[cur_router][dest_router];
  vector<int> min_ports;
  for (size_t i = 0; i < gHBMNetAdj[cur_router].size(); i++) {
    int nb   = gHBMNetAdj[cur_router][i].first;
    int port = gHBMNetAdj[cur_router][i].second;
    if (gHBMNetDistFabric[nb][dest_router] == cur_dist - 1)
      min_ports.push_back(port);
  }
  assert(!min_ports.empty());
  outputs->AddRange(min_ports[RandomInt(min_ports.size() - 1)], 1, gNumVCs - 1, 1);
}

// ============================================================
//  MINIMAL ADAPTIVE ROUTING — Duato's protocol
//
//  Two-level decision:
//  1) Group minimal ports by next-hop router, compute average
//     credit usage per port for each direction.
//  2) Pick the direction with lowest average congestion.
//  3) Within that direction, pick a random port (all equivalent).
// ============================================================
void hbmnet_min_adaptive( const Router *r, const Flit *f, int in_channel,
                                   OutputSet *outputs, bool inject )
{
  int cur_router, dest_router;
  if (routing_preamble(r, f, outputs, inject, cur_router, dest_router))
    return;


  // Group minimal ports by next-hop router
  int cur_dist = gHBMNetDistFabric[cur_router][dest_router];
  map<int, vector<int>> dir_ports;  // next-hop router -> list of ports
  map<int, int> dir_credit;         // next-hop router -> total credit used
  for (size_t i = 0; i < gHBMNetAdj[cur_router].size(); i++) {
    int nb   = gHBMNetAdj[cur_router][i].first;
    int port = gHBMNetAdj[cur_router][i].second;
    if (gHBMNetDistFabric[nb][dest_router] == cur_dist - 1) {
      dir_ports[nb].push_back(port);
      dir_credit[nb] += r->GetUsedCredit(port);
    }
  }
  assert(!dir_ports.empty());

  // Pick direction with lowest average credit per port
  int best_nb = -1;
  double best_avg = numeric_limits<double>::max();
  for (map<int, vector<int>>::iterator it = dir_ports.begin(); it != dir_ports.end(); ++it) {
    double avg = (double)dir_credit[it->first] / (double)it->second.size();
    if (avg < best_avg) {
      best_avg = avg;
      best_nb = it->first;
    }
  }

  // Random port within chosen direction
  vector<int> &ports = dir_ports[best_nb];
  outputs->AddRange(ports[RandomInt(ports.size() - 1)], 1, gNumVCs - 1, 1);
}

// ============================================================
//  UGAL ROUTING — Universal Globally-Adaptive Load-balanced
//
//  At injection, compares minimal vs non-minimal (via random
//  intermediate) using avg_credit_per_port * hop_count.
//  Phase 0: route toward intermediate. Phase 1: route toward dest.
//  Per-hop: pick direction with lowest avg credit, random port within.
//  Config: ugal_threshold biases toward minimal routing.
// ============================================================

// Helper: compute average credit per port for a given direction (next-hop router)
// among minimal ports toward target_router.
static double avg_credit_toward(const Router *r, int cur_router, int target_router)
{
  int target_dist = gHBMNetDistFabric[cur_router][target_router];
  int total_credit = 0;
  int port_count = 0;
  for (size_t i = 0; i < gHBMNetAdj[cur_router].size(); i++) {
    int nb   = gHBMNetAdj[cur_router][i].first;
    int port = gHBMNetAdj[cur_router][i].second;
    if (gHBMNetDistFabric[nb][target_router] == target_dist - 1) {
      total_credit += r->GetUsedCredit(port);
      port_count++;
    }
  }
  if (port_count == 0) return 0.0;
  return (double)total_credit / (double)port_count;
}

void hbmnet_ugal( const Router *r, const Flit *f, int in_channel,
                           OutputSet *outputs, bool inject )
{
  int cur_router, dest_router;
  if (routing_preamble(r, f, outputs, inject, cur_router, dest_router))
    return;


  int K = gHBMNetK;

  // ===== UGAL decision at source router =====
  if (in_channel < gHBMNetRouterConc[cur_router]) {

    if (is_crossbar_router(cur_router)) {
      f->ph = 1;
    } else {
      // -- Minimal path cost: avg credit per port * hop distance --
      int min_dist = gHBMNetDistFabric[cur_router][dest_router];
      double min_avg = avg_credit_toward(r, cur_router, dest_router);

      // -- Non-minimal intermediate selection --
      int intm_router = -1;
      double nonmin_cost = 0.0;

      vector<int> intm_candidates;
      for (int h = 0; h < K; h++) {
        int intm_r = hbm_router_id(h);
        if (intm_r != cur_router && intm_r != dest_router)
          intm_candidates.push_back(intm_r);
      }

      if (intm_candidates.empty()) {
        f->ph = 1;
      } else {
        if (gHBMNetUGALIntmSelect == 0) {
          // Random intermediate
          intm_router = intm_candidates[RandomInt(intm_candidates.size() - 1)];
        } else {
          // Least-cost intermediate
          double best_cost = numeric_limits<double>::max();
          for (size_t c = 0; c < intm_candidates.size(); c++) {
            int intm_r = intm_candidates[c];
            int nt = gHBMNetDistFabric[cur_router][intm_r]
                   + gHBMNetDistFabric[intm_r][dest_router];
            double nc = avg_credit_toward(r, cur_router, intm_r);
            double cost = nt * nc;
            if (cost < best_cost) {
              best_cost = cost;
              intm_router = intm_r;
            }
          }
        }

        int nonmin_total = gHBMNetDistFabric[cur_router][intm_router]
                         + gHBMNetDistFabric[intm_router][dest_router];
        double nonmin_avg = avg_credit_toward(r, cur_router, intm_router);
        nonmin_cost = nonmin_total * nonmin_avg;

        if (min_dist * min_avg <= nonmin_cost + gHBMNetUGALThreshold) {
          f->ph = 1;  // minimal
          ++gUGALMinDecisions;
        } else {
          f->ph = 0;  // non-minimal
          f->intm = intm_router;
          ++gUGALNonMinDecisions;
        }
      }
    }
  }

  // ===== Phase transition: reached intermediate → switch to minimal =====
  if (f->ph == 0 && cur_router == f->intm) {
    f->ph = 1;
  }

  // ===== Per-hop: direction with lowest avg credit, random port within =====
  int target_router = (f->ph == 0) ? f->intm : dest_router;
  int target_dist = gHBMNetDistFabric[cur_router][target_router];
  map<int, vector<int>> dir_ports;
  map<int, int> dir_credit;
  for (size_t i = 0; i < gHBMNetAdj[cur_router].size(); i++) {
    int nb   = gHBMNetAdj[cur_router][i].first;
    int port = gHBMNetAdj[cur_router][i].second;
    if (gHBMNetDistFabric[nb][target_router] == target_dist - 1) {
      dir_ports[nb].push_back(port);
      dir_credit[nb] += r->GetUsedCredit(port);
    }
  }
  assert(!dir_ports.empty());

  int best_nb = -1;
  double best_avg = numeric_limits<double>::max();
  for (map<int, vector<int>>::iterator it = dir_ports.begin(); it != dir_ports.end(); ++it) {
    double avg = (double)dir_credit[it->first] / (double)it->second.size();
    if (avg < best_avg) {
      best_avg = avg;
      best_nb = it->first;
    }
  }

  vector<int> &ports = dir_ports[best_nb];
  outputs->AddRange(ports[RandomInt(ports.size() - 1)], 1, gNumVCs - 1, 1);
}

// ============================================================
//  VALIANT ROUTING — Oblivious non-minimal routing
//
//  Always routes via a random intermediate HBM router (no congestion
//  check). Phase 0: minimal toward intermediate. Phase 1: minimal
//  toward destination. Per-hop uses random minimal port selection.
//  Crossbar sources skip the detour (no benefit).
// ============================================================
void hbmnet_valiant( const Router *r, const Flit *f, int in_channel,
                              OutputSet *outputs, bool inject )
{
  int cur_router, dest_router;
  if (routing_preamble(r, f, outputs, inject, cur_router, dest_router))
    return;


  int K = gHBMNetK;

  // ===== Valiant decision at source router =====
  if (in_channel < gHBMNetRouterConc[cur_router]) {
    if (is_crossbar_router(cur_router)) {
      f->ph = 1;  // crossbar can't usefully detour
    } else {
      vector<int> intm_candidates;
      for (int h = 0; h < K; h++) {
        int intm_r = hbm_router_id(h);
        if (intm_r != cur_router && intm_r != dest_router)
          intm_candidates.push_back(intm_r);
      }
      if (intm_candidates.empty()) {
        f->ph = 1;
      } else {
        f->ph = 0;  // always non-minimal
        f->intm = intm_candidates[RandomInt(intm_candidates.size() - 1)];
      }
    }
  }

  // Phase transition: reached intermediate → switch to minimal
  if (f->ph == 0 && cur_router == f->intm) {
    f->ph = 1;
  }

  // Per-hop: random minimal toward target
  int target_router = (f->ph == 0) ? f->intm : dest_router;
  int target_dist = gHBMNetDistFabric[cur_router][target_router];
  vector<int> val_ports;
  for (size_t i = 0; i < gHBMNetAdj[cur_router].size(); i++) {
    int nb   = gHBMNetAdj[cur_router][i].first;
    int port = gHBMNetAdj[cur_router][i].second;
    if (gHBMNetDistFabric[nb][target_router] == target_dist - 1)
      val_ports.push_back(port);
  }
  assert(!val_ports.empty());
  outputs->AddRange(val_ports[RandomInt(val_ports.size() - 1)], 1, gNumVCs - 1, 1);
}

// ============================================================
//  Hybrid routing: probabilistic baseline + adaptive
//  With probability baseline_ratio → baseline tree routing
//  Otherwise → hybrid_routing function (fabric/min_adaptive/ugal)
//  Decision is made once at injection and stored in f->intm.
// ============================================================
static const int HYBRID_BASELINE_MARKER = -100;

void hbmnet_hybrid( const Router *r, const Flit *f, int in_channel,
                             OutputSet *outputs, bool inject )
{
  // At injection, decide routing mode for this flit
  if (inject) {
    if (RandomInt(9999) < (int)(gHBMNetBaselineRatio * 10000.0)) {
      f->intm = HYBRID_BASELINE_MARKER;
    }
    // Handle injection output (same for all routing modes)
    int dummy1, dummy2;
    routing_preamble(r, f, outputs, inject, dummy1, dummy2);
    return;
  }

  // Subsequent hops: dispatch based on per-flit routing decision
  if (f->intm == HYBRID_BASELINE_MARKER) {
    hbmnet_baseline(r, f, in_channel, outputs, inject);
  } else {
    switch (gHBMNetHybridRouting) {
      case 1:  hbmnet_min_adaptive(r, f, in_channel, outputs, inject); break;
      case 2:  hbmnet_ugal(r, f, in_channel, outputs, inject); break;
      case 3:  hbmnet_valiant(r, f, in_channel, outputs, inject); break;
      case 4:  hbmnet_min_oblivious(r, f, in_channel, outputs, inject); break;
      default: hbmnet_baseline(r, f, in_channel, outputs, inject); break;
    }
  }
}
