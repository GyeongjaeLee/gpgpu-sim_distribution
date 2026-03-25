#include "booksim.hpp"
#include <vector>
#include <sstream>
#include <cassert>
#include <iostream>
#include <queue>
#include <map>
#include <limits>

#include "hbmnet_accelsim.hpp"
#include "misc_utils.hpp"
#include "random_utils.hpp"

// ============================================================
//  Global tables for routing
// ============================================================
vector<vector<AccelAdjEntry>> gHBMNetAccelAdj;
vector<vector<int>> gHBMNetAccelDistHit;
vector<vector<int>> gHBMNetAccelDistMissBaseline;
vector<vector<int>> gHBMNetAccelDistMissFabric;
vector<int> gHBMNetAccelRouterConc;
vector<int> gHBMNetAccelRouterFirstNode;
int gHBMNetAccelK = 0;

// Reverse mapping: node -> router, node -> port
static vector<int> gAccelNodeToRouter;
static vector<int> gAccelNodeToPort;

// Topology parameters cached for routing
static int gAccelNumSMs;
static int gAccelNumL2;
static int gAccelL;          // l2_per_hbm
static int gAccelInterleave;
static int gAccelUGALThreshold;
static int gAccelUGALIntmSelect;
static double gAccelHitRate;     // baseline_ratio repurposed as L2 hit rate
static int gAccelHybridRouting;  // 0=baseline,1=min_adaptive,2=ugal,3=valiant,4=min_oblivious

// UGAL routing decision counters
long long gAccelUGALMinDecisions    = 0;
long long gAccelUGALNonMinDecisions = 0;

// Escape VC usage counters
long long gAccelEscapeVCEjects = 0;
long long gAccelTotalEjects    = 0;

void hbmnet_accelsim_reset_stats() {
  gAccelUGALMinDecisions    = 0;
  gAccelUGALNonMinDecisions = 0;
  gAccelEscapeVCEjects      = 0;
  gAccelTotalEjects         = 0;
}

int hbmnet_accelsim_node_to_router(int node) { return gAccelNodeToRouter[node]; }
int hbmnet_accelsim_node_to_port(int node)   { return gAccelNodeToPort[node]; }

// ============================================================
//  Router ID helpers
//    Router 0,1    : Crossbar 0,1
//    Router 2..K+1 : MC 0..K-1
//    Router K+2..2K+1: HBM 0..K-1
// ============================================================
static inline int accel_crossbar_id(int p) { return p; }
static inline int accel_mc_id(int mc_idx)  { return mc_idx + 2; }
static inline int accel_hbm_id(int hbm_idx){ return hbm_idx + 2 + gHBMNetAccelK; }
static inline bool accel_is_xbar(int r)    { return r < 2; }
static inline bool accel_is_mc(int r)      { return r >= 2 && r < 2 + gHBMNetAccelK; }
static inline bool accel_is_hbm(int r)     { return r >= 2 + gHBMNetAccelK; }
static inline int accel_mc_idx(int r)      { return r - 2; }
static inline int accel_hbm_idx(int r)     { return r - 2 - gHBMNetAccelK; }

// 2D layout helpers
static inline int accel_hbm_col(int h)  { return h / (gHBMNetAccelK / 2); }
static inline int accel_hbm_row(int h)  { return h % (gHBMNetAccelK / 2); }
static inline int accel_hbm_part(int h) {
  return (accel_hbm_row(h) < gHBMNetAccelK / 4) ? 0 : 1;
}

// ============================================================
//  Link type validity for routing
//
//  Hit flits:
//    At Xbar: XBAR_HBM, XBAR_XBAR
//    At HBM:  XBAR_HBM (reverse, toward Xbar)
//
//  Miss flits — physical (includes MC→HBM final hop, used by escape VC):
//    At Xbar: XBAR_MC, XBAR_XBAR
//    At MC:   MC_HBM, MC_MC, XBAR_MC (back to Xbar)
//
//  Miss flits — routing decisions (distance table, adaptive selection):
//    Operates on Xbar+MC logical network only (HBM excluded).
//    At Xbar: XBAR_MC, XBAR_XBAR
//    At MC:   MC_MC, XBAR_MC (back to Xbar)
//    MC→HBM is deterministic final hop, NOT part of adaptive routing.
// ============================================================
static inline bool accel_valid_hit_out(int r, AccelLinkType t) {
  if (accel_is_xbar(r))
    return t == ACCEL_LINK_XBAR_HBM || t == ACCEL_LINK_XBAR_XBAR;
  if (accel_is_hbm(r))
    return t == ACCEL_LINK_XBAR_HBM;  // HBM → Xbar direction
  return false;
}

// Physical miss links (escape VC, port filtering)
static inline bool accel_valid_miss_out(int r, AccelLinkType t) {
  if (accel_is_xbar(r))
    return t == ACCEL_LINK_XBAR_MC || t == ACCEL_LINK_XBAR_XBAR;
  if (accel_is_mc(r))
    return t == ACCEL_LINK_MC_HBM || t == ACCEL_LINK_MC_MC || t == ACCEL_LINK_XBAR_MC;
  if (accel_is_hbm(r))
    return t == ACCEL_LINK_MC_HBM;  // HBM → MC direction for L2 miss
  return false;
}

// Miss routing decisions (distance table, adaptive/minimal port selection)
// Only Xbar+MC network: excludes MC→HBM (deterministic final hop)
static inline bool accel_valid_miss_routing_out(int r, AccelLinkType t) {
  if (accel_is_xbar(r))
    return t == ACCEL_LINK_XBAR_MC || t == ACCEL_LINK_XBAR_XBAR;
  if (accel_is_mc(r))
    return t == ACCEL_LINK_MC_MC || t == ACCEL_LINK_XBAR_MC;
  return false;
}

// ============================================================
//  Constructor
// ============================================================
HBMNetAccelSim::HBMNetAccelSim( const Configuration &config, const string & name )
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
void HBMNetAccelSim::_ComputeSize( const Configuration &config )
{
  _num_sms        = config.GetInt("num_sms");
  _num_l2_slices  = config.GetInt("num_l2_slices");
  _num_hbm_stacks = config.GetInt("num_hbm_stacks");
  _l2_per_hbm     = _num_l2_slices / _num_hbm_stacks;
  _l2_interleave  = config.GetInt("l2_interleave");
  _is_fabric      = config.GetInt("is_fabric");

  assert(_num_sms % 2 == 0);
  assert(_num_l2_slices % _num_hbm_stacks == 0);
  assert(_num_hbm_stacks % 4 == 0);

  _xbar_xbar_latency   = config.GetInt("xbar_xbar_latency");
  _xbar_hbm_latency    = config.GetInt("xbar_hbm_latency");
  _xbar_mc_latency     = config.GetInt("xbar_mc_latency");
  _mc_hbm_latency      = config.GetInt("mc_hbm_latency");
  _mc_mc_latency       = config.GetInt("mc_mc_latency");
  _xbar_xbar_bandwidth = config.GetInt("xbar_xbar_bandwidth");
  _xbar_hbm_bandwidth  = config.GetInt("xbar_hbm_bandwidth");
  _xbar_mc_bandwidth   = config.GetInt("xbar_mc_bandwidth");
  _mc_hbm_bandwidth    = config.GetInt("mc_hbm_bandwidth");
  _mc_mc_bandwidth     = config.GetInt("mc_mc_bandwidth");

  // Cache for routing
  int K = _num_hbm_stacks;
  gHBMNetAccelK      = K;
  gAccelNumSMs       = _num_sms;
  gAccelNumL2        = _num_l2_slices;
  gAccelL            = _l2_per_hbm;
  gAccelInterleave   = _l2_interleave;
  gAccelUGALThreshold  = config.GetInt("ugal_threshold");
  gAccelUGALIntmSelect = config.GetInt("ugal_intm_select");
  gAccelHitRate        = config.GetFloat("baseline_ratio");  // repurposed as hit rate

  string hr = config.GetStr("hybrid_routing");
  if      (hr == "min_adaptive")  gAccelHybridRouting = 1;
  else if (hr == "ugal")          gAccelHybridRouting = 2;
  else if (hr == "valiant")       gAccelHybridRouting = 3;
  else if (hr == "min_oblivious") gAccelHybridRouting = 4;
  else                            gAccelHybridRouting = 0;

  // Total routers: 2 Xbars + K MCs + K HBMs
  _size = 2 + 2 * K;

  // Per-router concentration and first node
  gHBMNetAccelRouterConc.resize(_size);
  gHBMNetAccelRouterFirstNode.resize(_size);

  // Xbar 0: SM 0..N/2-1
  gHBMNetAccelRouterConc[0] = _num_sms / 2;
  gHBMNetAccelRouterFirstNode[0] = 0;
  // Xbar 1: SM N/2..N-1
  gHBMNetAccelRouterConc[1] = _num_sms / 2;
  gHBMNetAccelRouterFirstNode[1] = _num_sms / 2;

  // MC routers: no nodes attached (pure transit)
  for (int m = 0; m < K; m++) {
    gHBMNetAccelRouterConc[accel_mc_id(m)] = 0;
    gHBMNetAccelRouterFirstNode[accel_mc_id(m)] = _num_sms + _num_l2_slices; // sentinel
  }

  // HBM routers: L2 slice nodes
  for (int h = 0; h < K; h++) {
    gHBMNetAccelRouterConc[accel_hbm_id(h)] = _l2_per_hbm;
    gHBMNetAccelRouterFirstNode[accel_hbm_id(h)] = _num_sms + h * _l2_per_hbm;
  }

  _nodes = _num_sms + _num_l2_slices;

  // Node → router/port mapping
  gAccelNodeToRouter.resize(_nodes);
  gAccelNodeToPort.resize(_nodes);
  for (int rtr = 0; rtr < _size; rtr++) {
    int first = gHBMNetAccelRouterFirstNode[rtr];
    int conc  = gHBMNetAccelRouterConc[rtr];
    for (int p = 0; p < conc; p++) {
      if (first + p < _nodes) {
        gAccelNodeToRouter[first + p] = rtr;
        gAccelNodeToPort[first + p]   = p;
      }
    }
  }

  // Count channels
  int xbar_xbar_ch = _xbar_xbar_bandwidth * 2;
  int xbar_hbm_ch  = K * _xbar_hbm_bandwidth * 2;  // hit path
  int xbar_mc_ch   = K * _xbar_mc_bandwidth * 2;    // miss path entry
  int mc_hbm_ch    = K * _mc_hbm_bandwidth * 2;     // miss path local
  int mc_mc_ch     = 0;
  if (_is_fabric) {
    int vert_pairs = 2 * (K / 2 - 1);
    mc_mc_ch = vert_pairs * _mc_mc_bandwidth * 2;
  }
  _channels = xbar_xbar_ch + xbar_hbm_ch + xbar_mc_ch + mc_hbm_ch + mc_mc_ch;

  cout << "HBMNetAccelSim Config:"
       << " num_sms=" << _num_sms
       << " num_l2=" << _num_l2_slices
       << " K=" << K
       << " L=" << _l2_per_hbm
       << " is_fabric=" << _is_fabric
       << " hit_rate=" << gAccelHitRate << endl;
  cout << "  xbar_xbar_bw=" << _xbar_xbar_bandwidth
       << " xbar_hbm_bw=" << _xbar_hbm_bandwidth
       << " xbar_mc_bw=" << _xbar_mc_bandwidth
       << " mc_hbm_bw=" << _mc_hbm_bandwidth
       << " mc_mc_bw=" << _mc_mc_bandwidth << endl;
  cout << "  routers=" << _size
       << " (2 Xbar + " << K << " MC + " << K << " HBM)"
       << " nodes=" << _nodes
       << " channels=" << _channels << endl;
}

// ============================================================
//  _BuildNet
// ============================================================
void HBMNetAccelSim::_BuildNet( const Configuration &config )
{
  ostringstream name;
  int chan_idx = 0;
  int K = _num_hbm_stacks;

  // Compute router degrees
  vector<int> degree(_size, 0);

  // Xbar routers: concentration + hit links + miss links + inter-partition
  for (int p = 0; p < 2; p++) {
    degree[p] = gHBMNetAccelRouterConc[p];  // SM eject ports
    for (int h = 0; h < K; h++) {
      if (accel_hbm_part(h) == p) {
        degree[p] += _xbar_hbm_bandwidth;  // hit path to HBM h
        degree[p] += _xbar_mc_bandwidth;   // miss path to MC h
      }
    }
    degree[p] += _xbar_xbar_bandwidth;  // inter-partition
  }

  // MC routers: 0 concentration + Xbar link + local HBM link + fabric
  for (int m = 0; m < K; m++) {
    int rid = accel_mc_id(m);
    degree[rid] = 0;
    degree[rid] += _xbar_mc_bandwidth;   // link to partition Xbar
    degree[rid] += _mc_hbm_bandwidth;    // link to local HBM
    if (_is_fabric) {
      int row = accel_hbm_row(m);
      int rows_per_col = K / 2;
      if (row > 0) degree[rid] += _mc_mc_bandwidth;
      if (row < rows_per_col - 1) degree[rid] += _mc_mc_bandwidth;
    }
  }

  // HBM routers: concentration + Xbar hit link + MC miss link
  for (int h = 0; h < K; h++) {
    int rid = accel_hbm_id(h);
    degree[rid] = gHBMNetAccelRouterConc[rid];  // L2 eject ports
    degree[rid] += _xbar_hbm_bandwidth;  // hit path from Xbar
    degree[rid] += _mc_hbm_bandwidth;    // miss path from local MC
  }

  // Create routers
  for (int rtr = 0; rtr < _size; rtr++) {
    if (accel_is_xbar(rtr))
      name << "router_xbar" << rtr;
    else if (accel_is_mc(rtr)) {
      int m = accel_mc_idx(rtr);
      name << "router_mc" << m
           << "_c" << accel_hbm_col(m)
           << "_r" << accel_hbm_row(m);
    } else {
      int h = accel_hbm_idx(rtr);
      name << "router_hbm" << h
           << "_c" << accel_hbm_col(h)
           << "_r" << accel_hbm_row(h);
    }
    _routers[rtr] = Router::NewRouter(config, this, name.str(), rtr,
                                       degree[rtr], degree[rtr]);
    _timed_modules.push_back(_routers[rtr]);
    name.str("");
  }

  // Inject/Eject (only for Xbar and HBM routers that have nodes)
  int inject_lat = config.GetInt("inject_latency");
  int eject_lat  = config.GetInt("eject_latency");
  for (int n = 0; n < _nodes; n++) {
    int router = gAccelNodeToRouter[n];
    _routers[router]->AddInputChannel(_inject[n], _inject_cred[n]);
    _routers[router]->AddOutputChannel(_eject[n], _eject_cred[n]);
    _inject[n]->SetLatency(inject_lat);
    _eject[n]->SetLatency(eject_lat);
  }

  // Build adjacency
  gHBMNetAccelAdj.clear();
  gHBMNetAccelAdj.resize(_size);
  vector<int> opc(_size);
  for (int rtr = 0; rtr < _size; rtr++)
    opc[rtr] = gHBMNetAccelRouterConc[rtr];  // network ports start after eject ports

  auto add_link = [&](int u, int v, int latency, AccelLinkType type) {
    assert(chan_idx < _channels);
    _chan[chan_idx]->SetLatency(latency);
    _chan_cred[chan_idx]->SetLatency(latency);
    _routers[u]->AddOutputChannel(_chan[chan_idx], _chan_cred[chan_idx]);
    _routers[v]->AddInputChannel(_chan[chan_idx], _chan_cred[chan_idx]);
    gHBMNetAccelAdj[u].push_back({v, opc[u]++, type});
    chan_idx++;
  };

  // 1) Xbar ↔ HBM (hit path)
  for (int h = 0; h < K; h++) {
    int part = accel_hbm_part(h);
    int xbar = accel_crossbar_id(part);
    int hbm  = accel_hbm_id(h);
    for (int b = 0; b < _xbar_hbm_bandwidth; b++) {
      add_link(xbar, hbm, _xbar_hbm_latency, ACCEL_LINK_XBAR_HBM);
      add_link(hbm, xbar, _xbar_hbm_latency, ACCEL_LINK_XBAR_HBM);
    }
  }

  // 2) Xbar ↔ MC (miss path entry)
  for (int m = 0; m < K; m++) {
    int part = accel_hbm_part(m);
    int xbar = accel_crossbar_id(part);
    int mc   = accel_mc_id(m);
    for (int b = 0; b < _xbar_mc_bandwidth; b++) {
      add_link(xbar, mc, _xbar_mc_latency, ACCEL_LINK_XBAR_MC);
      add_link(mc, xbar, _xbar_mc_latency, ACCEL_LINK_XBAR_MC);
    }
  }

  // 3) MC ↔ HBM (local, 1-to-1)
  for (int h = 0; h < K; h++) {
    int mc  = accel_mc_id(h);
    int hbm = accel_hbm_id(h);
    for (int b = 0; b < _mc_hbm_bandwidth; b++) {
      add_link(mc, hbm, _mc_hbm_latency, ACCEL_LINK_MC_HBM);
      add_link(hbm, mc, _mc_hbm_latency, ACCEL_LINK_MC_HBM);
    }
  }

  // 4) Xbar ↔ Xbar (inter-partition)
  for (int b = 0; b < _xbar_xbar_bandwidth; b++) {
    add_link(0, 1, _xbar_xbar_latency, ACCEL_LINK_XBAR_XBAR);
    add_link(1, 0, _xbar_xbar_latency, ACCEL_LINK_XBAR_XBAR);
  }

  // 5) MC ↔ MC (fabric, vertical only)
  if (_is_fabric) {
    int rows_per_col = K / 2;
    for (int col = 0; col < 2; col++) {
      for (int row = 0; row < rows_per_col - 1; row++) {
        int m1 = col * rows_per_col + row;
        int m2 = col * rows_per_col + row + 1;
        int mc1 = accel_mc_id(m1);
        int mc2 = accel_mc_id(m2);
        for (int b = 0; b < _mc_mc_bandwidth; b++) {
          add_link(mc1, mc2, _mc_mc_latency, ACCEL_LINK_MC_MC);
          add_link(mc2, mc1, _mc_mc_latency, ACCEL_LINK_MC_MC);
        }
      }
    }
  }

  assert(chan_idx == _channels);
}

// ============================================================
//  BFS helper — compute shortest distances on directed adjacency
// ============================================================
static void accel_bfs_distances(int size,
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
//  Builds three distance tables using directed BFS:
//    Hit:          links valid for hit flits
//    Miss baseline: links valid for miss flits, no MC-MC fabric
//    Miss fabric:   links valid for miss flits, with MC-MC fabric
// ============================================================
void HBMNetAccelSim::_BuildRoutingTable()
{
  int N = _size;

  // Build directed adjacencies based on link-type constraints
  vector<vector<pair<int,int>>> hitAdj(N);
  vector<vector<pair<int,int>>> missBaseAdj(N);
  vector<vector<pair<int,int>>> missFabAdj(N);

  for (int r = 0; r < N; r++) {
    for (size_t i = 0; i < gHBMNetAccelAdj[r].size(); i++) {
      const AccelAdjEntry &e = gHBMNetAccelAdj[r][i];

      if (accel_valid_hit_out(r, e.type))
        hitAdj[r].push_back(make_pair(e.neighbor, e.port));

      // Miss routing decisions use only Xbar+MC network (no MC→HBM)
      if (accel_valid_miss_routing_out(r, e.type)) {
        // Miss baseline: exclude MC-MC fabric
        if (e.type != ACCEL_LINK_MC_MC)
          missBaseAdj[r].push_back(make_pair(e.neighbor, e.port));
        // Miss fabric: include all valid miss routing links
        missFabAdj[r].push_back(make_pair(e.neighbor, e.port));
      }
    }
  }

  accel_bfs_distances(N, hitAdj, gHBMNetAccelDistHit);
  accel_bfs_distances(N, missBaseAdj, gHBMNetAccelDistMissBaseline);
  accel_bfs_distances(N, missFabAdj, gHBMNetAccelDistMissFabric);

  // Debug output
  auto print_table = [&](const string &title, const vector<vector<int>> &dist) {
    cout << "=== " << title << " ===" << endl;
    for (int s = 0; s < N; s++) {
      for (int d = 0; d < N; d++)
        cout << dist[s][d] << "\t";
      cout << endl;
    }
  };
  print_table("Hit Distance Table", gHBMNetAccelDistHit);
  print_table("Miss Baseline Distance Table", gHBMNetAccelDistMissBaseline);
  print_table("Miss Fabric Distance Table", gHBMNetAccelDistMissFabric);

  // Print adjacency
  cout << "=== AccelSim Adjacency ===" << endl;
  for (int r = 0; r < N; r++) {
    cout << "  Router " << r << " (";
    if (accel_is_xbar(r))
      cout << "Xbar" << r;
    else if (accel_is_mc(r))
      cout << "MC" << accel_mc_idx(r)
           << " c" << accel_hbm_col(accel_mc_idx(r))
           << "r" << accel_hbm_row(accel_mc_idx(r));
    else
      cout << "HBM" << accel_hbm_idx(r)
           << " c" << accel_hbm_col(accel_hbm_idx(r))
           << "r" << accel_hbm_row(accel_hbm_idx(r));
    cout << ", conc=" << gHBMNetAccelRouterConc[r] << "): ";
    map<int, int> nb_count;
    for (size_t i = 0; i < gHBMNetAccelAdj[r].size(); i++)
      nb_count[gHBMNetAccelAdj[r][i].neighbor]++;
    for (map<int,int>::iterator it = nb_count.begin(); it != nb_count.end(); ++it) {
      int nb = it->first;
      cout << "R" << nb << "(";
      if (accel_is_xbar(nb)) cout << "Xbar" << nb;
      else if (accel_is_mc(nb)) cout << "MC" << accel_mc_idx(nb);
      else cout << "HBM" << accel_hbm_idx(nb);
      cout << ")x" << it->second << " ";
    }
    cout << endl;
  }

  // L2 mapping
  cout << "=== L2 Slice Mapping (AccelSim) ===" << endl;
  cout << "  Mode: " << (gAccelInterleave ? "interleaved" : "sequential") << endl;
  for (int h = 0; h < gHBMNetAccelK; h++) {
    int rid = accel_hbm_id(h);
    cout << "  HBM" << h << " (Router " << rid
         << ", partition " << accel_hbm_part(h) << "): nodes ["
         << gHBMNetAccelRouterFirstNode[rid] << ", "
         << gHBMNetAccelRouterFirstNode[rid] + gHBMNetAccelRouterConc[rid] - 1 << "]"
         << endl;
  }
}

// ============================================================
//  RegisterRoutingFunctions
// ============================================================
void HBMNetAccelSim::RegisterRoutingFunctions()
{
  gRoutingFunctionMap["baseline_hbmnet_accelsim"]      = &hbmnet_accelsim_baseline;
  gRoutingFunctionMap["min_adaptive_hbmnet_accelsim"]  = &hbmnet_accelsim_min_adaptive;
  gRoutingFunctionMap["ugal_hbmnet_accelsim"]          = &hbmnet_accelsim_ugal;
  gRoutingFunctionMap["valiant_hbmnet_accelsim"]       = &hbmnet_accelsim_valiant;
  gRoutingFunctionMap["min_oblivious_hbmnet_accelsim"] = &hbmnet_accelsim_min_oblivious;
  gRoutingFunctionMap["hybrid_hbmnet_accelsim"]        = &hbmnet_accelsim_hybrid;
}

// ============================================================
//  Routing helpers
// ============================================================
static inline int accel_eject_port(int dest_node) {
  return hbmnet_accelsim_node_to_port(dest_node);
}

// Hit/miss markers stored in f->intm
static const int ACCELSIM_HIT_MARKER  = -100;
static const int ACCELSIM_MISS_MARKER = -200;

// ============================================================
//  Get the miss-routing target for a destination node.
//  SM→L2 (dest at HBM): target is paired MC router
//  L2→SM (dest at Xbar): target is dest Xbar router
// ============================================================
static inline int accel_miss_target(int dest_node) {
  int dest_router = hbmnet_accelsim_node_to_router(dest_node);
  if (accel_is_hbm(dest_router))
    return accel_mc_id(accel_hbm_idx(dest_router));
  return dest_router;  // Xbar
}

// ============================================================
//  Hit baseline next-hop (deterministic tree)
//  Xbar → HBM (direct) or Xbar → Xbar → HBM (cross-partition)
//  HBM → Xbar (reply traffic)
// ============================================================
static int accel_get_hit_baseline_next(int cur, int dest) {
  if (accel_is_xbar(cur)) {
    if (accel_is_xbar(dest)) return dest;
    if (accel_is_hbm(dest)) {
      int h = accel_hbm_idx(dest);
      if (accel_hbm_part(h) == cur) return dest;
      return accel_crossbar_id(1 - cur);
    }
    return accel_crossbar_id(1 - cur);
  }
  if (accel_is_hbm(cur)) {
    int h = accel_hbm_idx(cur);
    return accel_crossbar_id(accel_hbm_part(h));
  }
  assert(false);
  return -1;
}

// ============================================================
//  Miss baseline next-hop (deterministic, for escape VC)
//
//  Operates on the Xbar+MC logical network + deterministic end-hops.
//  miss_target is MC (SM→L2) or Xbar (L2→SM).
//
//  At HBM:     → local MC (deterministic first hop for L2 miss)
//  At Xbar:    → toward miss_target (direct or cross-partition)
//  At dest MC: → local HBM (deterministic final hop for SM→L2)
//  At intm MC: → fabric (same column) or back to partition Xbar
// ============================================================
static int accel_get_miss_baseline_next(int cur, int miss_target) {
  if (accel_is_hbm(cur)) {
    // HBM → local MC (deterministic first hop)
    return accel_mc_id(accel_hbm_idx(cur));
  }

  if (accel_is_xbar(cur)) {
    if (accel_is_mc(miss_target)) {
      int m = accel_mc_idx(miss_target);
      if (accel_hbm_part(m) == cur) return miss_target;
      return accel_crossbar_id(1 - cur);
    }
    // target is Xbar (L2→SM): cross-partition
    assert(cur != miss_target);  // should have been caught by eject
    return miss_target;
  }

  if (accel_is_mc(cur)) {
    if (accel_is_mc(miss_target) && cur == miss_target) {
      // Arrived at dest MC → deterministic to local HBM
      return accel_hbm_id(accel_mc_idx(cur));
    }
    // Baseline escape path always goes through partition Xbar
    // (baseline distance table excludes MC-MC fabric links)
    return accel_crossbar_id(accel_hbm_part(accel_mc_idx(cur)));
  }

  assert(false);
  return -1;
}

// ============================================================
//  Decide hit/miss at injection (hybrid routing only)
//  Both SM and L2 nodes use probabilistic hit/miss based on hit rate
// ============================================================
static void accel_decide_hit_miss(const Flit *f) {
  if (RandomInt(9999) < (int)(gAccelHitRate * 10000.0))
    f->intm = ACCELSIM_HIT_MARKER;
  else
    f->intm = ACCELSIM_MISS_MARKER;
}

// ============================================================
//  Collect minimal ports toward target
//  Miss routing uses accel_valid_miss_routing_out (Xbar+MC only)
// ============================================================
static void accel_collect_minimal_ports(int cur, int target,
                                        bool is_hit,
                                        const vector<vector<int>> &dist,
                                        vector<int> &ports)
{
  ports.clear();
  int cur_dist = dist[cur][target];
  if (cur_dist <= 0) return;
  for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
    if (is_hit ? accel_valid_hit_out(cur, e.type)
               : accel_valid_miss_routing_out(cur, e.type)) {
      if (dist[e.neighbor][target] == cur_dist - 1)
        ports.push_back(e.port);
    }
  }
}

// ============================================================
//  Collect minimal directions (grouped by next-hop router)
//  Also computes per-direction credit usage for adaptive routing
// ============================================================
static void accel_collect_minimal_dirs(const Router *r, int cur, int target,
                                       bool is_hit,
                                       const vector<vector<int>> &dist,
                                       map<int, vector<int>> &dir_ports,
                                       map<int, int> &dir_credit)
{
  dir_ports.clear();
  dir_credit.clear();
  int cur_dist = dist[cur][target];
  if (cur_dist <= 0) return;
  for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
    if (is_hit ? accel_valid_hit_out(cur, e.type)
               : accel_valid_miss_routing_out(cur, e.type)) {
      if (dist[e.neighbor][target] == cur_dist - 1) {
        dir_ports[e.neighbor].push_back(e.port);
        dir_credit[e.neighbor] += r->GetUsedCredit(e.port);
      }
    }
  }
}

// ============================================================
//  Average credit per port toward target (for UGAL cost)
// ============================================================
static double accel_avg_credit_toward(const Router *r, int cur, int target,
                                      bool is_hit,
                                      const vector<vector<int>> &dist)
{
  int cur_dist = dist[cur][target];
  if (cur_dist <= 0) return 0.0;
  int total_credit = 0, count = 0;
  for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
    if (is_hit ? accel_valid_hit_out(cur, e.type)
               : accel_valid_miss_routing_out(cur, e.type)) {
      if (dist[e.neighbor][target] == cur_dist - 1) {
        total_credit += r->GetUsedCredit(e.port);
        count++;
      }
    }
  }
  if (count == 0) return 0.0;
  return (double)total_credit / (double)count;
}

// ============================================================
//  Pick least-congested direction, random port within
// ============================================================
static int accel_pick_adaptive_port(const Router *r, int cur, int target,
                                    bool is_hit,
                                    const vector<vector<int>> &dist)
{
  map<int, vector<int>> dir_ports;
  map<int, int> dir_credit;
  accel_collect_minimal_dirs(r, cur, target, is_hit, dist, dir_ports, dir_credit);
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
  return ports[RandomInt(ports.size() - 1)];
}

// ============================================================
//  UGAL decision helper (miss path)
//  miss_target: MC (SM→L2) or Xbar (L2→SM)
// ============================================================
static void accel_ugal_decision(const Router *r, const Flit *f,
                                int cur, int miss_target) {
  int K = gHBMNetAccelK;
  int min_dist = gHBMNetAccelDistMissFabric[cur][miss_target];
  double min_avg = accel_avg_credit_toward(r, cur, miss_target, false,
                                           gHBMNetAccelDistMissFabric);

  vector<int> candidates;
  for (int m = 0; m < K; m++) {
    int mc = accel_mc_id(m);
    if (mc != cur && mc != miss_target
        && gHBMNetAccelDistMissFabric[cur][mc] > 0)
      candidates.push_back(mc);
  }

  if (candidates.empty()) {
    f->intm = ACCELSIM_MISS_MARKER;
    f->ph = 1;
    return;
  }

  int intm;
  if (gAccelUGALIntmSelect == 0) {
    intm = candidates[RandomInt(candidates.size() - 1)];
  } else {
    double best_cost = numeric_limits<double>::max();
    intm = candidates[0];
    for (size_t c = 0; c < candidates.size(); c++) {
      int mc = candidates[c];
      int d = gHBMNetAccelDistMissFabric[cur][mc]
            + gHBMNetAccelDistMissFabric[mc][miss_target];
      double cr = accel_avg_credit_toward(r, cur, mc, false,
                                          gHBMNetAccelDistMissFabric);
      if (d * cr < best_cost) {
        best_cost = d * cr;
        intm = mc;
      }
    }
  }

  int nonmin_dist = gHBMNetAccelDistMissFabric[cur][intm]
                  + gHBMNetAccelDistMissFabric[intm][miss_target];
  double nonmin_avg = accel_avg_credit_toward(r, cur, intm, false,
                                              gHBMNetAccelDistMissFabric);

  if (min_dist * min_avg <= nonmin_dist * nonmin_avg + gAccelUGALThreshold) {
    f->intm = ACCELSIM_MISS_MARKER;
    f->ph = 1;
    ++gAccelUGALMinDecisions;
  } else {
    f->intm = intm;
    f->ph = 0;
    ++gAccelUGALNonMinDecisions;
  }
}

// ============================================================
//  Valiant decision helper (miss path)
//  miss_target: MC (SM→L2) or Xbar (L2→SM)
// ============================================================
static void accel_valiant_decision(const Flit *f, int cur, int miss_target) {
  int K = gHBMNetAccelK;
  vector<int> candidates;
  for (int m = 0; m < K; m++) {
    int mc = accel_mc_id(m);
    if (mc != cur && mc != miss_target
        && gHBMNetAccelDistMissFabric[cur][mc] > 0)
      candidates.push_back(mc);
  }
  if (candidates.empty()) {
    f->intm = ACCELSIM_MISS_MARKER;
    f->ph = 1;
  } else {
    f->intm = candidates[RandomInt(candidates.size() - 1)];
    f->ph = 0;
  }
}

// ============================================================
//  Deterministic MC↔HBM hop helper
//  Used for: HBM→MC (L2 miss first hop) and MC→HBM (miss final hop)
// ============================================================
static void accel_add_mc_hbm_ports(int cur, OutputSet *outputs) {
  vector<int> ports;
  for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
    if (e.type == ACCEL_LINK_MC_HBM)
      ports.push_back(e.port);
  }
  assert(!ports.empty());
  outputs->AddRange(ports[RandomInt(ports.size() - 1)], 0, gNumVCs - 1);
}

// ============================================================
//  Common routing preamble (non-inject)
//
//  1. Eject at destination router
//  2. Deterministic HBM→MC / MC→HBM hops for miss flits
//  3. Escape VC (VC 0, Duato's protocol)
//
//  Returns true if routing is complete (eject or deterministic hop).
//  Sets miss_target for caller's use in data VC routing.
// ============================================================
static bool accelsim_routing_preamble(const Router *r, const Flit *f,
                                      OutputSet *outputs,
                                      int &cur_router, int &dest_router,
                                      int &miss_target)
{
  outputs->Clear();

  cur_router  = r->GetID();
  dest_router = hbmnet_accelsim_node_to_router(f->dest);

  // Eject at destination router
  if (cur_router == dest_router) {
    ++gAccelTotalEjects;
    if (f->vc == 0) ++gAccelEscapeVCEjects;
    outputs->AddRange(accel_eject_port(f->dest), 0, gNumVCs - 1);
    return true;
  }

  bool is_hit = (f->intm == ACCELSIM_HIT_MARKER);
  miss_target = is_hit ? -1 : accel_miss_target(f->dest);

  if (!is_hit) {
    // Deterministic HBM→MC hop (L2 miss injection first hop)
    if (accel_is_hbm(cur_router)) {
      accel_add_mc_hbm_ports(cur_router, outputs);
      return true;
    }
    // Deterministic MC→HBM hop (at dest MC, SM→L2 miss final hop)
    if (accel_is_mc(cur_router) && accel_is_mc(miss_target)
        && cur_router == miss_target) {
      accel_add_mc_hbm_ports(cur_router, outputs);
      return true;
    }
  }

  assert(gNumVCs >= 2);

  // Escape VC (VC 0): deterministic baseline routing
  int esc_next;
  if (is_hit)
    esc_next = accel_get_hit_baseline_next(cur_router, dest_router);
  else
    esc_next = accel_get_miss_baseline_next(cur_router, miss_target);

  vector<int> esc_ports;
  for (size_t i = 0; i < gHBMNetAccelAdj[cur_router].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur_router][i];
    if (e.neighbor == esc_next) {
      if (is_hit ? accel_valid_hit_out(cur_router, e.type)
                 : accel_valid_miss_out(cur_router, e.type))
        esc_ports.push_back(e.port);
    }
  }
  assert(!esc_ports.empty());
  outputs->AddRange(esc_ports[RandomInt(esc_ports.size() - 1)], 0, 0, 0);

  return false;
}

// ============================================================
//  Helper: add hit data ports (baseline tree, VCs 1..V-1)
// ============================================================
static void accel_add_hit_data_ports(int cur, int dest, OutputSet *outputs) {
  int next = accel_get_hit_baseline_next(cur, dest);
  vector<int> ports;
  for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
    if (e.neighbor == next && accel_valid_hit_out(cur, e.type))
      ports.push_back(e.port);
  }
  assert(!ports.empty());
  outputs->AddRange(ports[RandomInt(ports.size() - 1)], 1, gNumVCs - 1, 1);
}

// ============================================================
//  BASELINE ROUTING (non-hybrid: all miss)
//
//  All traffic takes the miss path: Xbar→MC→HBM or HBM→MC→Xbar
//  Deterministic baseline routing on Xbar+MC network
// ============================================================
void hbmnet_accelsim_baseline( const Router *r, const Flit *f, int in_channel,
                                OutputSet *outputs, bool inject )
{
  if (inject) {
    outputs->Clear();
    f->intm = ACCELSIM_MISS_MARKER;
    outputs->AddRange(-1, 1, gNumVCs - 1);
    return;
  }

  int cur, dest, miss_target;
  if (accelsim_routing_preamble(r, f, outputs, cur, dest, miss_target)) return;

  int next = accel_get_miss_baseline_next(cur, miss_target);
  vector<int> ports;
  for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
    if (e.neighbor == next && accel_valid_miss_out(cur, e.type))
      ports.push_back(e.port);
  }
  assert(!ports.empty());
  outputs->AddRange(ports[RandomInt(ports.size() - 1)], 1, gNumVCs - 1, 1);
}

// ============================================================
//  MINIMAL OBLIVIOUS ROUTING (non-hybrid: all miss)
//
//  Random minimal port using miss fabric distance table
// ============================================================
void hbmnet_accelsim_min_oblivious( const Router *r, const Flit *f, int in_channel,
                                     OutputSet *outputs, bool inject )
{
  if (inject) {
    outputs->Clear();
    f->intm = ACCELSIM_MISS_MARKER;
    outputs->AddRange(-1, 1, gNumVCs - 1);
    return;
  }

  int cur, dest, miss_target;
  if (accelsim_routing_preamble(r, f, outputs, cur, dest, miss_target)) return;

  vector<int> ports;
  accel_collect_minimal_ports(cur, miss_target, false,
                              gHBMNetAccelDistMissFabric, ports);
  assert(!ports.empty());
  outputs->AddRange(ports[RandomInt(ports.size() - 1)], 1, gNumVCs - 1, 1);
}

// ============================================================
//  MINIMAL ADAPTIVE ROUTING (non-hybrid: all miss)
//
//  Least-congested minimal direction, random port within
// ============================================================
void hbmnet_accelsim_min_adaptive( const Router *r, const Flit *f, int in_channel,
                                    OutputSet *outputs, bool inject )
{
  if (inject) {
    outputs->Clear();
    f->intm = ACCELSIM_MISS_MARKER;
    outputs->AddRange(-1, 1, gNumVCs - 1);
    return;
  }

  int cur, dest, miss_target;
  if (accelsim_routing_preamble(r, f, outputs, cur, dest, miss_target)) return;

  int port = accel_pick_adaptive_port(r, cur, miss_target, false,
                                      gHBMNetAccelDistMissFabric);
  outputs->AddRange(port, 1, gNumVCs - 1, 1);
}

// ============================================================
//  UGAL ROUTING (non-hybrid: all miss)
//
//  UGAL among MC/Xbar routers on the Xbar+MC network
//  Decision deferred to first routing hop (injection has no router)
// ============================================================
void hbmnet_accelsim_ugal( const Router *r, const Flit *f, int in_channel,
                            OutputSet *outputs, bool inject )
{
  if (inject) {
    outputs->Clear();
    f->intm = ACCELSIM_MISS_MARKER;
    f->ph = 0;
    outputs->AddRange(-1, 1, gNumVCs - 1);
    return;
  }

  int cur, dest, miss_target;
  if (accelsim_routing_preamble(r, f, outputs, cur, dest, miss_target)) return;

  // UGAL decision at first routing hop (f->intm == MISS_MARKER && f->ph == 0)
  bool has_intm = (f->intm != ACCELSIM_MISS_MARKER && f->intm >= 0);
  if (!has_intm && f->intm == ACCELSIM_MISS_MARKER && f->ph == 0) {
    accel_ugal_decision(r, f, cur, miss_target);
    has_intm = (f->intm != ACCELSIM_MISS_MARKER && f->intm >= 0);
  }

  // Phase transition: arrived at intermediate
  if (has_intm && f->ph == 0 && cur == f->intm)
    f->ph = 1;

  int target = (has_intm && f->ph == 0) ? f->intm : miss_target;
  int port = accel_pick_adaptive_port(r, cur, target, false,
                                      gHBMNetAccelDistMissFabric);
  outputs->AddRange(port, 1, gNumVCs - 1, 1);
}

// ============================================================
//  VALIANT ROUTING (non-hybrid: all miss)
//
//  Always route via random intermediate MC
//  Phase 0: minimal toward intermediate
//  Phase 1: minimal toward miss_target
// ============================================================
void hbmnet_accelsim_valiant( const Router *r, const Flit *f, int in_channel,
                               OutputSet *outputs, bool inject )
{
  if (inject) {
    outputs->Clear();
    f->intm = ACCELSIM_MISS_MARKER;
    f->ph = 0;
    outputs->AddRange(-1, 1, gNumVCs - 1);
    return;
  }

  int cur, dest, miss_target;
  if (accelsim_routing_preamble(r, f, outputs, cur, dest, miss_target)) return;

  // Valiant decision at first routing hop
  bool has_intm = (f->intm != ACCELSIM_MISS_MARKER && f->intm >= 0);
  if (!has_intm && f->intm == ACCELSIM_MISS_MARKER && f->ph == 0) {
    accel_valiant_decision(f, cur, miss_target);
    has_intm = (f->intm != ACCELSIM_MISS_MARKER && f->intm >= 0);
  }

  // Phase transition
  if (has_intm && f->ph == 0 && cur == f->intm)
    f->ph = 1;

  int target = (has_intm && f->ph == 0) ? f->intm : miss_target;
  vector<int> ports;
  accel_collect_minimal_ports(cur, target, false,
                              gHBMNetAccelDistMissFabric, ports);
  assert(!ports.empty());
  outputs->AddRange(ports[RandomInt(ports.size() - 1)], 1, gNumVCs - 1, 1);
}

// ============================================================
//  HYBRID ROUTING
//
//  Only routing function that uses hit/miss distinction.
//  At injection: probabilistic hit/miss via hit_rate
//  Hit flits: deterministic baseline tree on Xbar↔HBM links
//  Miss flits: sub-routing function selected by hybrid_routing config
//    0=baseline, 1=min_adaptive, 2=ugal, 3=valiant, 4=min_oblivious
// ============================================================
void hbmnet_accelsim_hybrid( const Router *r, const Flit *f, int in_channel,
                              OutputSet *outputs, bool inject )
{
  if (inject) {
    outputs->Clear();
    accel_decide_hit_miss(f);
    f->ph = 0;
    outputs->AddRange(-1, 1, gNumVCs - 1);
    return;
  }

  int cur, dest, miss_target;
  if (accelsim_routing_preamble(r, f, outputs, cur, dest, miss_target)) return;

  // Hit flits: deterministic baseline tree
  if (f->intm == ACCELSIM_HIT_MARKER) {
    accel_add_hit_data_ports(cur, dest, outputs);
    return;
  }

  // Miss flits: UGAL/Valiant decision at first routing hop
  bool has_intm = (f->intm != ACCELSIM_MISS_MARKER && f->intm >= 0);
  if (!has_intm && f->intm == ACCELSIM_MISS_MARKER && f->ph == 0) {
    switch (gAccelHybridRouting) {
      case 2:  accel_ugal_decision(r, f, cur, miss_target); break;
      case 3:  accel_valiant_decision(f, cur, miss_target); break;
      default: f->ph = 1; break;
    }
    has_intm = (f->intm != ACCELSIM_MISS_MARKER && f->intm >= 0);
  }

  // Phase transition for UGAL/Valiant intermediates
  if (has_intm && f->ph == 0 && cur == f->intm)
    f->ph = 1;

  int target = (has_intm && f->ph == 0) ? f->intm : miss_target;

  switch (gAccelHybridRouting) {
    case 1:  // min_adaptive
    case 2:  // ugal (per-hop is adaptive)
    {
      int port = accel_pick_adaptive_port(r, cur, target, false,
                                          gHBMNetAccelDistMissFabric);
      outputs->AddRange(port, 1, gNumVCs - 1, 1);
      break;
    }
    case 3:  // valiant (per-hop is random minimal)
    case 4:  // min_oblivious
    {
      vector<int> ports;
      accel_collect_minimal_ports(cur, target, false,
                                  gHBMNetAccelDistMissFabric, ports);
      assert(!ports.empty());
      outputs->AddRange(ports[RandomInt(ports.size() - 1)], 1, gNumVCs - 1, 1);
      break;
    }
    default: // baseline
    {
      int next = accel_get_miss_baseline_next(cur, miss_target);
      vector<int> ports;
      for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
        const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
        if (e.neighbor == next && accel_valid_miss_out(cur, e.type))
          ports.push_back(e.port);
      }
      assert(!ports.empty());
      outputs->AddRange(ports[RandomInt(ports.size() - 1)], 1, gNumVCs - 1, 1);
      break;
    }
  }
}
