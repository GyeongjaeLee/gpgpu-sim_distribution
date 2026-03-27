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
int gHBMNetAccelP = 0;
int gHBMNetAccelHPS = 0;

// Reverse mapping: node -> router, node -> port
static vector<int> gAccelNodeToRouter;
static vector<int> gAccelNodeToPort;

// Reverse input port map: [router][input_port] -> source router (-1 for inject ports)
static vector<vector<int>> gAccelInputPortSrc;

// Link latency cache indexed by AccelLinkType
static int gAccelLinkLatency[5] = {0, 0, 0, 0, 0};

// Near-min adaptive config
static int gAccelNearMinK = 1;
static double gAccelNearMinPenalty = 1.0;  // penalty multiplier for near-min cost

// Topology parameters cached for routing
static int gAccelNumSMs;
static int gAccelNumL2;
static int gAccelL;          // l2_per_hbm
static int gAccelInterleave;
static int gAccelUGALThreshold;
static int gAccelUGALIntmSelect;
static double gAccelHitRate;     // baseline_ratio repurposed as L2 hit rate
static int gAccelHybridRouting;  // 0=baseline,1=min_adaptive,2=ugal,3=valiant,4=min_oblivious,5=near_min_adaptive

// UGAL routing decision counters
long long gAccelUGALMinDecisions    = 0;
long long gAccelUGALNonMinDecisions = 0;

// Escape VC usage counters
long long gAccelEscapeVCEjects = 0;
long long gAccelTotalEjects    = 0;

// Near-min adaptive decision counters
long long gAccelNearMinMinDecisions    = 0;
long long gAccelNearMinNonMinDecisions = 0;

// Near-min path-level tracking
long long gAccelNearMinPathsUsed  = 0;  // packets that took near-min at least once
long long gAccelTotalMissPackets  = 0;  // total miss packets

// Near-min direction tracking: [src_router][dst_router] -> count of near-min decisions
static map<pair<int,int>, long long> gAccelNearMinDirCount;

// Per-link-type traversal counters (flit traversals, incremented at routing time)
// Index by AccelLinkType: 0=XBAR_XBAR, 1=XBAR_HBM, 2=XBAR_MC, 3=MC_HBM, 4=MC_MC
long long gAccelLinkTypeTraversals[5] = {0, 0, 0, 0, 0};
// Per-router-pair directional traversals for detailed analysis
static map<pair<int,int>, long long> gAccelLinkDirTraversals;
// Per-router-pair cumulative saturation (for computing average at print time)
static map<pair<int,int>, double> gAccelLinkDirSatSum;

// Cached buffer capacity per output port: num_vcs * vc_buf_size
static int gAccelPortBufCapacity = 0;

void hbmnet_accelsim_reset_stats() {
  gAccelUGALMinDecisions    = 0;
  gAccelUGALNonMinDecisions = 0;
  gAccelEscapeVCEjects      = 0;
  gAccelTotalEjects         = 0;
  gAccelNearMinMinDecisions    = 0;
  gAccelNearMinNonMinDecisions = 0;
  gAccelNearMinPathsUsed    = 0;
  gAccelTotalMissPackets    = 0;
  gAccelNearMinDirCount.clear();
  for (int i = 0; i < 5; i++) gAccelLinkTypeTraversals[i] = 0;
  gAccelLinkDirTraversals.clear();
  gAccelLinkDirSatSum.clear();
}

int hbmnet_accelsim_node_to_router(int node) { return gAccelNodeToRouter[node]; }
int hbmnet_accelsim_node_to_port(int node)   { return gAccelNodeToPort[node]; }

// Helper to get router name string (forward-declared, uses ID helpers below)
static string accel_router_name(int r);

// ============================================================
//  Router ID helpers
//    Router 0..P-1       : Crossbar 0..P-1
//    Router P..P+K-1     : MC 0..K-1
//    Router P+K..P+2K-1  : HBM 0..K-1
// ============================================================
static inline int accel_crossbar_id(int p) { return p; }
static inline int accel_mc_id(int mc_idx)  { return mc_idx + gHBMNetAccelP; }
static inline int accel_hbm_id(int hbm_idx){ return hbm_idx + gHBMNetAccelP + gHBMNetAccelK; }
static inline bool accel_is_xbar(int r)    { return r < gHBMNetAccelP; }
static inline bool accel_is_mc(int r)      { return r >= gHBMNetAccelP && r < gHBMNetAccelP + gHBMNetAccelK; }
static inline bool accel_is_hbm(int r)     { return r >= gHBMNetAccelP + gHBMNetAccelK; }
static inline int accel_mc_idx(int r)      { return r - gHBMNetAccelP; }
static inline int accel_hbm_idx(int r)     { return r - gHBMNetAccelP - gHBMNetAccelK; }

// 2D layout helpers
//   Column 0: indices 0..K/2-1,  Column 1: indices K/2..K-1
//   Each Xbar p owns rows p*HPS..(p+1)*HPS-1 in each column
static inline int accel_hbm_col(int h)  { return h / (gHBMNetAccelK / 2); }
static inline int accel_hbm_row(int h)  { return h % (gHBMNetAccelK / 2); }
static inline int accel_hbm_part(int h) {
  return accel_hbm_row(h) / gHBMNetAccelHPS;
}

// Router name helper for stats output
static string accel_router_name(int r) {
  ostringstream oss;
  if (accel_is_xbar(r))
    oss << "Xbar" << r;
  else if (accel_is_mc(r)) {
    int m = accel_mc_idx(r);
    oss << "MC" << m << "(c" << accel_hbm_col(m) << "r" << accel_hbm_row(m) << ")";
  } else {
    int h = accel_hbm_idx(r);
    oss << "HBM" << h << "(c" << accel_hbm_col(h) << "r" << accel_hbm_row(h) << ")";
  }
  return oss.str();
}

// ============================================================
//  Print link utilization and near-min direction stats
// ============================================================
void hbmnet_accelsim_print_link_stats() {
  static const char* link_type_names[] = {
    "XBAR_XBAR", "XBAR_HBM", "XBAR_MC", "MC_HBM", "MC_MC"
  };

  long long total = 0;
  for (int i = 0; i < 5; i++) total += gAccelLinkTypeTraversals[i];

  cout << "=== Link Type Utilization ===" << endl;
  for (int i = 0; i < 5; i++) {
    double pct = (total > 0) ? 100.0 * (double)gAccelLinkTypeTraversals[i] / (double)total : 0.0;
    cout << "  " << link_type_names[i] << ": " << gAccelLinkTypeTraversals[i]
         << " (" << pct << "%)" << endl;
  }
  cout << "  Total link traversals: " << total << endl;

  // Per-direction traversals (grouped by link type) with saturation
  // Also accumulate per-type aggregate saturation (weighted average)
  double type_sat_sum[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  long long type_count[5] = {0, 0, 0, 0, 0};

  // First pass: accumulate per-type stats
  for (map<pair<int,int>, long long>::iterator it = gAccelLinkDirTraversals.begin();
       it != gAccelLinkDirTraversals.end(); ++it) {
    int src = it->first.first;
    int dst = it->first.second;
    AccelLinkType lt = ACCEL_LINK_XBAR_XBAR;
    for (size_t j = 0; j < gHBMNetAccelAdj[src].size(); j++) {
      if (gHBMNetAccelAdj[src][j].neighbor == dst) {
        lt = gHBMNetAccelAdj[src][j].type;
        break;
      }
    }
    long long count = it->second;
    type_count[(int)lt] += count;
    map<pair<int,int>, double>::iterator sit = gAccelLinkDirSatSum.find(it->first);
    if (sit != gAccelLinkDirSatSum.end())
      type_sat_sum[(int)lt] += sit->second;
  }

  // Print per-direction detail
  cout << "=== Per-Direction Miss Link Traversals ===" << endl;
  for (int t = 0; t < 5; t++) {
    // Skip XBAR_HBM (hit path only) and MC_HBM (deterministic, not adaptive)
    if (t == ACCEL_LINK_XBAR_HBM || t == ACCEL_LINK_MC_HBM) continue;
    bool header_printed = false;
    for (map<pair<int,int>, long long>::iterator it = gAccelLinkDirTraversals.begin();
         it != gAccelLinkDirTraversals.end(); ++it) {
      int src = it->first.first;
      int dst = it->first.second;
      AccelLinkType lt = ACCEL_LINK_XBAR_XBAR;
      for (size_t j = 0; j < gHBMNetAccelAdj[src].size(); j++) {
        if (gHBMNetAccelAdj[src][j].neighbor == dst) {
          lt = gHBMNetAccelAdj[src][j].type;
          break;
        }
      }
      if ((int)lt != t) continue;
      if (!header_printed) {
        double type_avg_sat = (type_count[t] > 0)
            ? type_sat_sum[t] / (double)type_count[t] * 100.0 : 0.0;
        cout << "  [" << link_type_names[t] << "] avg_sat=" << type_avg_sat << "%" << endl;
        header_printed = true;
      }
      long long count = it->second;
      double dir_pct = (total > 0) ? 100.0 * (double)count / (double)total : 0.0;
      double avg_sat = 0.0;
      map<pair<int,int>, double>::iterator sit = gAccelLinkDirSatSum.find(it->first);
      if (sit != gAccelLinkDirSatSum.end() && count > 0)
        avg_sat = sit->second / (double)count;
      cout << "    " << accel_router_name(src) << " -> " << accel_router_name(dst)
           << ": count=" << count
           << " (" << dir_pct << "%)"
           << " avg_sat=" << (avg_sat * 100.0) << "%" << endl;
    }
  }

  // Machine-parseable aggregate line for run_all_moe.py
  double sat_xx = (type_count[ACCEL_LINK_XBAR_XBAR] > 0)
      ? type_sat_sum[ACCEL_LINK_XBAR_XBAR] / (double)type_count[ACCEL_LINK_XBAR_XBAR] * 100.0 : 0.0;
  double sat_xmc = (type_count[ACCEL_LINK_XBAR_MC] > 0)
      ? type_sat_sum[ACCEL_LINK_XBAR_MC] / (double)type_count[ACCEL_LINK_XBAR_MC] * 100.0 : 0.0;
  double sat_mc = (type_count[ACCEL_LINK_MC_MC] > 0)
      ? type_sat_sum[ACCEL_LINK_MC_MC] / (double)type_count[ACCEL_LINK_MC_MC] * 100.0 : 0.0;
  cout << "Miss link avg saturation:"
       << " XBAR_XBAR=" << sat_xx << "%"
       << " XBAR_MC=" << sat_xmc << "%"
       << " MC_MC=" << sat_mc << "%" << endl;

  // Near-min adaptive decision stats
  long long nm_total = gAccelNearMinMinDecisions + gAccelNearMinNonMinDecisions;
  if (nm_total > 0) {
    double nm_pct = 100.0 * (double)gAccelNearMinNonMinDecisions / (double)nm_total;
    cout << "=== Near-Min Adaptive Decisions ===" << endl;
    cout << "  Near-min decisions: total=" << nm_total
         << " min=" << gAccelNearMinMinDecisions
         << " non-min=" << gAccelNearMinNonMinDecisions
         << " near-min ratio=" << nm_pct << "%" << endl;

    // Path-level near-min usage
    if (gAccelTotalMissPackets > 0) {
      double path_pct = 100.0 * (double)gAccelNearMinPathsUsed / (double)gAccelTotalMissPackets;
      cout << "  Near-min path usage: " << gAccelNearMinPathsUsed
           << " / " << gAccelTotalMissPackets
           << " (" << path_pct << "%)" << endl;
    }

    if (!gAccelNearMinDirCount.empty()) {
      cout << "  Near-min direction breakdown:" << endl;
      for (map<pair<int,int>, long long>::iterator it = gAccelNearMinDirCount.begin();
           it != gAccelNearMinDirCount.end(); ++it) {
        int src = it->first.first;
        int dst = it->first.second;
        AccelLinkType lt = ACCEL_LINK_XBAR_XBAR;
        for (size_t j = 0; j < gHBMNetAccelAdj[src].size(); j++) {
          if (gHBMNetAccelAdj[src][j].neighbor == dst) {
            lt = gHBMNetAccelAdj[src][j].type;
            break;
          }
        }
        // Also show avg saturation when near-min was chosen for this direction
        double avg_sat = 0.0;
        map<pair<int,int>, double>::iterator sit = gAccelLinkDirSatSum.find(it->first);
        long long dir_traversals = 0;
        map<pair<int,int>, long long>::iterator dit = gAccelLinkDirTraversals.find(it->first);
        if (dit != gAccelLinkDirTraversals.end()) dir_traversals = dit->second;
        if (sit != gAccelLinkDirSatSum.end() && dir_traversals > 0)
          avg_sat = sit->second / (double)dir_traversals;
        cout << "    " << accel_router_name(src) << " -> " << accel_router_name(dst)
             << " [" << link_type_names[(int)lt] << "]: " << it->second
             << " avg_sat=" << (avg_sat * 100.0) << "%" << endl;
      }
    }
  }
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
  _num_xbars      = config.GetInt("num_xbars");
  _hbm_per_side   = config.GetInt("hbm_per_side");
  _num_hbm_stacks = _num_xbars * _hbm_per_side * 2;
  _l2_per_hbm     = _num_l2_slices / _num_hbm_stacks;
  _l2_interleave  = config.GetInt("l2_interleave");
  _is_fabric      = config.GetInt("is_fabric");

  int P = _num_xbars;
  int K = _num_hbm_stacks;

  assert(P >= 1);
  assert(_hbm_per_side >= 1);
  assert(_num_sms % P == 0 && "num_sms must be a multiple of num_xbars");
  assert(_num_l2_slices % K == 0 && "num_l2_slices must be a multiple of num_hbm_stacks");

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
  gHBMNetAccelK      = K;
  gHBMNetAccelP      = P;
  gHBMNetAccelHPS    = _hbm_per_side;
  gAccelNumSMs       = _num_sms;
  gAccelNumL2        = _num_l2_slices;
  gAccelL            = _l2_per_hbm;
  gAccelInterleave   = _l2_interleave;
  gAccelUGALThreshold  = config.GetInt("ugal_threshold");
  gAccelUGALIntmSelect = config.GetInt("ugal_intm_select");
  gAccelHitRate        = config.GetFloat("baseline_ratio");  // repurposed as hit rate
  gAccelNearMinK       = config.GetInt("near_min_k");
  gAccelNearMinPenalty = config.GetFloat("near_min_penalty");

  // Cache link latencies by type for near-min penalty calculation
  gAccelLinkLatency[ACCEL_LINK_XBAR_XBAR] = _xbar_xbar_latency;
  gAccelLinkLatency[ACCEL_LINK_XBAR_HBM]  = _xbar_hbm_latency;
  gAccelLinkLatency[ACCEL_LINK_XBAR_MC]   = _xbar_mc_latency;
  gAccelLinkLatency[ACCEL_LINK_MC_HBM]    = _mc_hbm_latency;
  gAccelLinkLatency[ACCEL_LINK_MC_MC]     = _mc_mc_latency;

  // Cache buffer capacity for saturation calculation
  gAccelPortBufCapacity = config.GetInt("num_vcs") * config.GetInt("vc_buf_size");

  string hr = config.GetStr("hybrid_routing");
  if      (hr == "min_adaptive")       gAccelHybridRouting = 1;
  else if (hr == "ugal")               gAccelHybridRouting = 2;
  else if (hr == "valiant")            gAccelHybridRouting = 3;
  else if (hr == "min_oblivious")      gAccelHybridRouting = 4;
  else if (hr == "near_min_adaptive")  gAccelHybridRouting = 5;
  else                                 gAccelHybridRouting = 0;

  // Total routers: P Xbars + K MCs + K HBMs
  _size = P + 2 * K;

  // Per-router concentration and first node
  gHBMNetAccelRouterConc.resize(_size);
  gHBMNetAccelRouterFirstNode.resize(_size);

  // Xbar p: SM nodes [p*N/P .. (p+1)*N/P - 1]
  int sms_per_xbar = _num_sms / P;
  for (int p = 0; p < P; p++) {
    gHBMNetAccelRouterConc[p] = sms_per_xbar;
    gHBMNetAccelRouterFirstNode[p] = p * sms_per_xbar;
  }

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
  // Xbar chain: P-1 pairs, each pair has bandwidth * 2 (bidirectional)
  int xbar_xbar_ch = (P - 1) * _xbar_xbar_bandwidth * 2;
  int xbar_hbm_ch  = K * _xbar_hbm_bandwidth * 2;  // hit path
  int xbar_mc_ch   = K * _xbar_mc_bandwidth * 2;    // miss path entry
  int mc_hbm_ch    = K * _mc_hbm_bandwidth * 2;     // miss path local
  int mc_mc_ch     = 0;
  if (_is_fabric) {
    // Each column has K/2 rows, so K/2 - 1 adjacent pairs per column, 2 columns
    int vert_pairs = 2 * (K / 2 - 1);
    mc_mc_ch = vert_pairs * _mc_mc_bandwidth * 2;
  }
  _channels = xbar_xbar_ch + xbar_hbm_ch + xbar_mc_ch + mc_hbm_ch + mc_mc_ch;

  cout << "HBMNetAccelSim Config:"
       << " num_sms=" << _num_sms
       << " num_l2=" << _num_l2_slices
       << " P=" << P
       << " H=" << _hbm_per_side
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
       << " (" << P << " Xbar + " << K << " MC + " << K << " HBM)"
       << " nodes=" << _nodes
       << " channels=" << _channels << endl;
  double hbm_sp = config.GetFloat("hbm_internal_speedup");
  cout << "  internal_speedup=" << config.GetFloat("internal_speedup");
  if (hbm_sp > 0.0)
    cout << "  hbm_internal_speedup=" << hbm_sp << " (HBM routers only)";
  cout << endl;
}

// ============================================================
//  _BuildNet
// ============================================================
void HBMNetAccelSim::_BuildNet( const Configuration &config )
{
  ostringstream name;
  int chan_idx = 0;
  int K = _num_hbm_stacks;

  int P = _num_xbars;

  // Compute router degrees
  vector<int> degree(_size, 0);

  // Xbar routers: concentration + hit links + miss links + chain links
  for (int p = 0; p < P; p++) {
    degree[p] = gHBMNetAccelRouterConc[p];  // SM eject ports
    for (int h = 0; h < K; h++) {
      if (accel_hbm_part(h) == p) {
        degree[p] += _xbar_hbm_bandwidth;  // hit path to HBM h
        degree[p] += _xbar_mc_bandwidth;   // miss path to MC h
      }
    }
    // Chain links: left neighbor (p-1) and right neighbor (p+1)
    if (p > 0)     degree[p] += _xbar_xbar_bandwidth;
    if (p < P - 1) degree[p] += _xbar_xbar_bandwidth;
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

  // Optional per-type internal speedup override for HBM routers.
  // 0.0 (default) = use the global internal_speedup for all routers.
  double hbm_speedup = config.GetFloat("hbm_internal_speedup");

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

    if (accel_is_hbm(rtr) && hbm_speedup > 0.0) {
      // Give HBM routers a dedicated internal speedup so the crossbar can
      // service all HBM→MC output ports each cycle, making mc_hbm_bandwidth
      // the true bottleneck (TSV model) rather than the router crossbar.
      Configuration hbm_config = config;
      hbm_config.Assign("internal_speedup", hbm_speedup);
      _routers[rtr] = Router::NewRouter(hbm_config, this, name.str(), rtr,
                                         degree[rtr], degree[rtr]);
    } else {
      _routers[rtr] = Router::NewRouter(config, this, name.str(), rtr,
                                         degree[rtr], degree[rtr]);
    }
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

  // Reverse input port map: [router][input_port] -> source router
  // Inject ports (0..conc-1) have source = -1
  gAccelInputPortSrc.clear();
  gAccelInputPortSrc.resize(_size);
  for (int rtr = 0; rtr < _size; rtr++)
    gAccelInputPortSrc[rtr].assign(degree[rtr], -1);

  // Input port counters (network input ports start after inject ports)
  vector<int> ipc(_size);
  for (int rtr = 0; rtr < _size; rtr++)
    ipc[rtr] = gHBMNetAccelRouterConc[rtr];

  auto add_link = [&](int u, int v, int latency, AccelLinkType type) {
    assert(chan_idx < _channels);
    _chan[chan_idx]->SetLatency(latency);
    _chan_cred[chan_idx]->SetLatency(latency);
    _routers[u]->AddOutputChannel(_chan[chan_idx], _chan_cred[chan_idx]);
    _routers[v]->AddInputChannel(_chan[chan_idx], _chan_cred[chan_idx]);
    gHBMNetAccelAdj[u].push_back({v, opc[u]++, type});
    gAccelInputPortSrc[v][ipc[v]++] = u;  // record source router for U-turn detection
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

  // 4) Xbar ↔ Xbar (linear chain: Xbar_p ↔ Xbar_{p+1})
  for (int p = 0; p < P - 1; p++) {
    for (int b = 0; b < _xbar_xbar_bandwidth; b++) {
      add_link(accel_crossbar_id(p), accel_crossbar_id(p + 1),
               _xbar_xbar_latency, ACCEL_LINK_XBAR_XBAR);
      add_link(accel_crossbar_id(p + 1), accel_crossbar_id(p),
               _xbar_xbar_latency, ACCEL_LINK_XBAR_XBAR);
    }
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
  gRoutingFunctionMap["baseline_hbmnet_accelsim"] = &hbmnet_accelsim_baseline;
  gRoutingFunctionMap["hybrid_hbmnet_accelsim"]   = &hbmnet_accelsim_hybrid;
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
//  Xbar → HBM (direct) or Xbar chain → HBM (cross-partition)
//  HBM → Xbar (reply traffic)
//
//  With P Xbars in a chain, at Xbar_p:
//    if dest HBM belongs to partition p → go direct to HBM
//    else → step toward target partition via chain
// ============================================================
static int accel_get_hit_baseline_next(int cur, int dest) {
  if (accel_is_xbar(cur)) {
    int cur_p = cur;  // Xbar ID == partition index
    if (accel_is_xbar(dest)) {
      // Route toward dest Xbar via chain
      return accel_crossbar_id(cur_p + (dest > cur_p ? 1 : -1));
    }
    if (accel_is_hbm(dest)) {
      int dest_p = accel_hbm_part(accel_hbm_idx(dest));
      if (dest_p == cur_p) return dest;  // direct to local HBM
      // Step toward target partition
      return accel_crossbar_id(cur_p + (dest_p > cur_p ? 1 : -1));
    }
    assert(false);
    return -1;
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
//  At Xbar:    → toward miss_target via chain (direct if same partition)
//  At dest MC: → local HBM (deterministic final hop for SM→L2)
//  At intm MC: → back to local partition Xbar (baseline = no fabric)
// ============================================================
static int accel_get_miss_baseline_next(int cur, int miss_target) {
  if (accel_is_hbm(cur)) {
    // HBM → local MC (deterministic first hop)
    return accel_mc_id(accel_hbm_idx(cur));
  }

  if (accel_is_xbar(cur)) {
    int cur_p = cur;
    if (accel_is_mc(miss_target)) {
      int target_p = accel_hbm_part(accel_mc_idx(miss_target));
      if (target_p == cur_p) return miss_target;  // direct to local MC
      // Step toward target partition via Xbar chain
      return accel_crossbar_id(cur_p + (target_p > cur_p ? 1 : -1));
    }
    // target is Xbar (L2→SM): route through chain
    assert(accel_is_xbar(miss_target));
    int target_p = miss_target;
    assert(cur_p != target_p);  // should have been caught by eject
    return accel_crossbar_id(cur_p + (target_p > cur_p ? 1 : -1));
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
    // Track near-min path usage at ejection
    if (f->nm_used)
      ++gAccelNearMinPathsUsed;
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

// Forward declaration: link traversal tracking (defined later, after near-min)
static void accel_track_link_traversal(const Router *r, int cur, int chosen_port);

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
    ++gAccelTotalMissPackets;
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
  int chosen = ports[RandomInt(ports.size() - 1)];
  accel_track_link_traversal(r, cur, chosen);
  outputs->AddRange(chosen, 1, gNumVCs - 1, 1);
}


// ============================================================
//  Link traversal tracking helper
//  Called at each miss routing decision to count link-type utilization
//  and record per-direction saturation (used_credits / capacity).
//
//  Saturation = average across all parallel ports in the same direction
//  of (GetUsedCredit(port) / gAccelPortBufCapacity).
// ============================================================
static void accel_track_link_traversal(const Router *r, int cur, int chosen_port) {
  // Find chosen port's neighbor and link type
  int neighbor = -1;
  AccelLinkType lt = ACCEL_LINK_XBAR_XBAR;
  for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
    if (e.port == chosen_port) {
      neighbor = e.neighbor;
      lt = e.type;
      break;
    }
  }
  if (neighbor < 0) return;

  gAccelLinkTypeTraversals[lt]++;
  pair<int,int> dir_key = make_pair(cur, neighbor);
  gAccelLinkDirTraversals[dir_key]++;

  // Compute saturation: average used_credit/capacity across all ports in this direction
  if (gAccelPortBufCapacity > 0) {
    int total_used = 0, port_count = 0;
    for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
      const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
      if (e.neighbor == neighbor && e.type == lt) {
        total_used += r->GetUsedCredit(e.port);
        port_count++;
      }
    }
    if (port_count > 0) {
      double sat = (double)total_used / (double)(port_count * gAccelPortBufCapacity);
      gAccelLinkDirSatSum[dir_key] += sat;
    }
  }
}

// ============================================================
//  NEAR-MINIMAL ADAPTIVE ROUTING (non-hybrid: all miss)
//
//  Extension of min_adaptive that allows up to k extra hops
//  beyond minimum distance, using MC-MC fabric links.
//
//  At each hop:
//    - Min directions (dist decreases by 1): cost = avg_saturation
//    - Near-min directions (dist stays same, k budget > 0):
//        cost = avg_saturation + latency_penalty * saturation_factor
//      where latency_penalty = link_latency / max_link_latency (normalized)
//      and saturation_factor amplifies the penalty when congested.
//    - Pick least-cost direction overall
//    - U-turn prevention: exclude the router the flit came from
//    - k budget tracked in f->nm_budget (separate from f->ph)
//    - f->nm_used tracks whether near-min was ever taken on this path
//
//  Cost function rationale:
//    Min cost    = saturation (lower is better: less congested)
//    NearMin cost = saturation + normalized_latency * (1 + saturation)
//    The (1+saturation) factor means near-min is more expensive when
//    the detour direction is already congested, preventing piling up.
//    The normalized latency penalty reflects the actual delay cost of
//    the extra hop (MC-MC links are expensive: 160 cycles).
// ============================================================
static int accel_pick_near_min_port(const Router *r, const Flit *f,
                                     int in_channel,
                                     int cur, int target,
                                     const vector<vector<int>> &dist)
{
  int cur_dist = dist[cur][target];
  assert(cur_dist > 0);

  // Determine previous router for U-turn prevention
  int prev_router = -1;
  if (in_channel >= 0 && in_channel < (int)gAccelInputPortSrc[cur].size())
    prev_router = gAccelInputPortSrc[cur][in_channel];

  int k_remaining = f->nm_budget;

  // Max link latency for normalization (cached)
  static int max_lat = 0;
  if (max_lat == 0) {
    for (int i = 0; i < 5; i++)
      if (gAccelLinkLatency[i] > max_lat) max_lat = gAccelLinkLatency[i];
    if (max_lat == 0) max_lat = 1;
  }

  // Collect min and near-min directions grouped by neighbor
  vector<int> dir_nb;
  vector<vector<int>> dir_ports;
  vector<int> dir_credit;
  vector<bool> dir_is_nearmin;
  vector<int> dir_latency;

  for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
    if (!accel_valid_miss_routing_out(cur, e.type)) continue;
    if (e.neighbor == prev_router && prev_router >= 0) continue;

    int nb_dist = dist[e.neighbor][target];
    if (nb_dist < 0) continue;

    bool is_min = (nb_dist == cur_dist - 1);
    bool is_nearmin = (nb_dist == cur_dist && k_remaining > 0);
    if (!is_min && !is_nearmin) continue;

    // Find existing direction entry for this neighbor
    int idx = -1;
    for (size_t j = 0; j < dir_nb.size(); j++) {
      if (dir_nb[j] == e.neighbor) { idx = (int)j; break; }
    }
    if (idx < 0) {
      idx = (int)dir_nb.size();
      dir_nb.push_back(e.neighbor);
      dir_ports.push_back(vector<int>());
      dir_credit.push_back(0);
      dir_is_nearmin.push_back(is_nearmin);
      dir_latency.push_back(gAccelLinkLatency[e.type]);
    }
    dir_ports[idx].push_back(e.port);
    dir_credit[idx] += r->GetUsedCredit(e.port);
  }

  assert(!dir_nb.empty());

  // Pick best direction: lowest cost
  int best_idx = -1;
  double best_cost = numeric_limits<double>::max();

  for (size_t i = 0; i < dir_nb.size(); i++) {
    // Saturation: avg used credit / buffer capacity per port
    double avg_credit = (double)dir_credit[i] / (double)dir_ports[i].size();
    double saturation = (gAccelPortBufCapacity > 0)
        ? avg_credit / (double)gAccelPortBufCapacity : avg_credit;
    double cost;
    if (dir_is_nearmin[i]) {
      // Near-min cost: penalized by latency, scaled by near_min_penalty
      // penalty=1.0 → full latency penalty (conservative, original behavior)
      // penalty=0.0 → no penalty, pure congestion comparison
      double lat_penalty = (double)dir_latency[i] / (double)max_lat;
      cost = saturation + gAccelNearMinPenalty * lat_penalty * (1.0 + saturation);
    } else {
      // Min cost: boosted by (1 - penalty) to make near-min relatively cheaper
      // penalty=1.0 → no boost (min cost = sat, conservative)
      // penalty=0.0 → full boost (min cost = sat + lat*(1+sat), aggressive)
      cost = saturation;
    }
    if (cost < best_cost) {
      best_cost = cost;
      best_idx = (int)i;
    }
  }

  assert(best_idx >= 0);

  if (dir_is_nearmin[best_idx]) {
    f->nm_budget = k_remaining - 1;
    f->nm_used = true;
    ++gAccelNearMinNonMinDecisions;
    gAccelNearMinDirCount[make_pair(cur, dir_nb[best_idx])]++;
  } else {
    ++gAccelNearMinMinDecisions;
  }

  vector<int> &ports = dir_ports[best_idx];
  int chosen = ports[RandomInt(ports.size() - 1)];

  // Track link utilization and saturation
  accel_track_link_traversal(r, cur, chosen);

  return chosen;
}


// ============================================================
//  HYBRID ROUTING
//
//  Only routing function that uses hit/miss distinction.
//  At injection: probabilistic hit/miss via hit_rate.
//  Hit flits: deterministic baseline tree on Xbar<->HBM links.
//  Miss flits: sub-routing function selected by hybrid_routing config:
//    0=baseline, 1=min_adaptive, 2=ugal, 3=valiant, 4=min_oblivious,
//    5=near_min_adaptive
//
//  Field usage (completely separated):
//    UGAL/Valiant: f->ph = phase (0=toward intermediate, 1=toward final)
//                  f->intm = intermediate router ID or MISS_MARKER
//    Near-min:     f->nm_budget = remaining k budget
//                  f->nm_used = true if near-min hop ever taken
//                  f->ph = unused (0)
//    Others:       f->ph = 0, f->intm = MISS_MARKER
// ============================================================
void hbmnet_accelsim_hybrid( const Router *r, const Flit *f, int in_channel,
                              OutputSet *outputs, bool inject )
{
  if (inject) {
    outputs->Clear();
    accel_decide_hit_miss(f);

    if (f->intm == ACCELSIM_MISS_MARKER) {
      ++gAccelTotalMissPackets;
    }

    // Initialize per-routing-type fields
    f->ph = 0;
    f->nm_budget = 0;
    f->nm_used = false;
    if (gAccelHybridRouting == 5 && f->intm == ACCELSIM_MISS_MARKER) {
      f->nm_budget = gAccelNearMinK;  // near-min k budget
    }
    // UGAL/Valiant: ph=0 means "decision not yet made"
    // (decision happens at first routing hop, not injection)

    outputs->AddRange(-1, 1, gNumVCs - 1);
    return;
  }

  int cur, dest, miss_target;
  if (accelsim_routing_preamble(r, f, outputs, cur, dest, miss_target)) return;

  // ── Hit flits: deterministic baseline tree ──
  if (f->intm == ACCELSIM_HIT_MARKER) {
    accel_add_hit_data_ports(cur, dest, outputs);
    return;
  }

  // ── Miss flits: sub-routing ──
  int port = -1;

  switch (gAccelHybridRouting) {
    case 0:  // baseline (deterministic)
    {
      int next = accel_get_miss_baseline_next(cur, miss_target);
      vector<int> ports;
      for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
        const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
        if (e.neighbor == next && accel_valid_miss_out(cur, e.type))
          ports.push_back(e.port);
      }
      assert(!ports.empty());
      port = ports[RandomInt(ports.size() - 1)];
      break;
    }

    case 1:  // min_adaptive
    {
      port = accel_pick_adaptive_port(r, cur, miss_target, false,
                                      gHBMNetAccelDistMissFabric);
      break;
    }

    case 2:  // ugal
    {
      // UGAL decision at first routing hop: ph==0 means decision pending
      bool has_intm = (f->intm != ACCELSIM_MISS_MARKER && f->intm >= 0);
      if (!has_intm && f->ph == 0) {
        accel_ugal_decision(r, f, cur, miss_target);
        has_intm = (f->intm != ACCELSIM_MISS_MARKER && f->intm >= 0);
      }
      // Phase transition: arrived at intermediate
      if (has_intm && f->ph == 0 && cur == f->intm)
        f->ph = 1;

      int target = (has_intm && f->ph == 0) ? f->intm : miss_target;
      port = accel_pick_adaptive_port(r, cur, target, false,
                                      gHBMNetAccelDistMissFabric);
      break;
    }

    case 3:  // valiant
    {
      bool has_intm = (f->intm != ACCELSIM_MISS_MARKER && f->intm >= 0);
      if (!has_intm && f->ph == 0) {
        accel_valiant_decision(f, cur, miss_target);
        has_intm = (f->intm != ACCELSIM_MISS_MARKER && f->intm >= 0);
      }
      if (has_intm && f->ph == 0 && cur == f->intm)
        f->ph = 1;

      int target = (has_intm && f->ph == 0) ? f->intm : miss_target;
      vector<int> ports;
      accel_collect_minimal_ports(cur, target, false,
                                  gHBMNetAccelDistMissFabric, ports);
      assert(!ports.empty());
      port = ports[RandomInt(ports.size() - 1)];
      break;
    }

    case 4:  // min_oblivious
    {
      vector<int> ports;
      accel_collect_minimal_ports(cur, miss_target, false,
                                  gHBMNetAccelDistMissFabric, ports);
      assert(!ports.empty());
      port = ports[RandomInt(ports.size() - 1)];
      break;
    }

    case 5:  // near_min_adaptive (uses nm_budget, not ph)
    {
      port = accel_pick_near_min_port(r, f, in_channel, cur, miss_target,
                                      gHBMNetAccelDistMissFabric);
      // accel_pick_near_min_port already calls accel_track_link_traversal
      outputs->AddRange(port, 1, gNumVCs - 1, 1);
      return;  // early return: tracking already done inside
    }

    default:
      assert(false && "Unknown hybrid_routing type");
      return;
  }

  assert(port >= 0);
  // Track link traversal for non-near-min miss routing
  accel_track_link_traversal(r, cur, port);
  outputs->AddRange(port, 1, gNumVCs - 1, 1);
}
