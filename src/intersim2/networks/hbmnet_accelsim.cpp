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
int gHBMNetAccelSMPerXbar = 0;

// Reverse mapping: node -> router, node -> port
static vector<int> gAccelNodeToRouter;
static vector<int> gAccelNodeToPort;

// Reverse input port map: [router][input_port] -> source router (-1 for inject ports)
static vector<vector<int>> gAccelInputPortSrc;

// Link latency cache indexed by AccelLinkType
static int gAccelLinkLatency[5] = {0, 0, 0, 0, 0};

static int gAccelNearMinK = 1;
static double gAccelNearMinPenalty = 1.0;  // penalty multiplier for near-min cost
static int gAccelNearMinStrict = 0;      // 0=relieved(no non-min in opp partitions), 1=strict(no non-min in dest partition)
static int gAccelNearMinRemoteOnly = 0;  // 0=allow non-min everywhere, 1=only for remote (cross-partition) access

// Topology parameters cached for routing
static int gAccelNumSMs;
static int gAccelNumL2;
static int gAccelL;          // l2_per_hbm
static int gAccelInterleave;
static double gAccelHitRate;     // baseline_ratio repurposed as L2 hit rate
static int gAccelHybridRouting;  // 0=baseline,1=min_adaptive,4=min_oblivious,
                                 // 5=near_min_adaptive,6=near_min_random,7=fixed_min
                                 // Set via hybrid_routing= config (routing_function must be hybrid)

// Escape VC usage counters
long long gAccelEscapeVCEjects = 0;
long long gAccelTotalEjects    = 0;

// Near-min adaptive decision counters
long long gAccelNearMinMinDecisions    = 0;
long long gAccelNearMinNonMinDecisions = 0;
long long gAccelNearMinPlus2Decisions  = 0;  // +2 hop (backward) decisions

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
// Per-router-pair cumulative saturation, sample count, and max saturation
static map<pair<int,int>, double> gAccelLinkDirSatSum;
static map<pair<int,int>, long long> gAccelLinkDirSatCount;
static map<pair<int,int>, double> gAccelLinkDirSatMax;

// Windowed max link utilization (packets / cycle, window = ACCEL_UTIL_WINDOW cycles)
static const int ACCEL_UTIL_WINDOW = 200;
static long long gAccelElapsedCycles = 0;           // cycles since last reset_stats()
static long long gAccelLastWindowCycle = 0;          // cycle at last window flush
static map<pair<int,int>, long long> gAccelLinkDirWindowCount;  // count at window start
static long long gAccelLinkTypeWindowCount[5] = {0, 0, 0, 0, 0};
static map<pair<int,int>, double> gAccelLinkDirUtilMax;
static double gAccelLinkTypeUtilMax[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

// Cached buffer capacity per output port for data VCs: (num_vcs - 1) * vc_buf_size
static int gAccelPortBufCapacity = 0;

void hbmnet_accelsim_reset_stats() {
  gAccelEscapeVCEjects      = 0;
  gAccelTotalEjects         = 0;
  gAccelNearMinMinDecisions    = 0;
  gAccelNearMinNonMinDecisions = 0;
  gAccelNearMinPlus2Decisions  = 0;
  gAccelNearMinPathsUsed    = 0;
  gAccelTotalMissPackets    = 0;
  gAccelNearMinDirCount.clear();
  for (int i = 0; i < 5; i++) gAccelLinkTypeTraversals[i] = 0;
  gAccelLinkDirTraversals.clear();
  gAccelLinkDirSatSum.clear();
  gAccelLinkDirSatCount.clear();
  gAccelLinkDirSatMax.clear();
  gAccelLinkDirUtilMax.clear();
  for (int i = 0; i < 5; i++) gAccelLinkTypeUtilMax[i] = 0.0;
  gAccelElapsedCycles = 0;
  gAccelLastWindowCycle = 0;
  gAccelLinkDirWindowCount.clear();
  for (int i = 0; i < 5; i++) gAccelLinkTypeWindowCount[i] = 0;
}

// ============================================================
//  Per-cycle tick for windowed max link utilization tracking.
//  Must be called once per simulation cycle after stats reset.
// ============================================================
void hbmnet_accelsim_tick() {
  gAccelElapsedCycles++;
  long long window_elapsed = gAccelElapsedCycles - gAccelLastWindowCycle;
  if (window_elapsed < ACCEL_UTIL_WINDOW) return;

  double w = (double)window_elapsed;

  // Update max util for each direction
  for (map<pair<int,int>, long long>::iterator it = gAccelLinkDirTraversals.begin();
       it != gAccelLinkDirTraversals.end(); ++it) {
    map<pair<int,int>, long long>::iterator wit = gAccelLinkDirWindowCount.find(it->first);
    long long start_count = (wit != gAccelLinkDirWindowCount.end()) ? wit->second : 0;
    double util = (double)(it->second - start_count) / w;
    map<pair<int,int>, double>::iterator mit = gAccelLinkDirUtilMax.find(it->first);
    if (mit == gAccelLinkDirUtilMax.end() || util > mit->second)
      gAccelLinkDirUtilMax[it->first] = util;
  }

  // Update max util for each link type
  for (int t = 0; t < 5; t++) {
    double util = (double)(gAccelLinkTypeTraversals[t] - gAccelLinkTypeWindowCount[t]) / w;
    if (util > gAccelLinkTypeUtilMax[t]) gAccelLinkTypeUtilMax[t] = util;
  }

  // Advance window start
  gAccelLastWindowCycle = gAccelElapsedCycles;
  for (map<pair<int,int>, long long>::iterator it = gAccelLinkDirTraversals.begin();
       it != gAccelLinkDirTraversals.end(); ++it) {
    gAccelLinkDirWindowCount[it->first] = it->second;
  }
  for (int t = 0; t < 5; t++) {
    gAccelLinkTypeWindowCount[t] = gAccelLinkTypeTraversals[t];
  }
}

// ============================================================
//  Called from FlitChannel::Send() when a flit actually traverses a link.
//  This counts real link usage (post allocation), not routing decisions.
// ============================================================
void hbmnet_accelsim_count_link_traversal(int src_router, int dst_router) {
  // Find link type from adjacency table
  AccelLinkType lt = ACCEL_LINK_XBAR_XBAR;
  bool found = false;
  for (size_t i = 0; i < gHBMNetAccelAdj[src_router].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[src_router][i];
    if (e.neighbor == dst_router) {
      lt = e.type;
      found = true;
      break;
    }
  }
  if (!found) return;

  gAccelLinkTypeTraversals[lt]++;
  pair<int,int> dir_key = make_pair(src_router, dst_router);
  gAccelLinkDirTraversals[dir_key]++;
}

int hbmnet_accelsim_node_to_router(int node) { return gAccelNodeToRouter[node]; }
int hbmnet_accelsim_node_to_port(int node)   { return gAccelNodeToPort[node]; }

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
  double elapsed = (double)gAccelElapsedCycles;

  // Flush the final partial window so max util covers remaining cycles
  if (gAccelElapsedCycles > gAccelLastWindowCycle) {
    double w = (double)(gAccelElapsedCycles - gAccelLastWindowCycle);
    for (map<pair<int,int>, long long>::iterator it = gAccelLinkDirTraversals.begin();
         it != gAccelLinkDirTraversals.end(); ++it) {
      map<pair<int,int>, long long>::iterator wit = gAccelLinkDirWindowCount.find(it->first);
      long long start_count = (wit != gAccelLinkDirWindowCount.end()) ? wit->second : 0;
      double util = (double)(it->second - start_count) / w;
      map<pair<int,int>, double>::iterator mit = gAccelLinkDirUtilMax.find(it->first);
      if (mit == gAccelLinkDirUtilMax.end() || util > mit->second)
        gAccelLinkDirUtilMax[it->first] = util;
    }
    for (int t = 0; t < 5; t++) {
      double util = (double)(gAccelLinkTypeTraversals[t] - gAccelLinkTypeWindowCount[t]) / w;
      if (util > gAccelLinkTypeUtilMax[t]) gAccelLinkTypeUtilMax[t] = util;
    }
  }

  // === Link Type Utilization ===
  cout << "=== Link Type Utilization ===" << endl;
  cout << "  Elapsed cycles: " << gAccelElapsedCycles
       << "  Window size: " << ACCEL_UTIL_WINDOW << endl;
  for (int i = 0; i < 5; i++) {
    double avg_util = (elapsed > 0) ? (double)gAccelLinkTypeTraversals[i] / elapsed : 0.0;
    cout << "  " << link_type_names[i] << ": count=" << gAccelLinkTypeTraversals[i]
         << " avg_util=" << avg_util << " pkt/cyc"
         << " max_util=" << gAccelLinkTypeUtilMax[i] << " pkt/cyc" << endl;
  }
  cout << "  Total link traversals: " << total << endl;

  // Per-direction traversals (grouped by link type) with saturation
  // Accumulate per-type aggregate saturation (weighted average) and max
  double type_sat_sum[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  long long type_sat_count[5] = {0, 0, 0, 0, 0};
  double type_sat_max[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  long long type_count[5] = {0, 0, 0, 0, 0};

  for (map<pair<int,int>, long long>::iterator it = gAccelLinkDirTraversals.begin();
       it != gAccelLinkDirTraversals.end(); ++it) {
    int src = it->first.first;
    int dst = it->first.second;
    AccelLinkType lt = ACCEL_LINK_XBAR_XBAR;
    for (size_t j = 0; j < gHBMNetAccelAdj[src].size(); j++) {
      if (gHBMNetAccelAdj[src][j].neighbor == dst) { lt = gHBMNetAccelAdj[src][j].type; break; }
    }
    long long count = it->second;
    type_count[(int)lt] += count;
    map<pair<int,int>, double>::iterator sit = gAccelLinkDirSatSum.find(it->first);
    if (sit != gAccelLinkDirSatSum.end()) type_sat_sum[(int)lt] += sit->second;
    map<pair<int,int>, long long>::iterator scit = gAccelLinkDirSatCount.find(it->first);
    if (scit != gAccelLinkDirSatCount.end()) type_sat_count[(int)lt] += scit->second;
    map<pair<int,int>, double>::iterator mit = gAccelLinkDirSatMax.find(it->first);
    if (mit != gAccelLinkDirSatMax.end() && mit->second > type_sat_max[(int)lt])
      type_sat_max[(int)lt] = mit->second;
  }

  // Print per-direction detail
  cout << "=== Per-Direction Miss Link Traversals ===" << endl;
  for (int t = 0; t < 5; t++) {
    if (t == ACCEL_LINK_XBAR_HBM || t == ACCEL_LINK_MC_HBM) continue;
    bool header_printed = false;
    for (map<pair<int,int>, long long>::iterator it = gAccelLinkDirTraversals.begin();
         it != gAccelLinkDirTraversals.end(); ++it) {
      int src = it->first.first;
      int dst = it->first.second;
      AccelLinkType lt = ACCEL_LINK_XBAR_XBAR;
      for (size_t j = 0; j < gHBMNetAccelAdj[src].size(); j++) {
        if (gHBMNetAccelAdj[src][j].neighbor == dst) { lt = gHBMNetAccelAdj[src][j].type; break; }
      }
      if ((int)lt != t) continue;
      if (!header_printed) {
        double type_avg_sat = (type_sat_count[t] > 0)
            ? type_sat_sum[t] / (double)type_sat_count[t] * 100.0 : 0.0;
        cout << "  [" << link_type_names[t] << "]"
             << " avg_sat=" << type_avg_sat << "%"
             << " max_sat=" << (type_sat_max[t] * 100.0) << "%" << endl;
        header_printed = true;
      }
      long long count = it->second;
      double avg_util = (elapsed > 0) ? (double)count / elapsed : 0.0;
      map<pair<int,int>, double>::iterator uit = gAccelLinkDirUtilMax.find(it->first);
      double max_util = (uit != gAccelLinkDirUtilMax.end()) ? uit->second : 0.0;
      double avg_sat = 0.0;
      map<pair<int,int>, double>::iterator sit = gAccelLinkDirSatSum.find(it->first);
      map<pair<int,int>, long long>::iterator scit = gAccelLinkDirSatCount.find(it->first);
      long long sat_n = (scit != gAccelLinkDirSatCount.end()) ? scit->second : 0;
      if (sit != gAccelLinkDirSatSum.end() && sat_n > 0)
        avg_sat = sit->second / (double)sat_n;
      double max_sat = 0.0;
      map<pair<int,int>, double>::iterator msit = gAccelLinkDirSatMax.find(it->first);
      if (msit != gAccelLinkDirSatMax.end()) max_sat = msit->second;
      cout << "    " << accel_router_name(src) << " -> " << accel_router_name(dst)
           << ": count=" << count
           << " avg_util=" << avg_util << " pkt/cyc"
           << " max_util=" << max_util << " pkt/cyc"
           << " avg_sat=" << (avg_sat * 100.0) << "%"
           << " max_sat=" << (max_sat * 100.0) << "%" << endl;
    }
  }

  // Machine-parseable aggregate line for run_all_moe.py
  double sat_xx = (type_sat_count[ACCEL_LINK_XBAR_XBAR] > 0)
      ? type_sat_sum[ACCEL_LINK_XBAR_XBAR] / (double)type_sat_count[ACCEL_LINK_XBAR_XBAR] * 100.0 : 0.0;
  double sat_xmc = (type_sat_count[ACCEL_LINK_XBAR_MC] > 0)
      ? type_sat_sum[ACCEL_LINK_XBAR_MC] / (double)type_sat_count[ACCEL_LINK_XBAR_MC] * 100.0 : 0.0;
  double sat_mc = (type_sat_count[ACCEL_LINK_MC_MC] > 0)
      ? type_sat_sum[ACCEL_LINK_MC_MC] / (double)type_sat_count[ACCEL_LINK_MC_MC] * 100.0 : 0.0;
  cout << "Miss link avg saturation:"
       << " XBAR_XBAR=" << sat_xx << "%"
       << " XBAR_MC=" << sat_xmc << "%"
       << " MC_MC=" << sat_mc << "%" << endl;
  cout << "Miss link max saturation:"
       << " XBAR_XBAR=" << (type_sat_max[ACCEL_LINK_XBAR_XBAR] * 100.0) << "%"
       << " XBAR_MC=" << (type_sat_max[ACCEL_LINK_XBAR_MC] * 100.0) << "%"
       << " MC_MC=" << (type_sat_max[ACCEL_LINK_MC_MC] * 100.0) << "%" << endl;

  // Near-min adaptive decision stats
  long long nm_total = gAccelNearMinMinDecisions + gAccelNearMinNonMinDecisions;
  if (nm_total > 0) {
    double nm_pct = 100.0 * (double)gAccelNearMinNonMinDecisions / (double)nm_total;
    cout << "=== Near-Min Adaptive Decisions ===" << endl;
    cout << "  Near-min decisions: total=" << nm_total
         << " min=" << gAccelNearMinMinDecisions
         << " non-min=" << gAccelNearMinNonMinDecisions
         << " near-min ratio=" << nm_pct << "%" << endl;
    cout << "  +2 hop decisions: " << gAccelNearMinPlus2Decisions << endl;

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
          if (gHBMNetAccelAdj[src][j].neighbor == dst) { lt = gHBMNetAccelAdj[src][j].type; break; }
        }
        double avg_sat = 0.0;
        map<pair<int,int>, double>::iterator sit = gAccelLinkDirSatSum.find(it->first);
        map<pair<int,int>, long long>::iterator scit2 = gAccelLinkDirSatCount.find(it->first);
        long long sat_n2 = (scit2 != gAccelLinkDirSatCount.end()) ? scit2->second : 0;
        if (sit != gAccelLinkDirSatSum.end() && sat_n2 > 0)
          avg_sat = sit->second / (double)sat_n2;
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
  _num_xbars      = config.GetInt("num_xbars");
  _hbm_per_side   = config.GetInt("hbm_per_side");
  _num_hbm_stacks = _num_xbars * _hbm_per_side * 2;

  _sm_per_xbar    = config.GetInt("sm_per_xbar");
  _num_sms        = _sm_per_xbar * _num_xbars;

  _l2_per_hbm     = config.GetInt("l2_per_hbm");
  _num_l2_slices  = _num_hbm_stacks * _l2_per_hbm;
  _l2_interleave  = config.GetInt("l2_interleave");
  _is_fabric      = config.GetInt("is_fabric");

  int P = _num_xbars;
  int K = _num_hbm_stacks;

  assert(P >= 1);
  assert(_hbm_per_side >= 1);
  assert(_sm_per_xbar >= 1);

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
  gHBMNetAccelK         = K;
  gHBMNetAccelP         = P;
  gHBMNetAccelHPS       = _hbm_per_side;
  gHBMNetAccelSMPerXbar = _sm_per_xbar;
  gAccelNumSMs       = _num_sms;
  gAccelNumL2        = _num_l2_slices;
  gAccelL            = _l2_per_hbm;
  gAccelInterleave   = _l2_interleave;
  gAccelHitRate        = config.GetFloat("baseline_ratio");  // repurposed as hit rate
  gAccelNearMinK       = config.GetInt("near_min_k");
  gAccelNearMinPenalty = config.GetFloat("near_min_penalty");
  gAccelNearMinStrict     = config.GetInt("near_min_strict");
  gAccelNearMinRemoteOnly = config.GetInt("near_min_remote_only");

  // Cache link latencies by type for near-min penalty calculation
  gAccelLinkLatency[ACCEL_LINK_XBAR_XBAR] = _xbar_xbar_latency;
  gAccelLinkLatency[ACCEL_LINK_XBAR_HBM]  = _xbar_hbm_latency;
  gAccelLinkLatency[ACCEL_LINK_XBAR_MC]   = _xbar_mc_latency;
  gAccelLinkLatency[ACCEL_LINK_MC_HBM]    = _mc_hbm_latency;
  gAccelLinkLatency[ACCEL_LINK_MC_MC]     = _mc_mc_latency;

  // Cache buffer capacity for saturation calculation (data VCs only, excluding escape VC 0)
  gAccelPortBufCapacity = (config.GetInt("num_vcs") - 1) * config.GetInt("vc_buf_size");

  string hr = config.GetStr("hybrid_routing");
  if      (hr == "min_adaptive")       gAccelHybridRouting = 1;
  else if (hr == "min_oblivious")      gAccelHybridRouting = 4;
  else if (hr == "near_min_adaptive")  gAccelHybridRouting = 5;
  else if (hr == "near_min_random")    gAccelHybridRouting = 6;
  else if (hr == "fixed_min")          gAccelHybridRouting = 7;
  else /* "baseline" or unknown */     gAccelHybridRouting = 0;

  // Total routers: P Xbars + K MCs + K HBMs
  _size = P + 2 * K;

  // Per-router concentration and first node
  gHBMNetAccelRouterConc.resize(_size);
  gHBMNetAccelRouterFirstNode.resize(_size);

  // Xbar p: SM nodes [p*sm_per_xbar .. (p+1)*sm_per_xbar - 1]
  for (int p = 0; p < P; p++) {
    gHBMNetAccelRouterConc[p] = _sm_per_xbar;
    gHBMNetAccelRouterFirstNode[p] = p * _sm_per_xbar;
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
       << " sm_per_xbar=" << _sm_per_xbar
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

  // Create routers
  // Optional per-type internal speedup override for HBM routers.
  // 0.0 (default) = use the global internal_speedup for all routers.
  double hbm_speedup = config.GetFloat("hbm_internal_speedup");

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

// Return the partition (Xbar index) for any router type.
//   Xbar: partition = Xbar ID
//   MC / HBM: partition derived from HBM-index row → row / HPS
static inline int accel_router_part(int router) {
  if (accel_is_xbar(router)) return router;
  if (accel_is_mc(router))  return accel_hbm_part(accel_mc_idx(router));
  if (accel_is_hbm(router)) return accel_hbm_part(accel_hbm_idx(router));
  return 0;
}

// Return the destination partition from a miss_target router
//   (MC router → its partition; Xbar router → its ID)
static inline int accel_miss_target_part(int miss_target) {
  if (accel_is_mc(miss_target))
    return accel_hbm_part(accel_mc_idx(miss_target));
  return miss_target;  // Xbar ID == partition
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
//  Data-VC credit helper: sum occupancy of VCs 1..num_vcs-1 for a port.
//  Excludes escape VC (VC 0) which is tracked separately.
// ============================================================
static inline int accel_data_vc_credit(const Router *r, int port) {
  int used = 0;
  for (int vc = 1; vc < gNumVCs; vc++)
    used += r->GetUsedCreditVC(port, vc);
  return used;
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
        dir_credit[e.neighbor] += accel_data_vc_credit(r, e.port);
      }
    }
  }
}

// ============================================================
//  Credit-cycling port selection helper
//
//  Given a list of parallel ports for a chosen direction, start from a
//  random port and cycle through (round-robin) until finding one with
//  available credit on data VCs (VC 1..N-1).  If all ports are exhausted,
//  the last tried port is returned and the preamble's escape VC (pri=0)
//  acts as fallback via Duato's protocol.
// ============================================================
static int accel_pick_port_credit_cycle(const Router *r, vector<int> &ports) {
  int num_ports = (int)ports.size();
  int start = RandomInt(num_ports - 1);
  int chosen = ports[start];
  for (int i = 1; i < num_ports; i++) {
    if (accel_data_vc_credit(r, chosen) < gAccelPortBufCapacity) break;
    chosen = ports[(start + i) % num_ports];
  }
  return chosen;
}

// ============================================================
//  Pick random direction, random port within (oblivious)
// ============================================================
static int accel_pick_oblivious_port(const Router *r, int cur, int target,
                                     bool is_hit,
                                     const vector<vector<int>> &dist)
{
  map<int, vector<int>> dir_ports;
  map<int, int> dir_credit;
  accel_collect_minimal_dirs(r, cur, target, is_hit, dist, dir_ports, dir_credit);
  assert(!dir_ports.empty());

  map<int, vector<int>>::iterator it = dir_ports.begin();
  advance(it, RandomInt(dir_ports.size() - 1));
  return accel_pick_port_credit_cycle(r, it->second);
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

  return accel_pick_port_credit_cycle(r, dir_ports[best_nb]);
}

// ============================================================
//  Deterministic MC↔HBM hop helper
//  Used for: HBM→MC (L2 miss first hop) and MC→HBM (miss final hop)
// ============================================================
// MC→HBM (eject hop): all VCs are fine, no priority needed.
static void accel_add_mc_to_hbm_ports(int cur, OutputSet *outputs) {
  vector<int> ports;
  for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
    if (e.type == ACCEL_LINK_MC_HBM)
      ports.push_back(e.port);
  }
  assert(!ports.empty());
  outputs->AddRange(ports[RandomInt(ports.size() - 1)], 0, gNumVCs - 1);
}

// HBM→MC (inject/first hop): prefer data VCs (pri=1); escape VC 0 (pri=0)
// only as last resort so the flit does not start its fabric journey on VC 0.
static void accel_add_hbm_to_mc_ports(int cur, OutputSet *outputs) {
  vector<int> ports;
  for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
    if (e.type == ACCEL_LINK_MC_HBM)
      ports.push_back(e.port);
  }
  assert(!ports.empty());
  int port = ports[RandomInt(ports.size() - 1)];
  if (gNumVCs > 1)
    outputs->AddRange(port, 1, gNumVCs - 1, 1);
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

  // 1. Eject at destination router
  if (cur_router == dest_router) {
    ++gAccelTotalEjects;
    if (f->nm_used)
      ++gAccelNearMinPathsUsed;
      
    if (accel_is_xbar(cur_router) && f->vc == 0) {
        ++gAccelEscapeVCEjects;
    }
    
    outputs->AddRange(accel_eject_port(f->dest), 0, gNumVCs - 1);
    return true;
  }

  bool is_hit = (f->intm == ACCELSIM_HIT_MARKER);
  miss_target = is_hit ? -1 : accel_miss_target(f->dest);

  if (!is_hit) {
    // Deterministic HBM→MC hop (L2 miss injection first hop)
    if (accel_is_hbm(cur_router)) {
      accel_add_hbm_to_mc_ports(cur_router, outputs);
      return true;
    }
    // Deterministic MC→HBM hop (at dest MC, SM→L2 miss final hop)
    if (accel_is_mc(cur_router) && accel_is_mc(miss_target)
        && cur_router == miss_target) {
        
      if (f->vc == 0) {
          ++gAccelEscapeVCEjects;
      }
      
      accel_add_mc_to_hbm_ports(cur_router, outputs);
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
  
  if (f->vc == 0) {
      outputs->AddRange(esc_ports[RandomInt(esc_ports.size() - 1)], 0, 0, 0);
      return true;
  }

  // Add escape VC as low-priority fallback for data flits (Duato's protocol).
  // iq_router suppresses pri=0 when pri=1 data VCs are available.
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

// Forward declaration: output-side saturation tracking (defined after near-min)
static void accel_track_link_sat(const Router *r, int cur, int chosen_port);

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
  int chosen = accel_pick_port_credit_cycle(r, ports);
  accel_track_link_sat(r, cur, chosen);
  outputs->AddRange(chosen, 1, gNumVCs - 1, 1);
}


// ============================================================
//  Output-side saturation tracking helper (event-driven).
//  Called at each miss routing decision to record per-direction saturation.
//  Saturation = average data-VC credit across all parallel ports in the
//  same direction / gAccelPortBufCapacity, excluding escape VC 0.
//  Util counting is done separately via FlitChannel::Send().
// ============================================================
static void accel_track_link_sat(const Router *r, int cur, int chosen_port) {
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

  pair<int,int> dir_key = make_pair(cur, neighbor);

  if (gAccelPortBufCapacity > 0) {
    int total_used = 0, port_count = 0;
    for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
      const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
      if (e.neighbor == neighbor && e.type == lt) {
        total_used += accel_data_vc_credit(r, e.port);
        port_count++;
      }
    }
    if (port_count > 0) {
      double sat = (double)total_used / (double)(port_count * gAccelPortBufCapacity);
      gAccelLinkDirSatCount[dir_key]++;
      gAccelLinkDirSatSum[dir_key] += sat;
      map<pair<int,int>, double>::iterator msit = gAccelLinkDirSatMax.find(dir_key);
      if (msit == gAccelLinkDirSatMax.end() || sat > msit->second)
        gAccelLinkDirSatMax[dir_key] = sat;
    }
  }
}

// ============================================================
//  Helper: gate the near-min budget based on strict/relieved mode.
//
//  Goal of near-min: replace one Xbar-Xbar hop with extra MC-MC hops.
//  Xbar↔MC crossing happens exactly once; +k budget = k extra MC-MC steps.
//  e.g. B200 +1 hop:  Xbar0→MC0→MC1→MC2  (enters MC fabric one row early)
//       B200 +2 hop:  Xbar0→MC0→MC1→MC2→MC3
//
//  Strict:   no non-min once packet is inside dest partition (conservative).
//  Relieved: no non-min outside the natural routing range [lo, hi]
//            where lo=min(src_part,dest_part), hi=max(src_part,dest_part).
//            Blocks "overshoot" past dest but allows non-min within dest partition.
//  Applies to both Xbar→MC and MC→Xbar traffic directions.
// ============================================================
static int accel_nm_effective_budget(int k_remaining, int cur, int miss_target, int src_part) {
  if (k_remaining <= 0) return 0;

  int dest_part = accel_miss_target_part(miss_target);
  int cur_part  = accel_router_part(cur);

  // Strict mode
  if (gAccelNearMinStrict) {
    if (cur_part == dest_part) return 0;
  } 
  // Relieved mode
  else {
    int lo = (src_part < dest_part) ? src_part : dest_part;
    int hi = (src_part > dest_part) ? src_part : dest_part;
    if (cur_part < lo || cur_part > hi) return 0;
  }

  return k_remaining;
}

// ============================================================
//  NEAR-MINIMAL ROUTING — shared infrastructure
//
//  Goal: bypass one Xbar-Xbar hop by entering MC fabric one row early
//  and traversing extra MC-MC links instead.
//  e.g. instead of Xbar0→Xbar1→MC2, take Xbar0→MC0→MC1→MC2 (+1 hop).
//  Xbar↔MC crossing happens exactly once; budget = number of extra MC-MC hops.
//
//  Hop classification by dist to miss_target:
//    min:   nb_dist == cur_dist - 1             (no budget consumed)
//    +1:    nb_dist == cur_dist,   k_remaining >= 1  (1 budget)
//    +2:    nb_dist == cur_dist+1, k_remaining >= 2  (2 budget, 2x penalty)
//
//  Hard constraints (shared by adaptive and random):
//    - U-turn prevention (skip prev_router)
//    - At MC (SM→L2 target): skip XBAR-MC links — Xbar↔MC crosses once only
//    - Non-min Xbar→MC entry: must be same column as dest MC (no cross-column entry)
//    - Strict / relieved mode: via accel_nm_effective_budget
//    - Remote-only: nm_budget set to 0 at injection for local (same-partition) access
//    - Forward-escape filter: non-min hops must lead to forward-progressing escape
// ============================================================

// Direction candidate for near-min routing
struct NearMinDir {
  int neighbor;
  vector<int> ports;
  int total_credit;  // sum of used credit across parallel ports
  int hop_type;      // 0=min, 1=+1, 2=+2
  int latency;       // link latency (for adaptive penalty calculation)
};

// Collect candidate directions for near-min routing (shared between adaptive/random)
static void accel_collect_near_min_dirs(const Router *r, const Flit *f,
                                         int in_channel,
                                         int cur, int target,
                                         const vector<vector<int>> &dist,
                                         vector<NearMinDir> &dirs)
{
  dirs.clear();
  int cur_dist = dist[cur][target];
  assert(cur_dist > 0);

  int prev_router = -1;
  if (in_channel >= 0 && in_channel < (int)gAccelInputPortSrc[cur].size())
    prev_router = gAccelInputPortSrc[cur][in_channel];

  int src_part = f->ph;
  int k_remaining = accel_nm_effective_budget(f->nm_budget, cur, target, src_part);

  int cur_part  = accel_router_part(cur);
  int dest_part = accel_miss_target_part(target);

  bool force_mc_mc = false;
  if (accel_is_mc(cur) && accel_is_mc(target)) {
    if (accel_hbm_col(accel_mc_idx(cur)) == accel_hbm_col(accel_mc_idx(target))) {
      int cur_row = accel_hbm_row(accel_mc_idx(cur));
      int tgt_row = accel_hbm_row(accel_mc_idx(target));
      int mc_mc_dist = abs(cur_row - tgt_row);
      int global_dist = dist[cur][target];
      if (mc_mc_dist - global_dist <= k_remaining) {
        force_mc_mc = true;
      }
    }
  }

  for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
    if (!accel_valid_miss_routing_out(cur, e.type)) continue;

    // 1. U-turn prevention
    if (e.neighbor == prev_router && prev_router >= 0) continue;

    int nb_dist = dist[e.neighbor][target];
    if (nb_dist < 0) continue;

    // 2. Distance-based classification
    bool is_min   = (nb_dist == cur_dist - 1);
    bool is_plus1 = (nb_dist == cur_dist     && k_remaining >= 1);
    bool is_plus2 = (nb_dist == cur_dist + 1 && k_remaining >= 2);
    if (!is_min && !is_plus1 && !is_plus2) continue;

    int nb_part = accel_router_part(e.neighbor);

    // 3. Directionality constraint for MC-MC and Xbar-Xbar
    if (e.type == ACCEL_LINK_MC_MC || e.type == ACCEL_LINK_XBAR_XBAR) {
      if (cur_part < dest_part && nb_part < cur_part) continue;
      if (cur_part > dest_part && nb_part > cur_part) continue;
      if (cur_part == dest_part && nb_part != dest_part) continue;
    }

    // Prevent backward MC-MC hops within the same column
    if (e.type == ACCEL_LINK_MC_MC && accel_is_mc(target)) {
      if (accel_hbm_col(accel_mc_idx(cur)) == accel_hbm_col(accel_mc_idx(target))) {
        int cur_row = accel_hbm_row(accel_mc_idx(cur));
        int tgt_row = accel_hbm_row(accel_mc_idx(target));
        int nb_row  = accel_hbm_row(accel_mc_idx(e.neighbor));
        if (tgt_row > cur_row && nb_row < cur_row) continue;
        if (tgt_row < cur_row && nb_row > cur_row) continue;
      }
    }

    // 4. Xbar-MC link constraint
    if (e.type == ACCEL_LINK_XBAR_MC) {
      if (accel_is_mc(cur)) {
        if (force_mc_mc) continue;

        bool is_l2_sm = accel_is_xbar(target);
        bool is_cross_col = accel_is_mc(target) &&
            (accel_hbm_col(accel_mc_idx(cur)) != accel_hbm_col(accel_mc_idx(target)));
        bool is_same_col_out_of_budget = accel_is_mc(target) && !is_cross_col && !force_mc_mc;

        if (!is_l2_sm && !is_cross_col && !is_same_col_out_of_budget) continue;
      }
      else if (accel_is_xbar(cur)) {
        if (accel_is_xbar(target)) continue;
        if (accel_is_mc(target)) {
          if (accel_hbm_col(accel_mc_idx(e.neighbor)) != accel_hbm_col(accel_mc_idx(target))) continue;
          int tgt_row = accel_hbm_row(accel_mc_idx(target));
          int nb_row  = accel_hbm_row(accel_mc_idx(e.neighbor));
          int mc_mc_dist = abs(nb_row - tgt_row);
          int total_extra_cost = mc_mc_dist + 1 - cur_dist;
          if (total_extra_cost > k_remaining) continue;
        }
      }
    }

    // 5. Forward-escape filter for non-min hops.
    //    When a non-min hop is taken (early fabric entry, etc.), the escape path
    //    (VC 0 baseline) from the neighbor must still make forward progress on the
    //    baseline path compared to the current router.  If the escape from the
    //    neighbor would land at a point >= current baseline distance to target,
    //    it means a stuck packet would be forced to escape *backward* through the
    //    Xbar chain, creating extra load on the escape VC and risking deadlock.
    if (!is_min) {
      int esc_from_nb = accel_get_miss_baseline_next(e.neighbor, target);
      if (esc_from_nb >= 0 &&
          gHBMNetAccelDistMissBaseline[esc_from_nb][target] >=
          gHBMNetAccelDistMissBaseline[cur][target]) continue;
    }

    int hop_type = is_min ? 0 : (is_plus1 ? 1 : 2);

    // Group by neighbor router
    int idx = -1;
    for (size_t j = 0; j < dirs.size(); j++) {
      if (dirs[j].neighbor == e.neighbor) { idx = (int)j; break; }
    }
    if (idx < 0) {
      NearMinDir d;
      d.neighbor = e.neighbor;
      d.total_credit = 0;
      d.hop_type = hop_type;
      d.latency = gAccelLinkLatency[e.type];
      dirs.push_back(d);
      idx = (int)dirs.size() - 1;
    }
    dirs[idx].ports.push_back(e.port);
    dirs[idx].total_credit += accel_data_vc_credit(r, e.port);
  }
}

// Update near-min statistics after direction selection
static void accel_update_near_min_stats(const Flit *f, int cur,
                                         const NearMinDir &chosen) {
  if (chosen.hop_type == 0) {
    ++gAccelNearMinMinDecisions;
  } else if (chosen.hop_type == 1) {
    f->nm_budget = f->nm_budget - 1;
    f->nm_used = true;
    ++gAccelNearMinNonMinDecisions;
    gAccelNearMinDirCount[make_pair(cur, chosen.neighbor)]++;
  } else {
    f->nm_budget = f->nm_budget - 2;
    f->nm_used = true;
    ++gAccelNearMinNonMinDecisions;
    ++gAccelNearMinPlus2Decisions;
    gAccelNearMinDirCount[make_pair(cur, chosen.neighbor)]++;
  }
}

// ============================================================
//  NEAR-MINIMAL ADAPTIVE ROUTING
//
//  Cost-based direction selection among min + near-min candidates.
//  Credit-cycling among parallel ports in the chosen direction.
// ============================================================
static int accel_pick_near_min_port(const Router *r, const Flit *f,
                                     int in_channel,
                                     int cur, int target,
                                     const vector<vector<int>> &dist)
{
  vector<NearMinDir> dirs;
  accel_collect_near_min_dirs(r, f, in_channel, cur, target, dist, dirs);
  assert(!dirs.empty());

  static int max_lat = 0;
  if (max_lat == 0) {
    for (int i = 0; i < 5; i++)
      if (gAccelLinkLatency[i] > max_lat) max_lat = gAccelLinkLatency[i];
    if (max_lat == 0) max_lat = 1;
  }

  // Pick best direction: lowest cost
  int best_idx = -1;
  double best_cost = numeric_limits<double>::max();

  for (size_t i = 0; i < dirs.size(); i++) {
    double avg_credit = (double)dirs[i].total_credit / (double)dirs[i].ports.size();
    double saturation = (gAccelPortBufCapacity > 0)
        ? avg_credit / (double)gAccelPortBufCapacity : avg_credit;
    double cost;
    if (dirs[i].hop_type == 0) {
      cost = saturation;
    } else {
      double lat_penalty = (double)dirs[i].latency / (double)max_lat;
      double multiplier = (dirs[i].hop_type == 1) ? 1.0 : 2.0;
      cost = saturation + multiplier * gAccelNearMinPenalty * lat_penalty * (1.0 + saturation);
    }
    if (cost < best_cost ||
        (cost == best_cost && best_idx >= 0 && dirs[i].neighbor < dirs[best_idx].neighbor)) {
      best_cost = cost;
      best_idx = (int)i;
    }
  }

  assert(best_idx >= 0);

  accel_update_near_min_stats(f, cur, dirs[best_idx]);
  int chosen = accel_pick_port_credit_cycle(r, dirs[best_idx].ports);
  accel_track_link_sat(r, cur, chosen);
  return chosen;
}

// ============================================================
//  NEAR-MINIMAL RANDOM (OBLIVIOUS) ROUTING
//
//  Random direction selection among min + near-min candidates.
//  Credit-cycling among parallel ports in the chosen direction.
// ============================================================
static int accel_pick_near_min_random_port(const Router *r, const Flit *f,
                                            int in_channel,
                                            int cur, int target,
                                            const vector<vector<int>> &dist)
{
  vector<NearMinDir> dirs;
  accel_collect_near_min_dirs(r, f, in_channel, cur, target, dist, dirs);
  assert(!dirs.empty());

  // Random direction selection
  int chosen_idx = RandomInt((int)dirs.size() - 1);

  accel_update_near_min_stats(f, cur, dirs[chosen_idx]);
  int chosen = accel_pick_port_credit_cycle(r, dirs[chosen_idx].ports);
  accel_track_link_sat(r, cur, chosen);
  return chosen;
}


// ============================================================
//  FABRIC FIXED MIN ROUTING
//
//  Deterministic min routing that always prefers MC-MC fabric paths
//  over Xbar-Xbar chain when MC-MC min paths exist.
//
//  Priority at each hop:
//    1. MC-MC min (stay in fabric, cur=MC)
//    2. Xbar-MC min (enter fabric from Xbar)
//    3. Other min (Xbar-Xbar chain)
// ============================================================
static int accel_pick_fixed_min_port(const Router *r, int cur, int target,
                                     const vector<vector<int>> &dist)
{
  int cur_dist = dist[cur][target];
  assert(cur_dist > 0);

  vector<int> mc_mc_ports, xbar_mc_ports, other_ports;

  for (size_t i = 0; i < gHBMNetAccelAdj[cur].size(); i++) {
    const AccelAdjEntry &e = gHBMNetAccelAdj[cur][i];
    if (!accel_valid_miss_routing_out(cur, e.type)) continue;
    if (dist[e.neighbor][target] != cur_dist - 1) continue;  // not min

    if      (e.type == ACCEL_LINK_MC_MC)   mc_mc_ports.push_back(e.port);
    else if (e.type == ACCEL_LINK_XBAR_MC) xbar_mc_ports.push_back(e.port);
    else                                   other_ports.push_back(e.port);
  }

  // Prefer MC-MC > Xbar-MC (entering fabric) > Xbar-Xbar, with credit cycling
  if (!mc_mc_ports.empty())  return accel_pick_port_credit_cycle(r, mc_mc_ports);
  if (!xbar_mc_ports.empty()) return accel_pick_port_credit_cycle(r, xbar_mc_ports);
  assert(!other_ports.empty());
  return accel_pick_port_credit_cycle(r, other_ports);
}


// ============================================================
//  HYBRID ROUTING
//
//  Only routing function that uses hit/miss distinction.
//  At injection: probabilistic hit/miss via hit_rate.
//  Hit flits: deterministic baseline tree on Xbar<->HBM links.
//  Miss flits: sub-routing function selected by hybrid_routing config:
//    0=baseline, 1=min_adaptive, 4=min_oblivious,
//    5=near_min_adaptive, 6=near_min_random, 7=fixed_min
//
//  Field usage (completely separated):
//    Near-min:     f->nm_budget = remaining k budget
//                  f->nm_used = true if near-min hop ever taken
//                  f->ph = src_part (source partition, for strict/relieved mode)
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
    f->nm_budget = 0;
    f->nm_used = false;
    f->ph = 0;
    if (f->intm == ACCELSIM_MISS_MARKER &&
        (gAccelHybridRouting == 5 || gAccelHybridRouting == 6)) {
      // Store src_part in f->ph for strict/relieved mode gating
      int src_router = gAccelNodeToRouter[f->src];
      f->ph = accel_router_part(src_router);  // 0 for MC src (shouldn't happen)

      // Remote-only: disable non-min for same-partition (local) access
      int dest_part = accel_miss_target_part(accel_miss_target(f->dest));
      bool is_remote = (f->ph != dest_part);
      f->nm_budget = (!gAccelNearMinRemoteOnly || is_remote) ? gAccelNearMinK : 0;
    }
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
      port = accel_pick_port_credit_cycle(r, ports);
      break;
    }

    case 1:  // min_adaptive
    {
      port = accel_pick_adaptive_port(r, cur, miss_target, false,
                                      gHBMNetAccelDistMissFabric);
      break;
    }

    case 4:  // min_oblivious
    {
      port = accel_pick_oblivious_port(r, cur, miss_target, false,
                                       gHBMNetAccelDistMissFabric);
      break;
    }

    case 5:  // near_min_adaptive (uses nm_budget, not ph)
    {
      port = accel_pick_near_min_port(r, f, in_channel, cur, miss_target,
                                      gHBMNetAccelDistMissFabric);
      break;
    }

    case 6:  // near_min_random (same constraints as near_min_adaptive, random selection)
    {
      port = accel_pick_near_min_random_port(r, f, in_channel, cur, miss_target,
                                             gHBMNetAccelDistMissFabric);
      break;
    }

    case 7:  // fixed_min (prefer MC-MC fabric over Xbar-Xbar chain)
    {
      port = accel_pick_fixed_min_port(r, cur, miss_target,
                                       gHBMNetAccelDistMissFabric);
      break;
    }

    default:
      assert(false && "Unknown hybrid_routing type");
      return;
  }

  assert(port >= 0);
  // Track link traversal for non-near-min miss routing
  accel_track_link_sat(r, cur, port);
  outputs->AddRange(port, 1, gNumVCs - 1, 1);
}
