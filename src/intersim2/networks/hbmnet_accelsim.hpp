#ifndef _HBMNET_ACCELSIM_HPP_
#define _HBMNET_ACCELSIM_HPP_

#include "network.hpp"
#include "routefunc.hpp"
#include <vector>
#include <utility>

using namespace std;

// ============================================================
//  HBMNetAccelSim topology
//
//  Models L2-HBM path for AccelSim integration.
//  AccelSim structure: SM ─ BookSim ─ L2
//  This topology models the L2→HBM path with explicit MC routers.
//
//  Hit path  (L2 cache hit):  SM → Xbar → HBM → L2       (direct)
//  Miss path (L2 cache miss): SM → Xbar → MC  → HBM → L2 (via MC)
//
//  Node mapping:
//    Node 0..N-1       : SM nodes
//    Node N..N+M-1     : L2 slice nodes
//
//  Config parameters:
//    num_xbars (P)     : number of Xbar routers (e.g. 1, 2, 4)
//    hbm_per_side (H)  : HBM stacks per side per Xbar (e.g. 1, 2, 3)
//    K = P * H * 2     : total HBM stacks
//
//  Router layout:
//    Router 0..P-1       : Crossbar 0..P-1
//    Router P..P+K-1     : MC routers 0..K-1 (Memory Controllers)
//    Router P+K..P+2K-1  : HBM routers 0..K-1
//
//  Links:
//    Xbar ↔ HBM   : hit path (direct L2 cache hit)
//    Xbar ↔ MC    : miss path entry (L2 cache miss)
//    MC   ↔ HBM   : miss path local (1-to-1, MC_h ↔ HBM_h)
//    MC   ↔ MC    : fabric links (vertical, only when is_fabric=1)
//    Xbar ↔ Xbar  : linear chain (Xbar_p ↔ Xbar_{p+1})
//
//  Hit rate (baseline_ratio config) determines the fraction of flits
//  that take the direct Xbar→HBM path vs the Xbar→MC→HBM path.
//
//  Bandwidth constraint: xbar_hbm_bw == mc_hbm_bw for fair
//  processing of hit/miss flits at HBM routers.
//
//  2D physical layout:
//    Column 0: indices 0..K/2-1,  Column 1: indices K/2..K-1
//    Each Xbar p owns rows p*H..(p+1)*H-1 in each column.
//    MC fabric links are vertical (same column, adjacent rows).
// ============================================================

// Link type for adjacency entries
enum AccelLinkType {
  ACCEL_LINK_XBAR_XBAR = 0,
  ACCEL_LINK_XBAR_HBM  = 1,  // hit path
  ACCEL_LINK_XBAR_MC   = 2,  // miss path entry
  ACCEL_LINK_MC_HBM    = 3,  // miss path local
  ACCEL_LINK_MC_MC     = 4   // miss path fabric
};

struct AccelAdjEntry {
  int neighbor;
  int port;
  AccelLinkType type;
};

// Full adjacency with link types
extern vector<vector<AccelAdjEntry>> gHBMNetAccelAdj;

// Distance tables (directed BFS, respecting link-type constraints):
//   Hit:          Xbar↔HBM, Xbar↔Xbar
//   Miss baseline: Xbar→MC, MC→HBM, Xbar↔Xbar  (no MC↔MC)
//   Miss fabric:   Xbar→MC, MC→HBM, MC↔MC, Xbar↔Xbar
extern vector<vector<int>> gHBMNetAccelDistHit;
extern vector<vector<int>> gHBMNetAccelDistMissBaseline;
extern vector<vector<int>> gHBMNetAccelDistMissFabric;

// Per-router concentration and first node ID
extern vector<int> gHBMNetAccelRouterConc;
extern vector<int> gHBMNetAccelRouterFirstNode;

// Cached topology parameters for routing helpers
extern int gHBMNetAccelK;          // total HBM stacks
extern int gHBMNetAccelP;          // number of Xbars
extern int gHBMNetAccelHPS;        // hbm_per_side
extern int gHBMNetAccelSMPerXbar;  // SMs per Xbar (concentration)

// UGAL routing decision counters
extern long long gAccelUGALMinDecisions;
extern long long gAccelUGALNonMinDecisions;

// Escape VC usage counters
extern long long gAccelEscapeVCEjects;
extern long long gAccelTotalEjects;

// Near-min adaptive decision counters
extern long long gAccelNearMinMinDecisions;
extern long long gAccelNearMinNonMinDecisions;
extern long long gAccelNearMinPlus2Decisions;

// Near-min path-level tracking
extern long long gAccelNearMinPathsUsed;   // packets that took near-min at least once
extern long long gAccelTotalMissPackets;    // total miss packets (denominator)

// Per-link-type traversal counters (indexed by AccelLinkType)
extern long long gAccelLinkTypeTraversals[5];

// Reset all routing statistics
void hbmnet_accelsim_reset_stats();

// Print link utilization and near-min direction stats
void hbmnet_accelsim_print_link_stats();

class HBMNetAccelSim : public Network {

private:
  int _num_sms;
  int _num_l2_slices;
  int _num_xbars;         // P
  int _hbm_per_side;      // H (HBM stacks per side per Xbar)
  int _num_hbm_stacks;    // K = P * H * 2
  int _sm_per_xbar;       // SMs per Xbar (concentration); num_sms = _sm_per_xbar * P
  int _l2_per_hbm;        // L2 slices per HBM stack; num_l2_slices = _l2_per_hbm * K
  int _l2_interleave;
  int _is_fabric;

  int _xbar_xbar_latency;
  int _xbar_hbm_latency;   // hit path
  int _xbar_mc_latency;    // miss path entry
  int _mc_hbm_latency;     // miss path local
  int _mc_mc_latency;      // miss path fabric

  int _xbar_xbar_bandwidth;
  int _xbar_hbm_bandwidth;
  int _xbar_mc_bandwidth;
  int _mc_hbm_bandwidth;
  int _mc_mc_bandwidth;

  void _ComputeSize( const Configuration &config );
  void _BuildNet( const Configuration &config );
  void _BuildRoutingTable();

  int _HBMPartition(int hbm_idx) const {
    int row = hbm_idx % (_num_hbm_stacks / 2);
    return row / _hbm_per_side;
  }

public:
  HBMNetAccelSim( const Configuration &config, const string & name );

  static void RegisterRoutingFunctions();
};

// Routing functions
// is_fabric=0: baseline (deterministic miss path)
// is_fabric=1: hybrid  (hit/miss split; miss sub-routing via hybrid_routing config)
void hbmnet_accelsim_baseline( const Router *r, const Flit *f, int in_channel,
                                OutputSet *outputs, bool inject );

void hbmnet_accelsim_hybrid( const Router *r, const Flit *f, int in_channel,
                              OutputSet *outputs, bool inject );

// Node -> Router mapping
int hbmnet_accelsim_node_to_router(int node);
int hbmnet_accelsim_node_to_port(int node);

#endif
