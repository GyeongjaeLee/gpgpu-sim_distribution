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
//  Router layout:
//    Router 0           : Crossbar 0 (partition 0)
//    Router 1           : Crossbar 1 (partition 1)
//    Router 2..K+1      : MC routers 0..K-1 (Memory Controllers)
//    Router K+2..2K+1   : HBM routers 0..K-1
//
//  Links:
//    Xbar ↔ HBM   : hit path (direct L2 cache hit)
//    Xbar ↔ MC    : miss path entry (L2 cache miss)
//    MC   ↔ HBM   : miss path local (1-to-1, MC_h ↔ HBM_h)
//    MC   ↔ MC    : fabric links (vertical, only when is_fabric=1)
//    Xbar ↔ Xbar  : inter-partition
//
//  Hit rate (baseline_ratio config) determines the fraction of flits
//  that take the direct Xbar→HBM path vs the Xbar→MC→HBM path.
//
//  Bandwidth constraint: xbar_hbm_bw == mc_hbm_bw for fair
//  processing of hit/miss flits at HBM routers.
//
//  2D physical layout (same as HBMNet):
//    Column 0: indices 0..K/2-1,  Column 1: indices K/2..K-1
//    Partition 0: rows 0..K/4-1,  Partition 1: rows K/4..K/2-1
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

// Cached K (num_hbm_stacks) for routing helpers
extern int gHBMNetAccelK;

// UGAL routing decision counters
extern long long gAccelUGALMinDecisions;
extern long long gAccelUGALNonMinDecisions;

// Escape VC usage counters
extern long long gAccelEscapeVCEjects;
extern long long gAccelTotalEjects;

// Reset all routing statistics
void hbmnet_accelsim_reset_stats();

class HBMNetAccelSim : public Network {

private:
  int _num_sms;
  int _num_l2_slices;
  int _num_hbm_stacks;    // K
  int _l2_per_hbm;        // L = M / K
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
    return (row < _num_hbm_stacks / 4) ? 0 : 1;
  }

public:
  HBMNetAccelSim( const Configuration &config, const string & name );

  static void RegisterRoutingFunctions();
};

// Routing functions
void hbmnet_accelsim_baseline( const Router *r, const Flit *f, int in_channel,
                                OutputSet *outputs, bool inject );

void hbmnet_accelsim_min_adaptive( const Router *r, const Flit *f, int in_channel,
                                    OutputSet *outputs, bool inject );

void hbmnet_accelsim_min_oblivious( const Router *r, const Flit *f, int in_channel,
                                     OutputSet *outputs, bool inject );

void hbmnet_accelsim_ugal( const Router *r, const Flit *f, int in_channel,
                            OutputSet *outputs, bool inject );

void hbmnet_accelsim_valiant( const Router *r, const Flit *f, int in_channel,
                               OutputSet *outputs, bool inject );

void hbmnet_accelsim_hybrid( const Router *r, const Flit *f, int in_channel,
                              OutputSet *outputs, bool inject );

// Node -> Router mapping
int hbmnet_accelsim_node_to_router(int node);
int hbmnet_accelsim_node_to_port(int node);

#endif
