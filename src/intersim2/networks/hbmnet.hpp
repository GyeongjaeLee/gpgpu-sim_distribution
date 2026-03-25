#ifndef _HBMNET_HPP_
#define _HBMNET_HPP_

#include "network.hpp"
#include "routefunc.hpp"
#include <vector>
#include <utility>

using namespace std;

// ============================================================
//  HBMNet topology
//
//  Node mapping 
//    Node 0 .. N-1      : SM nodes
//    Node N .. N+M-1    : L2 slice nodes
//
//  Router layout:
//    Router 0            : Partitioned Crossbar 0 (top partition)
//    Router 1            : Partitioned Crossbar 1 (bottom partition)
//    Router 2 .. K+1     : Xbar HBM routers (K total)
//
//  2D physical layout (K HBMs in 2 columns x K/2 rows):
//    Column 0: HBM indices 0 .. K/2-1  (top to bottom)
//    Column 1: HBM indices K/2 .. K-1  (top to bottom)
//    Partition 0 (top):  rows 0 .. K/4-1
//    Partition 1 (bot):  rows K/4 .. K/2-1
//
//  Example K=8:
//    HBM0  HBM4    <- partition 0
//    HBM1  HBM5    <- partition 0
//    [Crossbar 0]
//    [Crossbar 1]
//    HBM2  HBM6    <- partition 1
//    HBM3  HBM7    <- partition 1
//
//  Fabric adds HBM-HBM links between vertically adjacent routers
//  on the same shoreline (no horizontal — GPU die is in between).
// ============================================================

// Distance tables:
//   gHBMNetDistBaseline[r][r] = min hops using only baseline links
//                              (Xbar-HBM, Xbar-Xbar; no HBM-HBM)
//   gHBMNetDistFabric[r][r]   = min hops using all links including HBM-HBM fabric
extern vector<vector<int>> gHBMNetDistBaseline;
extern vector<vector<int>> gHBMNetDistFabric;

// Adjacency table: gHBMNetAdj[router] = list of (neighbor_router, output_port)
// Contains ALL physical links (baseline + fabric).
extern vector<vector<pair<int,int>>> gHBMNetAdj;

// UGAL routing decision counters
extern long long gUGALMinDecisions;
extern long long gUGALNonMinDecisions;

// Escape VC usage: flits that arrived at destination on vc==0
extern long long gEscapeVCEjects;
extern long long gTotalEjects;

// Reset all routing statistics (call at the start of each simulation run).
void hbmnet_reset_stats();

// Per-router concentration
extern vector<int> gHBMNetRouterConc;
// Per-router first node ID
extern vector<int> gHBMNetRouterFirstNode;

class HBMNet : public Network {

private:
  int _num_sms;           // N: total SM count
  int _num_l2_slices;     // M: total L2 slice count
  int _num_hbm_stacks;    // K: number of Xbar HBM routers
  int _l2_per_hbm;        // L = M / K: L2 slices per HBM router
  int _l2_interleave;     // 0: sequential, 1: interleaved L2 mapping

  int _is_fabric;         // enable HBM-HBM fabric links

  int _xbar_xbar_latency;  // Crossbar <-> Crossbar
  int _xbar_hbm_latency;   // Crossbar <-> Xbar HBM
  int _hbm_hbm_latency;    // Xbar HBM <-> Xbar HBM (fabric)
  int _xbar_xbar_bandwidth; // inter-partition port count
  int _xbar_hbm_bandwidth;  // crossbar-HBM port count (per HBM)
  int _hbm_hbm_bandwidth;   // HBM-HBM fabric port count (per pair)

  int _xbar_xbar_channels;
  int _xbar_hbm_channels;
  int _hbm_hbm_channels;

  void _ComputeSize( const Configuration &config );
  void _BuildNet( const Configuration& config );
  void _BuildRoutingTable();

  // HBM index (0..K-1) -> 2D position
  int _HBMCol(int hbm_idx) const { return hbm_idx / (_num_hbm_stacks / 2); }
  int _HBMRow(int hbm_idx) const { return hbm_idx % (_num_hbm_stacks / 2); }
  int _HBMPartition(int hbm_idx) const {
    int row = _HBMRow(hbm_idx);
    return (row < _num_hbm_stacks / 4) ? 0 : 1;
  }

public:
  HBMNet( const Configuration &config, const string & name );

  static void RegisterRoutingFunctions();
};

// Routing functions
void hbmnet_baseline( const Router *r, const Flit *f, int in_channel,
                        OutputSet *outputs, bool inject );

void hbmnet_min_adaptive( const Router *r, const Flit *f, int in_channel,
                      OutputSet *outputs, bool inject );

void hbmnet_ugal( const Router *r, const Flit *f, int in_channel,
                      OutputSet *outputs, bool inject );

void hbmnet_valiant( const Router *r, const Flit *f, int in_channel,
                      OutputSet *outputs, bool inject );

void hbmnet_min_oblivious( const Router *r, const Flit *f, int in_channel,
                      OutputSet *outputs, bool inject );

void hbmnet_hybrid( const Router *r, const Flit *f, int in_channel,
                      OutputSet *outputs, bool inject );

// Node -> Router mapping
int hbmnet_node_to_router(int node);
int hbmnet_node_to_port(int node);

#endif
