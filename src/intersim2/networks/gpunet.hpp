#ifndef _GPUNET_HPP_
#define _GPUNET_HPP_
#include <cassert>
#include <vector>
#include <string>

#include "network.hpp"
#include "routefunc.hpp"

// Global tracking functions for GPUNet link utilization
void gpunet_tick();
void gpunet_count_link_traversal(int src_router, int dst_router);
void gpunet_print_link_stats();

// ================================================================
// GPUNet: Hierarchical GPU on-chip network topology
//
// Models the A100/H100-style partitioned crossbar network:
//
//   SM side (path: SM → ... → Crossbar):
//     SM nodes → Level 0 (e.g. TPC) → Level 1 (e.g. GPC) → Crossbar
//
//   L2 side (path: Crossbar → ... → L2 node):
//     Crossbar → Level 0 (e.g. L2Slice group) → L2 nodes
//
//   Partitions: multiple crossbar routers connected via interpartition links.
//
// Array index conventions (all arrays are indexed 0 .. level-1):
//   Index 0 is always the level CLOSEST to the leaf nodes.
//     SM side:  index 0 = TPC level (closest to SMs)
//     L2 side:  index 0 = L2Slice level (closest to L2 nodes)
// ================================================================

class GPUNet : public Network {

  // ---- SM-side hierarchy ----------------------------------------
  int _l;               // number of intermediate SM-side levels (0 = SM direct to crossbar)
  int _n_sm;            // total SM (shader) nodes across all partitions
  int _sm_per_partition;

  // _units[i]: concentration ratio at SM-side level i
  //   i=0: SMs per TPC router
  //   i=1: TPC routers per GPC router
  //   ...
  vector<int> _units;

  // _ports[i]: number of parallel physical ports on the upward link at SM-side level i
  //   i=0: ports on TPC → GPC (or crossbar if _l==1) link
  //   i=1: ports on GPC → Crossbar link
  //   ...
  vector<int> _ports;

  // _sm_total[i]: total number of routers at SM-side level i, per partition
  //   i=0: TPC count per partition
  //   i=1: GPC count per partition
  //   ...
  vector<int> _sm_total;

  // _latency[i]: wire latency (cycles) on SM-side link at level i
  //   i=0: TPC → GPC latency
  //   i=1: GPC → Crossbar latency
  //   ...
  vector<int> _latency;

  // ---- L2-side hierarchy ----------------------------------------
  int _l2_l;            // number of intermediate L2-side levels (0 = L2 direct to crossbar)
  int _n_l2;            // total L2 slice nodes across all partitions
  int _l2_per_partition;

  // _l2_units[i]: concentration ratio at L2-side level i
  //   i=0: L2 nodes per L2Slice router
  //   i=1: L2Slice routers per upper-level router
  //   ...
  vector<int> _l2_units;

  // _l2_ports[i]: number of parallel physical ports on the upward link at L2-side level i
  //   i=0: ports on L2Node → L2Slice (or crossbar if _l2_l==1) link
  //   i=1: ports on L2Slice → Crossbar link
  //   ...
  vector<int> _l2_ports;

  // _l2_total[i]: total number of routers at L2-side level i, per partition
  //   i=0: L2Slice count per partition
  //   ...
  vector<int> _l2_total;

  // _l2_latency[i]: wire latency (cycles) on L2-side link at level i
  //   i=0: L2Node → L2Slice latency
  //   ...
  vector<int> _l2_latency;

  // ---- Partition/Crossbar --------------------------------------
  bool _use_partition;  // true if network is partitioned
  int _n_partition;     // number of partitions
  int _interpartition_ports;
  int _interpartition_latency;

  int _inj_lat;         // injection latency
  int _ej_lat;          // ejection latency

  // ---- Router ID layout (per partition) ----
  // Within each partition block of size _routers_per_partition:
  //   _sm_offsets[i]  : start index of SM-side level-i routers
  //   _xbar_offset    : index of the crossbar router
  //   _l2_offsets[i]  : start index of L2-side level-i routers
  //
  // Full router ID: p * _routers_per_partition + local_offset
  vector<int> _sm_offsets;
  int         _xbar_offset;
  vector<int> _l2_offsets;
  int         _routers_per_partition;

  void _ComputeSize(const Configuration& config);
  void _BuildNet(const Configuration& config);

public:
  // Returns the number of children for parent parent_idx given total_children
  // and concentration ratio. Handles the uneven last group via ceil division.
  static int _GetNumChildren(int parent_idx, int total_children, int ratio);

  GPUNet(const Configuration& config, const string& name);
  static void RegisterRoutingFunctions();
};

void hierarchical_gpunet(const Router *r, const Flit *f, int in_channel,
                          OutputSet *outputs, bool inject);

#endif
