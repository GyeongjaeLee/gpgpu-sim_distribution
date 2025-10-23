#ifndef _GPUNET_HPP_
#define _GPUNET_HPP_
#include <cassert>
#include <vector>

#include "network.hpp"
#include "routefunc.hpp"

class GPUNet : public Network {
  
  // # of layer that requests traverse to get partition crossbars
  // ex) (_l == 3): Shader->TPC->GPC->Crossbar, (_l == 4): Shader->TPC->CPC->GPC->Crossbar
  int _l;

  int _n_shader;
  int _n_l2slice;
  int _l2slice_per_partition;
  
  // # of lower-level units connected to a single higher-level module.
  // l =          0,             1,      ...,     _l-1  
  // ex) [Shader per TPCs, TPCs per XPC, ..., GPCs per Crossbars]
  vector<int> _ratio;
  // ex) [TPCs, ..., GPCs, Crossbars]
  vector<int> _total_units;
  // _offsets for Router and Channel id of each layer
  vector<int> _offsets;
  // l =        0,     ...,    _l-2,     _l-1                _l
  // ex) [Shader->TPC, ..., CPC->GPC, GPC->Crossbar, Crossbar->L2Slice]
  // _speedups requires _l + 1 entries to indicate injection/ejection speedup which should be 1
  vector<int> _speedups;
    
  // A100, H100 supports partitioned GPUNet, not supported in V100
  bool _use_partition;
  // # of partitions, _p = 1 for non-partitioned GPUNet
  int _n_partition;
  int _interpartition_speedup;

  

  vector<pair<int, int> > _l2slice_coords;

  void _ComputeSize(const Configuration& config);
  void _BuildNet(const Configuration& config);

  // Set latency and bandwidth for a channel based on its layer
  void _SetupChannels();
  void _SetChannelProperties(FlitChannel* channel, CreditChannel* credit_channel,
                   int l, bool is_interpartition = false);
  int _GetWireLatency(int l, bool is_interpartition = false) const;
  int _GetChannelBandwidth(int l, bool is_interpartition = false) const;
  
  int _WireLatency(int l) const;
  int _FloorplanLatency(int src_x, int src_y, int dst_x, int dst_y) const;

public:

  GPUNet( const Configuration& config, const string & name );
  static void RegisterRoutingFunctions();
};

void hierarchical_gpunet( const Router *r, const Flit *f, int in_channel,
                   OutputSet *outputs, bool inject );

#endif
