#ifndef _GPUNET_MULTIPORT_HPP_
#define _GPUNET_MULTIPORT_HPP_

#include <cassert>
#include <vector>

#include "network.hpp"
#include "routefunc.hpp"

// GPUNet with multiple parallel interpartition ports instead of increased bandwidth
class GPUNetMultiport : public Network {
  
  int _n_shader;
  int _n_l2slice;
  int _l2slice_per_partition;
  
  // Partition configuration
  bool _use_partition;
  int _n_partition;
  
  // Number of parallel interpartition ports (replaces bandwidth increase)
  int _interpartition_ports;
  
  // Latencies
  int _inject_eject_latency;
  int _interpartition_latency;

  void _ComputeSize(const Configuration& config);
  void _BuildNet(const Configuration& config);
  void _SetupChannels();
  void _SetChannelProperties(FlitChannel* channel, CreditChannel* credit_channel, int latency);

public:
  GPUNetMultiport(const Configuration& config, const string & name);
  static void RegisterRoutingFunctions();
};

void hierarchical_gpunet_multiport(const Router *r, const Flit *f, int in_channel,
                                   OutputSet *outputs, bool inject);

#endif
