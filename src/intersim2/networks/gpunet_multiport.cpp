#include "booksim.hpp"
#include <vector>
#include <sstream>
#include <cmath>

#include "gpunet_multiport.hpp"
#include "misc_utils.hpp"
#include "globals.hpp"

// Global variables for routing function
int gMultiportPartitions;
int gMultiportPorts;
int gShaderPerPartition;
int gL2PerPartition;

GPUNetMultiport::GPUNetMultiport(const Configuration& config, const string & name)
: Network(config, name)
{
  _ComputeSize(config);
  _Alloc();
  _BuildNet(config);
}

void GPUNetMultiport::_ComputeSize(const Configuration& config)
{
  _n_shader = config.GetInt("shader");
  _n_l2slice = config.GetInt("l2slice");
  _nodes = _n_shader + _n_l2slice;
  
  _use_partition = (config.GetInt("use_partition") == 1);
  _n_partition = _use_partition ? config.GetInt("n_partition") : 1;
  _l2slice_per_partition = _n_l2slice / _n_partition;
  
  _interpartition_ports = config.GetInt("interpartition_ports");
  if (_interpartition_ports <= 0) {
    _interpartition_ports = 1;
  }
  
  _inject_eject_latency = config.GetInt("inject_eject_latency");
  _interpartition_latency = config.GetInt("interpartition_latency");
  
  _size = 2 * _n_partition;
  
  if (_use_partition) {
    _channels = 2 * _n_partition * (_n_partition - 1) * _interpartition_ports;
  } else {
    _channels = 0;
  }
  
  gMultiportPartitions = _n_partition;
  gMultiportPorts = _interpartition_ports;
  gShaderPerPartition = _n_shader / _n_partition;
  gL2PerPartition = _l2slice_per_partition;
  gNodes = _nodes;
}

void GPUNetMultiport::_BuildNet(const Configuration& config)
{
  ostringstream name;
  int shader_per_partition = _n_shader / _n_partition;
  
  // STEP 1: Create routers
  for (int p = 0; p < _n_partition; ++p) {
    int req_router_id = p;
    name.str("");
    name << "router_request_partition_" << p;
    
    int req_inputs = shader_per_partition;
    int req_outputs = _l2slice_per_partition;
    
    if (_use_partition) {
      req_inputs += (_n_partition - 1) * _interpartition_ports;
      req_outputs += (_n_partition - 1) * _interpartition_ports;
    }
    
    _routers[req_router_id] = Router::NewRouter(config, this, name.str(), 
                                                req_router_id, req_inputs, req_outputs, 1);
    _timed_modules.push_back(_routers[req_router_id]);
    
    int rep_router_id = _n_partition + p;
    name.str("");
    name << "router_reply_partition_" << p;
    
    _routers[rep_router_id] = Router::NewRouter(config, this, name.str(),
                                                rep_router_id, req_outputs, req_inputs, 1);
    _timed_modules.push_back(_routers[rep_router_id]);
  }
  
  // STEP 2: Connect injection channels
  for (int n = 0; n < _nodes; ++n) {
    if (n < _n_shader) {
      int partition = n / shader_per_partition;
      if (partition >= _n_partition) partition = _n_partition - 1;
      
      _routers[partition]->AddInputChannel(_inject[n], _inject_cred[n]);
    } else {
      int l2_index = n - _n_shader;
      int partition = l2_index / _l2slice_per_partition;
      if (partition >= _n_partition) partition = _n_partition - 1;
      
      _routers[_n_partition + partition]->AddInputChannel(_inject[n], _inject_cred[n]);
    }
  }
  
  // STEP 3: Connect ejection channels
  for (int n = 0; n < _nodes; ++n) {
    if (n < _n_shader) {
      int partition = n / shader_per_partition;
      if (partition >= _n_partition) partition = _n_partition - 1;
      
      _routers[_n_partition + partition]->AddOutputChannel(_eject[n], _eject_cred[n]);
    } else {
      int l2_index = n - _n_shader;
      int partition = l2_index / _l2slice_per_partition;
      if (partition >= _n_partition) partition = _n_partition - 1;
      
      _routers[partition]->AddOutputChannel(_eject[n], _eject_cred[n]);
    }
  }
  
  // STEP 4: Connect interpartition channels
  if (_use_partition) {
    int channel_idx = 0;
    
    for (int src_p = 0; src_p < _n_partition; ++src_p) {
      for (int dst_p = 0; dst_p < _n_partition; ++dst_p) {
        if (src_p == dst_p) continue;
        
        for (int port = 0; port < _interpartition_ports; ++port) {
          _routers[src_p]->AddOutputChannel(_chan[channel_idx], _chan_cred[channel_idx]);
          _routers[dst_p]->AddInputChannel(_chan[channel_idx], _chan_cred[channel_idx]);
          
          int reply_channel = _channels / 2 + channel_idx;
          _routers[_n_partition + dst_p]->AddOutputChannel(_chan[reply_channel], 
                                                           _chan_cred[reply_channel]);
          _routers[_n_partition + src_p]->AddInputChannel(_chan[reply_channel],
                                                          _chan_cred[reply_channel]);
          channel_idx++;
        }
      }
    }
  }
  
  _SetupChannels();
}

void GPUNetMultiport::_SetupChannels()
{
  for (int i = 0; i < _nodes; ++i) {
    _SetChannelProperties(_inject[i], _inject_cred[i], _inject_eject_latency);
    _SetChannelProperties(_eject[i], _eject_cred[i], _inject_eject_latency);
  }
  
  if (_use_partition) {
    for (int c = 0; c < _channels; ++c) {
      _SetChannelProperties(_chan[c], _chan_cred[c], _interpartition_latency);
    }
  }
}

void GPUNetMultiport::_SetChannelProperties(FlitChannel* channel, 
                                           CreditChannel* credit_channel, 
                                           int latency)
{
  channel->SetLatency(latency);
  channel->SetBandwidth(1);
  credit_channel->SetLatency(latency);
  credit_channel->SetBandwidth(1);
}

void GPUNetMultiport::RegisterRoutingFunctions()
{
  gRoutingFunctionMap["hierarchical_gpunet_multiport"] = &hierarchical_gpunet_multiport;
}

void hierarchical_gpunet_multiport(const Router *r, const Flit *f, int in_channel,
                                  OutputSet *outputs, bool inject)
{
  int vcBegin = 0, vcEnd = gNumVCs - 1;
  
  if (inject) {
    outputs->AddRange(-1, vcBegin, vcEnd);
    return;
  }
  
  int dest = f->dest;
  int router_id = r->GetID();
  
  int total_shaders = gShaderPerPartition * gMultiportPartitions;
  
  bool is_request_router = (router_id < gMultiportPartitions);
  int router_partition = is_request_router ? router_id : (router_id - gMultiportPartitions);
  
  int dest_partition = (dest < total_shaders) ?
                       (dest / gShaderPerPartition) :
                       ((dest - total_shaders) / gL2PerPartition);
  
  int out_port = -1;
  
  if (is_request_router) {
    if (router_partition == dest_partition) {
      int local_l2_index = (dest - total_shaders) % gL2PerPartition;
      out_port = local_l2_index;
    } else {
      int base_port = gL2PerPartition;
      int partition_offset = 0;
      
      for (int p = 0; p < gMultiportPartitions; ++p) {
        if (p == router_partition) continue;
        
        if (p == dest_partition) {
          int port_select = f->id % gMultiportPorts;
          out_port = base_port + partition_offset + port_select;
          break;
        }
        
        partition_offset += gMultiportPorts;
      }
    }
  } else {
    if (router_partition == dest_partition) {
      int local_shader_index = dest % gShaderPerPartition;
      out_port = local_shader_index;
    } else {
      int base_port = gShaderPerPartition;
      int partition_offset = 0;
      
      for (int p = 0; p < gMultiportPartitions; ++p) {
        if (p == router_partition) continue;
        
        if (p == dest_partition) {
          int port_select = f->id % gMultiportPorts;
          out_port = base_port + partition_offset + port_select;
          break;
        }
        
        partition_offset += gMultiportPorts;
      }
    }
  }
  
  if (f->watch) {
    *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
               << "Routing flit " << f->id
               << " -> out_port=" << out_port
               << endl;
  }
  
  outputs->Clear();
  outputs->AddRange(out_port, vcBegin, vcEnd);
}