#include "booksim.hpp"
#include <vector>
#include <sstream>
#include <cmath>

#include "gpunet.hpp"
#include "misc_utils.hpp"
#include "globals.hpp"

int gX; // # of partition crossbars
vector<int> gU; // units per layer

GPUNet::GPUNet( const Configuration& config, const string & name )
: Network ( config, name )
{
  _ComputeSize(config);
  _Alloc();
  _BuildNet(config);
}

void GPUNet::_ComputeSize( const Configuration& config )
{
  // Number of layers
  _l = config.GetInt("l");
  _use_partition = (config.GetInt("use_partition") == 1);
  
  // Nodes (shader and L2 slices)
  _n_shader = config.GetInt("shader");
  _n_l2slice = config.GetInt("l2slice");
  
  _nodes = _n_shader + _n_l2slice;

  _ratio = config.GetIntArray("units");
  if (_ratio.empty() || (_ratio.size() < size_t(_l))) {
    _ratio.resize(_l, 1);
  }
  
  _total_units.resize(_l);
  for (int l = 0; l < _l; ++l) {
    if (l == 0)
      _total_units[l] = _n_shader / _ratio[l];
    else
      _total_units[l] = _total_units[l - 1] / _ratio[l];
  }

  _offsets.resize(_l, 0);
  for (int l = 1; l < _l; ++l) {
    _offsets[l] = _total_units[l - 1] + _offsets[l - 1];
  }

  // Routers for shader-to-L2 Request and Reply Network
  _size = 0;
  for (int l = 0; l < _l; l++) {
    _size += 2 * _total_units[l];
  }

  // Channels for shader-to-L2 Network
  _channels = 0;
  _n_partition = 1;
  for (int l = 0; l < _l; l++) {
    if (l < _l - 1) {
      _channels += 2 * _total_units[l];
    } else {
      // Fully-connected partitioned crossbars
      _n_partition = _use_partition ? _total_units[l] : 1;
      _channels += 2 * _n_partition * (_n_partition - 1);
    }
  }

  _l2slice_per_partition = _n_l2slice / _n_partition; // L2 slices per partition

  _speedups = config.GetIntArray("speedups");
  if (_speedups.empty() || (_speedups.size() < size_t(_l + 1))) {
    _speedups.resize(_l + 1, 1);
  }

  _interpartition_speedup = config.GetInt("interpartition_speedup");

  gN = _l;
  gX = _n_partition;
  gU = _ratio;
}

void GPUNet::_BuildNet(const Configuration& config)
{
  ostringstream name;
  int c, id;

  // STEP 1: Create all routers first
  for (int l = 0; l < _l; ++l) {
    for (int addr = 0; addr < _total_units[l]; ++addr) {
      id = _offsets[l] + addr;
      
      int bottom_ports = (l < _l - 1) ? _ratio[l] : (_ratio[l] + (_n_partition - 1));
      int top_ports = (l < _l - 1) ? 1 : (_l2slice_per_partition + (_n_partition - 1));

      name.str("");
      name << "router_" << "request" << "_" << l << "_" << addr;
      _routers[id] = Router::NewRouter(config, this, name.str(), id, bottom_ports, top_ports, (l == _l - 1) ? _interpartition_speedup : _speedups[l + 1]);
      _timed_modules.push_back(_routers[id]);

      name.str("");
      name << "router_" << "reply" << "_" << l << "_" << addr;
      _routers[id + _size / 2] = Router::NewRouter(config, this, name.str(), id + _size / 2, top_ports, bottom_ports, (l == _l - 1) ? _interpartition_speedup : _speedups[l + 1]);
      _timed_modules.push_back(_routers[id + _size / 2]);
    }
  }
  
  // STEP 2: Connect shader->TPC (injection) and TPC->shader (ejection) first
  for (int addr = 0; addr < _total_units[0]; ++addr) {
    for (int port = 0; port < _ratio[0]; ++port) {
      id = _offsets[0] + addr;
      c = addr * _ratio[0] + port;  // shader node index

      // Request network: shader -> TPC (injection)
      _routers[id]->AddInputChannel(_inject[c], _inject_cred[c]);
      
      // Reply network: TPC -> shader (ejection)
      _routers[id + _size / 2]->AddOutputChannel(_eject[c], _eject_cred[c]);
    }
  }
  
  // STEP 3: Connect L2->Crossbar (injection) and Crossbar->L2 (ejection)
  for (int addr = 0; addr < _total_units[_l - 1]; ++addr) {
    for (int port = 0; port < _l2slice_per_partition; ++port) {
      id = _offsets[_l - 1] + addr;
      c = _n_shader + addr * _l2slice_per_partition + port;  // L2 node index
      
      // Request network: Crossbar -> L2 (ejection)
      _routers[id]->AddOutputChannel(_eject[c], _eject_cred[c]);
      
      // Reply network: L2 -> Crossbar (injection)
      _routers[id + _size / 2]->AddInputChannel(_inject[c], _inject_cred[c]);
    }
  }
  
  // STEP 4: Connect internal network channels
  // 4.1: Connect Request Network internal channels
  for (int l = 0; l < _l; ++l) {
    for (int addr = 0; addr < _total_units[l]; ++addr) {
      id = _offsets[l] + addr;
      
      // Connect bottom channels (from lower layer)
      if (l > 0) {
        for (int port = 0; port < _ratio[l]; ++port) {
          c = _offsets[l - 1] + addr * _ratio[l] + port;
          _routers[id]->AddInputChannel(_chan[c], _chan_cred[c]);
        }
      }

      // Connect top channels (to higher layer)
      if (l < _l - 1) {
        c = _offsets[l] + addr;
        _routers[id]->AddOutputChannel(_chan[c], _chan_cred[c]);
      }

      // Connect inter-partition channels for the last layer
      if (l == _l - 1) {
        for (int port = 0; port < (_n_partition - 1); ++port) {
          int src_partition = port;
          if (src_partition >= addr) {
            src_partition++;
          }

          int src_outport = addr;
          if (src_outport > src_partition) {
            src_outport--;
          }
          
          // For each partition router, output and input channels
          // are connected in sequential port order.
          // ports for inter-partition channels start after intra-partition channels
          // Output channel
          c = _offsets[l] + addr * (_n_partition - 1) + port;
          _routers[id]->AddOutputChannel(_chan[c], _chan_cred[c]);

#ifdef GPUNET_DEBUG
          cout << "Connecting inter-partition channel " << c
               << " as an output chnanel of partition " << addr
               << " through outport " << port << endl;
#endif
          
          // Input channel
          c = _offsets[l] + src_partition * (_n_partition - 1) + src_outport;
          _routers[id]->AddInputChannel(_chan[c], _chan_cred[c]);

#ifdef GPUNET_DEBUG
          cout << "Connecting inter-partition channel " << c
               << " as an input channel from partition " << src_partition
               << " using outport " << src_outport
               << " to partition " << addr
               << " through inport " << port << endl;
#endif

        }
      }
    }
  }
  
  // 4.2: Connect Reply Network internal channels
  for (int l = 0; l < _l; ++l) {
    for (int addr = 0; addr < _total_units[l]; ++addr) {
      id = _offsets[l] + addr + _size / 2;
      
      // Connect bottom channels (to lower layer)
      if (l > 0) {
        for (int port = 0; port < _ratio[l]; ++port) {
          c = _offsets[l - 1] + addr * _ratio[l] + port + _channels / 2;
          _routers[id]->AddOutputChannel(_chan[c], _chan_cred[c]);
        }
      }

      // Connect top channels (from higher layer)
      if (l < _l - 1) {
        c = _offsets[l] + addr + _channels / 2;
        _routers[id]->AddInputChannel(_chan[c], _chan_cred[c]);
      }
      
      // Connect inter-partition channels for the last layer
      if (l == _l - 1) {
        for (int port = 0; port < (_n_partition - 1); ++port) {
          int src_partition = port;
          if (src_partition >= addr) {
            src_partition++;
          }

          int src_outport = addr;
          if (src_outport > src_partition) {
            src_outport--;
          }

          // For each partition router, output and input channels
          // are connected in sequential port order.
          // ports for inter-partition channels start after intra-partition channels
          // Output channel
          c = _offsets[l] + addr * (_n_partition - 1) + port + _channels / 2;
          _routers[id]->AddOutputChannel(_chan[c], _chan_cred[c]);
          
#ifdef GPUNET_DEBUG
          cout << "Connecting inter-partition channel " << c
               << " as an output chnanel of partition " << addr
               << " through outport " << port << endl;
#endif
          // Input channel
          c = _offsets[l] + src_partition * (_n_partition - 1) + src_outport + _channels / 2;
          _routers[id]->AddInputChannel(_chan[c], _chan_cred[c]);

#ifdef GPUNET_DEBUG
          cout << "Connecting inter-partition channel " << c
               << " as an input channel from partition " << src_partition
               << " using outport " << src_outport
               << " to partition " << addr
               << " through inport " << port << endl;
#endif

        }
      }
    }
  }

  _SetupChannels();
}

// Set up all channel properties
void GPUNet::_SetupChannels()
{
  // Injection and Ejection channels
  for (int i = 0; i < _n_shader; i++) {
    _SetChannelProperties(_inject[i], _inject_cred[i], 0);
    _SetChannelProperties(_eject[i], _eject_cred[i], 0);
  }

  for (int i = _n_shader; i < _nodes; i++) {
    _SetChannelProperties(_inject[i], _inject_cred[i], _l);
    _SetChannelProperties(_eject[i], _eject_cred[i], _l);
  }
  
  for (int l = 1; l < _l; l++) {
    int start = _offsets[l - 1];
    int end = _offsets[l];
    
    // TPC <-> CPC, CPC <-> GPC, GPC <-> Crossbar channels
    for (int c = start; c < end; c++) {
      _SetChannelProperties(_chan[c], _chan_cred[c], l);
      _SetChannelProperties(_chan[c + _channels / 2], _chan_cred[c + _channels / 2], l);
    }
  }
  
  // interpartition channels for the last layer
  if (_use_partition) {
    int p_start = _offsets[_l - 1];
    int p_end = _offsets[_l - 1] + _n_partition * (_n_partition - 1);
    for (int c = p_start; c < p_end; c++) {
      _SetChannelProperties(_chan[c], _chan_cred[c], _l - 1, true);
      _SetChannelProperties(_chan[c + _channels / 2], _chan_cred[c + _channels / 2], _l - 1, true);
    }
  }
}

// Set channel latency and bandwidth based on layer properties
void GPUNet::_SetChannelProperties(FlitChannel* channel, CreditChannel* credit_channel, int l, bool is_interpartition)
{
  int latency = _GetWireLatency(l, is_interpartition);
  int bandwidth = _GetChannelBandwidth(l, is_interpartition);
  
  channel->SetLatency(latency);
  channel->SetBandwidth(bandwidth);
  credit_channel->SetLatency(latency);
  credit_channel->SetBandwidth(bandwidth);
}

int GPUNet::_GetWireLatency(int l, bool is_interpartition) const
{
  // Using an arithmetic progression for latencies.
  // S = n/2 * (2a + (n-1)d)
  // where, S = target total latency, n = _l+1, d = latency increment, a = base latency, 
  // Latency for channel to layer 'l' is a + l * d

  // Adjust latency_increment to make the total latency without partition
  // similar to the average of remote and local access latencies with partition.
  double s = _use_partition ? 75.0 : 75.0;
  int n = _l + 1;
  double d = 10.0;
  // Ensure d is less than 2S/(n(n-1)) to avoid negative local latency.
  assert(d < 2.0 * s / (n * (n - 1)));
  
  double a = s / n - (n - 1) * d / 2.0;
  
  int local_latency = round(a + l * d);
  const int remote_latency = 50;

  return is_interpartition ? remote_latency : local_latency;
}

int GPUNet::_GetChannelBandwidth(int l, bool is_interpartition) const
{
 return is_interpartition ? _interpartition_speedup : _speedups[l];
}


void GPUNet::RegisterRoutingFunctions()
{
  gRoutingFunctionMap["hierarchical_gpunet"] = &hierarchical_gpunet;
}

void hierarchical_gpunet(const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject)
{
  int vcBegin = 0, vcEnd = gNumVCs - 1;

  int out_port;
  
  if (inject) {
    outputs->AddRange(-1, vcBegin, vcEnd);
    return;
  }

  int src = f->src;
  int dest = f->dest;
  int hops = f->hops;

  // # of total shader nodes
  int shader = gX;
  for (int i = 0; i < gN; ++i) {
    shader *= gU[i];
  }
  int shader_per_partition = shader / gX;
  int l2slice = gNodes - shader;
  int l2slice_per_partition = l2slice / gX;

  assert((src < shader && dest >= shader) || (src >= shader && dest < shader));
  bool is_request = dest > src;

  int src_partition = is_request ? (src / shader_per_partition) : ((src - shader) / l2slice_per_partition);
  int dest_partition = is_request ? ((dest - shader) / l2slice_per_partition) : (dest / shader_per_partition);
  bool is_remote = (dest_partition != src_partition);

  // if remote partition access is required, add one additional hop
  // in the case of fully-connected partitioned network
  int total_hops = is_remote ? (gN + 1) : gN;
  int l = is_request ? hops : (total_hops - hops - 1);

  if (is_request) {
    // request network
    if (l < gN - 1) {
      out_port = 0;
    } else {
      // partition layer
      assert(r->NumOutputs() == l2slice_per_partition + (gX - 1));

      if (is_remote && (l == total_hops - 2)) {
        // inter-partition communication
        int dest_port = (dest_partition > src_partition) ? (dest_partition - 1) : dest_partition;
        out_port = l2slice_per_partition + dest_port;
      } else {
        // intra-partition communication or already crossed inter-partition
        out_port = (dest - shader) % l2slice_per_partition;
      }
    }
  } else {
    // reply network
    if (l < gN - 1) {
      int shader_cluster = 1;
      for (int i = 0; i < l; ++i) {
        shader_cluster *= gU[i];
      }
      out_port = (dest % (shader_cluster * gU[l])) / shader_cluster;
    } else {
      // partition layer
      assert(r->NumInputs() == l2slice_per_partition + (gX - 1));

      if (is_remote && (l == total_hops - 1)) {
        // inter-partition communication
        int dest_port = (dest_partition > src_partition) ? (dest_partition - 1) : dest_partition; 
        out_port = gU[gN - 1] + dest_port;
      } else {
        // intra-partition communication or already crossed inter-partition
        out_port = (dest % shader_per_partition) / (shader_per_partition / gU[gN - 1]);
      }
    }
  }

  if (f->watch) {
    *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
               << "Adding VC range ["
               << vcBegin << ","
               << vcEnd << "]"
               << " at output port " << out_port
               << " for flit " << f->id
               << " (input port " << in_channel
               << ", destination " << f->dest << ")"
               << "." << endl;
  }

  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}