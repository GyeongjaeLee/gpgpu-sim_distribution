#include <algorithm>
#include <cassert>
#include <cmath>
#include <sstream>
#include <vector>
#include "booksim.hpp"

#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include "globals.hpp"
#include "gpunet.hpp"
#include "misc_utils.hpp"
#include "random_utils.hpp"

using namespace std;

// ================================================================
// Routing globals — populated by _ComputeSize, read by routing function
// ================================================================

// SM-side hierarchy
static int gL;              // intermediate SM levels
static vector<int> gUnits;  // [l]: SMs/TPCs per parent at each level
static vector<int> gPorts;  // [l]: parallel ports on upward link at each level
static vector<int> gSmTotal;  // [l]: router count per level per partition
static vector<int>
    gSmSpan;  // [l]: leaf SM nodes covered per router at each level
static vector<int>
    gSmOffsets;  // [l]: router ID offset per level within partition

// L2-side hierarchy
static int gL2L;
static vector<int> gL2Units;    // [l2_l]
static vector<int> gL2Ports;    // [l2_l]
static vector<int> gL2Total;    // [l2_l]
static vector<int> gL2Span;     // [l2_l]
static vector<int> gL2Offsets;  // [l2_l]

// Crossbar port layout
//   Inputs:  SM-side ports [0 .. gXbarSmPorts-1]
//            interpartition [gXbarSmPorts .. gXbarSmPorts+gXbarInterPorts-1]
//   Outputs: L2-side ports [0 .. gXbarL2Ports-1]
//            interpartition [gXbarL2Ports .. gXbarL2Ports+gXbarInterPorts-1]
static int gXbarOffset;
static int gXbarSmPorts;
static int gXbarL2Ports;
static int gXbarInterPorts;

// Partition
static int gNPartition;
static int gInterpartitionPorts;

// Node counts
int gNSM;        // total SM nodes
int gSMPerPart;  // SM nodes per partition
int gL2PerPart;  // L2 nodes per partition
static int gRoutersPerPartition;

// Credit threshold for port selection
static int gMaxCredit;

// ================================================================
// selectPort: random port selection with credit-aware fallback
// ================================================================
static int selectPort(const Router* r, int base_port, int num_ports) {
  if (num_ports <= 1) return base_port;

  int start = RandomInt(num_ports - 1);
  int selected = base_port + start;

  for (int i = 0; i < num_ports; i++) {
    int port = base_port + (start + i) % num_ports;
    if (r->GetUsedCredit(port) < gMaxCredit) return port;
    selected = port;
  }
  return selected;  // all ports saturated; use last checked
}

// ================================================================
// GPUNet
// ================================================================

int GPUNet::_GetNumChildren(int parent_idx, int total_children, int ratio) {
  int start = parent_idx * ratio;
  if (start >= total_children) return 0;
  return min(ratio, total_children - start);
}

GPUNet::GPUNet(const Configuration& config, const string& name)
    : Network(config, name) {
  _ComputeSize(config);
  _Alloc();
  _BuildNet(config);
}

void GPUNet::_ComputeSize(const Configuration& config) {
  _l = config.GetInt("l");
  _l2_l = config.GetInt("l2_l");
  _use_partition = (config.GetInt("use_partition") == 1);
  _n_sm = config.GetInt("n_sm");
  _n_l2 = config.GetInt("n_l2");
  _nodes = _n_sm + _n_l2;

  _n_partition = _use_partition ? config.GetInt("n_partition") : 1;
  _sm_per_partition = _n_sm / _n_partition;
  _l2_per_partition = _n_l2 / _n_partition;

  // SM-side hierarchy
  _units = config.GetIntArray("units");
  _ports = config.GetIntArray("ports");
  if ((int)_units.size() < _l) _units.resize(_l, 1);
  if ((int)_ports.size() < _l) _ports.resize(_l, 1);

  _sm_total.resize(_l);
  for (int i = 0; i < _l; i++) {
    int children = (i == 0) ? _sm_per_partition : _sm_total[i - 1];
    _sm_total[i] = (children + _units[i] - 1) / _units[i];  // ceil
  }

  // L2-side hierarchy
  _l2_units = config.GetIntArray("l2_units");
  _l2_ports = config.GetIntArray("l2_ports");
  if ((int)_l2_units.size() < _l2_l) _l2_units.resize(_l2_l, 1);
  if ((int)_l2_ports.size() < _l2_l) _l2_ports.resize(_l2_l, 1);

  _l2_total.resize(_l2_l);
  for (int i = 0; i < _l2_l; i++) {
    int children = (i == 0) ? _l2_per_partition : _l2_total[i - 1];
    _l2_total[i] = (children + _l2_units[i] - 1) / _l2_units[i];
  }

  // Router ID offsets within one partition block
  int offset = 0;
  _sm_offsets.resize(_l);
  for (int i = 0; i < _l; i++) {
    _sm_offsets[i] = offset;
    offset += _sm_total[i];
  }
  _xbar_offset = offset++;
  _l2_offsets.resize(_l2_l);
  for (int i = 0; i < _l2_l; i++) {
    _l2_offsets[i] = offset;
    offset += _l2_total[i];
  }
  _routers_per_partition = offset;

  // Total routers: partitions * _routers_per_partition
  _size = _n_partition * _routers_per_partition;

  // Channel counts
  _interpartition_ports = (_use_partition && _n_partition > 1)
                              ? config.GetInt("interpartition_ports")
                              : 0;
  if (_interpartition_ports < 0) _interpartition_ports = 0;

  int sm_ch_per_part = 0;
  for (int i = 0; i < _l; i++) sm_ch_per_part += _sm_total[i] * _ports[i];
  int l2_ch_per_part = 0;
  for (int i = 0; i < _l2_l; i++) l2_ch_per_part += _l2_total[i] * _l2_ports[i];
  int inter_ch_total =
      (_n_partition * (_n_partition - 1) * _interpartition_ports) / 2;
  int one_dir =
      _n_partition * (sm_ch_per_part + l2_ch_per_part) + inter_ch_total;
  _channels = 2 * one_dir;

  // Latency arrays
  _interpartition_latency = config.GetInt("interpartition_latency");

  _inj_lat = config.GetInt("inj_lat");
  _ej_lat = config.GetInt("ej_lat");

  _latency = config.GetIntArray("latency");
  if ((int)_latency.size() < _l) _latency.resize(_l, 1);
  _l2_latency = config.GetIntArray("l2_latency");
  if ((int)_l2_latency.size() < _l2_l) _l2_latency.resize(_l2_l, 1);

  // Populate routing globals
  gL = _l;
  gL2L = _l2_l;
  gNPartition = _n_partition;
  gInterpartitionPorts = _interpartition_ports;
  gNSM = _n_sm;
  gSMPerPart = _sm_per_partition;
  gL2PerPart = _l2_per_partition;
  gRoutersPerPartition = _routers_per_partition;
  gUnits = _units;
  gSmTotal = _sm_total;
  gSmOffsets = _sm_offsets;
  gPorts = _ports;
  gL2Units = _l2_units;
  gL2Total = _l2_total;
  gL2Offsets = _l2_offsets;
  gL2Ports = _l2_ports;
  gXbarOffset = _xbar_offset;
  gXbarSmPorts =
      (_l > 0) ? (_sm_total[_l - 1] * _ports[_l - 1]) : _sm_per_partition;
  gXbarL2Ports = (_l2_l > 0) ? (_l2_total[_l2_l - 1] * _l2_ports[_l2_l - 1])
                             : _l2_per_partition;
  gXbarInterPorts = (_n_partition - 1) * _interpartition_ports;

  gSmSpan.resize(_l);
  for (int i = 0; i < _l; i++)
    gSmSpan[i] = (i == 0) ? _units[0] : gSmSpan[i - 1] * _units[i];

  gL2Span.resize(_l2_l);
  for (int i = 0; i < _l2_l; i++)
    gL2Span[i] = (i == 0) ? _l2_units[0] : gL2Span[i - 1] * _l2_units[i];

  gMaxCredit = config.GetInt("num_vcs") * config.GetInt("vc_buf_size");
  gN = _l;  // backward compatibility
}

void GPUNet::_BuildNet(const Configuration& config) {
  ostringstream name;

  int cid = 0;

  // add_link: add bidirectional channels between r1 and r2
  auto add_link = [&](int r1, int r2, int lat) {
    // r1 -> r2 (UP for child->parent, or interpartition)
    _chan[cid]->SetLatency(lat);
    _chan_cred[cid]->SetLatency(lat);
    _routers[r1]->AddOutputChannel(_chan[cid], _chan_cred[cid]);
    _routers[r2]->AddInputChannel(_chan[cid], _chan_cred[cid]);
    ++cid;

    // r2 -> r1 (DOWN for parent->child, or interpartition)
    _chan[cid]->SetLatency(lat);
    _chan_cred[cid]->SetLatency(lat);
    _routers[r2]->AddOutputChannel(_chan[cid], _chan_cred[cid]);
    _routers[r1]->AddInputChannel(_chan[cid], _chan_cred[cid]);
    ++cid;
  };

  // connect_tree: wire all links for one side of the hierarchy (SM or L2)
  auto connect_tree = [&](int levels, const vector<int>& lvl_total,
                          const vector<int>& lvl_units,
                          const vector<int>& lvl_ports,
                          const vector<int>& lvl_offsets,
                          const vector<int>& lvl_latency, int xbar_id,
                          int part_base) {
    for (int i = 0; i < levels; i++) {
      bool at_top = (i == levels - 1);
      for (int pj = 0; pj < (at_top ? 1 : lvl_total[i + 1]); pj++) {
        int parent_id =
            at_top ? xbar_id : (part_base + lvl_offsets[i + 1] + pj);

        int num_ch = at_top
                         ? lvl_total[i]
                         : _GetNumChildren(pj, lvl_total[i], lvl_units[i + 1]);

        for (int c = 0; c < num_ch; c++) {
          int child_j = at_top ? c : (pj * lvl_units[i + 1] + c);
          int child_id = part_base + lvl_offsets[i] + child_j;
          for (int k = 0; k < lvl_ports[i]; k++) {
            add_link(child_id, parent_id, lvl_latency[i]);
          }
        }
      }
    }
  };

  // ============================================================
  // STEP 1: Create all routers
  // ============================================================
  for (int p = 0; p < _n_partition; p++) {
    int part_base = p * _routers_per_partition;

    // SM-side routers (level 0 = TPC, closest to SMs)
    for (int i = 0; i < _l; i++) {
      for (int j = 0; j < _sm_total[i]; j++) {
        int num_ch =
            _GetNumChildren(j, (i == 0) ? _sm_per_partition : _sm_total[i - 1],
                            (i == 0) ? _units[0] : _units[i]);
        int down_ports = (i == 0) ? num_ch : (num_ch * _ports[i - 1]);
        int up_ports = _ports[i];
        int id = part_base + _sm_offsets[i] + j;

        name.str("");
        name << "sm_p" << p << "_l" << i << "_" << j;
        _routers[id] =
            Router::NewRouter(config, this, name.str(), id,
                              down_ports + up_ports, down_ports + up_ports);
        _timed_modules.push_back(_routers[id]);
      }
    }

    // Crossbar router
    {
      int xbar_sm =
          (_l > 0) ? (_sm_total[_l - 1] * _ports[_l - 1]) : _sm_per_partition;
      int xbar_l2 = (_l2_l > 0) ? (_l2_total[_l2_l - 1] * _l2_ports[_l2_l - 1])
                                : _l2_per_partition;
      int xbar_inter = (_n_partition - 1) * _interpartition_ports;
      int xbar_id = part_base + _xbar_offset;

      name.str("");
      name << "xbar_p" << p;
      int total_ports = xbar_sm + xbar_l2 + xbar_inter;
      _routers[xbar_id] = Router::NewRouter(config, this, name.str(), xbar_id,
                                            total_ports, total_ports);
      _timed_modules.push_back(_routers[xbar_id]);
    }

    // L2-side routers (level 0 = L2Slice, closest to L2 nodes)
    for (int i = 0; i < _l2_l; i++) {
      for (int j = 0; j < _l2_total[i]; j++) {
        int num_ch =
            _GetNumChildren(j, (i == 0) ? _l2_per_partition : _l2_total[i - 1],
                            (i == 0) ? _l2_units[0] : _l2_units[i]);
        int down_ports = (i == 0) ? num_ch : (num_ch * _l2_ports[i - 1]);
        int up_ports = _l2_ports[i];
        int id = part_base + _l2_offsets[i] + j;

        name.str("");
        name << "l2_p" << p << "_l" << i << "_" << j;
        _routers[id] =
            Router::NewRouter(config, this, name.str(), id,
                              down_ports + up_ports, down_ports + up_ports);
        _timed_modules.push_back(_routers[id]);
      }
    }
  }

  // ============================================================
  // STEP 2: Inject/eject channels
  // ============================================================

  for (int p = 0; p < _n_partition; p++) {
    int part_base = p * _routers_per_partition;

    // SM nodes
    for (int s = 0; s < _sm_per_partition; s++) {
      int sm_node = p * _sm_per_partition + s;
      int target = (_l > 0) ? (part_base + _sm_offsets[0] + s / _units[0])
                            : (part_base + _xbar_offset);

      _inject[sm_node]->SetLatency(_inj_lat);
      _inject_cred[sm_node]->SetLatency(_inj_lat);
      _eject[sm_node]->SetLatency(_ej_lat);
      _eject_cred[sm_node]->SetLatency(_ej_lat);

      _routers[target]->AddInputChannel(_inject[sm_node],
                                        _inject_cred[sm_node]);
      _routers[target]->AddOutputChannel(_eject[sm_node], _eject_cred[sm_node]);
    }

    // L2 nodes
    for (int s = 0; s < _l2_per_partition; s++) {
      int l2_node = _n_sm + p * _l2_per_partition + s;
      int target = (_l2_l > 0) ? (part_base + _l2_offsets[0] + s / _l2_units[0])
                               : (part_base + _xbar_offset);

      _inject[l2_node]->SetLatency(_inj_lat);
      _inject_cred[l2_node]->SetLatency(_inj_lat);
      _eject[l2_node]->SetLatency(_ej_lat);
      _eject_cred[l2_node]->SetLatency(_ej_lat);

      _routers[target]->AddInputChannel(_inject[l2_node],
                                        _inject_cred[l2_node]);
      _routers[target]->AddOutputChannel(_eject[l2_node], _eject_cred[l2_node]);
    }
  }

  // ============================================================
  // STEP 3 & 4: SM tree channels and L2 tree channels
  // ============================================================
  for (int p = 0; p < _n_partition; p++) {
    int part_base = p * _routers_per_partition;
    int xbar_id = part_base + _xbar_offset;

    // Connect SM-side hierarchy
    connect_tree(_l, _sm_total, _units, _ports, _sm_offsets, _latency, xbar_id,
                 part_base);

    // Connect L2-side hierarchy
    connect_tree(_l2_l, _l2_total, _l2_units, _l2_ports, _l2_offsets,
                 _l2_latency, xbar_id, part_base);
  }

  // ============================================================
  // STEP 5: Interpartition channels (crossbar <-> crossbar)
  // ============================================================
  if (_use_partition && _interpartition_ports > 0) {
    for (int sp = 0; sp < _n_partition; sp++) {
      int src_id = sp * _routers_per_partition + _xbar_offset;

      for (int dp = 0; dp < _n_partition; dp++) {
        if (sp == dp) continue;
        int dst_id = dp * _routers_per_partition + _xbar_offset;

        // add_link adds BOTH directions (src->dest and dest->src).
        if (sp < dp) {
          for (int k = 0; k < _interpartition_ports; k++) {
            add_link(src_id, dst_id, _interpartition_latency);
          }
        }
      }
    }
  }

  assert(cid == _channels);
}

// ================================================================
// Routing function registration
// ================================================================

void GPUNet::RegisterRoutingFunctions() {
  gRoutingFunctionMap["hierarchical_gpunet"] = &hierarchical_gpunet;
}

// ================================================================
// hierarchical_gpunet routing function
// ================================================================

void hierarchical_gpunet(const Router* r, const Flit* f, int in_channel,
                         OutputSet* outputs, bool inject) {
  int vcBegin = 0, vcEnd = gNumVCs - 1;

  if (inject) {
    outputs->AddRange(-1, vcBegin, vcEnd);
    return;
  }

  int dest = f->dest;
  int router_id = r->GetID();

  // Decode router position
  int cur_partition = router_id / gRoutersPerPartition;
  int local_offset = router_id % gRoutersPerPartition;

  // Identify SM and L2 endpoints and their partitions
  bool dest_is_sm = (dest < gNSM);
  int target_part =
      dest_is_sm ? (dest / gSMPerPart) : ((dest - gNSM) / gL2PerPart);

  // Identify which side/level the current router is in
  int side = -1, level = -1, addr = -1;

  for (int i = 0; i < gL && side < 0; i++) {
    if (local_offset >= gSmOffsets[i] &&
        local_offset < gSmOffsets[i] + gSmTotal[i]) {
      side = 0;
      level = i;
      addr = local_offset - gSmOffsets[i];
    }
  }
  if (side < 0 && local_offset == gXbarOffset) {
    side = 1;
  }
  for (int i = 0; i < gL2L && side < 0; i++) {
    if (local_offset >= gL2Offsets[i] &&
        local_offset < gL2Offsets[i] + gL2Total[i]) {
      side = 2;
      level = i;
      addr = local_offset - gL2Offsets[i];
    }
  }
  assert(side >= 0);

  int out_port = -1;

  if (side == 0) {
    // ===== SM-side router =====
    int num_ch = GPUNet::_GetNumChildren(
        addr, (level == 0) ? gSMPerPart : gSmTotal[level - 1], gUnits[level]);
    int num_down = (level == 0) ? num_ch : (num_ch * gPorts[level - 1]);
    int num_up = gPorts[level];

    if (!dest_is_sm) {
      // Dest is L2 -> go UP
      out_port = selectPort(r, num_down, num_up);
    } else {
      // Dest is SM -> go DOWN if descendant, else UP
      int local_sm = dest % gSMPerPart;
      int start_sm = addr * gSmSpan[level];
      if (target_part != cur_partition) {
        out_port = selectPort(r, num_down, num_up);
      } else if (local_sm >= start_sm && local_sm < start_sm + gSmSpan[level]) {
        if (level == 0) {
          out_port = local_sm - start_sm;
        } else {
          int child_abs = local_sm / gSmSpan[level - 1];
          int child_local = child_abs - addr * gUnits[level];
          out_port =
              selectPort(r, child_local * gPorts[level - 1], gPorts[level - 1]);
        }
      } else {
        out_port = selectPort(r, num_down, num_up);
      }
    }
  } else if (side == 1) {
    // ===== Crossbar =====
    if (target_part != cur_partition) {
      // Route to another partition
      int idx = (target_part < cur_partition) ? target_part : target_part - 1;
      out_port = selectPort(
          r, gXbarSmPorts + gXbarL2Ports + idx * gInterpartitionPorts,
          gInterpartitionPorts);
    } else {
      if (!dest_is_sm) {
        // Route to local L2 side
        int local_l2 = (dest - gNSM) % gL2PerPart;
        if (gL2L > 0) {
          int top_child = local_l2 / gL2Span[gL2L - 1];
          out_port =
              selectPort(r, gXbarSmPorts + top_child * gL2Ports[gL2L - 1],
                         gL2Ports[gL2L - 1]);
        } else {
          out_port = gXbarSmPorts + local_l2;
        }
      } else {
        // Route to local SM side
        int local_sm = dest % gSMPerPart;
        if (gL > 0) {
          int top_child = local_sm / gSmSpan[gL - 1];
          out_port = selectPort(r, top_child * gPorts[gL - 1], gPorts[gL - 1]);
        } else {
          out_port = local_sm;
        }
      }
    }
  } else if (side == 2) {
    // ===== L2-side router =====
    int num_ch = GPUNet::_GetNumChildren(
        addr, (level == 0) ? gL2PerPart : gL2Total[level - 1], gL2Units[level]);
    int num_down = (level == 0) ? num_ch : (num_ch * gL2Ports[level - 1]);
    int num_up = gL2Ports[level];

    if (dest_is_sm) {
      // Dest is SM -> go UP
      out_port = selectPort(r, num_down, num_up);
    } else {
      // Dest is L2 -> go DOWN if descendant, else UP
      int local_l2 = (dest - gNSM) % gL2PerPart;
      int start_l2 = addr * gL2Span[level];
      if (target_part != cur_partition) {
        out_port = selectPort(r, num_down, num_up);
      } else if (local_l2 >= start_l2 && local_l2 < start_l2 + gL2Span[level]) {
        if (level == 0) {
          out_port = local_l2 - start_l2;
        } else {
          int child_abs = local_l2 / gL2Span[level - 1];
          int child_local = child_abs - addr * gL2Units[level];
          out_port = selectPort(r, child_local * gL2Ports[level - 1],
                                gL2Ports[level - 1]);
        }
      } else {
        out_port = selectPort(r, num_down, num_up);
      }
    }
  }

  assert(out_port >= 0);

  if (f->watch) {
    *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
               << "Routing flit " << f->id << " (src=" << f->src
               << ", dest=" << dest << ")"
               << " -> out_port=" << out_port << endl;
  }

  outputs->Clear();
  outputs->AddRange(out_port, vcBegin, vcEnd);
}

// ================================================================
// GPUNet Link Utilization Statistics Tracking
// ================================================================

struct LinkStat {
  long long total = 0;
  long long window = 0;
  double peak = 0.0;
  int port_count = 0;
};

static std::map<string, LinkStat> gLinkStatsGlobal;
static std::map<string, LinkStat> gLinkStatsLocal;
static std::map<pair<int, int>, string> gLinkMapGlobal;
static std::map<pair<int, int>, string> gLinkMapLocal;

static long long gGPUNetCycles = 0;
static const int gGPUNetWindow = 200;

void gpunet_tick() {
  gGPUNetCycles++;
  if (gGPUNetCycles % gGPUNetWindow == 0) {
    // Calculate peak for global links
    for (auto& kv : gLinkStatsGlobal) {
      double util =
          (double)kv.second.window / (gGPUNetWindow * kv.second.port_count);
      if (util > kv.second.peak) kv.second.peak = util;
      kv.second.window = 0;
    }
    // Calculate peak for local links
    for (auto& kv : gLinkStatsLocal) {
      double util =
          (double)kv.second.window / (gGPUNetWindow * kv.second.port_count);
      if (util > kv.second.peak) kv.second.peak = util;
      kv.second.window = 0;
    }
  }
}

// Helper to deduce direction
static string GetLinkDirectionString(int src_r, int dst_r) {
  int src_part = src_r / gRoutersPerPartition;
  int dst_part = dst_r / gRoutersPerPartition;
  int src_off = src_r % gRoutersPerPartition;
  int dst_off = dst_r % gRoutersPerPartition;

  if (src_part != dst_part) return "Xbar -> Xbar (INT)";

  int src_side = -1, dst_side = -1;
  int src_level = -1, dst_level = -1;

  for (int i = 0; i < gL; i++) {
    if (src_off >= gSmOffsets[i] && src_off < gSmOffsets[i] + gSmTotal[i]) {
      src_side = 0;
      src_level = i;
    }
    if (dst_off >= gSmOffsets[i] && dst_off < gSmOffsets[i] + gSmTotal[i]) {
      dst_side = 0;
      dst_level = i;
    }
  }
  if (src_off == gXbarOffset) src_side = 1;
  if (dst_off == gXbarOffset) dst_side = 1;

  for (int i = 0; i < gL2L; i++) {
    if (src_off >= gL2Offsets[i] && src_off < gL2Offsets[i] + gL2Total[i]) {
      src_side = 2;
      src_level = i;
    }
    if (dst_off >= gL2Offsets[i] && dst_off < gL2Offsets[i] + gL2Total[i]) {
      dst_side = 2;
      dst_level = i;
    }
  }

  if (src_side == 0 && dst_side == 0)
    return (src_level < dst_level) ? "SM -> SM (UP)" : "SM -> SM (DOWN)";
  if (src_side == 0 && dst_side == 1) return "SM -> Xbar (UP)";
  if (src_side == 1 && dst_side == 0) return "Xbar -> SM (DOWN)";
  if (src_side == 1 && dst_side == 2) return "Xbar -> L2 (DOWN)";
  if (src_side == 2 && dst_side == 1) return "L2 -> Xbar (UP)";
  if (src_side == 2 && dst_side == 2)
    return (src_level < dst_level) ? "L2 -> L2 (UP)" : "L2 -> L2 (DOWN)";

  return "Unknown";
}

void gpunet_count_link_traversal(int src_router, int dst_router) {
  pair<int, int> p(src_router, dst_router);

  if (gLinkMapGlobal.find(p) == gLinkMapGlobal.end()) {
    string global_dir = GetLinkDirectionString(src_router, dst_router);
    // Since we don't know router names here without pointer, we rely on the
    // caller or just print IDs
    char local_name[128];
    snprintf(local_name, sizeof(local_name), "R%d -> R%d", src_router,
             dst_router);
    string local_dir = local_name;

    gLinkMapGlobal[p] = global_dir;
    gLinkMapLocal[p] = local_dir;

    // Increment port_count dynamically (each unique src->dst channel that is
    // traversed increments port_count if we assume each call for a NEW pair is
    // a distinct port... Wait, src->dst ID is unique per pair, but there might
    // be multiple ports. Actually, src_router and dst_router IDs uniquely
    // identify the PAIR of routers. Parallel ports have the SAME src and dst
    // router IDs. If we just count traversals between Router A and Router B,
    // and divide by gGPUNetCycles, we get flits/cycle between those two
    // routers. To get flits/cycle/port, we need to know the number of parallel
    // ports. This is complex to fetch dynamically. We can just print
    // "flits/cycle" and leave it as total bandwidth between the two entities!
    // This is much more accurate and less prone to division errors.
    if (gLinkStatsGlobal.find(global_dir) == gLinkStatsGlobal.end()) {
      gLinkStatsGlobal[global_dir].port_count = 1;  // Used as a dummy divisor
    }
    if (gLinkStatsLocal.find(local_dir) == gLinkStatsLocal.end()) {
      gLinkStatsLocal[local_dir].port_count = 1;
    }
  }

  string global_dir = gLinkMapGlobal[p];
  string local_dir = gLinkMapLocal[p];

  gLinkStatsGlobal[global_dir].total++;
  gLinkStatsGlobal[global_dir].window++;

  gLinkStatsLocal[local_dir].total++;
  gLinkStatsLocal[local_dir].window++;
}

void gpunet_print_link_stats() {
  cout << "===================================================================="
          "===="
       << endl;
  cout << " GPUNet Link Utilization Statistics (Total Cycles: " << gGPUNetCycles
       << ")" << endl;
  cout << "===================================================================="
          "===="
       << endl;

  cout << "--- GLOBAL DIRECTION STATS ---" << endl;
  for (auto const& kv : gLinkStatsGlobal) {
    double avg = (double)kv.second.total / gGPUNetCycles;
    cout << " [Global] " << left << setw(20) << kv.first << " | Avg: " << fixed
         << setprecision(4) << avg << " flits/cycle"
         << " | Peak (" << gGPUNetWindow << "c): " << kv.second.peak
         << " flits/cycle" << endl;
  }

  cout << "\n--- LOCAL ROUTER-TO-ROUTER STATS ---" << endl;
  for (auto const& kv : gLinkStatsLocal) {
    if (kv.second.total == 0) continue;  // skip unused
    double avg = (double)kv.second.total / gGPUNetCycles;
    cout << " [Local] " << left << setw(20) << kv.first << " | Avg: " << fixed
         << setprecision(4) << avg << " flits/cycle"
         << " | Peak (" << gGPUNetWindow << "c): " << kv.second.peak
         << " flits/cycle" << endl;
  }
  cout << "===================================================================="
          "===="
       << endl;
}
