# HBMNet AccelSim Experiment Guide

HBMNet 구조를 AccelSim(GPGPU-Sim) trace 기반 시뮬레이션으로 실험하고 결과를 분석하는 전체 파이프라인을 설명합니다.

---

## 디렉터리 구조

```
GPUNoC/
├── gpgpu-sim_distribution/          ← 이 파일이 위치한 곳
│   ├── experiments.csv              ← 실험할 구조·대역폭 정의
│   ├── hitrate.csv                  ← benchmark별 L2 hit rate
│   ├── avg_miss_rate.py             ← hit rate 측정 보조 스크립트
│   ├── experiments_loader.py        ← 공통 파라미터 로더 (모든 스크립트 공유)
│   ├── gen_accelsim_configs.py      ← [Step 1] 실험 config 파일 생성
│   ├── run_accelsim.py              ← [Step 2] 시뮬레이션 실행
│   ├── parse_accelsim_results.py    ← [Step 3] 로그 파싱 → CSV
│   ├── plot_accelsim_link_stats.py  ← [Step 4a] 링크 통계 3D bar 시각화
│   ├── plot_accelsim_perf.py        ← [Step 4b] 성능·L2 BW 비교 bar chart
│   └── configs/tested-cfgs/         ← 생성된 HW config 저장 위치
│
└── (AccelSim root)/
    ├── util/job_launching/
    │   ├── apps/define-all-apps.yml         ← benchmark 목록
    │   └── configs/define-standard-cfgs.yml ← config 목록 (자동 삽입됨)
    ├── hw_run/traces/device-0/12.8/         ← trace 파일 위치
    └── sim_run_12.8/                        ← 시뮬레이션 로그 출력 위치
```

---

## Step 0: 실험 파라미터 정의 (`experiments.csv`)

모든 structure/bandwidth 정보는 `experiments.csv` 한 파일에서 관리됩니다.
새로운 실험을 추가하려면 이 파일에 행을 추가하면 되고, 모든 스크립트에 자동으로 반영됩니다.

```
Structure   Core  HBM  SM_per_Xbar
B100_Local    1    4    74
H100          1    6    132
B100_Global   2    8    74
B100_Core_Rotate  2  12  74
Rubin_Ultra   4   16   112

Bandwidth         GPU-to-GPU  GPU-to-HBM  HBM-to-HBM  TSV   L2_per_HBM
B200+HBM3e        10          1           1            1     32
Rubin_Ultra+HBM4  10          2           2            2     64
Shoreline_1x      10          3.33        3.89         3.33  64
...
```

- **Structure**: Core = Xbar 수, HBM = 전체 HBM stack 수, SM_per_Xbar = Xbar당 SM 수
- **Bandwidth**: GPU-to-GPU, GPU-to-HBM 등 링크별 상대 대역폭 비율, L2_per_HBM = HBM stack당 L2 slice 수

---

## Step 1: Hit Rate 측정

`run_accelsim.py`는 시뮬레이션 전에 benchmark의 L2 hit rate를 `.icnt` config에 패치합니다.
이 값은 **기존 시뮬레이션 로그**에서 먼저 측정해야 합니다.

### 1-1. 기존 로그에서 hit rate 계산

`avg_miss_rate.py`를 로그 파일이 있는 경로로 복사하거나 경로를 직접 지정하여 실행합니다.

```bash
# 로그 파일을 직접 경로로 지정
python avg_miss_rate.py /path/to/sim_run_12.8/<app>/<args>/<config>-SASS/<app>*.o<N>

# 또는 현재 디렉터리에 로그가 있다면
python avg_miss_rate.py gpgpu-sim.log
```

출력 예시:

```
-------------------------------------------------------
Kernel     | L2 Accesses     | L2 Miss Rate
-------------------------------------------------------
Kernel 1   | 123456          | 0.2300
Kernel 2   | 234567          | 0.2100
...
-------------------------------------------------------

[결과] 전체 누적 L2 Accesses: 12345678
[결과] 실제 통합 L2 Miss Rate (Access 가중치): 0.2300 (23.00%)
```

### 1-2. hitrate.csv에 등록

`avg_miss_rate.py`가 출력하는 **Miss Rate**를 1에서 빼면 hit rate입니다.
결과를 `hitrate.csv`에 추가합니다 (탭 구분):

```
Benchmark                                   Hit_rate
rodinia-3.1:bfs-rodinia-3.1                 0.77
polybench:polybench-gemm                    0.99
GPU_Microbenchmark:mem_bw                   0
```

- Benchmark 이름은 AccelSim의 `define-all-apps.yml`의 `suite:app` 형식과 동일하게 작성
- Hit rate = 0 이면 baseline_ratio가 패치되지 않음 (miss-only 트래픽)
- 값이 없으면 `[WARN]`이 출력되고 패치 없이 실행됨

---

## Step 2: Benchmark 목록 확인

실험할 benchmark는 AccelSim의 `define-all-apps.yml`에서 확인합니다:

```bash
cat ../../util/job_launching/apps/define-all-apps.yml
```

형식 예시:

```yaml
rodinia-3.1:
  exec_dir: ...
  apps:
    bfs-rodinia-3.1:
      args:
        - ./data/graph1MW_6.txt
    nn-rodinia-3.1:
      args:
        - filelist_4 -r 5 -lat 30 -lng 90
```

`--benchmark` 인자로 사용할 때는 `suite:app` 형식으로 지정합니다:

- `rodinia-3.1:bfs-rodinia-3.1`
- `polybench:polybench-gemm`
- `Deepbench_nvidia_tencore:gemm_bench-tencore`
- `GPU_Microbenchmark:mem_bw`

---

## Step 3: Config 생성 (`gen_accelsim_configs.py`)

`experiments.csv`에 정의된 (structure × bandwidth × routing) 조합마다 `.icnt` 및 `gpgpusim.config` 파일을 생성합니다.
생성된 config는 `configs/tested-cfgs/SM100_<struct>_<bw>_<routing>/` 에 저장되고,
`define-standard-cfgs.yml`에도 자동 삽입됩니다.

```bash
# 특정 조합만 생성
python gen_accelsim_configs.py \
    --structure B100_Global \
    --bandwidth B200+HBM3e \
    --routing baseline min_adaptive near_min_adaptive \
    --near-min-p 0.0 1.0

# experiments.csv의 모든 조합 생성
python gen_accelsim_configs.py --all-configs --routing baseline near_min_adaptive

# 변경 사항만 미리 확인 (dry-run)
python gen_accelsim_configs.py --all-configs --routing baseline --dry-run
```

> **주의**: `run_accelsim.py` 실행 전에 반드시 먼저 실행해야 합니다.

---

## Step 4: 시뮬레이션 실행 (`run_accelsim.py`)

내부적으로 AccelSim의 `run_simulations.py`를 호출합니다.
각 (structure × bandwidth × routing) 조합에 대해 순차적으로 실행합니다.

```bash
# 단일 조합
python run_accelsim.py \
    --benchmark rodinia-3.1:bfs-rodinia-3.1 \
    --structure B100_Global \
    --bandwidth B200+HBM3e \
    --routing baseline min_adaptive near_min_adaptive \
    --near-min-p 0.0 1.0

# 여러 structure/bandwidth 동시 실행
python run_accelsim.py \
    --benchmark polybench:polybench-gemm \
    --structure B100_Local B100_Global H100 \
    --bandwidth B200+HBM3e Shoreline_1x \
    --routing baseline near_min_adaptive \
    --near-min-p 0.0

# experiments.csv 전체 sweep (dry-run으로 먼저 확인)
python run_accelsim.py \
    --benchmark GPU_Microbenchmark:mem_bw \
    --all-configs \
    --routing baseline \
    --dry-run
```

**job name 자동 생성 규칙** (AccelSim `-N` 인자):

```
{struct_abbrev}{bw_abbrev}{route_abbrev}_{bench_abbrev}
예) B100_Global + B200+HBM3e + near_min_p0.0 + bfs-rodinia-3.1  →  BGH3nm00_bfs
```

로그 출력 경로:

```
../../sim_run_12.8/<app>/<sanitized_args>/<config>-SASS/<app>*.o<N>
```

---

## Step 5: 로그 파싱 (`parse_accelsim_results.py`)

시뮬레이션 로그를 읽어 CSV 파일로 집계합니다.

```bash
# 단일 benchmark, 자동 prefix 생성 (bfs_BG_H3_bas+mina+nm00_*)
python parse_accelsim_results.py \
    --benchmark rodinia-3.1:bfs-rodinia-3.1 \
    --structure B100_Global \
    --bandwidth B200+HBM3e \
    --routing baseline min_adaptive near_min_adaptive \
    --near-min-p 0.0 1.0

# 명시적 output prefix 지정
python parse_accelsim_results.py \
    --benchmark polybench:polybench-gemm \
    --structure B100_Local B100_Global H100 \
    --bandwidth B200+HBM3e Shoreline_1x \
    --routing baseline near_min_adaptive \
    --near-min-p 0.0 \
    --output gemm_results

# 모든 config
python parse_accelsim_results.py \
    --benchmark GPU_Microbenchmark:mem_bw \
    --all-configs \
    --routing baseline min_adaptive
```

### 출력 파일


| 파일                         | 내용                                                                                |
| ------------------------------ | ------------------------------------------------------------------------------------- |
| `{prefix}_summary.csv`       | benchmark·struct·bw·routing별 집계 (cycles, L2 BW, hit rate, 링크 saturation 등) |
| `{prefix}_directions.csv`    | src→dst 방향별 링크 사용량·포화도                                                 |
| `{prefix}_nm_directions.csv` | Near-Min Adaptive 전용: non-min 방향 선택 breakdown                                 |

---

## Step 6a: 링크 통계 시각화 (`plot_accelsim_link_stats.py`)

src→dst 방향별 링크 사용량/포화도를 구조별로 3D bar chart로 그립니다.

```bash
# 단일 benchmark, baseline + near_min_adaptive 비교
python plot_accelsim_link_stats.py \
    bfs_BG_H3_bas+mina_directions.csv \
    --structure B100_Global \
    --bandwidth B200+HBM3e \
    --routing baseline near_min_adaptive \
    --metric both \
    -o bfs_link.png

# 여러 benchmark CSV 합치기
python plot_accelsim_link_stats.py \
    bfs_BG_H3_bas_directions.csv \
    gemm_BG_H3_bas_directions.csv \
    --routing baseline \
    --metric util \
    -o multi_link_util.png

# Near-Min direction breakdown (_nm_directions.csv, --near 자동 감지)
python plot_accelsim_link_stats.py \
    bfs_RU_H4_nm00_nm_directions.csv \
    --metric sat \
    -o bfs_nearmin_sat.png

# 모든 config 표시
python plot_accelsim_link_stats.py \
    results_directions.csv \
    --all-configs
```

**레이아웃**: 구조별 figure, 행 = bandwidth × benchmark, 열 = routing 조합
**matching summary 자동 탐색**: `_directions.csv` → `_summary.csv` (같은 prefix)
**--near 자동 감지**: 파일명에 `_nm_directions`가 포함되면 자동 활성화

---

## Step 6b: 성능·L2 BW 비교 (`plot_accelsim_perf.py`)

여러 benchmark의 `_summary.csv`를 읽어 정규화된 성능(speedup)과 L2 BW를 bar chart로 비교합니다.
항상 두 파일을 자동 저장합니다: `{prefix}_perf.png`, `{prefix}_bw.png`

```bash
# 기본 비교 (baseline 기준 정규화)
python plot_accelsim_perf.py \
    bfs_BG_H3_bas+mina.csv gemm_BG_H3_bas+mina.csv \
    --structure B100_Global H100 \
    --bandwidth B200+HBM3e \
    --routing baseline min_adaptive near_min_adaptive \
    --near-min-p 0.0 1.0 \
    -o comparison

# 단일 benchmark, 전체 sweep
python plot_accelsim_perf.py \
    mem_bw_summary.csv \
    --structure B100_Local B100_Global B100_Core_Rotate \
    --bandwidth B200+HBM3e Shoreline_1x Shoreline_2x \
    --routing baseline near_min_adaptive \
    --near-min-p 0.0 \
    -o mem_bw_sweep

# --structure/--bandwidth 생략 시 CSV에서 자동 탐지
python plot_accelsim_perf.py bfs_summary.csv \
    --routing baseline near_min_adaptive
```

**정규화 기준** (worst expected = ref = 1.0):

- struct: experiments.csv 순서 첫 번째
- bandwidth: experiments.csv 순서 첫 번째
- routing: baseline (없으면 첫 번째 선택)

**클럭 주파수** (cycles → latency 변환):


| Structure   | 클럭     |
| ------------- | ---------- |
| H100        | 1980 MHz |
| B100_*      | 1965 MHz |
| Rubin_Ultra | 2380 MHz |

---

## 전체 실행 예시 (bfs benchmark)

```bash
cd gpgpu-sim_distribution

# 1. Hit rate 측정 (기존 로그에서)
python avg_miss_rate.py ../../sim_run_12.8/bfs-rodinia-3.1/__data_graph1MW_6_txt/B200-SASS/bfs-rodinia-3.1*.o1
# → miss rate 0.23 → hit rate 0.77 → hitrate.csv에 기록

# 2. Config 생성
python gen_accelsim_configs.py \
    --structure B100_Global \
    --bandwidth B200+HBM3e \
    --routing baseline near_min_adaptive \
    --near-min-p 0.0 1.0

# 3. 시뮬레이션 실행
python run_accelsim.py \
    --benchmark rodinia-3.1:bfs-rodinia-3.1 \
    --structure B100_Global \
    --bandwidth B200+HBM3e \
    --routing baseline near_min_adaptive \
    --near-min-p 0.0 1.0

# 4. 로그 파싱
python parse_accelsim_results.py \
    --benchmark rodinia-3.1:bfs-rodinia-3.1 \
    --structure B100_Global \
    --bandwidth B200+HBM3e \
    --routing baseline near_min_adaptive \
    --near-min-p 0.0 1.0
# → bfs_BG_H3_bas+nm00+nm10_summary.csv 등 생성

# 5a. 링크 통계 플롯
python plot_accelsim_link_stats.py \
    bfs_BG_H3_bas+nm00+nm10_directions.csv \
    --metric both \
    -o bfs_link.png

# 5b. 성능 플롯
python plot_accelsim_perf.py \
    bfs_BG_H3_bas+nm00+nm10_summary.csv \
    --routing baseline near_min_adaptive \
    --near-min-p 0.0 1.0 \
    -o bfs_perf
# → bfs_perf_perf.png, bfs_perf_bw.png 저장
```

---

## 새 실험 추가 방법

### 새 Structure 추가

`experiments.csv`의 Structure 섹션에 행 추가:

```
Structure       Core  HBM  SM_per_Xbar
...
MyNewChip        3    12   80           ← 추가
```

### 새 Bandwidth 설정 추가

`experiments.csv`의 Bandwidth 섹션에 행 추가:

```
Bandwidth         GPU-to-GPU  GPU-to-HBM  HBM-to-HBM  TSV   L2_per_HBM
...
MyBW_config       12          2.5         2.5          2.5   64          ← 추가
```

추가 후 `gen_accelsim_configs.py`를 다시 실행하면 새 조합의 config가 자동 생성됩니다.

### 새 Benchmark 추가

1. `define-all-apps.yml`에 benchmark가 등록되어 있는지 확인
2. 기존 baseline config로 **한 번 시뮬레이션** 실행
3. `avg_miss_rate.py`로 hit rate 측정
4. `hitrate.csv`에 `suite:app` 형식으로 등록
5. `run_accelsim.py`로 전체 실험 실행

---

## 주요 파라미터 약어 표


| Structure        | 약어 | Bandwidth        | 약어 |
| ------------------ | ------ | ------------------ | ------ |
| B100_Local       | BL   | B200+HBM3e       | H3   |
| H100             | H1   | Rubin_Ultra+HBM4 | H4   |
| B100_Global      | BG   | Shoreline_1x     | S1   |
| B100_Core_Rotate | BCR  | Shoreline_1.5x   | S15  |
| Rubin_Ultra      | RU   | Shoreline_2x     | S2   |


| Routing                   | 약어 |
| --------------------------- | ------ |
| baseline                  | bas  |
| min_oblivious             | mino |
| min_adaptive              | mina |
| near_min_adaptive (p=0.0) | nm00 |
| near_min_adaptive (p=1.0) | nm10 |
| ugal                      | ug   |
| valiant                   | val  |
