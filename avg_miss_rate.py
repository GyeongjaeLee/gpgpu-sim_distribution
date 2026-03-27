import re
import sys

def calculate_weighted_miss_rate(log_file_path):
    # 정규표현식 패턴 설정
    # 주의: 로그 포맷에 따라 공백이 다를 수 있으므로 \s* 로 유연하게 매칭합니다.
    cycle_pattern = re.compile(r'gpu_sim_cycle\s*=\s*(\d+)')
    miss_rate_pattern = re.compile(r'L2_total_cache_miss_rate\s*=\s*([0-9.]+)')

    cycles = []
    miss_rates = []

    try:
        with open(log_file_path, 'r') as file:
            for line in file:
                # gpu_sim_cycle 찾기
                match_cycle = cycle_pattern.search(line)
                if match_cycle:
                    cycles.append(int(match_cycle.group(1)))
                
                # L2_total_cache_miss_rate 찾기
                match_miss = miss_rate_pattern.search(line)
                if match_miss:
                    miss_rates.append(float(match_miss.group(1)))

    except FileNotFoundError:
        print(f"오류: '{log_file_path}' 파일을 찾을 수 없습니다.")
        return

    # 추출된 데이터 개수 확인 (커널 개수와 일치해야 함)
    if len(cycles) == 0 or len(miss_rates) == 0:
        print("로그 파일에서 데이터를 찾지 못했습니다. 패턴을 확인해주세요.")
        return
    
    if len(cycles) != len(miss_rates):
        print(f"경고: 추출된 Cycle 개수({len(cycles)})와 Miss Rate 개수({len(miss_rates)})가 다릅니다!")
        print("로그가 중간에 끊겼거나 포맷이 일관되지 않을 수 있습니다. 짝이 맞는 곳까지만 계산합니다.")
    
    # 두 리스트 중 더 짧은 길이에 맞춰 계산 (안전장치)
    min_len = min(len(cycles), len(miss_rates))
    
    total_weighted_miss_rate = 0.0
    total_cycles = 0

    print("-" * 50)
    print(f"{'Kernel':<10} | {'Cycles':<15} | {'L2 Miss Rate':<15}")
    print("-" * 50)

    # 가중 평균 계산
    for i in range(min_len):
        c = cycles[i]
        m = miss_rates[i]
        
        print(f"Kernel {i+1:<3} | {c:<15} | {m:.4f}")
        
        total_weighted_miss_rate += (c * m)
        total_cycles += c

    print("-" * 50)

    if total_cycles > 0:
        final_miss_rate = total_weighted_miss_rate / total_cycles
        print(f"\n[결과] 전체 누적 사이클: {total_cycles}")
        print(f"[결과] Cycle 가중 평균 L2 Miss Rate: {final_miss_rate:.4f} ({final_miss_rate * 100:.2f}%)")
    else:
        print("계산할 사이클이 없습니다.")

if __name__ == "__main__":
    # 터미널에서 실행할 때 로그 파일 경로를 인자로 받을 수 있도록 처리
    if len(sys.argv) > 1:
        log_path = sys.argv[1]
    else:
        log_path = "gpgpu-sim.log" # 기본 파일명 (필요시 수정하세요)
        
    calculate_weighted_miss_rate(log_path)