import re
import sys

def calculate_access_weighted_miss_rate(log_file_path):
    # 1. 정규표현식 수정: cycle 대신 L2_total_cache_accesses를 찾습니다.
    access_pattern = re.compile(r'L2_total_cache_accesses\s*=\s*(\d+)')
    miss_rate_pattern = re.compile(r'L2_total_cache_miss_rate\s*=\s*([0-9.]+)')

    accesses = []
    miss_rates = []

    try:
        with open(log_file_path, 'r') as file:
            for line in file:
                # L2_total_cache_accesses 찾기
                match_access = access_pattern.search(line)
                if match_access:
                    accesses.append(int(match_access.group(1)))
                
                # L2_total_cache_miss_rate 찾기
                match_miss = miss_rate_pattern.search(line)
                if match_miss:
                    miss_rates.append(float(match_miss.group(1)))

    except FileNotFoundError:
        print(f"오류: '{log_file_path}' 파일을 찾을 수 없습니다.")
        return

    # 데이터 개수 확인
    if len(accesses) == 0 or len(miss_rates) == 0:
        print("로그 파일에서 데이터를 찾지 못했습니다. 로그에 L2_total_cache_accesses가 있는지 확인해주세요.")
        return
    
    if len(accesses) != len(miss_rates):
        print(f"경고: 추출된 Access 개수({len(accesses)})와 Miss Rate 개수({len(miss_rates)})가 다릅니다!")
    
    min_len = min(len(accesses), len(miss_rates))
    
    total_calculated_misses = 0.0
    total_accesses = 0

    print("-" * 55)
    print(f"{'Kernel':<10} | {'L2 Accesses':<15} | {'L2 Miss Rate':<15}")
    print("-" * 55)

    # 2. 로직 수정: Access 기반 가중 평균 계산
    for i in range(min_len):
        a = accesses[i]
        m = miss_rates[i]
        
        print(f"Kernel {i+1:<8} | {a:<15} | {m:.4f}")
        
        # (접근 횟수 * 미스율) = 해당 커널에서 발생한 실제 미스 횟수 추정치
        total_calculated_misses += (a * m)
        total_accesses += a

    print("-" * 55)

    if total_accesses > 0:
        # 최종 미스율 = (전체 미스 횟수 총합) / (전체 접근 횟수 총합)
        final_miss_rate = total_calculated_misses / total_accesses
        print(f"\n[결과] 전체 누적 L2 Accesses: {total_accesses}")
        print(f"[결과] 실제 통합 L2 Miss Rate (Access 가중치): {final_miss_rate:.4f} ({final_miss_rate * 100:.2f}%)")
    else:
        print("계산할 Access 데이터가 없습니다.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        log_path = sys.argv[1]
    else:
        log_path = "gpgpu-sim.log"
        
    calculate_access_weighted_miss_rate(log_path)