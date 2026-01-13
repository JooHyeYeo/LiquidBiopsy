#%% conda activate aligner

import subprocess
import tempfile
import os
from collections import defaultdict
import pandas as pd


class Bowtie2OffTargetAnalyzer:
    def __init__(self, index_path):
        self.index_path = index_path
        # 인덱스 파일 존재 확인
        if not os.path.exists(f"{index_path}.1.bt2"):
            raise FileNotFoundError(f"Bowtie2 index not found: {index_path}.1.bt2")

    def _expand_n(self, pam):
        if 'N' not in pam:
            return [pam]
        results = ['']
        for base in pam:
            if base == 'N':
                results = [r + b for r in results for b in 'ATGC']
            else:
                results = [r + base for r in results]
        return results

    def search_pam(self, guides, pam, mismatch=3, threads=2):
        expanded_pams = self._expand_n(pam)
        print(f"\n🔍 PAM: {pam} → {expanded_pams}")
        print(f"   Guides: {len(guides)}, Mismatch: {mismatch}")

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            for i, guide in enumerate(guides):
                guide_upper = guide.upper()
                for ep in expanded_pams:
                    query = guide_upper + ep
                    f.write(f">{i}|{guide_upper}|{ep}\n{query}\n")
            fasta_file = f.name

        sam_file = tempfile.mktemp(suffix='.sam')

        cmd = [
            'bowtie2', '-x', self.index_path,
            '-f', fasta_file, '-k', '10000',
            '--very-sensitive',
            '--score-min', 'L,-24,0',
            '--threads', str(threads),
            '-S', sam_file, '--no-hd'
        ]

        # 에러 메시지 출력되도록 수정
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"❌ Bowtie2 Error:\n{result.stderr}")
            raise RuntimeError(f"Bowtie2 failed with exit code {result.returncode}")

        hit_counts = defaultdict(int)
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                if fields[2] == '*' or int(fields[1]) & 4:
                    continue

                nm = 0
                for field in fields[11:]:
                    if field.startswith('NM:i:'):
                        nm = int(field.split(':')[2])
                        break

                if nm <= mismatch:
                    parts = fields[0].split('|')
                    if len(parts) >= 2:
                        hit_counts[parts[1]] += 1

        os.unlink(fasta_file)
        os.unlink(sam_file)

        print(f"   ✅ {sum(hit_counts.values())} hits (mismatch <= {mismatch})")
        return dict(hit_counts)

    def analyze_all(self, guides, pam_list=['NAG', 'NGA', 'NGG'], mismatch=3):
        return {pam: self.search_pam(guides, pam, mismatch) for pam in pam_list}


print("✅ 함수 정의 완료! (NM 태그 필터링 추가)")

# 출력 폴더 생성
os.makedirs('/extdata3/YEO/LIB_DESIGN/20250109_LC_TSO500/Target_Candi_Off', exist_ok=True)

# 인덱스 경로 - 폴더 안에 있으면 이렇게
INDEX_PATH = "/extdata3/YEO/Bowtie/GRCh38_noalt_as/GRCh38_noalt_as"
# 또는 폴더 바로 밑에 있으면
# INDEX_PATH = "/extdata3/YEO/Bowtie/GRCh38_noalt_as"

path = '/extdata3/YEO/LIB_DESIGN/20250109_LC_TSO500/Target_Candi_merge'
lst = os.listdir(path)

for i, file in enumerate(lst):
    print(f"\n{'='*50}")
    print(f"📁 [{i+1}/{len(lst)}] {file}")
    
    bar = pd.read_csv(f'{path}/{file}', sep='\t')
    bar = bar.rename(columns = {'WT_sgRNA' : 'WT_guide'})

    analyzer = Bowtie2OffTargetAnalyzer(INDEX_PATH)
    guides = bar['WT_guide'].tolist()

    print(f"📊 총 {len(guides)}개 guide 분석 시작")
    results = analyzer.analyze_all(guides, ['NAG', 'NGA', 'NGG'], mismatch=0)

    for pam in ['NAG', 'NGA', 'NGG']:
        bar[f'offtarget_{pam}'] = bar['WT_guide'].apply(lambda x: results[pam].get(x.upper(), 0))

    bar['offtarget_total'] = bar['offtarget_NAG'] + bar['offtarget_NGA'] + bar['offtarget_NGG']
    name = file.split('.txt')[0]

    bar.to_csv(f'/extdata3/YEO/LIB_DESIGN/20250109_LC_TSO500/Target_Candi_Off/{name}_0MM.csv', index=False)
    print(f"✅ 저장 완료: {name}_0MM.csv")


for i, file in enumerate(lst):
    print(f"\n{'='*50}")
    print(f"📁 [{i+1}/{len(lst)}] {file}")
    
    bar = pd.read_csv(f'{path}/{file}', sep='\t')
    bar = bar.rename(columns = {'WT_sgRNA' : 'WT_guide'})

    analyzer = Bowtie2OffTargetAnalyzer(INDEX_PATH)
    guides = bar['WT_guide'].tolist()

    print(f"📊 총 {len(guides)}개 guide 분석 시작")
    results = analyzer.analyze_all(guides, ['NAG', 'NGA', 'NGG'], mismatch=1)

    for pam in ['NAG', 'NGA', 'NGG']:
        bar[f'offtarget_{pam}'] = bar['WT_guide'].apply(lambda x: results[pam].get(x.upper(), 0))

    bar['offtarget_total'] = bar['offtarget_NAG'] + bar['offtarget_NGA'] + bar['offtarget_NGG']
    name = file.split('.txt')[0]

    bar.to_csv(f'/extdata3/YEO/LIB_DESIGN/20250109_LC_TSO500/Target_Candi_Off/{name}_1MM.csv', index=False)
    print(f"✅ 저장 완료: {name}_1MM.csv")



for i, file in enumerate(lst):
    print(f"\n{'='*50}")
    print(f"📁 [{i+1}/{len(lst)}] {file}")
    
    bar = pd.read_csv(f'{path}/{file}', sep='\t')
    bar = bar.rename(columns = {'WT_sgRNA' : 'WT_guide'})

    analyzer = Bowtie2OffTargetAnalyzer(INDEX_PATH)
    guides = bar['WT_guide'].tolist()

    print(f"📊 총 {len(guides)}개 guide 분석 시작")
    results = analyzer.analyze_all(guides, ['NAG', 'NGA', 'NGG'], mismatch=2)

    for pam in ['NAG', 'NGA', 'NGG']:
        bar[f'offtarget_{pam}'] = bar['WT_guide'].apply(lambda x: results[pam].get(x.upper(), 0))

    bar['offtarget_total'] = bar['offtarget_NAG'] + bar['offtarget_NGA'] + bar['offtarget_NGG']
    name = file.split('.txt')[0]

    bar.to_csv(f'/extdata3/YEO/LIB_DESIGN/20250109_LC_TSO500/Target_Candi_Off/{name}_2MM.csv', index=False)
    print(f"✅ 저장 완료: {name}_2MM.csv")



for i, file in enumerate(lst):
    print(f"\n{'='*50}")
    print(f"📁 [{i+1}/{len(lst)}] {file}")
    
    bar = pd.read_csv(f'{path}/{file}', sep='\t')
    bar = bar.rename(columns = {'WT_sgRNA' : 'WT_guide'})

    analyzer = Bowtie2OffTargetAnalyzer(INDEX_PATH)
    guides = bar['WT_guide'].tolist()

    print(f"📊 총 {len(guides)}개 guide 분석 시작")
    results = analyzer.analyze_all(guides, ['NAG', 'NGA', 'NGG'], mismatch=3)

    for pam in ['NAG', 'NGA', 'NGG']:
        bar[f'offtarget_{pam}'] = bar['WT_guide'].apply(lambda x: results[pam].get(x.upper(), 0))

    bar['offtarget_total'] = bar['offtarget_NAG'] + bar['offtarget_NGA'] + bar['offtarget_NGG']
    name = file.split('.txt')[0]

    bar.to_csv(f'/extdata3/YEO/LIB_DESIGN/20250109_LC_TSO500/Target_Candi_Off/{name}_3MM.csv', index=False)
    print(f"✅ 저장 완료: {name}_3MM.csv")
