# https://colab.research.google.com/drive/1XACU_t9t-dZNzeGIJYotdxHuaojVcSLN


!apt-get install -y bowtie2

!wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
!unzip GRCh38_noalt_as.zip

from google.colab import files
import pandas as pd

uploaded = files.upload()




filename = list(uploaded.keys())[0]
bar = pd.read_csv(filename, sep='\t')

print(f"Shape: {bar.shape}")
bar['WT_guide'].head()





import subprocess
import tempfile
import os
from collections import defaultdict

class Bowtie2OffTargetAnalyzer:
    def __init__(self, index_path):
        self.index_path = index_path

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

        # 넉넉하게 검색 (나중에 mismatch로 필터링)
        cmd = [
            'bowtie2', '-x', self.index_path,
            '-f', fasta_file, '-k', '10000',
            '--very-sensitive',
            '--score-min', 'L,-24,0',  # 최대 4 mismatch까지 일단 검색
            '--threads', str(threads),
            '-S', sam_file, '--no-hd'
        ]

        subprocess.run(cmd, check=True, capture_output=True, text=True)

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

                # NM 태그에서 실제 mismatch 수 확인
                nm = 0
                for field in fields[11:]:
                    if field.startswith('NM:i:'):
                        nm = int(field.split(':')[2])
                        break

                # mismatch 필터링 ⭐
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





INDEX_PATH = "/content/GRCh38_noalt_as/GRCh38_noalt_as"

analyzer = Bowtie2OffTargetAnalyzer(INDEX_PATH)
guides = bar['WT_guide'].tolist()

print(f"📊 총 {len(guides)}개 guide 분석 시작")
results = analyzer.analyze_all(guides, ['NAG', 'NGA', 'NGG'], mismatch=2)  # 0, 1, 3 각각 테스트

# 결과 추가
for pam in ['NAG', 'NGA', 'NGG']:
    bar[f'offtarget_{pam}'] = bar['WT_guide'].apply(lambda x: results[pam].get(x.upper(), 0))

bar['offtarget_total'] = bar['offtarget_NAG'] + bar['offtarget_NGA'] + bar['offtarget_NGG']

print("✅ 완료!")
bar[['WT_guide', 'offtarget_NAG', 'offtarget_NGA', 'offtarget_NGG', 'offtarget_total']].head(10)






bar.to_csv('bar_with_offtargets.csv', index=False)
files.download('bar_with_offtargets.csv')



