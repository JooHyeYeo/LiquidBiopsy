#%% conda activate aligner
import subprocess
import tempfile
import os
from collections import defaultdict
import pandas as pd
from pyfaidx import Fasta  # pysam 대신 pyfaidx 사용


class Bowtie2OffTargetAnalyzer:
    def __init__(self, index_path, ref_fasta_path):
        self.index_path = index_path
        self.ref_fasta_path = ref_fasta_path
        
        if not os.path.exists(f"{index_path}.1.bt2"):
            raise FileNotFoundError(f"Bowtie2 index not found: {index_path}.1.bt2")
        if not os.path.exists(ref_fasta_path):
            raise FileNotFoundError(f"Reference fasta not found: {ref_fasta_path}")
        
        self.ref_fasta = Fasta(ref_fasta_path)

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

    def _reverse_complement(self, seq):
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))

    def _get_flanking_seq(self, chrom, pos, strand, seq_len=23, flank_left=4, flank_right=8):
        """
        Reference genome에서 flanking 서열 포함하여 추출
        pos: 1-based position (SAM format, forward strand 기준 leftmost)
        반환: 4bp + 23bp + 8bp = 35bp
        """
        try:
            if strand == '+':
                start = pos - 1 - flank_left
                end = pos - 1 + seq_len + flank_right
            else:
                start = pos - 1 - flank_right
                end = pos - 1 + seq_len + flank_left
            
            # 범위 체크
            chrom_len = len(self.ref_fasta[chrom])
            if start < 0:
                start = 0
            if end > chrom_len:
                end = chrom_len
            
            seq = str(self.ref_fasta[chrom][start:end]).upper()
            
            if strand == '-':
                seq = self._reverse_complement(seq)
            
            return seq
        except Exception as e:
            print(f"⚠️ Failed to fetch {chrom}:{pos} - {e}")
            return ''

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

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"❌ Bowtie2 Error:\n{result.stderr}")
            raise RuntimeError(f"Bowtie2 failed with exit code {result.returncode}")

        hit_counts = defaultdict(int)
        hit_details = defaultdict(list)
        
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 10:
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
                        guide_seq = parts[1]
                        hit_counts[guide_seq] += 1
                        
                        flag = int(fields[1])
                        strand = '-' if flag & 16 else '+'
                        chrom = fields[2]
                        pos = int(fields[3])
                        aligned_seq = fields[9]
                        pam_used = parts[2] if len(parts) >= 3 else pam
                        
                        off_target_35bp = self._get_flanking_seq(chrom, pos, strand)
                        
                        hit_details[guide_seq].append({
                            'chrom': chrom,
                            'pos': pos,
                            'strand': strand,
                            'aligned_seq': aligned_seq,
                            'aligned_guide': aligned_seq[:20],
                            'aligned_pam': aligned_seq[20:],
                            'Off_Target': off_target_35bp,
                            'mismatch': nm,
                            'pam_used': pam_used
                        })

        os.unlink(fasta_file)
        os.unlink(sam_file)

        print(f"   ✅ {sum(hit_counts.values())} hits (mismatch <= {mismatch})")
        return dict(hit_counts), dict(hit_details)

    def analyze_all(self, guides, pam_list=['NAG', 'NGA', 'NGG'], mismatch=3):
        all_counts = {}
        all_details = {}
        for pam in pam_list:
            counts, details = self.search_pam(guides, pam, mismatch)
            all_counts[pam] = counts
            all_details[pam] = details
        return all_counts, all_details
    
    def close(self):
        pass  # pyfaidx는 자동으로 닫힘


print("✅ 함수 정의 완료! (pyfaidx 사용)")

# 출력 폴더 생성
os.makedirs('/extdata3/YEO/LIB_DESIGN/20250109_EC_WGS/Target_Candi_Off', exist_ok=True)
os.makedirs('/extdata3/YEO/LIB_DESIGN/20250109_EC_WGS/Target_Candi_Off/details', exist_ok=True)

INDEX_PATH = "/extdata3/YEO/Bowtie/GRCh38_noalt_as/GRCh38_noalt_as"
REF_FASTA_PATH = "/extdata3/YEO/Bowtie/GRCh38_noalt_as/hg38.fa"

path = '/extdata3/YEO/LIB_DESIGN/20250109_EC_WGS/Target_Candi'
lst = os.listdir(path)

OFFTARGET_THRESHOLD = 100

for i, file in enumerate(lst):
    print(f"\n{'='*50}")
    print(f"📁 [{i+1}/{len(lst)}] {file}")
    
    bar = pd.read_csv(f'{path}/{file}', sep='\t')

    analyzer = Bowtie2OffTargetAnalyzer(INDEX_PATH, REF_FASTA_PATH)
    guides = bar['WT_guide'].tolist()

    print(f"📊 총 {len(guides)}개 guide 분석 시작")
    results_counts, results_details = analyzer.analyze_all(guides, ['NAG', 'NGA', 'NGG'], mismatch=1)

    for pam in ['NAG', 'NGA', 'NGG']:
        bar[f'offtarget_{pam}'] = bar['WT_guide'].apply(lambda x: results_counts[pam].get(x.upper(), 0))

    bar['offtarget_total'] = bar['offtarget_NAG'] + bar['offtarget_NGA'] + bar['offtarget_NGG']
    
    detail_rows = []
    for idx, row in bar.iterrows():
        guide = row['WT_guide'].upper()
        total_count = row['offtarget_total']
        
        if total_count <= OFFTARGET_THRESHOLD and total_count > 0:
            for pam in ['NAG', 'NGA', 'NGG']:
                if guide in results_details[pam]:
                    for hit in results_details[pam][guide]:
                        detail_rows.append({
                            'WT_guide': guide,
                            'offtarget_total': total_count,
                            'PAM_type': pam,
                            'chrom': hit['chrom'],
                            'position': hit['pos'],
                            'strand': hit['strand'],
                            'aligned_seq': hit['aligned_seq'],
                            'aligned_guide': hit['aligned_guide'],
                            'aligned_pam': hit['aligned_pam'],
                            'mismatch': hit['mismatch'],
                            'pam_used': hit['pam_used'],
                            'Off_Target': hit['Off_Target']
                        })
    
    name = file.split('.txt')[0]
    
    bar.to_csv(f'/extdata3/YEO/LIB_DESIGN/20250109_EC_WGS/Target_Candi_Off/{name}_1MM.csv', index=False)
    print(f"✅ 저장 완료: {name}_1MM.csv")
    
    if detail_rows:
        detail_df = pd.DataFrame(detail_rows)
        detail_df.to_csv(f'/extdata3/YEO/LIB_DESIGN/20250109_EC_WGS/Target_Candi_Off/details/{name}_1MM_details.csv', index=False)
        
        n_guides_with_details = detail_df['WT_guide'].nunique()
        print(f"📍 Off-target 상세 정보 저장: {name}_1MM_details.csv")
        print(f"   (off-target ≤ {OFFTARGET_THRESHOLD}인 가이드 {n_guides_with_details}개, 총 {len(detail_rows)}개 위치)")
    else:
        print(f"⚠️ Off-target ≤ {OFFTARGET_THRESHOLD}인 가이드 없음")
    
    analyzer.close()
