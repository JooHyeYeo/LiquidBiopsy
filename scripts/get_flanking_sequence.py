#%%
import os

# --- 작업 디렉토리 ---
new_directory = "/Users/jhyeo/Dropbox/2_LIQUID_BIOPSY/CODE_FLASH_sgRNA_design"
os.chdir(new_directory)

#%%
def reverse_complement(sSeq):
    dict_sBases = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
        'N': 'N', '.': '.', '*': '*',
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a'
    }
    return ''.join(dict_sBases.get(b, b) for b in sSeq)[::-1]

def load_chr_fasta(chrom):
    """
    hg38/chr{chrom}.fa 파일을 읽어서
    전체 염기서열을 하나의 문자열로 반환
    """
    fasta_path = f"hg38/chr{chrom}.fa"
    seq = []
    with open(fasta_path, "r") as f:
        f.readline()  # >chrX 헤더 제거
        for line in f:
            seq.append(line.strip())
    return ''.join(seq)

def get_flanking_seq(chrom, pos, flank=100):
    """
    chrom : "X"
    pos   : hg38 1-based coordinate
    flank : ±bp
    """
    genome = load_chr_fasta(chrom)

    start = max(0, pos - flank - 1)  # 1-based → 0-based
    end   = pos + flank              # slicing end-exclusive

    return genome[start:end]

#%%
# ---- 입력 ----
chrom = "X"
pos   = 101_497_425
flank = 100

seq = get_flanking_seq(chrom, pos, flank)

print(f">chr{chrom}:{pos-flank}-{pos+flank}")
print(seq)
print("Length:", len(seq))   # 201이어야 정상
