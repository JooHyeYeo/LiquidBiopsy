#%%
import pandas as pd
import numpy as np
from scipy.stats import beta as beta_dist
from pathlib import Path


#%%
# -----------------------------
# MRD용 기본 파라미터
# -----------------------------
# === 기준값만 변경 ===
THETA_T = 0.15   # Tumor_VAF >= 15%
THETA_N = 0.05   # Normal_VAF <= 5%

POSTERIOR_CUTOFF = 0.95  # MRD면 유지, 덜 엄격히면 0.95


ALPHA_PRIOR = 1.0
BETA_PRIOR  = 1.0

MIN_DP_T = 30
MIN_ALT_T = 4
MIN_DP_N = 30

# -----------------------------
# 안전한 파일 로더: 인코딩/구분자 자동 시도
# -----------------------------
def read_table_robust(path: Path) -> pd.DataFrame:
    encodings_to_try = ["utf-8-sig", "utf-16", "cp949", "euc-kr", "latin1"]
    seps_to_try = ["\t", ",", r"\s+"]  # 탭, 콤마, 공백(다중)

    last_err = None
    for enc in encodings_to_try:
        for sep in seps_to_try:
            try:
                if sep == r"\s+":
                    df = pd.read_csv(path, sep=sep, engine="python", encoding=enc)
                else:
                    df = pd.read_csv(path, sep=sep, encoding=enc)
                # 컬럼이 너무 적으면(구분자 실패) 다음 시도
                if df.shape[1] < 5:
                    continue
                return df
            except Exception as e:
                last_err = e
                continue
    raise RuntimeError(f"Failed to read {path.name}. Last error: {last_err}")

# -----------------------------
# posterior 계산 함수
# -----------------------------
def posterior_prob_ge(alt, ref, theta, alpha0=1.0, beta0=1.0):
    alt = np.asarray(alt, dtype=float)
    ref = np.asarray(ref, dtype=float)
    dp = alt + ref
    a_post = alpha0 + alt
    b_post = beta0 + ref
    out = np.full_like(dp, np.nan, dtype=float)
    valid = (dp > 0) & np.isfinite(a_post) & np.isfinite(b_post)
    out[valid] = 1.0 - beta_dist.cdf(theta, a_post[valid], b_post[valid])
    return out

def posterior_prob_le(alt, ref, theta, alpha0=1.0, beta0=1.0):
    alt = np.asarray(alt, dtype=float)
    ref = np.asarray(ref, dtype=float)
    dp = alt + ref
    a_post = alpha0 + alt
    b_post = beta0 + ref
    out = np.full_like(dp, np.nan, dtype=float)
    valid = (dp > 0) & np.isfinite(a_post) & np.isfinite(b_post)
    out[valid] = beta_dist.cdf(theta, a_post[valid], b_post[valid])
    return out

# -----------------------------
# 한 파일 처리
# -----------------------------
def process_one_file(
    in_path: Path,
    out_dir: Path,
    normal_ref_col="Normal_Ref_DP",
    normal_alt_col="Normal_Alt_DP",
    tumor_ref_col="Tumor_Ref_DP",
    tumor_alt_col="Tumor_Alt_DP",
):
    df = read_table_robust(in_path)

    # 혹시 컬럼명 앞뒤 공백/보이지 않는 문자(BOM 등) 제거
    df.columns = [str(c).strip().replace("\ufeff", "") for c in df.columns]

    # 필요한 컬럼 존재 확인
    needed = [normal_ref_col, normal_alt_col, tumor_ref_col, tumor_alt_col]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}. Available: {list(df.columns)[:30]}")

    # 숫자 변환
    for c in needed:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df["Normal_DP"] = df[normal_ref_col] + df[normal_alt_col]
    df["Tumor_DP"]  = df[tumor_ref_col]  + df[tumor_alt_col]

    df["P_tumor_ge_theta"] = posterior_prob_ge(
        alt=df[tumor_alt_col].to_numpy(),
        ref=df[tumor_ref_col].to_numpy(),
        theta=THETA_T,
        alpha0=ALPHA_PRIOR,
        beta0=BETA_PRIOR,
    )
    df["P_normal_le_theta"] = posterior_prob_le(
        alt=df[normal_alt_col].to_numpy(),
        ref=df[normal_ref_col].to_numpy(),
        theta=THETA_N,
        alpha0=ALPHA_PRIOR,
        beta0=BETA_PRIOR,
    )
    df["SomaticPosterior"] = df["P_tumor_ge_theta"] * df["P_normal_le_theta"]

    df["MRD_hardfilter_pass"] = (
        (df["Tumor_DP"] >= MIN_DP_T) &
        (df[tumor_alt_col] >= MIN_ALT_T) &
        (df["Normal_DP"] >= MIN_DP_N)
    )
    df["MRD_pass"] = df["MRD_hardfilter_pass"] & (df["SomaticPosterior"] >= POSTERIOR_CUTOFF)

    sample = in_path.stem  # EC0156_WGS
    out_all = out_dir / f"{sample}_with_posterior.tsv"
    out_mrd = out_dir / f"{sample}_MRD_candidates.tsv"

    df.to_csv(out_all, sep="\t", index=False)
    df[df["MRD_pass"]].to_csv(out_mrd, sep="\t", index=False)

    print(f"[DONE] {sample}: total={len(df)}  MRD_pass={int(df['MRD_pass'].sum())}")

#%%
# -----------------------------
# 폴더 내 일괄 처리
# -----------------------------
if __name__ == "__main__":
    base_dir = Path("LIB_DESIGN/20250910_EC_Samples/WGS/Raw")  # txt 있는 폴더에서 실행
    out_dir = base_dir / "posterior_results"
    out_dir.mkdir(exist_ok=True)

    for fp in sorted(base_dir.glob("*_WGS*.txt")):
        try:
            process_one_file(fp, out_dir)
        except Exception as e:
            print(f"[ERROR] {fp.name}: {e}")

# %%
