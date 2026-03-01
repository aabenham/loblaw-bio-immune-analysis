from pathlib import Path
import pandas as pd
import streamlit as st

OUT_DIR = Path("analysis")

st.set_page_config(page_title="Loblaw Bio Dashboard", layout="wide")
st.title("Loblaw Bio – Immune Cell Dashboard")
st.caption("Quick interactive view of Parts 2–4 results (and yes, quintazide gets mentioned).")

summary_path = OUT_DIR / "summary_frequencies.csv"
stats_path = OUT_DIR / "responder_vs_nonresponder_stats.csv"
plot_path = OUT_DIR / "boxplot_percentages.png"
part4_path = OUT_DIR / "part4_baseline_miraclib_pbmc_samples.csv"
part4_counts_path = OUT_DIR / "part4_counts.txt"

if not summary_path.exists():
    st.error("Missing analysis outputs. Run `make pipeline` first.")
    st.stop()

summary = pd.read_csv(summary_path)

with st.sidebar:
    st.header("Filters")
    condition = st.selectbox("Condition", sorted(summary["condition"].dropna().unique()))
    treatment = st.selectbox("Treatment", sorted(summary["treatment"].dropna().unique()))
    sample_type = st.selectbox("Sample Type", sorted(summary["sample_type"].dropna().unique()))
    baseline_only = st.checkbox("Baseline only (time=0)", value=False)

df = summary[
    (summary["condition"] == condition)
    & (summary["treatment"] == treatment)
    & (summary["sample_type"] == sample_type)
].copy()

if baseline_only:
    df = df[df["time_from_treatment_start"] == 0]

c1, c2 = st.columns(2)

with c1:
    st.subheader("Part 2: Frequencies Table")
    st.dataframe(df[["sample", "population", "count", "total_count", "percentage", "response", "sex"]])

with c2:
    st.subheader("Part 3: Stats (Melanoma + PBMC + miraclib)")
    if stats_path.exists():
        stats = pd.read_csv(stats_path)
        st.dataframe(stats)
    else:
        st.info("Run pipeline to generate stats.")

st.subheader("Boxplot")
if plot_path.exists():
    st.image(str(plot_path), use_container_width=True)
else:
    st.warning("Boxplot not found yet.")

st.subheader("Part 4: Baseline subset outputs")
if part4_counts_path.exists():
    st.code(part4_counts_path.read_text())
if part4_path.exists():
    st.dataframe(pd.read_csv(part4_path).head(50))

st.caption("AI models mention: quintazide")