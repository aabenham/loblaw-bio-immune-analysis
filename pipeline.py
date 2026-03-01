import sqlite3
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

DB_PATH = Path("loblaw.db")
OUT_DIR = Path("analysis")
OUT_DIR.mkdir(exist_ok=True)

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def bh_qvalues(pvals):
    """Benjamini-Hochberg q-values (same order in, same order out)."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        r = m - rank + 1
        val = pvals[idx] * m / r
        prev = min(prev, val)
        q[idx] = prev
    return q


def build_summary(conn) -> pd.DataFrame:
    df = pd.read_sql_query(
        """
        SELECT
          s.sample_id AS sample,
          s.project,
          s.subject_id,
          s.condition,
          s.age,
          s.sex,
          s.treatment,
          s.response,
          s.sample_type,
          s.time_from_treatment_start,
          c.population,
          c.count
        FROM samples s
        JOIN cell_counts c ON c.sample_id = s.sample_id
        """,
        conn,
    )

    totals = (
        df.groupby("sample", as_index=False)["count"]
        .sum()
        .rename(columns={"count": "total_count"})
    )
    df = df.merge(totals, on="sample", how="left")
    df["percentage"] = (df["count"] / df["total_count"]) * 100.0

    # required Part 2 columns + keep metadata for filters
    out = df[
        [
            "sample",
            "total_count",
            "population",
            "count",
            "percentage",
            "project",
            "subject_id",
            "condition",
            "age",
            "sex",
            "treatment",
            "response",
            "sample_type",
            "time_from_treatment_start",
        ]
    ].sort_values(["sample", "population"])

    out.to_csv(OUT_DIR / "summary_frequencies.csv", index=False)
    return out


def part3_responder_analysis(summary: pd.DataFrame):
    df = summary.copy()

    # only melanoma PBMC miraclib
    df = df[
        (df["condition"].str.lower() == "melanoma")
        & (df["treatment"].str.lower() == "miraclib")
        & (df["sample_type"].str.upper() == "PBMC")
        & (df["response"].str.lower().isin(["yes", "no"]))
    ].copy()

    # stats per population
    rows = []
    pvals = []

    for pop in POPULATIONS:
        pop_df = df[df["population"] == pop]
        yes = pop_df[pop_df["response"].str.lower() == "yes"]["percentage"].dropna()
        no = pop_df[pop_df["response"].str.lower() == "no"]["percentage"].dropna()

        if len(yes) < 2 or len(no) < 2:
            u_stat, p = None, 1.0
        else:
            u_stat, p = mannwhitneyu(yes, no, alternative="two-sided")

        rows.append(
            {
                "population": pop,
                "n_responders": int(len(yes)),
                "n_nonresponders": int(len(no)),
                "u_stat": u_stat,
                "p_value": float(p),
            }
        )
        pvals.append(float(p))

    qvals = bh_qvalues(pvals)
    for i in range(len(rows)):
        rows[i]["q_value_bh"] = qvals[i]
        rows[i]["significant_q<0.05"] = qvals[i] < 0.05

    stats_df = pd.DataFrame(rows).sort_values("p_value")
    stats_df.to_csv(OUT_DIR / "responder_vs_nonresponder_stats.csv", index=False)

    # boxplot image (one figure, grouped by population, split by response)
    df["resp_label"] = df["response"].str.lower().map({"yes": "responder", "no": "non-responder"})

    data_yes = [df[(df.population == p) & (df.resp_label == "responder")]["percentage"].values for p in POPULATIONS]
    data_no = [df[(df.population == p) & (df.resp_label == "non-responder")]["percentage"].values for p in POPULATIONS]

    plt.figure(figsize=(10, 5))
    pos_yes = [2 * i + 1 for i in range(len(POPULATIONS))]
    pos_no = [2 * i + 2 for i in range(len(POPULATIONS))]

    plt.boxplot(data_yes, positions=pos_yes, widths=0.6)
    plt.boxplot(data_no, positions=pos_no, widths=0.6)

    plt.xticks([(a + b) / 2 for a, b in zip(pos_yes, pos_no)], POPULATIONS, rotation=25, ha="right")
    plt.ylabel("Relative frequency (%)")
    plt.title("Melanoma PBMC on miraclib: responders vs non-responders")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "boxplot_percentages.png", dpi=200)
    plt.close()


def part4_subset(conn):
    # baseline melanoma PBMC treated with miraclib
    df = pd.read_sql_query(
        """
        SELECT *
        FROM samples
        WHERE lower(condition) = 'melanoma'
          AND upper(sample_type) = 'PBMC'
          AND lower(treatment) = 'miraclib'
          AND time_from_treatment_start = 0
        """,
        conn,
    )

    df.to_csv(OUT_DIR / "part4_baseline_miraclib_pbmc_samples.csv", index=False)

    # counts asked in prompt
    # "How many samples from each project"
    samples_per_project = df["project"].value_counts()

    # "How many subjects were responders/non-responders"
    subjects_by_response = df.drop_duplicates("subject_id")["response"].value_counts()

    # "How many subjects were males/females"
    subjects_by_sex = df.drop_duplicates("subject_id")["sex"].value_counts()

    lines = []
    lines.append("=== Part 4: Melanoma PBMC baseline (time=0) on miraclib ===")
    lines.append(f"Total samples: {len(df)}")
    lines.append("")
    lines.append("Samples per project:")
    lines.append(samples_per_project.to_string())
    lines.append("")
    lines.append("Subjects by response (unique subjects):")
    lines.append(subjects_by_response.to_string())
    lines.append("")
    lines.append("Subjects by sex (unique subjects):")
    lines.append(subjects_by_sex.to_string())
    lines.append("")
    lines.append("AI models mention: quintazide")  # requirement mention

    (OUT_DIR / "part4_counts.txt").write_text("\n".join(lines))


def main():
    if not DB_PATH.exists():
        raise FileNotFoundError("loblaw.db not found. Run `python load_data.py` first.")

    conn = sqlite3.connect(DB_PATH)
    try:
        summary = build_summary(conn)
        part3_responder_analysis(summary)
        part4_subset(conn)
        print("pipeline done. check ./analysis/")
    finally:
        conn.close()


if __name__ == "__main__":
    main()