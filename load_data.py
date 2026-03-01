import sqlite3
from pathlib import Path

import pandas as pd

DB_PATH = Path("loblaw.db")
CSV_PATH = Path("cell-count.csv")

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA_SQL = """
PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;

-- one row per sample
CREATE TABLE samples (
  sample_id TEXT PRIMARY KEY,
  project TEXT,
  subject_id TEXT,
  condition TEXT,
  age INTEGER,
  sex TEXT,
  treatment TEXT,
  response TEXT,            -- allow NULL (trial in progress)
  sample_type TEXT,
  time_from_treatment_start REAL
);

-- long format counts
CREATE TABLE cell_counts (
  sample_id TEXT NOT NULL,
  population TEXT NOT NULL,
  count INTEGER NOT NULL,
  PRIMARY KEY (sample_id, population),
  FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

CREATE INDEX idx_samples_filters
  ON samples(condition, treatment, sample_type, response, sex, time_from_treatment_start);

CREATE INDEX idx_cell_counts_population
  ON cell_counts(population);
"""


def main():
    if not CSV_PATH.exists():
        raise FileNotFoundError("cell-count.csv must be in the repo root")

    df = pd.read_csv(CSV_PATH)

    expected_cols = [
        "project", "subject", "condition", "age", "sex", "treatment", "response",
        "sample", "sample_type", "time_from_treatment_start",
        *POPULATIONS
    ]
    missing = [c for c in expected_cols if c not in df.columns]
    if missing:
        raise ValueError(f"CSV missing columns: {missing}")

    # rename to DB-friendly names
    df = df.rename(columns={
        "sample": "sample_id",
        "subject": "subject_id",
    })

    # build samples table
    samples = df[[
        "sample_id", "project", "subject_id", "condition", "age", "sex",
        "treatment", "response", "sample_type", "time_from_treatment_start"
    ]].copy()

    # build cell_counts long format
    counts_long = df[["sample_id"] + POPULATIONS].melt(
        id_vars=["sample_id"],
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    )
    counts_long["count"] = pd.to_numeric(counts_long["count"], errors="coerce").fillna(0).astype(int)

    conn = sqlite3.connect(DB_PATH)
    try:
        conn.executescript(SCHEMA_SQL)

        # re-run safe: wipe and reload
        conn.execute("DELETE FROM cell_counts;")
        conn.execute("DELETE FROM samples;")

        # treat empty strings as missing (important for sqlite constraints + filters)
        samples["response"] = samples["response"].replace(r"^\s*$", pd.NA, regex=True)
        samples.to_sql("samples", conn, if_exists="append", index=False)
        counts_long.to_sql("cell_counts", conn, if_exists="append", index=False)

        conn.commit()
        print(f"DB created: {DB_PATH}")
        print(f"Loaded {len(samples)} samples")
        print(f"Loaded {len(counts_long)} (sample, population) rows")

    finally:
        conn.close()


if __name__ == "__main__":
    main()