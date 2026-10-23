"""htn_utils.py - shared helpers and code sets for the htn-bp-followup-mimic analyses.

Every notebook starts with:

    import sys; sys.path.insert(0, f"{REPO_DIR}/src")
    from htn_utils import *
    con = connect(MIMIC_HOSP, MIMIC_ED)
"""

import json
from pathlib import Path


# --- DuckDB connection -------------------------------------------------------

def connect(hosp_dir, ed_dir=None, threads=4):
    """DuckDB connection with every MIMIC-IV CSV registered as a lazy view.

    hosp tables keep their names (admissions, omr, ...);
    ED tables get an  ed_  prefix (ed_edstays, ...).
    """
    import duckdb
    con = duckdb.connect()
    con.execute(f"PRAGMA threads={threads}")
    for f in sorted(Path(hosp_dir).glob("*.csv.gz")):
        name = f.name.replace(".csv.gz", "")
        con.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_csv_auto({str(f)!r})")
    if ed_dir is not None:
        for f in sorted(Path(ed_dir).glob("*.csv.gz")):
            name = "ed_" + f.name.replace(".csv.gz", "")
            con.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_csv_auto({str(f)!r})")
    views = [r[0] for r in con.execute("SHOW TABLES").fetchall()]
    print(f"{len(views)} views registered")
    return con


# --- file IO -----------------------------------------------------------------
# Patient-level tables go to DATA_DIR (outside the repo, never committed).
# Aggregate summaries (counts, rules) go to the repo's results/ folder.

def save_df(df, data_dir, name):
    path = Path(data_dir) / f"{name}.parquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    print(f"saved {path.name}  ({len(df):,} rows, {path.stat().st_size / 1e6:.1f} MB)")
    return path


def load_df(data_dir, name):
    import pandas as pd
    df = pd.read_parquet(Path(data_dir) / f"{name}.parquet")
    print(f"{name}: {len(df):,} rows x {df.shape[1]} cols")
    return df


def save_json(obj, results_dir, name):
    path = Path(results_dir) / f"{name}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str))
    print(f"saved {path.name}")
    return path


# --- checks against the paper ------------------------------------------------

def check(label, value, expected):
    """Print value next to the number reported in the paper."""
    ok = value == expected
    print(f"{'OK      ' if ok else 'MISMATCH'} {label}: {value:,} (paper: {expected:,})")
    return ok


# --- DuckDB SQL helpers ------------------------------------------------------
# DuckDB has no LIKE ANY(array), so prefix and substring matches are expanded
# into chained ORs.

def sql_prefix(col, prefixes):
    return "(" + " OR ".join(f"{col} LIKE {p + '%'!r}" for p in prefixes) + ")"


def sql_contains(col, terms):
    return "(" + " OR ".join(f"LOWER({col}) LIKE {'%' + t + '%'!r}" for t in terms) + ")"


# --- code sets (Table S1) ----------------------------------------------------

HTN_9 = ("401", "402", "403", "404")
HTN_10 = ("I10", "I11", "I12", "I13")

# secondary hypertension with a surgically correctable cause -> exclusion
CORRECTABLE_9 = ("4401", "2270", "1940", "2550", "25510", "25511", "25512",
                 "25513", "25514", "7471")
CORRECTABLE_10 = ("I701", "E275", "C74", "D35", "E260", "E24", "Q251")

MI_9, MI_10 = ("410",), ("I21", "I22")
STROKE_9, STROKE_10 = ("430", "431", "432", "433", "434", "436"), ("I60", "I61", "I62", "I63")
HF_9, HF_10 = ("428",), ("I50",)

CV_9 = MI_9 + STROKE_9 + HF_9
CV_10 = MI_10 + STROKE_10 + HF_10

ANTIHTN = [
    # ACE inhibitors
    "lisinopril", "enalapril", "ramipril", "captopril", "benazepril", "quinapril", "fosinopril",
    # angiotensin receptor blockers
    "losartan", "valsartan", "olmesartan", "irbesartan", "candesartan", "telmisartan",
    # beta-blockers
    "metoprolol", "atenolol", "carvedilol", "labetalol", "propranolol", "bisoprolol", "nebivolol",
    # calcium channel blockers
    "amlodipine", "nifedipine", "diltiazem", "verapamil", "nicardipine", "felodipine",
    # thiazide-type diuretics
    "hydrochlorothiazide", "chlorthalidone", "indapamide", "metolazone",
    # loop diuretics
    "furosemide", "torsemide", "bumetanide",
    # potassium-sparing diuretics
    "spironolactone", "eplerenone", "triamterene", "amiloride",
    # central agents and direct vasodilators
    "clonidine", "methyldopa", "hydralazine", "minoxidil",
    # alpha-blockers
    "doxazosin", "terazosin", "prazosin",
]


# --- study constants ---------------------------------------------------------

HOSP_DATA_END_YEAR = 2020   # latest real year of hospital and ED records
DEATH_EXCL_DAYS = 90        # exclude death within 90 days of index
LANDMARK_DAY = 730          # day-730 landmark
DEATH_CAPTURE_DAYS = 365    # deaths captured up to 1 year after last hospital/ED contact
FOLLOWUP_END_DAY = 1460     # follow-up ends by day 1,460 at the latest
