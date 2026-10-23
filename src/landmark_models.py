"""landmark_models.py - shared covariate matrix, splines, and pooling for notebooks 06 and 06b.

Used by notebooks 06 and 06b.
"""

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

ORDER = ["Minimal Contact", "Declining Contact", "Rising Contact",
         "Sustained Moderate Contact", "Sustained High Contact"]

# covariates entered as restricted cubic splines: 4 knots, and 3 knots (only 3 distinct percentiles)
SPL4 = ["age_index", "sbp_baseline", "egfr", "bmi", "log_enc_prior", "n_antihtn_drugs"]
SPL3 = ["elixhauser_vw", "log_inp_prior"]
MISS = ["egfr", "bmi", "hba1c"]


# --- splines -----------------------------------------------------------------

def rcs(x, t):
    """Harrell restricted cubic spline basis: [x, nonlinear terms ...]."""
    x = np.asarray(x, float); t = np.asarray(t, float); k = len(t)
    p = lambda u: np.maximum(u, 0) ** 3
    cols = [x]
    for j in range(k - 2):
        cols.append((p(x - t[j]) - p(x - t[k - 2]) * (t[k - 1] - t[j]) / (t[k - 1] - t[k - 2])
                     + p(x - t[k - 1]) * (t[k - 2] - t[j]) / (t[k - 1] - t[k - 2])) / (t[k - 1] - t[0]) ** 2)
    return np.column_stack(cols)


def make_knots(A, spline_vars, gap_col="gap_days"):
    """Knots at the 5th, 35th, 65th, and 95th percentiles of observed values.

    Variables with fewer than 4 distinct knots get a 3-knot version under the key '<var>_3'.
    """
    KN = {}
    for v in spline_vars + [gap_col]:
        k = np.unique(np.nanpercentile(A[v], [5, 35, 65, 95]).round(4))
        KN[v] = k.tolist() if len(k) == 4 else None
    for v in SPL3:
        KN[v + "_3"] = np.unique(np.nanpercentile(A[v], [5, 35, 65, 95]).round(4)).tolist()
    return KN


# --- design matrices ---------------------------------------------------------

def build_X(F, KN, linear=False):
    """Covariate matrix for the outcome models. linear=True gives the all-linear sensitivity version."""
    out = {}
    for v in SPL4 + SPL3:
        t = KN[v] if KN.get(v) else KN.get(v + "_3")
        if linear or t is None:
            out[v] = F[v].astype(float).values
        else:
            B = rcs(F[v], t)
            for j in range(B.shape[1]):
                out[v if j == 0 else f"{v}_s{j}"] = B[:, j]
    for v in ["dbp_baseline", "hba1c", "log_idx_los"]:
        out[v] = F[v].astype(float).values
    out["any_ed_prior"]       = (F.n_ed_prior > 0).astype(float).values
    out["log_ed_prior"]       = F.log_ed_prior.astype(float).values
    out["idx_admitted"]       = F.has_idx_adm.astype(float).values
    out["idx_emergency"]      = (F.idx_adm_type == "emergency").astype(float).values
    out["idx_ed_route"]       = F.idx_ed.astype(float).values
    out["idx_disch_facility"] = F.idx_disch.isin(["snf_rehab", "other"]).astype(float).values
    out["female"]             = (F.gender == "F").astype(float).values
    X = pd.DataFrame(out, index=F.index)
    X = pd.concat([X,
                   pd.get_dummies(F.race, prefix="race", dtype=float).drop(columns="race_White"),
                   pd.get_dummies(F.ins_grp, prefix="ins", dtype=float).drop(columns="ins_Other")], axis=1)
    return X.loc[:, X.std() > 0]


def class_dummies(F):
    return (pd.get_dummies(pd.Categorical(F.label, categories=ORDER), dtype=float)
              .drop(columns="Minimal Contact").set_index(F.index))


def gap_spline(F, KN):
    B = rcs(F.gap_days, KN["gap_days"])
    return pd.DataFrame(B, columns=["gap_days"] + [f"gap_days_s{j}" for j in range(1, B.shape[1])], index=F.index)


def with_indicators(F, cols=MISS):
    """Missing-indicator method: median-fill and add an indicator per column."""
    F = F.copy(); ind = pd.DataFrame(index=F.index)
    for c in cols:
        ind[f"{c}_miss"] = F[c].isna().astype(float)
        F[c] = F[c].fillna(F[c].median())
    return F, ind


# --- pooling over imputations (Rubin's rules) --------------------------------

def pool(res):
    """res = list of (coef, vcov, names), one per imputation -> (table, Wald, p)."""
    Q = np.array([r[0] for r in res]); U = np.array([r[1] for r in res]); m = len(res); k = res[0][2]
    qb = Q.mean(0); Ub = U.mean(0)
    B = np.atleast_2d(np.cov(Q, rowvar=False, ddof=1)) if len(k) > 1 else np.array([[Q.var(ddof=1)]])
    Tt = Ub + (1 + 1 / m) * B; se = np.sqrt(np.diag(Tt))
    wald = float(qb @ np.linalg.solve(Tt, qb)); p = chi2.sf(wald, len(k))
    fmi = np.diag((1 + 1 / m) * B) / np.diag(Tt)
    tab = pd.DataFrame({"HR": np.exp(qb), "lo": np.exp(qb - 1.96 * se), "hi": np.exp(qb + 1.96 * se),
                        "p": 2 * (1 - norm.cdf(np.abs(qb / se))), "FMI": fmi}, index=k)
    return tab, wald, p


def pool_full(res):
    """Pooled coefficients and total variance (for spline curves)."""
    Q = np.array([r[0] for r in res]); U = np.array([r[1] for r in res]); m = len(res)
    qb = Q.mean(0); B = np.atleast_2d(np.cov(Q, rowvar=False, ddof=1))
    return qb, U.mean(0) + (1 + 1 / m) * B


def pool_subset(res, keep):
    """Pool only some coefficients from a joint model (e.g. class given gap)."""
    sub = []
    for b, V, k in res:
        ix = [k.index(x) for x in keep]
        sub.append((b[ix], V[np.ix_(ix, ix)], keep))
    return pool(sub)


def single(r):
    """One fit (no imputation) -> same format as pool()."""
    b, V, k = r; se = np.sqrt(np.diag(V))
    w = float(b @ np.linalg.solve(V, b))
    tab = pd.DataFrame({"HR": np.exp(b), "lo": np.exp(b - 1.96 * se), "hi": np.exp(b + 1.96 * se)}, index=k)
    return tab, w, chi2.sf(w, len(k))


def fmt_hr(tab, k):
    return f"{tab.loc[k, 'HR']:.2f} ({tab.loc[k, 'lo']:.2f}-{tab.loc[k, 'hi']:.2f})"


# --- power ---------------------------------------------------------------------

def schoenfeld_mde(events, exposed_share, alpha=0.05, power=0.80):
    """Smallest hazard ratio detectable (Schoenfeld 1983): returns (HR > 1, HR < 1)."""
    z = norm.ppf(1 - alpha / 2) + norm.ppf(power)
    b = z / np.sqrt(events * exposed_share * (1 - exposed_share))
    return float(np.exp(b)), float(np.exp(-b))
