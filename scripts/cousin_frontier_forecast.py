"""
Frozen forecast for A080840(19) (cousin primes below 10^n, upper member counted)
and verdict on the a(19) term added 2026-06-02.

Recipe identical to twin_frontier_forecast.py, applied to A080840 with a(19) HELD OUT.
Additional gate: an independent numpy sieve must reproduce a(8) and a(9) exactly
(the n = 8 run also gives a runtime estimate for n = 9; ~1 GB RAM for the 1e9 sieve).
The Hardy-Littlewood singular series for the d = 4 pair family equals that of twins (2*C2).
Requires: mpmath, numpy.
"""
import json, math, time
import numpy as np
from mpmath import mp, mpf, li, log as mlog, sqrt as msqrt

mp.dps = 50
C2 = mpf('0.66016181584686957392781211001455577843262336028473')
LI2_OFFSET = li(2) - 2/mlog(2)
def Li2(x):
    x = mpf(x)
    return (li(x) - x/mlog(x)) - LI2_OFFSET
def HL(x):
    return 2*C2*Li2(x)

# --- data: OEIS A080840 (values as of 2026-07-02, transcribed from oeis.org;
# a revision-history snapshot is preserved in research/raw/) ---
A080840 = {1: 1, 2: 8, 3: 41, 4: 203, 5: 1216, 6: 8144, 7: 58622, 8: 440258,
           9: 3424680, 10: 27409999, 11: 224373161, 12: 1870585459,
           13: 15834656003, 14: 135779962760, 15: 1177207270204,
           16: 10304191320777, 17: 90948823579815, 18: 808675898548206}
A19_COUSIN = 7237516937438708  # added 2026-06-02 -- HELD OUT of calibration

out = {}

# ---------- GATE: own cousin sieve on n=8, n=9 (exact-match) ----------
def cousin_count_upper_below(N):
    """count pairs (p, p+4) both prime with upper member p+4 < N"""
    sieve = np.ones(N, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if sieve[p]:
            sieve[p*p::p] = False
    return int(np.count_nonzero(sieve[4:] & sieve[:-4]))

gate = {}
t0 = time.time()
c8 = cousin_count_upper_below(10**8)
t8 = time.time() - t0
gate["n8"] = {"own": c8, "oeis": A080840[8], "match": c8 == A080840[8], "time_s": round(t8, 2)}
budget_ok = t8 * 12 < 300  # skip n=9 if the n=8 timing suggests it would exceed ~5 min (n=9 ~ 10x n=8 + overhead)
gate["timing_pilot_n9_budget_ok"] = budget_ok
if gate["n8"]["match"] and budget_ok:
    t0 = time.time()
    c9 = cousin_count_upper_below(10**9)
    gate["n9"] = {"own": c9, "oeis": A080840[9], "match": c9 == A080840[9],
                  "time_s": round(time.time() - t0, 2)}
out["gate_own_sieve"] = gate
if not (gate.get("n8", {}).get("match") and gate.get("n9", {}).get("match")):
    print(json.dumps(out, indent=2)); raise SystemExit("GATE FAILED -- EARLY STOP")

# ---------- CALIBRATE on n<=18 ----------
tbl = []
for n in sorted(A080840):
    x = mpf(10)**n
    S = HL(x)
    d = mpf(A080840[n]) - S
    tbl.append({"n": n, "a_n": A080840[n], "HL": float(S), "delta": float(d),
                "rel_bias_pct": float(100*d/S), "w": float(d/msqrt(S))})
out["residual_table_n13_18"] = [r for r in tbl if r["n"] >= 13]
b = {r["n"]: mpf(r["rel_bias_pct"])/100 for r in tbl}
w = {r["n"]: r["w"] for r in tbl}
ratios = sorted(float(b[n+1]/b[n]) for n in range(14, 18))
gamma = (ratios[1] + ratios[2]) / 2
ws = [w[n] for n in range(13, 19)]
rms_w = math.sqrt(sum(v*v for v in ws)/len(ws))
out["calibration"] = {"w_13_18": ws, "rms_w": rms_w, "gamma_decade": gamma,
                      "w_range_all_n": [min(r["w"] for r in tbl), max(r["w"] for r in tbl)]}

# ---------- FORECAST (frozen from n<=18) ----------
x = mpf(10)**19
S = HL(x)
b_hat = b[18] * mpf(gamma)
P1 = S*(1 + b_hat)
P2 = mpf(A080840[18]) * S/HL(mpf(10)**18)
x_hat = (P1 + P2)/2
sig_spread = abs(P1 - P2)/2
sig_fluct = mpf(rms_w)*msqrt(S)
sigma = msqrt(sig_spread**2 + sig_fluct**2)
out["FORECAST_pi4_1e19"] = {
    "HL": float(S), "P1": float(P1), "P2": float(P2),
    "median": int(round(x_hat)), "sigma": float(sigma),
    "range90": [int(round(x_hat - mpf('1.645')*sigma)), int(round(x_hat + mpf('1.645')*sigma))],
    "range3sigma": [int(round(x_hat - 3*sigma)), int(round(x_hat + 3*sigma))]}

# ---------- VERDICT on the held-out a(19) ----------
innov = (mpf(A19_COUSIN) - x_hat) / sigma
out["VERDICT_a19_cousin"] = {
    "a19": A19_COUSIN,
    "innovation_sigma": float(innov),
    "abs_delta_vs_median": int(A19_COUSIN - int(round(x_hat))),
    "call": "OK (<3 sigma)" if abs(innov) < 3 else "CONFLICT (>=3 sigma)"}

print(json.dumps(out, indent=2))
