"""
Frozen forecast for A007508(19) and A007508(20), and verdict on the disputed a(19).

Stages: validation gates -> calibration (n <= 18 only) -> forecast registration -> verdict.
Recipe (fixed before any comparison):
  P1 = HL(x)*(1 + b_hat), b_hat = b_18 * gamma, gamma = median(b_{n+1}/b_n, n = 14..17)
  P2 = a_18 * HL(x)/HL(1e18)
  x_hat = (P1+P2)/2; sigma = sqrt(sigma_spread^2 + sigma_fluct^2),
  sigma_spread = |P1-P2|/2; sigma_fluct = RMS(w_n, n = 13..18)*sqrt(HL(x)), w_n = delta_n/sqrt(HL)
  90% range = x_hat +/- 1.645*sigma; failure criterion: an independently verified value outside 3 sigma.
The a(19) term added 2026-06-02 is HELD OUT of all calibration. Requires: mpmath.
"""
import json
from mpmath import mp, mpf, li, log as mlog, sqrt as msqrt

mp.dps = 50
C2 = mpf('0.66016181584686957392781211001455577843262336028473')

# --- exact Li2(x) = int_2^x dt/ln^2 t = (li(x) - x/ln x) - (li(2) - 2/ln 2) ---
LI2_OFFSET = li(2) - 2/mlog(2)
def Li2(x):
    x = mpf(x)
    return (li(x) - x/mlog(x)) - LI2_OFFSET

def HL(x):
    return 2*C2*Li2(x)

# --- independent trapezoid Li2 (for the instrument gate) ---
import math
def Li2_trap(N, steps=400000):
    a, b = math.log(2.0), math.log(N)
    h = (b - a) / steps
    s = 0.0
    prev_f = math.exp(a) / (a * a)
    for k in range(1, steps + 1):
        t = a + k * h
        f = math.exp(t) / (t * t)
        s += 0.5 * (prev_f + f) * h
        prev_f = f
    return s

# --- data: OEIS A007508 (fetched 2026-07-02; raw copies in research/raw/) ---
A007508 = {1: 2, 2: 8, 3: 35, 4: 205, 5: 1224, 6: 8169, 7: 58980, 8: 440312,
           9: 3424506, 10: 27412679, 11: 224376048, 12: 1870585220,
           13: 15834664872, 14: 135780321665, 15: 1177209242304,
           16: 10304195697298, 17: 90948839353159, 18: 808675888577436}
A19_DISPUTED = 7237516880334496  # added 2026-06-02 -- HELD OUT of calibration
SIEVE_CHECK = {9: 3424506, 10: 27412679, 11: 224376048, 12: 1870585220}  # independent sieve counts

out = {}

# ---------- STAGE 1: GATES ----------
gateA = all(A007508[n] == SIEVE_CHECK[n] for n in SIEVE_CHECK)
gateB_pairs = []
for n in (6, 9, 12):
    exact = Li2(mpf(10)**n)
    trap = Li2_trap(10**n)
    rel = abs(mpf(trap) - exact)/exact
    gateB_pairs.append({"n": n, "Li2_exact": float(exact), "Li2_trap": trap, "rel_diff": float(rel)})
gateB = all(p["rel_diff"] < 1e-6 for p in gateB_pairs)
anchor_8248 = float(HL(mpf(10)**6))
gateB = gateB and abs(anchor_8248 - 8248.0) < 1.0
out["gates"] = {"A_data_sieve_vs_OEIS": "4/4 PASS" if gateA else "FAIL",
                "B_instrument": {"pairs": gateB_pairs, "HL(1e6)": anchor_8248,
                                 "verdict": "PASS" if gateB else "FAIL"}}
if not (gateA and gateB):
    print(json.dumps(out, indent=2)); raise SystemExit("GATE FAILED -- EARLY STOP")

# ---------- STAGE 2: CALIBRATE on n<=18 ----------
tbl = []
for n in sorted(A007508):
    x = mpf(10)**n
    S = HL(x)
    d = mpf(A007508[n]) - S
    tbl.append({"n": n, "a_n": A007508[n], "HL": float(S),
                "delta": float(d), "rel_bias_pct": float(100*d/S),
                "w": float(d/msqrt(S))})
out["residual_table"] = tbl
b = {r["n"]: mpf(r["rel_bias_pct"])/100 for r in tbl}
w = {r["n"]: r["w"] for r in tbl}
ratios = sorted(float(b[n+1]/b[n]) for n in range(14, 18))
gamma = (ratios[1] + ratios[2]) / 2
ws = [w[n] for n in range(13, 19)]
import statistics
rms_w = math.sqrt(sum(v*v for v in ws)/len(ws))
out["calibration"] = {"w_13_18": ws, "mean_w": statistics.mean(ws), "std_w": statistics.pstdev(ws),
                      "rms_w": rms_w, "gamma_decade": gamma, "ratios_14_17": ratios}

# ---------- STAGE 3: REGISTER FORECASTS (frozen from n<=18) ----------
def forecast(n_target, b18, a18):
    x = mpf(10)**n_target
    S = HL(x)
    b_hat = b18 * mpf(gamma)**(n_target - 18)
    P1 = S*(1 + b_hat)
    P2 = mpf(a18) * S/HL(mpf(10)**18)
    x_hat = (P1 + P2)/2
    sig_spread = abs(P1 - P2)/2
    sig_fluct = mpf(rms_w)*msqrt(S)
    sigma = msqrt(sig_spread**2 + sig_fluct**2)
    return {"target": f"pi2(10^{n_target})", "HL": float(S), "b_hat_pct": float(100*b_hat),
            "P1": float(P1), "P2": float(P2),
            "median": int(round(x_hat)), "sigma": float(sigma),
            "sigma_rel_pct": float(100*sigma/x_hat),
            "range90": [int(round(x_hat - mpf('1.645')*sigma)), int(round(x_hat + mpf('1.645')*sigma))],
            "range3sigma": [int(round(x_hat - 3*sigma)), int(round(x_hat + 3*sigma))]}

f19 = forecast(19, b[18], A007508[18])
f20 = forecast(20, b[18], A007508[18])
out["FORECAST_pi2_1e19"] = f19
out["FORECAST_pi2_1e20"] = f20

# ---------- STAGE 4: VERDICT on the held-out a(19) ----------
innov = (mpf(A19_DISPUTED) - mpf(f19["median"])) / mpf(f19["sigma"])
out["VERDICT_a19"] = {"a19_disputed": A19_DISPUTED,
                      "innovation_sigma": float(innov),
                      "abs_delta_vs_median": int(A19_DISPUTED - f19["median"]),
                      "call": "OK (<3 sigma)" if abs(innov) < 3 else "CONFLICT (>=3 sigma)"}

print(json.dumps(out, indent=2))
