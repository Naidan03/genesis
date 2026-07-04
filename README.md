# Genesis

Pre-registered forecasts and independent numerical verification for prime-counting frontiers. Every prediction is registered before measurement, together with its σ and an explicit failure criterion (default 3σ).

## The a(19) anomaly (A007508 / A080840)

For [A007508](https://oeis.org/A007508) (twin prime pairs below 10^n) and [A080840](https://oeis.org/A080840) (cousin primes below 10^n; only the upper member of each pair is counted), define normalized residuals

w(n) = (a(n) − 2C₂·Li₂(10ⁿ)) / √(2C₂·Li₂(10ⁿ)),  Li₂(x) = ∫₂ˣ dt/ln²t,

where C₂ = 0.6601618158… is the twin prime constant (the same singular series applies to the d = 4 family). Over the eleven-decade calibration range n = 8..18, both sequences satisfy |w| ≤ 1.1 (and |w| ≤ 1.75 over all verified terms, n = 1..18).

The terms added on 2026-06-02 break this trend in a single decade: w ≈ −13.9 for A007508(19) and w ≈ −13.2 for A080840(19). Both families fall short of the Hardy–Littlewood prediction by ≈1.1×10⁹, with deficits differing by only 5.7×10⁷ — the signature of a single systematic error in one sieving campaign rather than a simultaneous breakdown of the prediction in two independent families. Three independent methods agree on the size of the conflict: a Simpson-quadrature z-proxy (−13.2σ) and two held-out forecasts calibrated only on n ≤ 18 (twin: −16.9σ; cousin: −17.28σ).

We do not claim a(19) is wrong — we report a measured conflict with falsifiable predictions for the true values. The A007508 entry states (H. Pfoertner, Jun 17 2026): "An independent verification of a(19) is required." This repository responds to that request.

## Reproduce

```
pip install mpmath numpy
python scripts/twin_frontier_forecast.py
python scripts/cousin_frontier_forecast.py
python scripts/frontier_check.py
```

`twin_frontier_forecast.py` runs the data and instrument validation gates, freezes the twin forecast on n ≤ 18, and reports the verdict on A007508(19) (−16.9σ). `cousin_frontier_forecast.py` validates an independent sieve against known terms (≈1 GB RAM), freezes the cousin forecast, and reports the verdict on A080840(19) (−17.28σ). `frontier_check.py` is an independent Simpson-quadrature cross-check of the residuals of both families.

All calibration uses n ≤ 18 only; the disputed a(19) values are held out. Data provenance: OEIS pages for A007508, A080840 and A152051 fetched on 2026-07-02; raw HTML/text responses with HTTP headers and SHA256 checksums are preserved in `research/raw/`. Windowed spot checks use the [primesieve](https://github.com/kimwalisch/primesieve) 12.14 command-line binary (not bundled): the twin pair count over [10¹⁹ − 10⁹, 10¹⁹] is 690,217, against a Hardy–Littlewood expectation of ≈689,800 — reproducible by anyone with the primesieve CLI.

## Repository map

- `scripts/` — all analysis code
- `research/raw/` — data provenance (raw fetches, HTTP headers, SHA256 checksums)

## License

MIT — see [LICENSE](LICENSE).
