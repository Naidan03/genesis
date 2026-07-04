# Numerical verification of the A007508/A080840 a(19) anomaly

For [A007508](https://oeis.org/A007508) (twin prime pairs below 10^n) and [A080840](https://oeis.org/A080840) (cousin primes below 10^n; only the upper member of each pair is counted), define normalized residuals

w(n) = (a(n) − 2C₂·Li₂(10ⁿ)) / √(2C₂·Li₂(10ⁿ)),  Li₂(x) = ∫₂ˣ dt/ln²t,

where C₂ = 0.6601618158… is the twin prime constant (the same singular series applies to the d = 4 family). Over the range n = 8..18, both sequences satisfy |w| ≤ 1.13 (and |w| ≤ 1.75 over all verified terms, n = 1..18).

The terms added on 2026-06-02 break this trend in a single decade: w ≈ −13.9 for A007508(19) and w ≈ −13.2 for A080840(19) (`scripts/frontier_check.py`). The two families fall short of the Hardy–Littlewood prediction by 1.18×10⁹ and 1.13×10⁹ respectively, with deficits differing by only 5.7×10⁷ — consistent with a single systematic error in one sieving campaign rather than a simultaneous breakdown of the prediction in two independent families. The extrapolation scripts below flag the same conflict through an independent hold-out test.

We do not claim a(19) is wrong — we report a measured conflict with the Hardy–Littlewood prediction. The A007508 entry states (H. Pfoertner, Jun 17 2026): "An independent verification of a(19) is required." This repository responds to that request.

## Reproduce

```
pip install mpmath numpy
python scripts/twin_frontier_forecast.py
python scripts/cousin_frontier_forecast.py
python scripts/frontier_check.py
```

`twin_frontier_forecast.py` runs consistency checks on the input data and the Li₂ quadrature, calibrates a Hardy–Littlewood extrapolation on n ≤ 18 only, and reports how far the disputed A007508(19) term falls from it. `cousin_frontier_forecast.py` validates an independent sieve against known terms (≈1 GB RAM) and applies the same analysis to A080840(19). `frontier_check.py` is an independent Simpson-quadrature cross-check of the residuals of both families.

All calibration uses n ≤ 18 only; the disputed a(19) values are held out. Data provenance: the OEIS A007508 page was fetched on 2026-07-02, and the raw HTML/text responses (with HTTP headers for the curl fetch and SHA256 checksums) are preserved in `research/raw/`, together with revision-history snapshots for A007508 and A080840; the A080840 terms used by the scripts were transcribed from oeis.org. A windowed spot check uses the [primesieve](https://github.com/kimwalisch/primesieve) 12.14 command-line binary (not bundled): the twin pair count over [10¹⁹ − 10⁹, 10¹⁹] is 690,217, against a Hardy–Littlewood expectation of ≈689,800 — reproducible by anyone with the primesieve CLI.

## Repository map

- `scripts/` — all analysis code
- `research/raw/` — data provenance (raw fetches, HTTP headers, SHA256 checksums)

## License

MIT — see [LICENSE](LICENSE).
