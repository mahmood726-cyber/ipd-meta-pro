# IPD Meta-Analysis Pro

IPD Meta-Analysis Pro is a browser-based application for individual participant data meta-analysis. The distributable app remains a single HTML file, but development now runs through a verified module manifest so core statistical code can be edited and tested in smaller pieces.

## Repository layout

- `ipd-meta-pro.html`: distributable single-file app.
- `dev/modules/`: editable split source for the main app bundle.
- `dev/modules/manifest.json`: assembly order for the generated HTML.
- `dev/modules/export_schema_module.js`: export/report schema helper, now bundled into the generated HTML.
- `dev/modules/02_00_vendor_xlsx.js` and `dev/modules/02_00_vendor_jspdf.js`: vendored runtime libraries inlined into the distributable build.
- `dev/modules/02_21a_embedded-validation-manifest.js`: build-generated validation manifest that embeds build-matched benchmark summaries into the distributable for offline evidence views.
- `dev/build.py`: split, rebuild, and verify the module manifest against the distributable HTML.
- `dev/build-scripts/`: validation, browser automation, and R parity gates.

## Development workflow

Use the module manifest as the source of truth for day-to-day changes:

```bash
python dev/build.py build
python dev/build.py verify
python dev/build.py bootstrap-from-html --force
```

- Edit the relevant files in `dev/modules/` rather than patching the large inline script in `ipd-meta-pro.html`.
- Run `build` after module edits to regenerate the distributable HTML.
- Run `verify` before commit or release to ensure the manifest and `ipd-meta-pro.html` are in sync.
- Startup-critical runtime dependencies are vendored into the manifest build so the generated HTML does not need CDN access to boot.
- The build also regenerates an embedded validation manifest from compatible benchmark artifacts in `dev/benchmarks/`, allowing the single-file app to show build-specific validation scorecards in `file://` mode.
- `bootstrap-from-html --force` is recovery-only. It exists for emergency reconstruction from a trusted built artifact, not as a normal development step.
- Legacy scripts that mutate `ipd-meta-pro.html` directly are retired; the only supported write path for the app artifact is `python dev/build.py build`.

## Validation

Fast local browser-facing gate:

```bash
python dev/build-scripts/release_checklist.py --skip-selenium
```

Full local browser gate:

```bash
python dev/build-scripts/release_checklist.py
```

Standalone browser harness for benchmark-backed unit coverage:

```bash
python dev/build-scripts/browser_test_runner.py
```

DOM-free core stats/meta validation:

```bash
node dev/tests/core_stats_meta_test.js
```

Focused real-user-flow smoke suite:

```bash
python dev/build-scripts/user_flow_smoke_test.py
```

Full parity gate against R reference implementations:

```bash
python dev/build-scripts/ipd_parity_gate.py --with-frontier-gap --with-loop7
```

`release_checklist.py` also validates that the core, browser harness, smoke, and parity scripts emit fresh PASS JSON artifacts in `dev/benchmarks/`, so stale benchmark files no longer count as a green gate.

The CI workflow uses `requirements-ci.txt` for Python browser dependencies and `dev/build-scripts/install_r_parity_dependencies.R` for the R benchmark stack.

## Reproducibility

- Manuscript: `F1000_Software_Tool_Article.md`
- Review checklist: `F1000_Submission_Checklist_RealReview.md`
- Release checklist: `RELEASE_CHECKLIST.md`
- Validation workflow: `.github/workflows/validation.yml`

## Methods

The IPD pooling engine is JavaScript, assembled from `dev/modules/` into the single-file `ipd-meta-pro.html` distributable. The numerically-relevant modules:

- **`02_06_stats.js` / `02_01_math-utils.js`** — core descriptive statistics and numerical helpers (Welford running variance, log-sum-exp, beta CDF for Clopper–Pearson intervals at `qbeta(α/2, x, n−x+1)`).
- **`02_08_meta-analysis.js`** — fixed-effect and random-effects pooling on the log-effect scale, with DerSimonian–Laird, REML, and Paule–Mandel τ² estimators and the HKSJ variance correction with the `max(1, Q/(k−1))` floor.
- **`02_09_survival-analysis.js`** — Kaplan–Meier estimation with arm-level IPD, log-rank, and Cox proportional-hazards with Efron tie handling.
- **`02_11_bayesian-mcmc.js`** — light Bayesian path with seeded Metropolis-within-Gibbs sampling; Rhat / ESS diagnostics expose convergence directly in the UI.
- **`02_12_publication-bias.js`** — Egger / radial / Peters tests, trim-and-fill (sensitivity only), conditional PET → PEESE.
- **`02_10_seeded-rng.js`** — Mulberry32 seeded PRNG so any bootstrap, MCMC, or simulation step is bit-reproducible given the same seed.

R-parity gate: `dev/build-scripts/ipd_parity_gate.py` re-runs the same IPD inputs through `metafor`, `meta`, `survival`, and bundled Bayesian R oracles, then enforces pass-rate floors per metric. The gate is part of the release checklist.

## Limitations

- **Single HTML file means single browser tab.** Memory bounds, not numerics, define the largest IPD network this can pool. For very large IPD studies, slice by outcome or trial subset.
- **DL is exposed as an alternative; REML is the default for k ≥ 5.** At k < 5 the engine reports REML τ² with a low-k warning rather than auto-falling-back, because Paule–Mandel and REML both have low-k pathologies; the user is asked to inspect and choose.
- **Cox PH ties use Efron, not exact.** Efron is the default for IPD with frequent ties; exact partial-likelihood is exposed via a flag but is slow for n > a few hundred per arm.
- **R-parity gate needs a working R install.** The gate is the only step that requires R; everything else runs in pure JS / Python. Missing R packages are reported as `unavailable` rather than masked as PASS.
- **Bayesian MCMC is single-chain by default in the UI.** Rhat cannot be computed at 1 chain; the UI surfaces ESS only and warns that production posterior reporting should rerun with ≥ 2 chains via the `02_11_bayesian-mcmc.js` config.
- **No automatic GRADE certainty rating.** GRADE prompts are templated in the UI but the certainty grade is the analyst's call, not the engine's.

## Conclusions

Use IPD Meta-Analysis Pro when (a) the analyst wants IPD-level pooling in a single offline-capable browser tab with an explicit R-parity release gate, (b) the methods set is frequentist NMA / survival / publication-bias plus a light Bayesian path, and (c) reproducibility of every release against R oracles is part of the publication contract. For full posterior inference, very-large-scale IPD federation, or specialised TSA designs not in the Lan–DeMets family, drive a Stan / R pipeline downstream and import the summaries back into this app's reporting layer.
