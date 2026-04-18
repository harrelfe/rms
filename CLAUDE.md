# CLAUDE.md — rms Package Development Context

This file provides persistent context for Claude Code sessions on the `rms` package
(Regression Modeling Strategies) by Frank E. Harrell Jr.

- CRAN: https://cran.r-project.org/package=rms
- GitHub: https://github.com/harrelfe/rms
- Current version: 8.1-1

---

## Package Overview

`rms` provides regression modeling, testing, estimation, validation, graphics, prediction,
and typesetting. It is especially designed for binary/ordinal logistic regression, Cox
regression (via `cph`), accelerated failure time models, OLS, the Buckley-James model,
GLS, GLMs, and quantile regression.

All model fitters store enhanced design attributes in the fit object (via `Design()` and
`DesignAssign()`), enabling a uniform downstream API: `predict`, `anova`, `summary`,
`plot`, `nomogram`, `validate`, `calibrate`, etc.

---

## Repository Structure

```
rms/
├── R source: mixed extensions — .r, .s, .R  (legacy S-Plus heritage; all are R code)
├── Fortran:  lrmll.f90, ormll.f90, mlmats.f, robcovf.f90
├── C:        init.c  (routine registration only)
├── NAMESPACE
├── DESCRIPTION
└── inst/tests/   (test scripts; run all after changes)
```

**File extension conventions:**
- `.s` — legacy S-Plus files, most of the core; treat as plain R
- `.r` / `.R` — newer R files
- All are sourced normally by R; no functional difference at runtime

---

## Development Environment

- macOS, zsh shell
- R ≥ 4.4.0
- Workflow: edit source → `R CMD check` → run all files in `inst/tests/` → recompile
  the three books that use `rms` (Regression Modeling Strategies, BBR, RMS course notes)

---

## Project Goals (Priority Order)

### Tier 1 — High value, low risk (do these first)
1. **`strat`/`strata` synonymy** — ~dozen lines across `cph.s` and `predictrms.s`;
   allow users to write either `strat(x)` or `strata(x)` in formulas
2. **`survfit.cph` no-newdata path** — replace the call to `survival:::survfit.coxph`
   with a proxy-object approach: construct a bare `coxph` object carrying just the
   fields `survival::survfit.coxph` needs, then call the public `survfit()` generic.
   This eliminates dependence on a private `survival` API. The proxy approach works
   for all outcome types `cph` supports including counting-process data.
3. **Defensive wrapping of remaining private `survival` API calls** — audit all
   `survival:::` triple-colon calls and either replace with public equivalents or
   wrap in `tryCatch` with clear error messages
4. **NAMESPACE audit** — identify and remove unnecessary exports

### Tier 2 — Medium value, self-contained
- `as.standard()` / `rms_compat` wrapper class with `coef`, `vcov`, and `predict`
  methods targeting the `insight` package API. House in a new `compat.r` file.
- Key mechanism: remap `fit$Design$mmcolnames` to coefficient names.
- Fix known `vcov.orm` dimensionality mismatch first: `coef()` returns all intercepts
  but `vcov()` returns only the mid-intercept — audit internal callers relying on
  `intercepts='mid'` before changing behavior.

### Tier 3 — Documentation
- Improve roxygen2 coverage for: `Design()`, `predictrms()`, `DesignAssign()`,
  `%ia%`, `strat`/`strata` synonymy

### Tier 4 — Mechanical modernization (low-risk, incremental)
- Replace `T`/`F` with `TRUE`/`FALSE` throughout
- Replace `1:length(x)` / `1:nrow(x)` patterns with `seq_along()` / `seq_len()`
- Consolidate `ia_operator.s` (currently duplicated logic)
- Reduce NAMESPACE export surface (Tier 1 audit feeds this)

### Tier 5 — Deferred
- `options(Design.attr)` global side-channel pattern — revisit only after Tiers 1–3

---

## Hard Constraints — Never Violate These

1. **No tidyverse.** Preserve base R style throughout. Do not introduce `dplyr`,
   `tidyr`, `purrr`, `tibble`, or any tidyverse-adjacent package.

2. **Do not refactor Fortran unless explicitly requested.** `lrmll.f90`, `ormll.f90`,
   `mlmats.f`, `robcovf.f90` are working, performance-critical code. Treat as
   read-only.

3. **Changes must be conservative and surgical.** This is a mature, widely-used
   package with many reverse dependencies. Prefer minimal diffs. Flag anything that
   could affect CRAN compliance or break reverse dependencies.

4. **Preserve backward compatibility.** Existing user code must continue to work.
   New behavior should be additive or opt-in where possible.

5. **Do not use `:::` to call internal functions of other packages** (especially
   `survival`) except as a last resort. This is one of the primary things we are
   *fixing*, not introducing.

---

## Critical Architectural Facts

### Short variable naming is load-bearing infrastructure
Names like `age`, `age'`, `age''` are hash keys into `parms`, `limits`, `values`,
and `datadist`. This is not cosmetic — a package-wide rename is impractical and would
break every user's saved fit objects. Do not propose renaming these.

### `Design()` / `DesignAssign()` / `predictrms()`
These three functions form the core of the rms object system. `Design()` processes
the model formula and stores metadata; `DesignAssign()` attaches it to the fit;
`predictrms()` is the workhorse predict method used by almost every model type.
Changes to any of these ripple widely — be especially careful.

### `%ia%` is strictly binary
The `%ia%` interaction operator does not extend to three-way interactions. Its
`nonlinear` tracking breaks when the left argument is itself a `%ia%` result.
For three-way interactions, users should use `*` — `Design()` handles these correctly.

### The `cph`/`coxph` gap
The gap between `cph` and `survival::coxph` is at the `Design()`/formula-processing
level, not the fitter level. Replacing `coxph.fit` with `coxph()` internally unlocks
nothing useful. The proxy-object strategy for `survfit.cph` remains sound.

### `insight` integration target
`insight` currently covers `lrm`, `orm`, `psm` only shallowly via `coef()` wrapping.
The `as.standard()` / `rms_compat` wrapper (Tier 2) is the rms-side fix. After that
is stable, a small PR to `insight` adding proper `cph` support (using existing `lrm`
methods as a template) is planned.

---

## Key Files by Concern

| Concern | Files |
|---|---|
| Cox regression | `cph.s`, `survfit_cph.s`, `survest_cph.s`, `residuals_cph.s`, `validate_cph.s`, `calibrate_cph.s` |
| Logistic regression | `lrm.s`, `lrm_fit.r`, `predict_lrm.s`, `residuals_lrm.s`, `validate_lrm.s` |
| Ordinal regression | `orm.s`, `orm_fit.s`, `adapt_orm.r`, `survest_orm.r`, `survplot_orm.r` |
| Core object system | `rms.s`, `rmsMisc.s`, `predictrms.s` |
| Formula/design processing | `rms_trans.s`, `ia_operator.s`, `datadist.s` |
| `survival` interface | `survfit_cph.s`, `npsurv.s`, `survplot_rms.s` |
| Validation/calibration | `predab_resample.s`, `validate_*.s`, `calibrate_*.s` |
| Planned new file | `compat.r` — `as.standard()` / `rms_compat` wrapper for `insight` |

---

## Testing

After any change:
1. `R CMD check` from the package root — must pass with 0 errors, 0 warnings
2. Run all scripts in `inst/tests/`
3. For changes touching `cph`, `lrm`, `orm`, or `predictrms`: also spot-check against
   the RMS course notes examples if feasible

---

## Dependencies Notes

- **`survival`** — primary integration target; reduce private API calls (Tier 1)
- **`Hmisc`** — sister package by the same author; `rms` depends on it heavily;
  changes to `Hmisc` can affect `rms`
- **`insight` / `easystats`** — interoperability target (Tier 2)
- **`marginaleffects`** — compatibility goal; depends on `insight` support being solid
- **`data.table`** — opportunistic use welcome where it offers clear benefit; do not
  impose it broadly
