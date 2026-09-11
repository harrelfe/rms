# Session Summary: Correlation Structures for Longitudinal Ordinal Data (MOST models)

**Purpose of this document:** a complete technical record of one long working session, sufficient to resume this line of work in a new chat or in Claude Code without re-deriving anything. Read top to bottom for the narrative; the boxed formulas and file list are the load-bearing reference material.

**Overarching goal:** specify a longitudinal ordinal state transition model (MOST) that (a) handles absorbing states, (b) fits realistic raw-data correlation patterns, and (c) stays computationally compatible with `orm`'s sparse-Hessian machinery (i.e., avoid full joint/copula likelihoods). The session built up the theory needed to translate a fitted marginal correlation structure into a well-behaved lag-1 (Markov-1) transition model with random effects.

---

## 1. Two custom `nlme::corStruct` classes were built, tested, and delivered

Both follow the same pattern: only the constructor is `@export`ed; internal S3 methods (`corMatrix`, `Initialize`, `coef`, `coef<-`, `print`) use `@exportS3Method`+`@noRd` so they register correctly in `NAMESPACE` without cluttering the public API. Both were verified via `roxygen2::roxygenise()` on a throwaway package, parameter-recovery simulation, `intervals()`, and `lme()` integration.

### `corLinTime.r`
```
corr(t1,t2) = tanh(k + a*|t1-t2| + b*(t1+t2)/2)
```
- 3 params: `k` (intercept), `a` (lag slope), `b` (non-isotropy slope), all on the tanh-linear-predictor scale (already unconstrained, no separate link needed).
- **Superseded by `corFloorExp`** because its asymptote is −1 (or +1) as lag→∞ — implausible for most applications — and it has no literal floor.
- Kept for reference; not the recommended structure going forward.

### `corFloorExp.r` (the recommended corStruct)
```
corr(t1,t2) = k + (1-k)*exp(-exp(a+b*m)*L),   L=|t1-t2|, m=(t1+t2)/2
```
- `k = plogis(k*)` internally (logit link) → automatically in `(0,1)`, a literal non-negative floor.
- Baseline decay rate `exp(a)` at `m=0`; non-isotropy via `b` (rate shifts with pair midpoint).
- `corr=1` exactly at `L=0`; `corr→k` exactly as `L→∞` — both hold by construction, not by fit.
- **Isotropic case (`b=0`) is provably PD for any point set** (convex combination of the all-ones matrix and a standard exponential/OU kernel). Non-isotropic case uses a nearest-PD numerical guard (rarely triggers at realistic parameter magnitudes — empirically ~1/3000 trials at `b~N(0,0.1²)`).
- File also fixes a real bug worth remembering: **non-ASCII em/en-dashes in roxygen comments broke `readLines()` and cascaded into `roxygen2` misparsing block boundaries** (manifesting as `Block must have a @name` and orphaned `@export` tags). Fixed by converting to plain ASCII; don't reintroduce smart punctuation in these files.

### `corExpNonIso` (defined inline, not yet a standalone file — offered twice, not yet built)
```
corr(t1,t2) = exp(-exp(a+b*m)*L)
```
`corFloorExp` with the floor forced to `k=0` — used when an explicit random intercept supplies the floor separately (see §2).

---

## 2. Combining a random intercept with a corStruct — the *right* way to get compound-symmetry + AR(1)

For a **traditional (no-lag) model**, `random=~1|subject` combined with `correlation=corExpNonIso(...)` (or any stationary decay corStruct) reproduces `corFloorExp`'s marginal correlation pattern with high fidelity — this is literally the same decomposition `corFloorExp` is built from (`\sigma_u=\sqrt{k}\sigma`, decay handled by the corStruct). At the *true* parameters this matches exactly (accuracy=1.0); a single finite-sample `lme()` fit showed `r=0.998` between implied and target correlations, with the residual gap being ordinary estimation noise (shrinks with more subjects), not model inadequacy.

**Important negative result:** using this *same* `w(t)`-weighted-random-intercept idea (§4) *without* a lag term does **not** produce compound-symmetry+AR(1) — it produces a qualitatively different "single-factor with time-varying loading" structure, where correlation at a fixed gap depends on *absolute position relative to baseline*, not on the gap itself (verified numerically: gap-5 correlation was 0.317 near baseline vs. 0.138 far from baseline, whereas true AR(1)+RE gives the same 0.683 regardless of position). If you want compound-symmetry+AR(1) with no lags, use the `random=~1|subject`+corStruct combination above, not `w(t)`.

---

## 3. Correlation-structure theory (useful facts, mostly about `corFloorExp`)

- **Mean-time vs. min-time parametrization of the lag term are algebraically identical** once `|t1-t2|` is already in the model (`min(t1,t2) = mean - lag/2`), verified: identical log-likelihood, identical fitted correlations, coefficients related by `a_min = a_mean + b_mean/2`, `b_min=b_mean`. Choice is about which mechanistic story you want to tell (symmetric drift vs. anchored-at-earlier-visit), not about fit.
- **`corFloorExp`'s lag-only decay does not resemble AR(1)** except in narrow, near-saturated parameter regimes. Structural differences: AR(1) has `\rho(0)=1` exactly and asymptotes to 0 (never negative); `corFloorExp`'s earlier tanh-based lag part (superseded) asymptoted to −1; even the exponential decay in `corFloorExp` itself is qualitatively different in shape from AR(1) unless parameters put you in a narrow regime.
- **Absorbing states contaminate raw-data correlation estimates** if not handled: post-absorption pairs show trivially perfect correlation for reasons unrelated to genuine persistence, inflating apparent floor and distorting apparent decay shape. Must either restrict to pre-absorption pairs (landmark-style truncation) or explicitly condition on being alive at both times, *before* using raw-data `corFloorExp` fits as a target.
- **Raw ordinal (numeric-scored) correlation is attenuated relative to the latent-scale correlation** (discretization/thresholding attenuation, same phenomenon polychoric correlation corrects for). Verified via a Gaussian-copula ordinal simulation: latent corr tracked the `corFloorExp` target closely (0.693 vs 0.707 target at short lag), raw-ordinal corr was systematically lower (0.632). Matters for interpreting any raw-data `corFloorExp` fit as a "target" for a latent-scale model.
- **Not all observed non-isotropy in ordinal data reflects genuine non-stationary dependence.** A simulated ordinal Markov chain with a perfectly *constant* (truly isotropic) lag-1 coefficient still showed raw correlation swinging from −0.034 to +0.037 across follow-up, purely from ceiling/floor compression interacting with an ordinary mean trend (marginal variance crashed from 1.39 to 0.145 near a ceiling, then recovered). A large fitted `b` in a raw-data `corFloorExp` fit is therefore a **mixture** of genuine dependence non-stationarity and this "mechanical/compositional" artifact — the diagnostic workflow below (§6) is designed to disentangle them.
- **A semiparametric ordinal copula approach (latent process with `corFloorExp` correlation, `orm`-style semiparametric marginal) was considered and set aside**: it works mechanically (verified by simulation) but (a) needs custom Stan/JAGS code (no off-the-shelf package builds `corFloorExp` in), (b) requires multivariate-normal-CDF-type joint likelihoods that scale poorly and don't decompose per-observation, and (c) **destroys the sparse/banded Hessian structure that makes `orm` fast**, since a joint likelihood over a subject's whole sequence is not a sum of independent per-observation contributions the way a Markov-1 factorization is. This is why the session concluded that generalizing Markov-1 (§4–5) is the right direction rather than the copula route.

---

## 4. The linear MOST model and its exact correlation structure

### Model
$$y_i(t) = a + b\,y_i(t-1) + ct + dt\,y_i(t-1) + w(t)\,u_i + e_i(t)$$
$$u_i\sim N(0,s_u^2) \text{ i.i.d.}, \quad e_i(t)\sim N(0,s_e^2) \text{ i.i.d.}, \quad u_i\perp e_i(t)$$

`y_i(0)` (baseline) is **fixed/conditioned upon** — data, not generated by the `u_i`/`e` mechanism. `t=1` is the **first follow-up occasion** (not baseline) — this indexing distinction matters throughout. `\phi(t):=b+dt` is the effective lag-1 coefficient (`d=0` = isotropic/constant persistence).

Time itself can be on an **arbitrary scale**: what's special about "the first occasion" is its *ordinal position* (`j=1`, no earlier modeled occasion for `u_i` to have entered through), not any numeric value of elapsed time. Keep occasion index `j` and actual elapsed time `t_j` (used inside `c\cdot t_j`, `d\cdot t_j`) conceptually distinct — they coincide only for regularly-spaced (e.g. daily) data.

### General closed-form correlation (any `w(t)`, any `\phi(t)`)
With `P(s,t):=\prod_{r=s+1}^t \phi(r)` (empty product = 1):
$$A(t)=\sum_{s=1}^t P(s,t)\Big|_{w\equiv1}, \quad\text{more generally } A(t)=\phi(t)A(t-1)+w(t),\ A(1)=w(1)$$
$$V(t) = \text{Var}(y(t)), \qquad \mathrm{Corr}(t_1,t_2) = \frac{s_u^2A(t_1)A(t_2) + P(t_1,t_2)\bigl(V(t_1)-s_u^2A(t_1)^2\bigr)}{\sqrt{V(t_1)V(t_2)}}$$

### Four nested special cases (all verified against Monte Carlo)
| Case | Formula |
|---|---|
| No RE, `d` general | `Corr(t1,t2) = P(t1,t2)*sqrt(V(t1)/V(t2))`, `V(t)=s_e^2 \sum_s P(s,t)^2` |
| `d=0`, RE present | `Corr = sqrt(k(t1)k(t2)) + (1-k(t1))b^h sqrt(V(t1)/V(t2))`, `A(t)=(1-b^t)/(1-b)`, `k(t)=s_u^2A(t)^2/V(t)` |
| `d=0`, no RE | `Corr = b^h sqrt(V(t1)/V(t2))` (classic finite-sample AR(1), converges to stationary `b^h` only as `t1\to\infty`) |
| General (both) | the boxed formula above |

`r := s_u^2/s_e^2` is the identified variance ratio; for `d=0`, `k_{\text{stationary}} = r/(r+(1-b)/(1+b))`, invertible: `r = k(1-b)/[(1-k)(1+b)]`. **The stationary `k+(1-k)b^h` form is exact only as `t_1\to\infty`** — near baseline there's a real, quantified deviation (e.g. −0.168 at `t_1=1` vs. the stationary approximation, for `b=0.6`, decaying geometrically at rate `~b^{t_1}`, negligible by `t_1\approx20`).

---

## 5. The `u_i`-accumulation problem and its resolution — the core technical result of the session

### The problem
A **flat** random intercept (`w(t)\equiv1`, i.e. what any standard `lme4`/`nlme::lme()` random-intercepts specification gives you) combined with `y(t-1)` as a covariate causes `u_i`'s effective loading to **accumulate**: `A(t)=\phi(t)A(t-1)+1`, climbing toward `1/(1-b)` (constant `\phi`) — this is not a software limitation, it's an algebraic fact about the recursion (verified: `A(10)=4.46` for `b=0.8` starting from `A(1)=1`). This causes variance to inflate sharply for `b` near 1 (decomposition showed the `r\cdot A(t)^2` term, not ordinary noise buildup, drives ~94% of the inflation at `b=0.8,r=2,t=10`).

**Consequence for matching `corFloorExp`:** this flat-RE-plus-lag class has a **hard ceiling** on reproducing a *constant* floor whenever persistence is non-isotropic (`d\neq0` in this model / non-isotropy in the target) — richer functional forms for the lag term (splines, up to 4-df) did not move this ceiling at all (~0.35–0.4 accuracy regardless). This is a structural fact, not a specification problem.

**Important qualification:** for *regularly-spaced* data with *constant* persistence (`d=0`), the ceiling **disappears entirely** if the model correctly reflects a mature/stationary baseline (not "fresh start") — the full `(k,a)` plane is then exactly reachable (verified exactly across a 9-point grid). The ceiling is specifically a non-isotropy (`d\neq0`) + fresh-baseline problem.

### The fix: weight `u_i` by `w(t)` instead of a flat 1

Progression of attempts (documented because each failure is informative):

1. **`w(t)=1-\phi(t)` for `t>1`, `w(1)=1`** (piecewise) — exact fix, `A(t)\equiv1` for *any* `\phi(t)`, constant or time-varying (verified). Ties `w` to `b,d`.
2. **Trying to avoid the piecewise definition by forcing `\phi(1)=0` via `d=-b`** — technically removes the special case, but forces `\phi(t)=-b(t-1)`, unboundedly negative, causing the *separate* noise-variance component to explode to millions by `t=10`. **Rejected.**
3. **Fully decoupled linear `w(t)=p+qt`, free `(p,q)`** — has a genuine **identifiability problem**: uniform rescaling `(p,q)\to(cp,cq)` with `s_u\to s_u/c` gives an identical model (verified to 8 decimals) — the overall scale of `(p,q)` is unidentified against `s_u` unless anchored (e.g. fix `p=1` or keep `w(1)=1`).
4. **Any *linear*-in-`(t-t_1)` compensating term, anchored or not** (`1+b(t-t_1)`, `1-s(t-t_1)`, `1-st`) — **all diverge** for large `t`, in whichever direction the sign dictates (verified repeatedly, including the case `s=b` which matches the exact fix at `t=1,2` then overshoots and diverges to large negative values by `t=12`). **General lesson: a bounded/contracting recursion cannot be stabilized by an unbounded (linear) forcing term, regardless of sign, magnitude, or how many free parameters it has.**
5. **`w(t) = 1-\lambda\phi(t)`** (one new parameter `\lambda`, tied to whatever `\phi(t)` is) — bounded for any `\lambda`, works identically well whether `d=0` or not (verified across constant and drifting `\phi`, out to `t=40`), recovers the exact fix at `\lambda=1`. Still tied to `\phi(t)`.
6. **Final, fully general recommendation** — smooth, saturating, **zero parameters shared with `a,b,c,d`**:
$$\boxed{w(t) = w_\infty + (1-w_\infty)\,\rho^{\,t-t_1}, \qquad t\ge t_1}$$
   - `w(t_1)=1` falls out automatically (`\rho^0=1`) — no piecewise definition needed at all.
   - **Guaranteed bounded for *any* real `w_\infty`** (verified: `k(t)` stays in `[0,1)` for `w_\infty` ranging from −3 to +20 — this is structural: `k(t)=A(t)^2s_u^2/V(t)` is a ratio of a non-negative quantity to itself-plus-something-non-negative, automatically bounded regardless of `A(t)`'s sign or size). **Do not restrict `w_\infty` to `[0,1]`** — it costs nothing to leave it free and doing so allows the data to reveal patterns (growing influence `w_\infty>1`, sign-flipping `w_\infty<0`) that an artificial restriction would rule out.
   - Verified bounded over long horizons (`t=50`) *provided `\phi(t)` itself stays bounded* — `\phi(t)=b+dt` used raw/unbounded will eventually exceed 1 regardless of `w(t)`'s behavior; this is a separate, already-known requirement (transform `t` via `\sqrt{t}`, a spline, or an `\exp()`-link before it enters `\phi(t)`, mirroring the design principle already used in `corFloorExp`).
   - Nests both extremes: `\rho\to0` recovers the flat-constant version; `w_\infty=1` recovers the original naive/unweighted model.
   - **Robust to omitting `d` (quantified):** simulated from a true `d\neq0` process, fit three versions — correct (`b,d` free): accuracy 1.00 (exact); `d` forced to 0 but `w(t)` free: accuracy 0.74; `d` forced to 0 **and** `w(t)` forced flat: accuracy 0.04. Nearly all the robustness comes from `w(t)`'s flexibility specifically, not from the Markov structure generally — `w(t)`'s free parameters partially absorb the misspecification (fitted `b` inflated from 0.5 true to 0.66, `\rho` collapsed toward 0) rather than the model degrading catastrophically.

**Practical estimation note (applies to every version above except the fully naive one):** because `w(t)` depends on parameters being estimated (either `b,d` directly, or `w_\infty,\rho` jointly with everything else), this is a genuinely nonlinear-in-parameters mixed model — not fittable via plain `lme4`/`nlme::lme()`. Use `nlme::nlme()` or an iterative/profile procedure (fix `w(t)` from provisional parameter estimates, refit, update, repeat).

---

## 6. Diagnostic workflow for real data (outlined, not yet built)

1. Fit `corFloorExp` on **absorption-adjusted**, numeric-scored (or rank/normal-score) raw data → empirical target `(k,a,b)`.
2. Simulate from the fitted lag-1 MOST + `w(t)`-weighted-RE model, using the actual visit-time design, many replicates.
3. Fit `corFloorExp` to each simulated replicate (or average) → distribution of implied `(k,a,b)` for comparison, not just a point estimate.
4. Compare: if empirical and simulated `(k,a,b)` agree, lag-1+`w(t)` is adequate. If not, the *direction* of mismatch is diagnostic:
   - Empirical `b`≠0 not reproduced → check whether it's genuine dependence non-stationarity or the mechanical/compositional artifact (§3) before concluding you need `d\neq0` in the transition model.
   - Empirical floor higher than achievable → possible between-subject heterogeneity in *volatility* (not just level) — consider a random-intercept + random-scale (location-scale) mixed model.
   - Decay shape at moderate lags not reproduced by single lag-1 term → consider a smoothly-weighted window of recent past states rather than jumping to lag-2 (lag-2 requires a pre-baseline measurement, which studies essentially never have).

**Not yet built:** the actual diagnostic tooling (absorption-adjusted target-fitting function; simulate-from-MOST-and-compare machinery). Offered multiple times, not yet requested.

---

## 7. File inventory

| File | Status | Contents |
|---|---|---|
| `corLinTime.r` | delivered, tested | tanh-linear corStruct (superseded by corFloorExp for most purposes) |
| `corFloorExp.r` | delivered, tested, bug-fixed | floor+exponential-decay corStruct (recommended) |
| `test_corLinTime.R`, `test_corFloorExp.R` | delivered | validation scripts (recovery, `intervals()`, `lme()`, AIC comparisons) |
| `corExpNonIso` | **not yet a standalone file** | inline-only; = corFloorExp with k=0; offered as a deliverable, not yet built |

---

## 8. Open items / natural next steps

1. Build `corExpNonIso.r` as a standalone, documented, tested file (parallel treatment to the other two) if still wanted.
2. Build the diagnostic tooling from §6 (absorption handling, simulate-from-MOST, compare-to-`corFloorExp`).
3. Translate the transition-model dataset construction (lag value, gap, absolute time, `w(t)`-weighted random effect) into an `orm`-compatible dataset/fitting procedure — the concrete implementation step that everything above was building toward. Needs the profile/iterative fitting approach since `w(t)`'s parameters aren't estimable via plain linear mixed model machinery.
4. Consider whether the "smoothly-weighted window of recent states" idea (§6, as a lag-2-avoiding way to relax the strict first-order-memory assumption) is worth developing further.
5. If pursuing the location-scale (patient-varying volatility) generalization flagged in §6/§3, that would need its own derivation — not yet started.
