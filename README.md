# TCMP and HWRA

**Temporally Correlated Mixture Prior** and **Heavy-Tail-Weighted Risk Allocation**

MATLAB implementations of two core components for robust GNSS/INS estimation
and integrity monitoring. Both components use the Gaussian/Student-t mixture
responsibilities produced by a variational Bayesian (VB) filter:

- **TCMP** combines a Markov prediction of channel consistency with current
  innovation evidence to initialize the mixture prior.
- **HWRA** uses the converged responsibilities to redistribute monitored
  fault-hypothesis probabilities while preserving their total mass.

This repository contains the core formulas. The navigation filter, VB
iterations, satellite-history management, fault detection/exclusion, and
protection-level (PL) solver are supplied by the application.

## Files

| File | Purpose |
|---|---|
| [TCMP.m](TCMP.m) | Markov prediction and adaptive Beta-prior initialization. |
| [HWRA.m](HWRA.m) | Heavy-tail-weighted fault-hypothesis probability allocation. |
| [README.md](README.md) | Usage, interfaces, and integration guidance. |

## Requirements

- MATLAB; no additional toolboxes are required by these two functions.
- Normalized innovations and channel responsibilities from an external filter.
- A fault-hypothesis matrix and its original probability assignment for HWRA.

The datasets and frameworks listed in [Notes](#notes) provide context for the
GNSS/INS experiments; they are not runtime dependencies of these functions.


## TCMP: Temporally Correlated Mixture Prior

```matlab
[a0,b0,aM,K] = TCMP(u,q_prev,gap,pGG,pTT,nu);
```

### Inputs

Vector inputs must contain the same number of channels, `N`. Row and column
vectors are accepted and converted to columns internally.

| Input | Description |
|---|---|
| `u` | Finite normalized innovations: residual divided by the square root of the nominal innovation variance. Signed or absolute values are accepted. |
| `q_prev` | Last converged Gaussian responsibility for each satellite, in `[0,1]`. Use `NaN` when no history is available. |
| `gap` | Positive integer number of Markov steps since each stored responsibility. `NaN` is allowed for channels without history. |
| `pGG` | Scalar probability of remaining in the Gaussian state for one step; strictly between 0 and 1. |
| `pTT` | Scalar probability of remaining in the heavy-tail state for one step; strictly between 0 and 1. |
| `nu` | Positive scalar degrees of freedom of the unit-scale Student-t component. |

Current experimental settings are `pGG=0.95`, `pTT=0.90`, and `nu=5`, with
one Markov step equal to one second. These settings are passed explicitly;
the function does not estimate or calibrate them.

### Outputs

All outputs are `N`-by-1 column vectors.

| Output | Description |
|---|---|
| `a0`, `b0` | Initial Beta shape parameters for the existing VB update, with `a0+b0=1`. |
| `aM` | Markov-predicted Gaussian probability; new channels use the stationary probability. |
| `K` | Adaptive weight applied to current innovation evidence. |

The prediction uses `piG=(1-pTT)/(2-pGG-pTT)` and `phi=pGG+pTT-1` to obtain
`aM=piG+phi^gap*(q_prev-piG)` for each known channel. Current evidence is
`c=BF/(1+BF)`, where `BF` compares the standard Gaussian and unit-scale
Student-t densities at `u`. Their fusion gives
`a0=aM+K.*(c-aM)` and `b0=1-a0`, with
`K=P./(P+R)`, `P=aM.*(1-aM)`, and `R=c.*(1-c)`.

The implementation evaluates evidence in log space and keeps Beta shapes
strictly positive. `P` and `R` are Bernoulli variance proxies, not independently
calibrated Kalman-filter error variances.

## HWRA: Heavy-Tail-Weighted Risk Allocation

```matlab
[p_bar,w] = HWRA(q,keep,p_prior);
```

### Inputs

Let `H` denote the number of hypotheses, including the no-fault hypothesis H0.

| Input | Shape | Description |
|---|---|---|
| `q` | `N` elements | Converged VB **Gaussian** responsibilities in `[0,1]`, aligned with the columns of `keep`. Heavy-tail probabilities are computed internally as `rho=1-q`. |
| `keep` | `H`-by-`N` | Binary matrix: `1` means retained/nonfaulty; `0` means excluded/faulty. Row 1 is H0 and contains only ones. |
| `p_prior` | `H` elements | Nonnegative original hypothesis probabilities, with H0 first. Their sum must not exceed 1, apart from numerical tolerance. |

Any remaining mass may be reserved for unmonitored hypotheses. For example,
`p_prior` in the quick start leaves a mass of `0.0001` outside the supplied set.

### Outputs

Both outputs are `H`-by-1 column vectors.

| Output | Description |
|---|---|
| `p_bar` | Reallocated hypothesis probabilities for the PL solver. H0 and the total monitored fault probability are preserved. |
| `w` | Conditional weights within the monitored fault event. `w(1)=0`; fault weights sum to one when the monitored fault mass is positive. |

Each fault-set score multiplies `rho` over excluded channels and `q` over
retained channels. Normalization excludes H0 and uses only positive-prior
fault rows. Identical physical fault sets are scored once, then their weight
is split equally among duplicate labels. The normalized weights are multiplied
by the original monitored fault mass, `sum(p_prior(2:end))`.

Zero-prior rows remain zero. If the monitored fault mass is zero, `p_bar`
equals `p_prior` and `w` is zero.

## Integration

1. Align each satellite with its stored posterior and elapsed step count.
2. Call TCMP and pass `a0`/`b0` to the existing VB iterations.
3. After convergence, pass `q` to HWRA using the current hypothesis matrix.
4. Calculate PL with `p_bar` and store the full-set converged `q` for the next epoch.

TCMP is stateless: the caller maintains the satellite cache and commits the
full-set posterior once per epoch. Store the converged `q`, rather than `a0`.
HWRA reallocates hypothesis probability mass; it does not change the PL risk
budget or estimate calibrated physical fault rates. The H0 PL coefficient,
fixed at `2` in the current experiments, is managed by the external PL solver
and is distinct from `p_prior(1)`.

## Notes

Dataset: [IPNL-POLYU](https://github.com/IPNL-POLYU)  
GNSS/INS datasets for performance evaluation under different scenarios.

GNSS/INS Tightly Coupled Framework: [TightlyCoupledINSGNSS](https://github.com/benzenemo/TightlyCoupledINSGNSS.git)  
Provides tightly-coupled integration of GNSS and INS, improving robustness in challenging environments.

ARAIM (via MAAST): [Stanford GPS Lab – MAAST](https://gps.stanford.edu/resources/software-tools/maast)  
Used for simulation, integrity monitoring, and performance testing.
