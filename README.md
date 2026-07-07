# basil

Survival prediction for severe limb ischaemia in R, using the prognostic model
from the BASIL trial (Bypass versus Angioplasty in Severe Ischaemia of the Leg).

The model is an accelerated failure time Weibull regression. Given a patient's
baseline characteristics it returns the predicted probability of survival at one
or more time points.

## Install

```r
# install.packages("remotes")
remotes::install_github("chrislongros/basil")
```

## Use

```r
library(basil)

basil(age = 80, tissue_loss = TRUE, bmi = "30+", creatinine = "medium",
      bollinger = "<5", smoking = "ex", mi_angina = TRUE, stroke_tia = TRUE,
      ankle_pressures_n = 0, max_ankle_pressure = 0)
#>      t0.5        t1        t2
#> 0.8505195 0.7693470 0.6539923
```

The return value is survival probabilities at `times` (default 0.5, 1 and 2
years), with the linear predictor attached as attribute `"lp"`.

Two variables are entered as categories rather than raw values. `bollinger`, the
below-knee Bollinger angiogram score, is one (`"<5"`, `">=5"`, `"unknown"`). The
other is `creatinine`, whose categories follow the model's serum creatinine cut
points:

| Category | Serum creatinine | `creatinine =` |
|----------|------------------|----------------|
| Low | < 88 µmol/L (< 1.0 mg/dL) | `"low"` |
| Medium (reference) | 88–115 µmol/L (1.0–1.3 mg/dL) | `"medium"` |
| High | > 115 µmol/L (> 1.3 mg/dL) | `"high"` |
| Missing | not measured | `"missing"` |

`basil_creatinine_category()` maps a raw value to the right category, in either
unit:

```r
basil_creatinine_category(103)                 # "medium"
basil_creatinine_category(1.5, units = "mgdl") # "high"
```

Age is clamped to 40–95 and maximum ankle pressure to 0–200, as in the deployed
model. See `?basil` for the full argument list.

## Source

- Bradbury AW, Adam DJ, Bell J, Forbes JF, Fowkes FGR, Gillespie I, Ruckley CV,
  Raab GM (2010). Bypass versus Angioplasty in Severe Ischaemia of the Leg
  (BASIL) trial: a survival prediction model to facilitate clinical decision
  making. *Journal of Vascular Surgery* 51(5 Suppl):52S–68S.

## Scope

This package is for research and audit. It is not a substitute for clinical
judgement, and its output must not be used as the sole basis for any clinical
decision.
Confirm the coefficients against the original model before use.

## License

MIT.
