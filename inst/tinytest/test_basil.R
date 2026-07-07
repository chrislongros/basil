library(basil)

## --- the BASIL spreadsheet's worked example --------------------------------
## age 80, tissue loss, BMI 30+, medium creatinine, Bollinger <5, ex-smoker,
## MI/angina, stroke/TIA, no measurable ankle pressures. Shrunk linear
## predictor 1.924469; survival 0.8505 / 0.7693 / 0.6540 at 0.5 / 1 / 2 yr.
ex <- basil(age = 80, tissue_loss = TRUE, bmi = "30+", creatinine = "medium",
            bollinger = "<5", smoking = "ex", mi_angina = TRUE, stroke_tia = TRUE,
            ankle_pressures_n = 0, max_ankle_pressure = 0)

expect_equal(attr(ex, "lp"), 1.924469, tolerance = 1e-5)
expect_equal(as.numeric(ex),
             c(0.8505195, 0.7693470, 0.6539923), tolerance = 1e-6)
expect_equal(names(ex), c("t0.5", "t1", "t2"))

## unshrunk linear predictor (raw sum of coefficients)
raw <- basil(age = 80, tissue_loss = TRUE, bmi = "30+", creatinine = "medium",
             bollinger = "<5", smoking = "ex", mi_angina = TRUE, stroke_tia = TRUE,
             ankle_pressures_n = 0, max_ankle_pressure = 0, shrunk = FALSE)
expect_equal(attr(raw, "lp"), 1.722760, tolerance = 1e-5)

## survival is monotonically non-increasing in time and bounded in (0, 1)
expect_true(all(diff(ex) <= 0))
expect_true(all(ex > 0 & ex < 1))

## --- clamping ---------------------------------------------------------------
## age below 40 and above 95 clamp to the endpoints
base <- list(tissue_loss = FALSE, bmi = "20-25", creatinine = "medium",
             bollinger = "<5", smoking = "non", mi_angina = FALSE,
             stroke_tia = FALSE, ankle_pressures_n = 0, max_ankle_pressure = 0)
lp <- function(age) attr(do.call(basil, c(list(age = age), base)), "lp")
expect_equal(lp(30), lp(40))
expect_equal(lp(120), lp(95))
expect_true(lp(40) > lp(95))

## max ankle pressure clamps to 0-200
mp <- function(p) attr(do.call(basil,
    c(modifyList(base, list(max_ankle_pressure = p)), list(age = 70))), "lp")
expect_equal(mp(-50), mp(0))
expect_equal(mp(300), mp(200))

## --- category handling ------------------------------------------------------
## reference categories (coefficient 0) leave the linear predictor unchanged
b0 <- do.call(basil, c(list(age = 70), base))
expect_equal(attr(b0, "lp"),
             attr(do.call(basil, c(list(age = 70),
                 modifyList(base, list(bmi = "20-25")))), "lp"))

## an unknown category level is rejected
expect_error(do.call(basil, c(list(age = 70),
    modifyList(base, list(bmi = "obese")))), "bmi")
expect_error(do.call(basil, c(list(age = 70),
    modifyList(base, list(creatinine = "normal")))), "creatinine")

## --- input validation -------------------------------------------------------
expect_error(do.call(basil, c(list(age = NA), base)), "numeric")
expect_error(do.call(basil,
    c(list(age = 70), modifyList(base, list(tissue_loss = "yes")))),
    "tissue_loss")
expect_error(do.call(basil,
    c(list(age = 70), modifyList(base, list(mi_angina = NA)))),
    "mi_angina")

## --- cohort interface -------------------------------------------------------
patients <- data.frame(
    age = c(80, 65), tissue_loss = c(TRUE, FALSE),
    bmi = c("30+", "20-25"), creatinine = c("medium", "low"),
    bollinger = c("<5", ">=5"), smoking = c("ex", "non"),
    mi_angina = c(TRUE, FALSE), stroke_tia = c(TRUE, FALSE),
    ankle_pressures_n = c(0, 2), max_ankle_pressure = c(0, 80),
    stringsAsFactors = FALSE)
scored <- basil_risk(patients)

## the first row matches the single-patient worked example
expect_equal(scored$lp[1], 1.924469, tolerance = 1e-5)
expect_equal(scored$surv_t2[1], 0.6539923, tolerance = 1e-6)
expect_equal(names(scored)[11:14], c("lp", "surv_t0.5", "surv_t1", "surv_t2"))

## a missing required column is reported
expect_error(basil_risk(patients[, -1]), "age")

## --- coefficients match the Raab prognostic spreadsheet ---------------------
## each unshrunk linear-predictor delta from a reference patient equals the
## spreadsheet coefficient exactly.
ref <- list(age = 70, tissue_loss = FALSE, bmi = "20-25", creatinine = "medium",
            bollinger = "<5", smoking = "non", mi_angina = FALSE,
            stroke_tia = FALSE, ankle_pressures_n = 0, max_ankle_pressure = 0)
rawlp <- function(...) attr(do.call(basil,
    c(modifyList(ref, list(...)), list(shrunk = FALSE))), "lp")
delta <- function(...) rawlp(...) - rawlp()

expect_equal(delta(tissue_loss = TRUE),       -0.80218, tolerance = 1e-8)
expect_equal(delta(bmi = "<20"),              -0.58388, tolerance = 1e-8)
expect_equal(delta(bmi = "25-30"),             0.02468, tolerance = 1e-8)
expect_equal(delta(bmi = "30+"),               0.77393, tolerance = 1e-8)
expect_equal(delta(bmi = "unknown"),          -0.21808, tolerance = 1e-8)
expect_equal(delta(creatinine = "low"),       -0.81762, tolerance = 1e-8)
expect_equal(delta(creatinine = "high"),      -0.75786, tolerance = 1e-8)
expect_equal(delta(creatinine = "missing"),   -0.97734, tolerance = 1e-8)
expect_equal(delta(bollinger = ">=5"),        -0.47981, tolerance = 1e-8)
expect_equal(delta(bollinger = "unknown"),    -0.05290, tolerance = 1e-8)
expect_equal(delta(smoking = "ex"),           -1.08946, tolerance = 1e-8)
expect_equal(delta(smoking = "current"),      -0.84268, tolerance = 1e-8)
expect_equal(delta(mi_angina = TRUE),         -0.74512, tolerance = 1e-8)
expect_equal(delta(stroke_tia = TRUE),        -0.56663, tolerance = 1e-8)
expect_equal(delta(ankle_pressures_n = 1),    -0.25982, tolerance = 1e-8)
expect_equal(delta(ankle_pressures_n = 2),     0.42329, tolerance = 1e-8)
expect_equal(delta(ankle_pressures_n = 3),     0.79429, tolerance = 1e-8)
expect_equal(delta(max_ankle_pressure = 100),  0.658,   tolerance = 1e-8)  # 100 * 0.00658
expect_equal(rawlp(age = 71) - rawlp(age = 70), -0.04932, tolerance = 1e-8)

## the spreadsheet's worked example, unshrunk survival
ex_unshrunk <- basil(age = 80, tissue_loss = TRUE, bmi = "30+",
    creatinine = "medium", bollinger = "<5", smoking = "ex", mi_angina = TRUE,
    stroke_tia = TRUE, ankle_pressures_n = 0, max_ankle_pressure = 0,
    shrunk = FALSE)
expect_equal(as.numeric(ex_unshrunk),
             c(0.8300289, 0.7395539, 0.6134708), tolerance = 1e-6)

## --- published worked examples (Bradbury and others, 2010, Table V) ---------
## patients A and B reproduce the paper's rounded 6-month/1-year/2-year
## survival percentages exactly. (Table V patients C and D do not reconcile
## with the paper's own coefficient table and are not asserted here.)
pctA <- round(100 * as.numeric(basil(age = 79, tissue_loss = TRUE,
    bmi = "20-25", creatinine = "low", bollinger = ">=5", smoking = "ex",
    mi_angina = TRUE, stroke_tia = FALSE, ankle_pressures_n = 0,
    max_ankle_pressure = 0)))
expect_equal(pctA, c(71, 57, 40))

pctB <- round(100 * as.numeric(basil(age = 80, tissue_loss = TRUE,
    bmi = "25-30", creatinine = "low", bollinger = "unknown", smoking = "ex",
    mi_angina = FALSE, stroke_tia = FALSE, ankle_pressures_n = 1,
    max_ankle_pressure = 60)))
expect_equal(pctB, c(84, 75, 63))

## --- creatinine category helper ---------------------------------------------
## cut points from the paper: <88 low, 88-115 medium, >115 high (umol/L)
expect_equal(basil_creatinine_category(c(70, 88, 100, 115, 130, NA)),
             c("low", "medium", "medium", "medium", "high", "missing"))

## mg/dL input is converted at 1 mg/dL = 88.4 umol/L (79.6 / 97.2 / 132.6)
expect_equal(basil_creatinine_category(c(0.9, 1.1, 1.5), units = "mgdl"),
             c("low", "medium", "high"))

## the helper output feeds straight into basil()
expect_equal(basil_creatinine_category(100), "medium")
expect_silent(basil(age = 70, tissue_loss = FALSE, bmi = "20-25",
    creatinine = basil_creatinine_category(100), bollinger = "<5",
    smoking = "non", mi_angina = FALSE, stroke_tia = FALSE,
    ankle_pressures_n = 0, max_ankle_pressure = 0))

expect_error(basil_creatinine_category("100"), "numeric")
