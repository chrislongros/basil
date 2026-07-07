.basil_beta <- function(x, name, table) {
    x <- as.character(x)
    if (length(x) != 1L || !x %in% names(table))
        stop(sprintf("`%s` must be one of: %s", name,
                     paste(names(table), collapse = ", ")), call. = FALSE)
    unname(table[x])
}

.basil_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x))
        stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
    x
}

#' BASIL predicted survival
#'
#' Predicts the probability of survival at one or more time points for a patient
#' with severe limb ischaemia, using the BASIL trial survival prediction model
#' (Bradbury and others, 2010). The model is an accelerated failure time Weibull
#' regression; coefficients are those of the official BASIL prognostic
#' spreadsheet.
#'
#' Two composite variables are entered as categories rather than raw values,
#' because their derivations are not fully specified in the spreadsheet:
#' \code{creatinine} as \code{"low"}, \code{"medium"}, \code{"high"} or
#' \code{"missing"}, and \code{bollinger} (the below-knee Bollinger angiogram
#' score) as \code{"<5"}, \code{">=5"} or \code{"unknown"}. Age is clamped to
#' 40-95 and maximum ankle pressure to 0-200, as in the deployed model.
#'
#' @param age Age in years (clamped to 40-95).
#' @param tissue_loss Logical; tissue loss present.
#' @param bmi Body mass index category: \code{"<20"}, \code{"20-25"},
#'   \code{"25-30"}, \code{"30+"} or \code{"unknown"}.
#' @param creatinine Serum creatinine category: \code{"low"}, \code{"medium"},
#'   \code{"high"} or \code{"missing"}.
#' @param bollinger Below-knee Bollinger angiogram score category: \code{"<5"},
#'   \code{">=5"} or \code{"unknown"}.
#' @param smoking Smoking status: \code{"non"}, \code{"ex"} or \code{"current"}.
#' @param mi_angina Logical; history of myocardial infarction or angina.
#' @param stroke_tia Logical; history of stroke or transient ischaemic attack.
#' @param ankle_pressures_n Number of below-knee ankle pressures measurable
#'   (0, 1, 2 or 3).
#' @param max_ankle_pressure Maximum ankle pressure obtained, mmHg (clamped to
#'   0-200).
#' @param times Numeric vector of times in years at which to return survival
#'   (default 0.5, 1 and 2).
#' @param shrunk Logical; if \code{TRUE} (default) apply the model's uniform
#'   shrinkage, as in the trial's primary output.
#' @return A named numeric vector of survival probabilities at \code{times},
#'   with the linear predictor attached as attribute \code{"lp"}.
#' @examples
#' # The BASIL spreadsheet's worked example
#' basil(age = 80, tissue_loss = TRUE, bmi = "30+", creatinine = "medium",
#'       bollinger = "<5", smoking = "ex", mi_angina = TRUE, stroke_tia = TRUE,
#'       ankle_pressures_n = 0, max_ankle_pressure = 0)
#' @export
basil <- function(age, tissue_loss, bmi, creatinine, bollinger, smoking,
                  mi_angina, stroke_tia, ankle_pressures_n, max_ankle_pressure,
                  times = c(0.5, 1, 2), shrunk = TRUE) {
    tissue_loss <- .basil_flag(tissue_loss, "tissue_loss")
    mi_angina   <- .basil_flag(mi_angina, "mi_angina")
    stroke_tia  <- .basil_flag(stroke_tia, "stroke_tia")
    if (!is.numeric(age) || !is.numeric(max_ankle_pressure) ||
        !is.numeric(times) || anyNA(c(age, max_ankle_pressure, times)))
        stop("`age`, `max_ankle_pressure` and `times` must be non-missing numeric.",
             call. = FALSE)

    age <- min(max(round(age), 40), 95)
    max_ankle_pressure <- min(max(max_ankle_pressure, 0), 200)

    bmi_b   <- .basil_beta(bmi, "bmi",
        c("<20" = -0.58388, "20-25" = 0, "25-30" = 0.02468,
          "30+" = 0.77393, "unknown" = -0.21808))
    creat_b <- .basil_beta(creatinine, "creatinine",
        c("low" = -0.81762, "medium" = 0, "high" = -0.75786, "missing" = -0.97734))
    boll_b  <- .basil_beta(bollinger, "bollinger",
        c("<5" = 0, ">=5" = -0.47981, "unknown" = -0.05290))
    smoke_b <- .basil_beta(smoking, "smoking",
        c("non" = 0, "ex" = -1.08946, "current" = -0.84268))
    ankle_b <- .basil_beta(ankle_pressures_n, "ankle_pressures_n",
        c("0" = 0, "1" = -0.25982, "2" = 0.42329, "3" = 0.79429))

    lp_raw <- 8.09782 +
        (if (tissue_loss) -0.80218 else 0) +
        bmi_b + creat_b + boll_b +
        (-0.04932 * age) +
        smoke_b +
        (if (mi_angina) -0.74512 else 0) +
        (if (stroke_tia) -0.56663 else 0) +
        ankle_b +
        (0.00658 * max_ankle_pressure)

    lp <- if (isTRUE(shrunk)) 2.529597 + 0.75 * (lp_raw - 2.529597) else lp_raw
    surv <- exp(-(exp(-lp) * times)^(1 / 1.437676))
    names(surv) <- paste0("t", times)
    attr(surv, "lp") <- lp
    surv
}

#' BASIL predicted survival for a cohort
#'
#' Applies \code{\link{basil}} to each row of a data frame, appending the linear
#' predictor and the predicted survival at each of \code{times} as columns.
#'
#' @param data A data frame with one row per patient and columns \code{age},
#'   \code{tissue_loss}, \code{bmi}, \code{creatinine}, \code{bollinger},
#'   \code{smoking}, \code{mi_angina}, \code{stroke_tia},
#'   \code{ankle_pressures_n} and \code{max_ankle_pressure}.
#' @param times Numeric vector of times in years (default 0.5, 1 and 2).
#' @param shrunk Logical; passed to \code{\link{basil}}.
#' @return \code{data} with an \code{lp} column and one \code{surv_t<time>}
#'   column for each entry in \code{times}.
#' @examples
#' patients <- data.frame(
#'     age = c(80, 65), tissue_loss = c(TRUE, FALSE),
#'     bmi = c("30+", "20-25"), creatinine = c("medium", "low"),
#'     bollinger = c("<5", ">=5"), smoking = c("ex", "non"),
#'     mi_angina = c(TRUE, FALSE), stroke_tia = c(TRUE, FALSE),
#'     ankle_pressures_n = c(0, 2), max_ankle_pressure = c(0, 80))
#' basil_risk(patients)
#' @export
basil_risk <- function(data, times = c(0.5, 1, 2), shrunk = TRUE) {
    required <- c("age", "tissue_loss", "bmi", "creatinine", "bollinger",
                  "smoking", "mi_angina", "stroke_tia", "ankle_pressures_n",
                  "max_ankle_pressure")
    absent <- setdiff(required, names(data))
    if (length(absent))
        stop("`data` is missing required columns: ",
             paste(absent, collapse = ", "), call. = FALSE)

    fits <- lapply(seq_len(nrow(data)), function(i)
        basil(data$age[i], data$tissue_loss[i], data$bmi[i], data$creatinine[i],
              data$bollinger[i], data$smoking[i], data$mi_angina[i],
              data$stroke_tia[i], data$ankle_pressures_n[i],
              data$max_ankle_pressure[i], times = times, shrunk = shrunk))

    out <- as.data.frame(data)
    out$lp <- vapply(fits, function(s) attr(s, "lp"), numeric(1))
    for (j in seq_along(times))
        out[[paste0("surv_t", times[j])]] <-
            vapply(fits, function(s) as.numeric(s[j]), numeric(1))
    out
}

#' Categorise a serum creatinine value for BASIL
#'
#' Maps a raw serum creatinine value to the category used by the
#' \code{creatinine} argument of \code{\link{basil}}. The cut points are those of
#' the BASIL model (Bradbury and others, 2010): below 88, 88-115, and above 115
#' \eqn{\mu}mol/L for low, medium and high, respectively.
#'
#' @param x Numeric vector of serum creatinine values. \code{NA} maps to
#'   \code{"missing"}.
#' @param units Units of \code{x}: \code{"umol"} (micromoles per litre, the
#'   default) or \code{"mgdl"} (milligrams per decilitre, converted to
#'   \eqn{\mu}mol/L at 1 mg/dL = 88.4 \eqn{\mu}mol/L before categorising).
#' @return A character vector of \code{"low"}, \code{"medium"}, \code{"high"} or
#'   \code{"missing"}, suitable for the \code{creatinine} argument of
#'   \code{\link{basil}}.
#' @examples
#' basil_creatinine_category(c(70, 100, 130, NA))
#' basil_creatinine_category(1.5, units = "mgdl")  # 1.5 mg/dL -> 132.6 umol/L -> "high"
#' @export
basil_creatinine_category <- function(x, units = c("umol", "mgdl")) {
    units <- match.arg(units)
    if (!is.numeric(x))
        stop("`x` must be numeric.", call. = FALSE)
    if (units == "mgdl")
        x <- x * 88.4

    category <- rep("missing", length(x))
    ok <- !is.na(x)
    category[ok & x < 88]               <- "low"
    category[ok & x >= 88 & x <= 115]   <- "medium"
    category[ok & x > 115]              <- "high"
    category
}
