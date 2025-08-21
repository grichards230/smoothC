#' find_C_smoothC()
#' @description
#' Calculates both the traditional concordance index (C) and the proposed smooth C-statistic
#' for a fitted survival model.
#'
#' @param object a fitted coxph-type model object.
#' @param weighted Boolean, default TRUE. Use IPCW weights to find C-indices?
#' @param dept.censor Boolean, default FALSE. Is censoring is dependent on predictors?
#'
#' @return c.vals List of two numeric values:
#'  * C : traditional concordance index value
#'  * smooth.C : proposed new c-statistic; coefficient of logit regression for standardized risk differences
#' @importFrom survival Surv coxph survfit concordance
find_C_smoothC = function(
        object,
        weighted = TRUE,
        dept.censor = FALSE
)
{
    # Extract event time and status and model predictors. Save to data frame
    predictor.vars = model.frame(object)[-1]
    p = ncol(predictor.vars)
    predictor.names = sapply(seq(1:p), function (x) paste0("Z", x))
    response.vars = object$y

    ## Error coding: Cannot take counting process or other types of Surv() objs
    if(ncol(response.vars) != 2) {
        return(simpleError(
            "Sorry, Smooth C currently only supports right-censored survival model. Please input model object in this format."
        ))
    }

    surv.data = data.frame(
        X = response.vars[, 1],
        delta = response.vars[, 2]
    )
    surv.data[, predictor.names] = predictor.vars

    # Calculate maximum follow-up time
    ### This prevents weights from being too large
    tmax = 0.5 * max(surv.data$X) ## 50% of max event time

    # Predict risk scores
    surv.data$risk = predict(object, type="lp")

    # Order by risk, descending so that index indicates risk score ranking
    surv.data = surv.data[order(surv.data$risk, decreasing = TRUE), ]

    # Convolve risk for all possible R_i - R_j
    # Extract entries above diagonal (i.e. no ties!)
    diff.mat = outer(surv.data$risk, surv.data$risk, "-")
    diff.idx = which(upper.tri(diff.mat), arr.ind = TRUE)

    # Filter indices for comparable pairs (i.e. earlier event is a failure)
    i.idx = diff.idx[, 1]
    j.idx = diff.idx[, 2]
    X_is = surv.data[i.idx, "X"]
    X_js = surv.data[j.idx, "X"]
    delta_is = surv.data[i.idx, "delta"]
    delta_js = surv.data[j.idx, "delta"]

    comp.pair.idx = which((X_is < X_js & delta_is == 1) |
                              (X_is > X_js & delta_js == 1))
    diff.idx = diff.idx[comp.pair.idx, ]

    # Define data frame of relevant values for comparable pairs
    # Note that i<j always => R_i - R_j > 0
    idx1 = diff.idx[, 1]
    idx2 = diff.idx[, 2]
    diff = diff.mat[diff.idx]
    diff = diff / sd(diff) #standardize risk differences
    X1 = surv.data$X[idx1]
    X2 = surv.data$X[idx2]
    X.min = pmin(X1, X2) #earliest event time of the pair. useful later
    Y = as.numeric(X1 < X2)
    pairs = (data.frame(Y, diff, X1, X2, X.min, idx1, idx2))

    # Critical: Sort pairs by indices in original dataset
    # This fixed my code. Not sure why though
    pairs = pairs[order(idx1, idx2), ]

    # Filter out:
    #pairs with equal risk (ties tell us nothing)
    #pairs beyond some observed point "tmax"
    pairs = pairs[which(pairs$diff > 0 & pairs$X.min < tmax), ]

    # Finally, find c-stat and smooth c-stat coefficients. Use weights if req.
    if(weighted == TRUE) {
        ## IPCW: Calculate weights
        ## Adapted from code from Nick Hartman
        if(dept.censor == TRUE) {
            ## Error coding II: Does not support factor variables
            if(any(!is.numeric(predictor.vars))) {
                return(simpleError(
                    "Sorry, Smooth C currently does not support dependent censoring for factor variables."
                ))
            }

            ## Find G per each observation, dependent on predictors
            ## weight = 1 / G1 * G2
            ##  = 1 / (S_0(tmin) ^ [exp(beta'*Z_i) + exp(beta'Z_j)]
            formula.obj.cens = as.formula(
                paste(
                    "Surv(X, 1-delta) ~",
                    paste(predictor.names, collapse = " + ")
                )
            )
            model.fit = coxph(formula.obj.cens, data = surv.data)

            ## Re-order to match survfit() function
            pairs = pairs[order(pairs$X.min), ]

            ## Find baseline survival estimate for each pair's first event time
            S_0 = summary(survfit(model.fit), times = pairs$X.min)$surv

            ## Find proportional hazard term per each pair
            Z_is = as.matrix(surv.data[pairs$idx1, predictor.names])
            Z_js = as.matrix(surv.data[pairs$idx2, predictor.names])
            beta_hat = coef(summary(model.fit))[, 1]
            ph.term_i=apply(Z_is, MARGIN = 1, FUN = function(x) {
                crossprod(beta_hat, x)
            })
            ph.term_j=apply(Z_js, MARGIN = 1, FUN = function(x) {
                crossprod(beta_hat, x)
            })

            ## Put it all together to define weights
            pairs$weights = 1 / (S_0 ^ (exp(ph.term_i)+exp(ph.term_j)))
        } else {
            ## Find G, for each pair directly
            ## weight = 1 / G^2; G corresponds to min(X1, X2)
            model.fit = survfit(Surv(X, 1-delta) ~ 1, data = surv.data)
            pairs = pairs[order(pairs$X.min),]
            #re-order to match survfit() function
            G = summary(model.fit, times = pairs$X.min)$surv
            pairs$weights = 1 / (G^2)
        }

        ## Failsafe: Cap weights at 1000
        pairs$weights = ifelse(pairs$weights < 1000, pairs$weights, 1000)

        # Finally, fit the models
        C = concordance(object, timewt = "n/G2")$concordance ####added weights
        suppressWarnings(
            smooth.C <- glm(
                Y ~ 0 + diff,
                data = pairs,
                family = "binomial",
                weights = pairs$weights)$coef
        ) ####added weights
        return(list(
            C = C,
            smooth.C = c(smooth.C, use.names = F)
        ))
    } else {
        C = concordance(object)$concordance
        smooth.C = glm(Y ~ 0 + diff,
                       data = pairs,
                       family = "binomial")$coef
        return(list(
            C = C,
            smooth.C = c(smooth.C, use.names = F)
        ))
    }
}
