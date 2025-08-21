#' compare_models()
#' @description
#' Allows the user to easily compare the gains in concordance, using both
#' traditional C-index and the proposed new Smooth C-index, between two nested
#' survival models. Returns concordance values and generates side-by-side plots
#' using these values as parameters to illustrate concordance gains.
#'
#' @param null.vars Vector of variable names for the null model
#' @param added.vars Vector of variable names to be added to the null model to yield the full model
#' @param time Character string. Name of the time-to-event variable
#' @param status Character string. Name of the event status variable
#' @param data Data frame containing time, status, and the variables specified in the null and full models
#' @param weighted Boolean, default TRUE
#' @param dept.censor Boolean, default FALSE
#' @return Data frame of traditional C-indices and proposed new Smooth C-indices
#' for both the null and full models.
#'
#' @importFrom sigmoid sigmoid
compare_models = function(null.vars, added.vars, time, status, data,
                          weighted = TRUE, dept.censor = FALSE) {
    ## Calculate and store concordance values
    null.form = paste0(
        "Surv(",
        time,
        ", ",
        status,
        ") ~ ",
        paste(null.vars, collapse = " + ")
    )
    cvals.null = find_C_smoothC(
        coxph(as.formula(null.form), data = data),
        weighted = weighted,
        dept.censor = dept.censor
    )
    full.form = paste0(c(null.form, added.vars), collapse = " + ")
    cvals.full = find_C_smoothC(
        coxph(as.formula(full.form), data = data),
        weighted = weighted,
        dept.censor = dept.censor
    )

    c.idx1 = cvals.null$C
    c.idx2 = cvals.full$C
    nu1 = cvals.null$smooth.C
    nu2 = cvals.full$smooth.C
    C.df = data.frame(
        null.model = c(c.idx1, nu1),
        full.model = c(c.idx2, nu2),
        row.names = c("C-index (trad)", "Smooth C-index")
    )

    ## Generate x & y values for C-index plot
    ## y-values are binary, based on sign of difference
    x.neg = c(-3, 0)
    x.pos = c(0, 3)
    y.neg1 = rep(1-c.idx1, 2)
    y.pos1 = rep(c.idx1, 2)
    y.neg2 = rep(1-c.idx2, 2)
    y.pos2 = rep(c.idx2, 2)

    ## Generate x-values for Smooth C plot
    x.val = seq(-3, 3, by = 0.01)

    ## Generate y-values for Smooth C plot using sigmoid function
    y.val1 = sigmoid(nu1 * x.val)
    y.val2 = sigmoid(nu2 * x.val)


    #### Plot concordance gains ####
    par(mfrow = c(1, 2))
    plot(
        y.neg1 ~ x.neg,
        type = "l",
        lwd = 5,
        ylab = "Probability T1 < T2",
        xlab = "Difference in Risk Scores (SD)",
        xlim = c(-3, 3),
        ylim = c(0, 1),
        main = "C-Index"
    )
    lines(y.pos1 ~ x.pos, type = "l", lwd = 5)
    lines(y.neg2 ~ x.neg, type = "l", lwd = 5, col = 2)
    lines(y.pos2 ~ x.pos, type = "l", lwd = 5, col = 2)
    plot(
        y.val1 ~ x.val,
        type = "l",
        lwd = 5,
        ylab = "Probability T1 < T2",
        xlab = "Difference in Risk Scores (SD)",
        xlim = c(-3, 3),
        ylim = c(0, 1),
        main = "Smooth C-Index"
    )
    lines(y.val2 ~ x.val, type = "l", lwd = 5, col = 2)

    ## Finally, return data frame of values
    return(C.df)
}
