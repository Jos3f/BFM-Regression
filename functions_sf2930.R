# Plot RSS, BIC, Cp, Adjusted R2, R2 for an object of the regsubsets class
plot.regsubsets <- function(regfit) {
    reg.summary <- summary(regfit)
    par(mfrow = c(2, 3))
    plot(reg.summary$rss,
         xlab = "Number of Variables ", ylab = "RSS", type = "l")
    plot(reg.summary$adjr2,
         xlab = "Number of Variables ", ylab = "Adjusted RSq", type = "l")
    points(which.max(reg.summary$adjr2),
           reg.summary$adjr2[which.max(reg.summary$adjr2)],
           col = "red", cex = 2, pch = 20)

    plot(reg.summary$cp, xlab = "Number of Variables ", ylab = "Cp", type = "l")
    points(which.min(reg.summary$cp),
           reg.summary$cp[which.min(reg.summary$cp)],
           col = "red", cex = 2, pch = 20)

    plot(reg.summary$bic,
         xlab = "Number of Variables ", ylab = "BIC", type = "l")
    points(which.min(reg.summary$bic),
           reg.summary$bic[which.min(reg.summary$bic)],
           col = "red", cex = 2, pch = 20)

    plot(reg.summary$rsq,
         xlab = "Number of Variables ", ylab = "R2", type = "l")
}

# Predict new values for an object of the regsubsets class
predict.regsubsets <- function(object, newdata, id, ...) {
    form <- as.formula(object$call[[2]])
    mat <- model.matrix(form, newdata)
    coefi <- coef(object, id = id)
    xvars <- names(coefi)
    mat[, xvars] %*% coefi
}