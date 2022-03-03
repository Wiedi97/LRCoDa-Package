#' @importFrom stats lm
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @import tidyr
#' @importFrom robCompositions pivotCoord
NULL


#' Classical and robust regression of non-compositional (real) response on compositional and/or non-compositional predictors
#'
#' Delivers appropriate inference for regression of y on compositional and/or non-compositional input X.
#'
#' Compositional explanatory variables should
#'
#' @aliases LRCoDa ilrregression robilrregression
#' @param y The response which should be non-compositional (Better explanation)
#' @param X The compositional and/or non-compositional predictors as a matrix, data.frame or numeric vector (Probably a refactoring done)
#' @param external Specify the columns which are not part of the composition
#' @param factor_column Specify the column which includes factor levels
#' @param method If robust, the fast MM-type linear regression is applied, while with method
#' \dQuote{classical}, the conventional least squares regression is applied.
#' @param method_pivot Method to choose the pivot coordinates (parameters pivotvar and norm have then no effect)
#' @return An object of class \sQuote{lmrob} or \sQuote{lm} and two summary
#' objects.
#' @author Roman Wiedemeier
#' @seealso \code{\link{lm}} and \code{\link{lmrob}}
#' @keywords models compositional data linear regression
#' @export
#' @examples
#'
#' ## How the content of sand of the agricultural and grazing land soils
#' ## in Germany depend on relative contributions of the main chemical trace elements,
#' ## their different soil types and the Annual mean temperature:
#' data("gemas")
#' gemas$COUNTRY <- as.factor(gemas$COUNTRY)
#' gemas_GER <- dplyr::filter(gemas, gemas$COUNTRY == 'GER')
#' y <- gemas_GER$sand
#' X <- dplyr::select(gemas_GER, c(MeanTemp, soilclass, Al:Zr))
#' LRCoDa(y, X, external = c('MeanTemp'), factor_column = 'soilclass', method='classical', method_pivot = 'orthonormal')
#' LRCoDa(y, X, external = c('MeanTemp'), factor_column = 'soilclass', method='robust', method_pivot = 'orthonormal')
LRCoDa <- function (y, X, external = NULL, factor_column = NULL, method = "robust", method_pivot = 'orthonormal') { # ltsReg mit lmrob ersetzen und dann sollte die Fehlermeldung verschwinden

  if (!is.null(external)) {
    external_col <- X %>% select(external)
    n_externals <- length(external_col)
  }
  if (!is.null(factor_column)) {
    factor_var <- X %>% select(factor_column)
    factor_col <- factor_var[, 1]
    n_levels <- length(unique(as.character(factor_col)))
  }
  ilrregression <- function(X, y, external, factor_column, method_pivot) {

    if (!is.null(factor_column) & !is.null(external)){
      X_selected <- X %>% select(-c(external, factor_column))
      ZV <- data.frame(Factor = factor_col, Externals = external_col, X = X_selected)
      d <- data.frame(y = y, X = ZV)
      lmcla <- lm(y ~ ., data = d)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:(ncol(X)-1-n_externals)) {
        X_prep <- X %>% select(-c(external, factor_column))   # Take X_selected
        Zj <- pivotCoord(cbind(X_prep[, j], X_prep[, -j]), method = method_pivot)
        ZVj <- data.frame(Factor = factor_col, Externals = external_col, Z = Zj)
        dj <- data.frame(y = y, Z = ZVj)
        res <- lm(y ~ ., data = dj)
        res.sum <- summary(res)
        if (j == 1) {
          ilr.sum$coefficients[1:(n_levels+n_externals+1), ] <- res.sum$coefficients[1:(n_levels+n_externals+1), ]
          ilr.sum$residuals <- res.sum$residuals
          ilr.sum$sigma <- res.sum$sigma
          ilr.sum$r.squared <- res.sum$r.squared
          ilr.sum$adj.r.squared <- res.sum$adj.r.squared
          ilr.sum$fstatistic <- res.sum$fstatistic
        }
        else {
          ilr.sum$coefficients[j + n_levels + n_externals, ] <- res.sum$coefficients[(n_levels+n_externals+1), ]

        }
      }
    }
    if (is.null(factor_column) & is.null(external)){
      d <- data.frame(y = y, X = X)
      lmcla <- lm(y ~ ., data = d)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:ncol(X)) {
        Zj <- pivotCoord(cbind(X[, j], X[, -j]), method = method_pivot)
        dj <- data.frame(y = y, Z = Zj)
        res <- lm(y ~ ., data = dj)
        res.sum <- summary(res)
        if (j == 1) {
          ilr.sum$coefficients[1:2, ] <- res.sum$coefficients[1:2, ]
          ilr.sum$residuals <- res.sum$residuals
          ilr.sum$sigma <- res.sum$sigma
          ilr.sum$r.squared <- res.sum$r.squared
          ilr.sum$adj.r.squared <- res.sum$adj.r.squared
          ilr.sum$fstatistic <- res.sum$fstatistic
        }
        else {
          ilr.sum$coefficients[j + 1, ] <- res.sum$coefficients[2, ]
        }
      }
    }
    if (!is.null(factor_column) & is.null(external)){
      X_selected <- X %>% select(-c(factor_column))
      ZV <- data.frame(Factor = factor_col, X = X_selected)
      d <- data.frame(y = y, X = ZV)
      lmcla <- lm(y ~ ., data = d)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:(ncol(X)-1)) {
        X_prep <- X %>% select(-c(factor_column))
        Zj <- pivotCoord(cbind(X_prep[, j], X_prep[, -j]), method = method_pivot)
        ZVj <- data.frame(Factor = factor_col, Z = Zj)
        dj <- data.frame(y = y, Z = ZVj)
        res <- lm(y ~ ., data = dj)
        res.sum <- summary(res)
        if (j == 1) {
          ilr.sum$coefficients[1:(n_levels+1), ] <- res.sum$coefficients[1:(n_levels+1), ]
          ilr.sum$residuals <- res.sum$residuals
          ilr.sum$sigma <- res.sum$sigma
          ilr.sum$r.squared <- res.sum$r.squared
          ilr.sum$adj.r.squared <- res.sum$adj.r.squared
          ilr.sum$fstatistic <- res.sum$fstatistic
        }
        else {
          ilr.sum$coefficients[j + n_levels, ] <- res.sum$coefficients[(n_levels+1), ]
        }
      }
    }
    if (is.null(factor_column) & !is.null(external)){
      X_selected <- X %>% select(-c(external))
      ZV <- data.frame(Externals = external_col, X = X_selected)
      d <- data.frame(y = y, X = ZV)
      lmcla <- lm(y ~ ., data = d)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:(ncol(X)-n_externals)) {
        X_prep <- X %>% select(-c(external))
        Zj <- pivotCoord(cbind(X_prep[, j], X_prep[, -j]), method = method_pivot)
        ZVj <- data.frame(Externals = external_col, Z = Zj)
        dj <- data.frame(y = y, Z = ZVj)
        res <- lm(y ~ ., data = dj)
        res.sum <- summary(res)
        if (j == 1) {
          ilr.sum$coefficients[1:(n_externals+2), ] <- res.sum$coefficients[1:(n_externals+2), ]
          ilr.sum$residuals <- res.sum$residuals
          ilr.sum$sigma <- res.sum$sigma
          ilr.sum$r.squared <- res.sum$r.squared
          ilr.sum$adj.r.squared <- res.sum$adj.r.squared
          ilr.sum$fstatistic <- res.sum$fstatistic
        }
        else {
          ilr.sum$coefficients[j + n_externals + 1, ] <- res.sum$coefficients[(n_externals+2), ]
        }
      }
    }
    list(lm = lmcla, lm = lmcla.sum, ilr = ilr.sum)
  }
  robilrregression <- function(X, y, external, factor_column, method_pivot) {

    if (!is.null(factor_column) & !is.null(external)){
      X_selected <- X %>% select(-c(external, factor_column))  ## maybe we can work with relocate from the dyplr Package
      ZV <- data.frame(Factor = factor_col, Externals = external_col, X = X_selected)  ## maybe we can work with relocate from the dyplr Package
      d <- data.frame(y = y, X = ZV)    ## Double indexing Main ELements with X prefix (X.X)
      lmcla <- robustbase::lmrob(y ~ ., data = d)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:(ncol(X)-1-n_externals)) {
        X_prep <- X %>% select(-c(external, factor_column))
        Zj <- pivotCoord(cbind(X_prep[, j], X_prep[, -j]), method = method_pivot)
        ZVj <- data.frame(Factor = factor_col, Externals = external_col, Z = Zj)
        dj <- data.frame(y = y, Z = ZVj)
        res <- robustbase::lmrob(y ~ ., data = dj)
        res.sum <- summary(res)
        if (j == 1) {
          ilr.sum$coefficients[1:(n_levels+n_externals+1), ] <- res.sum$coefficients[1:(n_levels+n_externals+1), ]
          ilr.sum$residuals <- res.sum$residuals
          ilr.sum$sigma <- res.sum$sigma
          ilr.sum$r.squared <- res.sum$r.squared
          ilr.sum$adj.r.squared <- res.sum$adj.r.squared
        }
        else {
          ilr.sum$coefficients[j + n_levels + n_externals, ] <- res.sum$coefficients[(n_levels+n_externals+1), ]
        }
      }
    }
    if (is.null(factor_column) & is.null(external)){
      d <- data.frame(y = y, X = X)
      lmcla <- robustbase::lmrob(y ~ ., data = d)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:ncol(X)) {
        Zj <- pivotCoord(cbind(X[, j], X[, -j]), method = method_pivot)
        dj <- data.frame(y = y, Z = Zj)
        res <- robustbase::lmrob(y ~ ., data = dj)
        res.sum <- summary(res)
        if (j == 1) {
          ilr.sum$coefficients[1:2, ] <- res.sum$coefficients[1:2, ]
          ilr.sum$residuals <- res.sum$residuals
          ilr.sum$sigma <- res.sum$sigma
          ilr.sum$r.squared <- res.sum$r.squared
          ilr.sum$adj.r.squared <- res.sum$adj.r.squared
        }
        else {
          ilr.sum$coefficients[j + 1, ] <- res.sum$coefficients[2, ]
        }
      }
    }
    if (!is.null(factor_column) & is.null(external)){
      X_selected <- X %>% select(-c(factor_column))
      ZV <- data.frame(Factor = factor_col, X = X_selected)
      d <- data.frame(y = y, X = ZV)
      lmcla <- robustbase::lmrob(y ~ ., data = d)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:(ncol(X)-1)) {
        X_prep <- X %>% select(-c(factor_column))
        Zj <- pivotCoord(cbind(X_prep[, j], X_prep[, -j]), method = method_pivot)
        ZVj <- data.frame(Factor = factor_col, Z = Zj)
        dj <- data.frame(y = y, Z = ZVj)
        res <- robustbase::lmrob(y ~ ., data = dj)
        res.sum <- summary(res)
        if (j == 1) {
          ilr.sum$coefficients[1:(n_levels+1), ] <- res.sum$coefficients[1:(n_levels+1), ]
          ilr.sum$residuals <- res.sum$residuals
          ilr.sum$sigma <- res.sum$sigma
          ilr.sum$r.squared <- res.sum$r.squared
          ilr.sum$adj.r.squared <- res.sum$adj.r.squared
        }
        else {
          ilr.sum$coefficients[j + n_levels, ] <- res.sum$coefficients[n_levels + 1, ]
        }
      }
    }
    if (is.null(factor_column) & !is.null(external)){
      X_selected <- X %>% select(-c(external))
      ZV <- data.frame(Externals = external_col, X = X_selected)
      d <- data.frame(y = y, X = ZV)
      lmcla <- robustbase::lmrob(y ~ ., data = d)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:(ncol(X)-n_externals)) {
        X_prep <- X %>% select(-c(external))
        Zj <- pivotCoord(cbind(X_prep[, j], X_prep[, -j]), method = method_pivot)
        ZVj <- data.frame(Externals = external_col, Z = Zj)
        dj <- data.frame(y = y, Z = ZVj)
        res <- robustbase::lmrob(y ~ ., data = dj)
        res.sum <- summary(res)
        if (j == 1) {
          ilr.sum$coefficients[1:(n_externals+2), ] <- res.sum$coefficients[1:(n_externals+2), ]
          ilr.sum$residuals <- res.sum$residuals
          ilr.sum$sigma <- res.sum$sigma
          ilr.sum$r.squared <- res.sum$r.squared
          ilr.sum$adj.r.squared <- res.sum$adj.r.squared
        }
        else {
          ilr.sum$coefficients[j + n_externals + 1, ] <- res.sum$coefficients[(n_externals+2), ]
        }
      }
    }
    list(lm = lmcla, lm = lmcla.sum, ilr = ilr.sum)
  }
  if (method == "classical") {
    reg <- ilrregression(X, y, external, factor_column, method_pivot)
  }
  else if (method == "robust") {
    reg <- robilrregression(X, y, external, factor_column, method_pivot)
  }
  return(reg)
}
