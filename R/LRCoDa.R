#' @importFrom stats lm
#' @import dplyr
#' @import tidyr
#' @importFrom robustbase lmrob
#' @importFrom robustbase lmrob.control
#' @importFrom robCompositions pivotCoord
#' @importFrom sjmisc is_empty
NULL


#' Classical and robust regression of non-compositional (real) response on compositional and non-compositional predictors combined
#'
#' Delivers appropriate inference for regression of y on compositional and non-compositional combined input X.
#'
#' Compositional explanatory variables should not be directly used in a linear regression model because any inference statistic
#' can become misleading. While various approaches for this problem were proposed, here an approach based on the pivot coordinates
#' is used. Further these compositional explanatory variables can be supplemented with external non-compositional data
#' and factor variables.
#'
#' @aliases LRCoDa ilrregression robilrregression
#' @param y The response which should be non-compositional
#' @param X The compositional and/or non-compositional predictors as a matrix, data.frame or numeric vector
#' @param external Specify the columns which are not part of the composition and not factors
#' @param factor_column Specify the column which includes factor levels
#' @param method If \dQuote{robust}, the fast MM-type linear regression is applied, while with method
#' \dQuote{classical}, the conventional least squares regression is applied.
#' @param pivot_norm if FALSE then the normalizing constant is not used, if TRUE sqrt((D-i)/(D-i+1))
#' is used (default). The user can also specify a self-defined constant
#' @param max_refinement_steps (for the fast-S algorithm): maximal number of refinement
#' steps for the “fully” iterated best candidates.
#' @return An object of class \sQuote{lmrob} or \sQuote{lm} and two summary
#' objects.
#' @author Roman Wiedemeier
#' @seealso \code{\link{lm}} and \code{\link{lmrob}}
#' @keywords models compositional data linear regression
#' @export
#' @examples
#'
#' ## How the content of sand of the agricultural
#' ## and grazing land soils in Germany depend on
#' ## relative contributions of the main chemical trace elements,
#' ## their different soil types and the Annual mean temperature:
#' data("gemas")
#' gemas$COUNTRY <- as.factor(gemas$COUNTRY)
#' gemas_GER <- dplyr::filter(gemas, gemas$COUNTRY == 'GER')
#' y <- gemas_GER$sand
#' X <- dplyr::select(gemas_GER, c(MeanTemp, soilclass, Al:Zr))
#' LRCoDa(y, X, external = c('MeanTemp'), factor_column = 'soilclass',
#' method='classical', pivot_norm = 'orthonormal')
#' LRCoDa(y, X, external = c('MeanTemp'), factor_column = 'soilclass',
#' method='robust', pivot_norm = 'orthonormal')
LRCoDa <- function (y, X, external = NULL, factor_column = NULL, method = "robust", pivot_norm = 'orthonormal', max_refinement_steps = 200) { # ltsReg mit lmrob ersetzen und dann sollte die Fehlermeldung verschwinden

  if (!is.null(external) & (typeof(external) != "character")) {
    stop("Invalid datatype for external")
  }
  if (!is.null(factor_column) & (typeof(factor_column) != "character")) {
    stop("Invalid datatype for factor_column")
  }
  if (is.null(factor_column) & (any(sapply(X, function(x) !is.numeric(x))))) {
    stop("X contains non-numerical variables but there is no factor_column defined")
  }
  if (length(factor_column) > 1) {
    stop("There are more than 1 factor variable defined")
  }
  if (any(is.na(y))){
    dat <- cbind(y, X)
    dat_missing <- dat %>% filter(is.na(y))
    n <- dim(dat_missing)[1]
    dat_new <- dat %>% filter(!is.na(y))

    X <- dat_new %>% select(-c(y))
    y <- dat_new %>% select(c(y))
  }
  if (!is.null(external)){
    external_col <- X %>% select(all_of(external))
    if (all(sapply(external_col, function(x) is.numeric(x)))){
      n_externals <- length(external_col)
    } else {
      stop("Datatype of all 'external' variables have to be numeric")
    }
  }
  if (!is.null(factor_column)) {
    factor_var <- X %>% select(all_of(factor_column))
    factor_col <- factor_var[, 1]
    if (is.factor(factor_col) | is.character(factor_col)){
      n_levels <- length(unique(as.character(factor_col)))
    } else {
      stop("Datatype of 'factor_column' has to be factor or character")
    }
    if (any(is_empty(unique(as.character(factor_col)), first.only = FALSE, all.na.empty = TRUE) == TRUE)){
      stop("Dataset contains levels with empty strings or missing values. Specify factor name or drop these observations.")
    }
  }

  ilrregression <- function(X, y, external, factor_column, pivot_norm) {

    if (!is.null(factor_column) & !is.null(external)){
      X_selected <- X %>% select(-all_of(c(external, factor_column)))
      ZV <- data.frame(cbind(y, X %>% relocate(all_of(c(factor_column, external))) %>% rename_with(.cols = colnames(X %>% select(all_of(factor_column))), function(x){paste0("Factor.", x)}) %>%
                               rename_with(.cols = colnames(X %>% select(all_of(external))), function(x){paste0("External.", x)}) %>%
                               rename_with(.cols = colnames(X %>% select(-all_of(c(external, factor_column)))), function(x){paste0("Internal.", x)})))
      lmcla <- lm(y ~ ., data = ZV)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:(ncol(X)-1-n_externals)) {
        Zj <- pivotCoord(cbind(X_selected[, j], X_selected[, -j]), method = pivot_norm)
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
      ZV <- data.frame(cbind(y, X %>% rename_with(.cols = everything(), function(x){paste0("Internal.", x)})))
      lmcla <- lm(y ~ ., data = ZV)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:ncol(X)) {
        Zj <- pivotCoord(cbind(X[, j], X[, -j]), method = pivot_norm)
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
      X_selected <- X %>% select(-all_of(c(factor_column)))
      ZV <- data.frame(cbind(y, X %>% relocate(all_of(c(factor_column))) %>% rename_with(.cols = colnames(X %>% select(all_of(factor_column))), function(x){paste0("Factor.", x)}) %>%
                               rename_with(.cols = colnames(X %>% select(-all_of(c(factor_column)))), function(x){paste0("Internal.", x)})))
      lmcla <- lm(y ~ ., data = ZV)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:(ncol(X)-1)) {
        Zj <- pivotCoord(cbind(X_selected[, j], X_selected[, -j]), method = pivot_norm)
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
      X_selected <- X %>% select(-all_of(c(external)))
      ZV <- data.frame(cbind(y, X %>% relocate(all_of(c(external))) %>% rename_with(.cols = colnames(X %>% select(all_of(external))), function(x){paste0("External.", x)}) %>%
                               rename_with(.cols = colnames(X %>% select(-all_of(c(external)))), function(x){paste0("Internal.", x)})))
      lmcla <- lm(y ~ ., data = ZV)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:(ncol(X)-n_externals)) {
        Zj <- pivotCoord(cbind(X_selected[, j], X_selected[, -j]), method = pivot_norm)
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

  robilrregression <- function(X, y, external, factor_column, pivot_norm) {
    cont_lmrob <- lmrob.control(fast.s.large.n = Inf, k.max = max_refinement_steps)

    if (!is.null(factor_column) & !is.null(external)){
      X_selected <- X %>% select(-all_of(c(external, factor_column)))   ## maybe we can work with relocate from the dyplr Package
      ZV <- data.frame(cbind(y, X %>% relocate(all_of(c(factor_column, external))) %>% rename_with(.cols = colnames(X %>% select(all_of(factor_column))), function(x){paste0("Factor.", x)}) %>%
                               rename_with(.cols = colnames(X %>% select(all_of(external))), function(x){paste0("External.", x)}) %>%
                               rename_with(.cols = colnames(X %>% select(-all_of(c(external, factor_column)))), function(x){paste0("Internal.", x)})))
      lmcla <- robustbase::lmrob(y ~ ., data = ZV, control = cont_lmrob)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:(ncol(X)-1-n_externals)) {
        Zj <- pivotCoord(cbind(X_selected[, j], X_selected[, -j]), method = pivot_norm)
        ZVj <- data.frame(Factor = factor_col, Externals = external_col, Z = Zj)
        dj <- data.frame(y = y, Z = ZVj)
        res <- robustbase::lmrob(y ~ ., data = dj, control = cont_lmrob)
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
      ZV <- data.frame(cbind(y, X %>% rename_with(.cols = everything(), function(x){paste0("Internal.", x)})))
      lmcla <- robustbase::lmrob(y ~ ., data = ZV, control = cont_lmrob)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:ncol(X)) {
        Zj <- pivotCoord(cbind(X[, j], X[, -j]), method = pivot_norm)
        dj <- data.frame(y = y, Z = Zj)
        res <- robustbase::lmrob(y ~ ., data = dj, control = cont_lmrob)
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
      X_selected <- X %>% select(-all_of(c(factor_column)))
      ZV <- data.frame(cbind(y, X %>% relocate(all_of(c(factor_column))) %>% rename_with(.cols = colnames(X %>% select(all_of(factor_column))), function(x){paste0("Factor.", x)}) %>%
                               rename_with(.cols = colnames(X %>% select(-all_of(c(factor_column)))), function(x){paste0("Internal.", x)})))
      lmcla <- robustbase::lmrob(y ~ ., data = ZV, control = cont_lmrob)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:(ncol(X)-1)) {
        Zj <- pivotCoord(cbind(X_selected[, j], X_selected[, -j]), method = pivot_norm)
        ZVj <- data.frame(Factor = factor_col, Z = Zj)
        dj <- data.frame(y = y, Z = ZVj)
        res <- robustbase::lmrob(y ~ ., data = dj, control = cont_lmrob)
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
      X_selected <- X %>% select(-all_of(c(external)))
      ZV <- data.frame(cbind(y, X %>% relocate(all_of(c(external))) %>% rename_with(.cols = colnames(X %>% select(all_of(external))), function(x){paste0("External.", x)}) %>%
                               rename_with(.cols = colnames(X %>% select(-all_of(c(external)))), function(x){paste0("Internal.", x)})))
      lmcla <- robustbase::lmrob(y ~ ., data = ZV, control = cont_lmrob)
      lmcla.sum <- summary(lmcla)
      ilr.sum <- lmcla.sum

      for (j in 1:(ncol(X)-n_externals)) {
        Zj <- pivotCoord(cbind(X_selected[, j], X_selected[, -j]), method = pivot_norm)
        ZVj <- data.frame(Externals = external_col, Z = Zj)
        dj <- data.frame(y = y, Z = ZVj)
        res <- robustbase::lmrob(y ~ ., data = dj, control = cont_lmrob)
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
    reg <- ilrregression(X, y, external, factor_column, pivot_norm)
  }
  else if (method == "robust") {
    reg <- robilrregression(X, y, external, factor_column, pivot_norm)
  }
  if (exists("dat_missing")) {
    message("There are ", n ," observations omitted due to missings in the target variable")
  }
  return(reg)
}
