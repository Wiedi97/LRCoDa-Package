lmCoDaX_Refactoring <- function (y, X, external = NULL, factor_column = NULL, method = "robust", method_pivot = 'orthonormal') { # ltsReg mit lmrob ersetzen und dann sollte die Fehlermeldung verschwinden
  
  if (!is.null(external)) {
    external_col <- X %>% select(external)
    n_externals <- length(external_col)  
  }
  if (!is.null(factor_column)) {
    factor_var <- X %>% select(factor_column)
    factor_col <- factor_var[, 1]
    n_levels <- length(levels(factor_col))  
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
      lmcla <- robustbase::ltsReg(y ~ ., data = d)  
      lmcla.sum <- summary(lmcla)  
      ilr.sum <- lmcla.sum
      
      for (j in 1:(ncol(X)-1-n_externals)) { 
        X_prep <- X %>% select(-c(external, factor_column))   
        Zj <- pivotCoord(cbind(X_prep[, j], X_prep[, -j]), method = method_pivot)    
        ZVj <- data.frame(Factor = factor_col, Externals = external_col, Z = Zj)  
        dj <- data.frame(y = y, Z = ZVj)
        res <- robustbase::ltsReg(y ~ ., data = dj)
        res.sum <- summary(res)
        if (j == 1) {
          ilr.sum$coefficients[1:(n_levels+n_externals+1), ] <- res.sum$coefficients[1:(n_levels+n_externals+1), ]  
          ilr.sum$residuals <- res.sum$residuals
          ilr.sum$sigma <- res.sum$sigma
          ilr.sum$r.squared <- res.sum$r.squared
          ilr.sum$adj.r.squared <- res.sum$adj.r.squared
          ilr.sum$fstatistic <- res.sum$fstatistic    ## bei lmrob gibt es keine F-Statistik mehr weil es ein Iterationsalgorithmus ist
        }
        else {
          ilr.sum$coefficients[j + n_levels + n_externals, ] <- res.sum$coefficients[(n_levels+n_externals+1), ]  
        }
      }
    }
    if (is.null(factor_column) & is.null(external)){
      d <- data.frame(y = y, X = X)    
      lmcla <- robustbase::ltsReg(y ~ ., data = d)  
      lmcla.sum <- summary(lmcla)  
      ilr.sum <- lmcla.sum 
      
      for (j in 1:ncol(X)) {
        Zj <- pivotCoord(cbind(X[, j], X[, -j]), method = method_pivot)
        dj <- data.frame(y = y, Z = Zj)
        res <- robustbase::ltsReg(y ~ ., data = dj)
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
      lmcla <- robustbase::ltsReg(y ~ ., data = d)  
      lmcla.sum <- summary(lmcla)  
      ilr.sum <- lmcla.sum
      
      for (j in 1:(ncol(X)-1)) {
        X_prep <- X %>% select(-c(factor_column)) 
        Zj <- pivotCoord(cbind(X_prep[, j], X_prep[, -j]), method = method_pivot)
        ZVj <- data.frame(Factor = factor_col, Z = Zj)  
        dj <- data.frame(y = y, Z = ZVj)  
        res <- robustbase::ltsReg(y ~ ., data = dj)
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
          ilr.sum$coefficients[j + n_levels, ] <- res.sum$coefficients[n_levels + 1, ]  
        }
      }
    }
    if (is.null(factor_column) & !is.null(external)){
      X_selected <- X %>% select(-c(external)) 
      ZV <- data.frame(Externals = external_col, X = X_selected)
      d <- data.frame(y = y, X = ZV)    
      lmcla <- robustbase::ltsReg(y ~ ., data = d)  
      lmcla.sum <- summary(lmcla)  
      ilr.sum <- lmcla.sum
      
      for (j in 1:(ncol(X)-n_externals)) {
        X_prep <- X %>% select(-c(external)) 
        Zj <- pivotCoord(cbind(X_prep[, j], X_prep[, -j]), method = method_pivot)
        ZVj <- data.frame(Externals = external_col, Z = Zj)  
        dj <- data.frame(y = y, Z = ZVj)  
        res <- robustbase::ltsReg(y ~ ., data = dj)
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
  if (method == "classical") {
    reg <- ilrregression(X, y, external, factor_column, method_pivot)
  }
  else if (method == "robust") {
    reg <- robilrregression(X, y, external, factor_column, method_pivot)
  }
  return(reg)
}  
