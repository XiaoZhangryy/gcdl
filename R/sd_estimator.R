#' Estimate Standard Deviation
#'
#' @export
#' @description Utilizes the RCV method to estimate standard deviation. This method is particularly useful in high-dimensional settings where traditional variance estimation methods may not be effective. For more details on the RCV method, see Fan et al. (2012).
#' 
#' @param x The design matrix, which is a \code{n} by \code{p} matrix, where \code{n} is the number of observations and \code{p} is the number of predictors.
#' @param y The response vector, with \code{n} elements corresponding to the observations in the design matrix.
#' @param L Graph Laplacian matrix, which is used to incorporate the graph structure into the estimation process.
#' @param nfolds The number of cross-validation folds. Default is 10.
#' 
#' @importFrom stats coef lm
#' @return The estimated standard deviation.
#' 
#' @references
#' Fan, J., Guo, S., & Hao, N. (2012). Variance estimation using refitted cross-validation in ultrahigh dimensional regression. Journal of the Royal Statistical Society Series B: Statistical Methodology, 74(1), 37-65.
#' 
#' @examples
#' set.seed(0)
#' data <- simu_data(200, 20, 9)
#' x <- data$x
#' y <- data$y
#' G <- data$G
#' lapmat <- Laplacian(G)
#' sd <- sd_estimator(x, y, lapmat, 10)
#' print(sd)
sd_estimator <- function(x, y, L, nfolds = 10) {
    n <- nrow(x)
    p <- ncol(x)
    ind <- sample(1:n, n / 2, replace = FALSE)
    x1 <- x[ind, ]
    y1 <- y[ind]
    x2 <- x[-ind, ]
    y2 <- y[-ind]
    
    # Function to fit model and compute sigmahat
    compute_sigmahat <- function(x_train, y_train, x_test, y_test, L, nfolds) {
        cv.fit <- cv.glmgraph(x_train, y_train, L,
                              penalty = "lasso", nfolds = nfolds,
                              lambda2 = c(0, 0.01 * 2^(seq(-7, 7, by = 0.5)))
        )
        selected_indices <- which(as.vector(coef(cv.fit, s = "lambda1.1se")[-1]) != 0)
        if (length(selected_indices) == 0) {
            return(sum(y_test^2) / (n / 2 - 1))
        } else {
            fit_lm <- lm(y_test ~ x_test[, selected_indices])
            return(sum((fit_lm$residuals)^2) / (n / 2 - 1 - length(selected_indices)))
        }
    }
    
    sigmahat1 <- compute_sigmahat(x1, y1, x2, y2, L, nfolds)
    sigmahat2 <- compute_sigmahat(x2, y2, x1, y1, L, nfolds)
    
    return(sqrt((sigmahat1 + sigmahat2) / 2))
}
