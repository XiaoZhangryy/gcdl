#' Calculate P-values based on Graph-Constrained Desparsified LASSO (GCDL) Method
#'
#' @export
#' @description Calculates P-values based on the graph-constrained desparsified LASSO (GCDL) method. This method incorporates the graph structure into the desparsified LASSO estimator, providing more accurate variable selection in high-dimensional settings.
#' 
#' @param x The design matrix, which is an \code{n} by \code{p} matrix, where \code{n} is the number of observations and \code{p} is the number of predictors.
#' @param y The response vector, with \code{n} elements corresponding to the observations in the design matrix.
#' @param G User-specified graph structure matrix. \code{G[i, j]} indicates the presence of an edge between nodes \code{i} and \code{j}.
#' @param nfolds The number of cross-validation folds. Default is 10.
#' @param Centering Logical. Indicator of whether the design matrix should be centered to column zero mean. Default is \code{TRUE}.
#' 
#' @importFrom glmgraph cv.glmgraph
#' @importFrom stats coef predict pnorm
#' @references 
#' Chen, L., Liu, H., Kocher, J. P. A., Li, H., & Chen, J. (2015). glmgraph: an R package for variable selection and predictive modeling of structured genomic data. Bioinformatics, 31(24), 3991-3993.
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{P_value} - Individual p-values for each parameter.
#'   \item \code{bhat} - The GCDL estimator.
#'   \item \code{betahat} - Initial estimate.
#'   \item \code{sigmahat} - The estimation of standard deviation obtained through the RCV method.
#'   \item \code{Se_bhat} - Individual standard deviation for each parameter.
#'   \item \code{Inv_Gram} - The approximate inverse of the matrix \eqn{x^\top x/n}.
#' }
#' 
#' @examples
#' set.seed(0)
#' data <- simu_data(200, 20, 9)
#' x <- data$x
#' y <- data$y
#' G <- data$G
#' res <- gcdl(x, y, G)
#' print(res)
gcdl <- function(x, y, G, nfolds = 10, Centering = TRUE) {
    n <- nrow(x)
    p <- ncol(x)
    
    # Centering the design matrix and response vector if needed
    if (Centering) {
        y <- y - mean(y)
        x <- apply(x, MARGIN = 2, FUN = function(X) (X - mean(X)))
    }
    
    # Calculate Laplacian matrix
    L <- Laplacian(G)
    
    # Fit the model using cross-validated glmgraph
    cv.fit <- cv.glmgraph(x, y, L,
                          penalty = "lasso", nfolds = nfolds,
                          lambda2 = c(0, 0.01 * 2^(seq(-7, 7, by = 0.5)))
    )
    
    # Initial estimate of beta
    betahatLaplacian <- as.vector(coef(cv.fit)[-1])
    
    # Calculate residuals
    residual.vector.Laplacian <- y - predict(cv.fit, x, s = "lambda1.min", type = "response")
    
    # Estimate standard deviation using RCV method
    sigmahat.Lap <- sd_estimator(x, y, L, nfolds)
    
    # Calculate Gram matrix and its inverse
    Gram <- t(x) %*% x / n
    inv_Gram_temp <- calculate_inv_gram(x, G)
    inv_Gram <- inv_Gram_temp$inverse_Gram
    Z <- inv_Gram_temp$Z
    
    # Calculate GCDL estimator
    bhat.Lap <- betahatLaplacian + (t(Z) %*% residual.vector.Laplacian) / (diag(t(Z) %*% x))
    
    # Calculate standard error for each parameter
    seLap <- (sigmahat.Lap * sqrt(colSums(Z^2))) / diag(t(Z) %*% x)
    
    # Calculate test statistics and p-values
    testLap <- (as.vector(bhat.Lap)) / (seLap)
    p_value_Lap <- 2 * (1 - pnorm(abs(testLap)))
    
    list(
        bhat = c(bhat.Lap), 
        P_value = p_value_Lap, 
        betahat = as.vector(coef(cv.fit)),
        Se_bhat = seLap, 
        sigmahat = sigmahat.Lap, 
        Inv_Gram = inv_Gram
    )
}
