#' Calculate Inverse Gram Matrix
#'
#' @export
#' @description Utilizes the graph structure to compute the approximate inverse of the Gram matrix \eqn{x^\top x/n}. This function performs nodewise regression using either ordinary least squares (OLS) or Lasso, depending on the degrees of freedom of the vertices.
#' 
#' @param x The design matrix, which is an \code{n} by \code{p} matrix, where \code{n} is the number of observations and \code{p} is the number of predictors.
#' @param G User-specified graph structure matrix. \code{G[i, j]} indicates the presence of an edge between nodes \code{i} and \code{j}.
#' @param k Integer. When the degrees of freedom of the vertex \code{j} are less than \code{k}, use ordinary least squares; otherwise, use Lasso. Default is \code{NULL}, which means \code{k} is set to \code{p}.
#' 
#' @importFrom glmnet cv.glmnet
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{inverse_Gram} - The approximate inverse of the matrix \eqn{x^\top x/n}.
#'   \item \code{Z} - The residuals of the nodewise regressions.
#' }
#' 
#' @examples
#' set.seed(0)
#' data <- simu_data(200, 20, 9)
#' x <- data$x
#' G <- data$G
#' inv_gram <- calculate_inv_gram(x, G)
#' print(inv_gram$inverse_Gram)
#' print(inv_gram$Z)
calculate_inv_gram <- function(x, G, k = NULL) {
    n <- nrow(x)
    p <- ncol(x)
    if (is.null(k)) k <- p
    Z <- matrix(0, n, p)
    diag(G) <- 0
    inverse_Gram <- matrix(0, p, p)
    for (j in 1:p) {
        D <- which(G[j, ] != 0)
        if (length(D) == 0) {
            inverse_Gram[, j] <- 0
            Z[, j] <- x[, j]
            inverse_Gram[j, j] <- n / sum(x[, j]^2)
        } else {
            if (length(D) >= 1 & length(D) <= k) {
                fit <- lm(x[, j] ~ x[, D] - 1)
                gamma <- fit$coefficients
                Z[, j] <- fit$residuals
                tao <- sum((fit$residuals)^2) / n
                inverse_Gram[j, j] <- 1 / tao
                inverse_Gram[, j][D] <- -gamma / tao
            }

            if (length(D) > k) {
                glmnetfit <- cv.glmnet(x[, D], x[, j])
                gamma <- as.vector(predict(glmnetfit, x[, D],
                    type = "coefficients",
                    s = "lambda.1se"
                ))[-1]
                prediction <- predict(glmnetfit, x[, D], s = "lambda.1se")
                Z[, j] <- x[, j] - prediction
                tao <- as.numeric((x[, j] %*%
                    (x[, j] - predict(glmnetfit, x[, D], s = "lambda.1se"))) / n)
                inverse_Gram[j, j] <- 1 / tao
                inverse_Gram[, j][D] <- -gamma / tao
            }
        }
    }
    list(inverse_Gram = inverse_Gram, Z = Z)
}
