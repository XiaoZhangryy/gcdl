#' Generate the Simulation Data
#'
#' @export
#' @description Generates simulation data based on a specified covariance structure. This function creates a design matrix, a response vector, and a graph structure matrix. 
#' 
#' @param n Integer. The sample size.
#' @param p Integer. The dimension of the covariates.
#' @param s0 Integer. The cardinality of the active set. Must be a multiple of 3.
#' @param Cov A \code{p} by \code{p} positive definite covariance matrix. Default is \code{NULL}. If not provided, a default covariance matrix with a specific structure is generated.
#' 
#' @importFrom mnormt rmnorm
#' @importFrom stats rnorm
#' @return A list containing:
#' \itemize{
#'   \item \code{x} - The design matrix, where each row is an observation vector.
#'   \item \code{y} - The response vector.
#'   \item \code{G} - The graph structure matrix.
#' }
#' 
#' @examples
#' set.seed(0)
#' data <- simu_data(200, 20, 9)
#' x <- data$x
#' y <- data$y
#' G <- data$G
simu_data <- function(n, p, s0, Cov = NULL) {
    # Check if s0 is a multiple of 3
    if (s0 %% 3 != 0) {
        stop(paste("The program stops running because", s0, "is not a multiple of 3."))
    }
    
    # Function to create a new covariance matrix
    Covariance_new <- function(p, s0, rho) {
        cluster <- s0 / 3
        cov <- matrix(0, 3, 3)
        COV <- matrix(0, p, p)
        for (i in 1:2) {
            for (j in (i + 1):3) {
                cov[i, j] <- rho
            }
        }
        Cov <- kronecker(diag(1, cluster), cov + t(cov))
        COV[1:s0, 1:s0] <- Cov
        diag(COV) <- rep(1, p)
        COV
    }

    # Function to generate data based on the covariance matrix
    Data <- function(n, p, s0, Cov) {
        x <- rmnorm(n, rep(0, p), Cov)
        beta <- rep(0, p)
        beta[1:s0] <- 3
        e <- rnorm(n)
        y <- x %*% beta + e
        list(x = x, y = y, beta = beta)
    }

    # If Cov is not provided, generate a default covariance matrix
    if (is.null(Cov)) Cov <- Covariance_new(p, s0, rho = 0.9)
    Theta <- solve(Cov)
    G <- ifelse(Theta != 0, 1, 0)
    da <- Data(n, p, s0, Cov)
    x <- da$x
    y <- da$y
    list(x = x, y = y, G = G)
}
