#' Calculate the Laplacian Matrix
#'
#' @export
#' @description This function computes the Laplacian matrix for a given graph structure matrix. 
#' 
#' @param G A square matrix representing the graph structure, where each entry \code{G[i, j]} 
#' indicates the weight of the edge between node \code{i} and node \code{j}. If \code{G[i, j]} is 
#' zero, it indicates that there is no edge between the nodes. Diagonal entries are ignored.
#' 
#' @return A square matrix representing the Laplacian matrix \code{L}. The Laplacian matrix is 
#' calculated as \code{L = D - A}, where \code{D} is the degree matrix and \code{A} is the adjacency matrix.
#' 
#' @details The adjacency matrix \code{A} is derived from the input graph structure matrix \code{G} by 
#' setting \code{A[i, j] = 1} if \code{G[i, j] != 0} and \code{A[i, j] = 0} otherwise. The degree matrix \code{D} 
#' is a diagonal matrix where each diagonal element \code{D[i, i]} is the sum of the corresponding row 
#' in the adjacency matrix \code{A}.
#' 
#' @examples
#' set.seed(0)
#' data <- simu_data(200, 20, 9)
#' G <- data$G
#' lapmat <- Laplacian(G)
#' print(lapmat)
#' 
#' @seealso \code{\link{simu_data}}
#'
#' @importFrom stats setNames
Laplacian <- function(G) {
    p <- ncol(G)
    A <- matrix(0, p, p)
    A <- ifelse(G != 0, 1, 0)
    diag(A) <- 0
    D <- diag(apply(abs(A), 1, sum))
    D - A
}
