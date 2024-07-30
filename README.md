# gcdl
The GCDL package provides tools for variable selection and uncertainty quantification in high-dimensional linear models while incorporating network graphical information. Genes typically function in networks and are correlated due to their functional connectivity. Traditional variable selection methods often fail to account for this correlation, leading to less accurate selections and uncertainty quantifications.

# Installation

    #install.packages("devtools")
    library(devtools)
    install_github("XiaoZhangryy/gcdl")

# Usage

   - [x] [hdtest-manual.pdf](https://github.com/XiaoZhangryy/gcdl/blob/master/inst/gcdl-manual.pdf) ---------- Details of the usage of the package.
# Example
    library(gcdl)

    set.seed(0)
    data <- simu_data(200, 20, 9)
    x <- data$x
    y <- data$y
    G <- data$G
    res <- gcdl(x, y, G)
    print(res)

# References
 Tan, X., Zhang, X., Cui, Y., and Liu, X. (2024). Uncertainty quantification in high-dimensional linear models incorporating graphical structures with applications to gene set analysis. Manuscript.

# Development
This R package is developed by Xiangyong Tan and Xiao Zhang (zhangxiao1994@cuhk.edu.cn).
