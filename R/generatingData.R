#' @title Generate data
#' @description Generates data using the mvrnorm function from the MASS package
#' according to user-specified means and effect sizes (ie. differences in means)
#' for variables.
#' @param num_samples_trt The number of samples in the group that has an effect
#' applied to the mean specified.
#' @param num_samples_control The number of samples in the group that maintains
#' the mean specified by the user.
#' @param means A vector specifying the mean for each feature in the dataset.
#' @param delta A vector specifying the difference in the mean for each feature
#' in the dataset.
#' @param cor_matrix The correlation matrix. If not specified an appropriately
#' sized identity matrix is generated. Default: NULL
#' @param num_reps Number of repetitions of the data that will be generated.
#' Default: 1000
#' @return A dataframe with (num_samples_trt+num_samples_control)*num_reps rows
#' and length(means) columns as well as a column for sample number, treatment
#' assignment, iteration number, and total number of samples n.
#' @examples
#' generate_data(50, 50, c(1, 2, 3), c(0.3, 0, 0.7))
#'
#' @rdname generate_data
#' @export
generate_data <- function(num_samples_trt, num_samples_control, means,
    delta, cor_matrix = NULL, num_reps = 1000) {
    # check cor_mat if null and if so generate an identity matrix
    if (is.null(cor_matrix)) {
        cor_matrix <- create_cor_matrix(length(means))
    }
    # check to make sure specified length for means is equal to
    # specified length for deltas
    if (length(means) != length(delta))
        stop("Error: Mean and delta vectors must be equal length.")

    all_reps_data <- data.frame
    n <- num_samples_trt + num_samples_control

    Mu_single <- matrix(nr = n, nc = length(means), rep(means, each = n))
    Mu_shiftermatrix <- matrix(nr = num_samples_trt, nc = length(delta),
        rep(delta, each = num_samples_trt))
    Zero_matrix <- matrix(nr = num_samples_control, nc = length(delta),
        rep(delta * 0, each = num_samples_control))
    Mu_shiftermatrix <- rbind(Mu_shiftermatrix, Zero_matrix)
    Mu <- Mu_single + Mu_shiftermatrix

    for (i in 1:num_reps) {
        #the means and effects are specified by the user while the errors are generated
        Epsilon <- mvrnorm(n, mu = rep(0, length(means)), Sigma = cor_matrix)
        Y <- Mu + Epsilon
        colnames(Y) <- colnames(Y, do.NULL = FALSE, prefix = "Gene_")

        treatment <- c(rep(1, num_samples_trt), rep(0, num_samples_control))
        reps <- c(rep(i, n))
        sample <- c(seq(1, n, 1))

        Y <- cbind(sample, Y, treatment, reps, n)

        if (i == 1) {
            all_reps_data <- Y
        } else {
            all_reps_data <- rbind(all_reps_data, Y)
        }
    }
    all_reps_data <- as.data.frame(all_reps_data)
    return(all_reps_data)
}

#' @title Create Custom Correlation Matrix
#' @description Creates a correlation matrix with a specified number of features
#' with strong, medium, and weak correlations. Background correlation can be set
#' by user with background= ; if no specifications given an identity matrix is
#' returned.
#' @param d A number. The dimensions of the matrix.
#' @param strong Number of strong correlations in the matrix. Default: 0
#' @param med Number of moderate correlations in the matrix. Default: 0
#' @param weak Number of weak correlations in the matrix. Default: 0
#' @param none Number of uncorrelated features. Default: 0
#' @param background Background correlation amongst all features. Default: 0
#' @return A matrix with the specified correlation structure.
#' @details Correlation matrix is built from upper left corner starting with
#' the first correlation strength specified and working iteratively through
#' correlation strength specifications along the diagonal.
#' @examples
#' create_cor_matrix(5)
#' create_cor_matrix(10, strong=3, weak=5)
#' create_core_matrix(20, strong=5, med=7, weak=5, background=0.1)
#' @rdname create_cor_matrix
#' @export
create_cor_matrix <- function(d, strong = 0, med = 0, weak = 0, none = 0,
    background = 0) {
    # future additions random parameter = TRUE then keep range, =FALSE then use
    # ranges;( Rob wants to fix )
    total <- strong + med + weak + none

    if (total > d)
        stop("Error: There are not enough dimensions in the matrix
         to represent that many correlations.")

    if (strong == 0 & med == 0 & weak == 0) {
        matrix <- diag(d)
    } else {

        if (strong != 0) {
            # generate enough correlation numbers to fill half a subblock
            strong_submatrix <- generate_submatrix(strong, 0.6, 0.9)
        } else strong_submatrix <- NULL

        if (med != 0) {
            med_submatrix <- generate_submatrix(med, 0.4, 0.6)
        } else med_submatrix <- NULL

        if (weak != 0) {
            weak_submatrix <- generate_submatrix(weak, 0.1, 0.4)
        } else weak_submatrix <- NULL

        if (d - total != 0 | none != 0) {
            none_submatrix <- diag(d - strong - med - weak)
        } else none_submatrix <- NULL

        matrix <- combine_matrices(strong_submatrix, med_submatrix,
            weak_submatrix, none_submatrix)
    }
    if (background != 0) {
        copy <- matrix
        matrix[copy == 0] <- background
        # fill it in with 'background'
    }


    if (!(is.positive.definite(matrix)))
        matrix <- make.positive.definite(matrix)
    return(matrix)
}

#' @title Generates submatrix with a given correlation structure.
#' @description Called from create_cor_matrix; Uses runif to generate random
#' numbers between a lower and upper bound (as specified in create_cor_matrix
#' depending on user's strength indication - strong, med, weak) for the strength
#' of the correlation.
#' @param dim Dimensions of the submatrix
#' @param lower_bound Lower boundary for generated random numbers
#' @param upper_bound Upper boundary for generated random numbers
#' @return A matrix with a strong, medium, or weak correlation structure.
#' @details Designed only to be called by the create_cor_matrix function.
#' @rdname generate_submatrix

generate_submatrix <- function(dim, lower_bound, upper_bound) {
    corr <- round(runif(round((dim * dim - dim)/2), lower_bound, upper_bound),
        2)
    submatrix <- diag(dim)
    submatrix[lower.tri(submatrix)] <- corr
    submatrix <- t(submatrix)
    submatrix[lower.tri(submatrix)] <- corr
    if (!(is.positive.definite(submatrix)))
        submatrix <- make.positive.definite(submatrix)

    return(submatrix)
}

#' @title Combines multiple smaller matrices along the diagonal.
#' @description Only called from create_cor_matrix.
#' @param s A Matrix with a strong correlation structure.
#' @param m A Matrix with a moderate correlation structure.
#' @param w A Matrix with a weak correlation structure.
#' @param n A Matrix with no correlated features. An identity matrix.
#' @return A Matrix that is a combination of all submatrices above.
#' @rdname combine_matrices
combine_matrices <- function(s, m, w, n) {
    if (!is.null(s) & !is.null(m) & !is.null(w) & !is.null(n)) {
        in_progress <- as.matrix(bdiag(s, m))
        in_progress2 <- as.matrix(bdiag(w, n))
        matrix <- as.matrix(bdiag(in_progress, in_progress2))
    } else if (!is.null(s) & !is.null(m) & !is.null(w)) {
        in_progress <- as.matrix(bdiag(s, m))
        matrix <- as.matrix(bdiag(in_progress, w))
    } else if (!is.null(s) & !is.null(m) & !is.null(n)) {
        in_progress <- as.matrix(bdiag(s, m))
        matrix <- as.matrix(bdiag(in_progress, n))
    } else if (!is.null(s) & !is.null(w) & !is.null(n)) {
        in_progress <- as.matrix(bdiag(s, w))
        matrix <- as.matrix(bdiag(in_progress, n))
    } else if (!is.null(m) & !is.null(w) & !is.null(n)) {
        in_progress <- as.matrix(bdiag(m, w))
        matrix <- as.matrix(bdiag(in_progress, n))
    } else if (!is.null(s) & !is.null(m)) {
        matrix <- as.matrix(bdiag(s, m))
    } else if (!is.null(s) & !is.null(w)) {
        matrix <- as.matrix(bdiag(s, w))
    } else if (!is.null(s) & !is.null(n)) {
        matrix <- as.matrix(bdiag(s, n))
    } else if (!is.null(m) & !is.null(n)) {
        matrix <- as.matrix(bdiag(m, n))
    } else if (!is.null(w) & !is.null(n)) {
        matrix <- as.matrix(bdiag(m, n))
    } else if (!is.null(m) & !is.null(w)) {
        matrix <- as.matrix(bdiag(m, w))
    }

    return(matrix)
}
