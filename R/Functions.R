#' Set alpha for simulation.
#'
#' @param alpha The alpha value you wish to use for tests of significance.
#' @examples
#' set.alpha(0.05)
#' @export
set.alpha <- local(function(alpha) {
    global_alpha <<- alpha
})

#' Simulate multiple scenarios
#'
#' Simulates scenarios with varying effect sizes or sample sizes; can accomodate
#' differing numbers of samples per treatment and control groups
#'
#' @param list_effect_sizes A list of effect_size objects.
#' @param list_treatment_groups A list of sample sizes for the group exhibiting effects.
#' @param list_control_groups A list of sample sizes for the group that does not exhibit effects.
#' @param means A vector which contains the mean for each individual feature.
#' @param cor_matrix A matrix which specifies the correlation structure amongst all features.
#' @param num_reps A number that describes how many replications of each scenario there will be.
#' @param method choose between "ofaat", "mv_glm", "lasso"
#' @param p_adjust Used for the internal implementation of p.adjust; takes the same arguments
#' @seealso \code{\link{analyze}} for details on method parameter
#' @seealso \code{\link{p.adjust}} for p value adjustments
#'
#' @return A dataframe which summarizes the results of each scenario that has been run
#'
#'  @section Notes:
#' List arguments must be passed as lists, even if only 1 item in the list.
#' If differing lists are not passed for list_treatment_groups and list_control_groups
#' arguments then default equal treatment/control group sizes assigned. However, at least one size
#' list must be passed.
#' Method must be specified.
#'
#' @examples
#' e <- effect_size(c(0.7, 0.2, 0.4))
#' list_e <- c(e)
#' group_sizes <- c(20,40,60,80,100)
#' simulate_scenarios(list_e, group_sizes, method="ofaat")
#' @export
simulate_scenarios <- function(list_effect_sizes = NULL, list_treatment_groups = NULL,
    list_control_groups = NULL, means = NULL, cor_matrix = NULL, num_reps = 1000,
    method = NULL, p_adjust = NULL) {
    # stores the results of all scenarios
    results_data <- data.frame()

    if (is.null(list_treatment_groups) & is.null(list_control_groups)) {
        stop("Sorry, but you must specify how many samples you would like
         to run the simulation with.")
    } else if (!is.null(list_treatment_groups) & is.null(list_control_groups)) {
        for (e in list_effect_sizes) {
            for (num1 in list_treatment_groups) {
                num2 = num1
                results <- single_scenario_simulation(num1, num2,
                  means = means, delta = e, cor_matrix = cor_matrix,
                  method = method)

                results_data <- rbind(results_data, results$overall_summary_table)
            }
        }
    } else {
        for (e in list_effect_sizes) {
            for (num1 in list_treatment_groups) {
                for (num2 in list_control_groups) {
                  results <- single_scenario_simulation(num1, num2,
                    means = means, delta = e, cor_matrix = cor_matrix,
                    method = method)

                  results_data <- rbind(results_data, results$overall_summary_table)
                }
            }
        }
    }
    return(results_data)
}

#' Simulate single scenario
#'
#' Simulates scenario with specified number of samples in treatment and control groups,
#' specified number of features OR effect size OR mean, a correlation structure between features, method
#' used for analysis, p value adjustment used in analysis, and the number of repetitions
#' that should be simulated for this scenario.
#'
#' @param num_samples_trt A number; specifies how many samples in group exhibiting effects.
#' @param num_samples_control A number; specifies how many samples in groups without effects.
#' @param num_genes A number; specifies how many features are being analyzed in this scenario.
#' @param means A vector which contains the mean for each individual feature.
#' @param delta An effect_size object containing the observed effect for each feature.
#' @param cor_matrix A matrix which specifies the correlation structure amongst all features.
#' @param num_reps A number that describes how many replications of each scenario there will be.
#' @param method choose between "ofaat", "mv_glm", "lasso"
#' @param p_adjust Used for the internal implementation of p.adjust; takes the same arguments
#' @seealso \code{\link{analyze}} for details on method parameter
#' @seealso \code{\link{p.adjust}} for p value adjustments
#' @return A list containing a summary table of overall simulation results,
#' an individual gene summary result table (with specified effect size and proportion),
#' and the total number of repetitions that the simulation carried out for this scenario.
#'
#' @section Notes:
#' One of num_samples_trt or num_samples_control must be specified.
#' One of num_genes, means, or delta must be specified.
#' Method must be specified.
#'
#' @section Warning:
#' global_alpha must be specified with set.alpha(alpha) for simulation to run to
#' completion.
#' @examples
#' single_scenario_simulation(num_samples_trt = 80, num_samples_control= 95, num_genes= 6, method="ofaat")
#' @export
single_scenario_simulation <- function(num_samples_trt = NULL, num_samples_control = NULL,
    num_genes = NULL, means = NULL, delta = NULL, cor_matrix = NULL,
    num_reps = 1000, method = NULL, p_adjust = NULL) {
    if (is.null(num_samples_trt) & is.null(num_samples_control)) {
        stop("Number of samples in simulation must be specified.")
    } else if (is.null(num_samples_trt)) {
        num_samples_trt <- num_samples_control
    } else if (is.null(num_samples_control)) {
        num_samples_control <- num_samples_trt
    }

    if (is.null(means) & is.null(delta) & is.null(num_genes)) {
        stop("Number of features in the simulation must be specified.")
    } else if (is.null(means) & is.null(delta) & !is.null(num_genes)) {
        means <- c(rep(0, num_genes))
        delta <- effect_size(rep(0, num_genes))
    } else if (is.null(means) & !is.null(delta)) {
        means <- c(rep(0, length(delta@effects)))
    } else if (!is.null(means) & is.null(delta)) {
        delta <- effect_size(rep(0, length(means)))
    }
    generated_data <- generate_data(num_samples_trt, num_samples_control,
        means, delta@effects, cor_matrix, num_reps = num_reps)

    # options for analyze function method field are 'ofaat', 'mv_glm',
    # and 'lasso' ofaat and mv_glm options support p value adjustments
    # see ?p.adjust for available adjustment options


    if (method == "mv_glm" & length(means) > num_samples_trt + num_samples_control) {
        stop("Warning: The analysis method specified is not reliable for a dataset of this shape.
           Please rerun simulation with lasso method.")
    }
    analyzed_data <- analyze(generated_data, method = method)

    # summarize information in all the reps; returns a list with the
    # summary table for the overall model and the summary for each
    # gene
    summarized_data <- summary(analyzed_data, delta, analysis_method = method)

    # add the number of iterations to summary list bind the specified
    # effect sizes to the summary for each gene
    summarized_data$overall_summary_table <- cbind(summarized_data$overall_summary_table,
        Treatment_group = num_samples_trt, Control_group = num_samples_control,
        delta = delta@id)
    summarized_data$number_repetitions <- num_reps

    return(summarized_data)
}





