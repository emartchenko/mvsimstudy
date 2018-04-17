#' mvsimstudy: A package for simulating the analysis of multivariate and
#' correlated data.
#'
#' The mvsimstudy package provides two overarching simulation functions -
#' simulate_single_scenario and simulate_scenarios. The first function carries
#' out the generation, analysis, and summarization of a single scenario,
#' returning a list that allows you to examine each individual feature analyzed
#' as well as a summary for the model overall. The simulate_scenarios function
#' simulates multiple scenarios and returns an overall summary for each of them.
#'
#' Along with these functions there are also generate_data, analyze, and summary
#' functions that are accessible if you would like to track intermediate states.
#'
#' @section Notes:
#' For the simulation to run to completion the level of significance must be
#' set with set.alpha() . This must be done once before package use commences.
#' Please see documentation for more details.
#'
#' @examples
#' A single scenario, uncorrelated data analyzed with 'one feature at a time'
#' method
#' \dontrun{
#' if(interactive()){
#' set.alpha(0.05)
#' e <- effect_size(c(0.3, 0, 0.7))
#' single_scenario_simulation(num_samples_trt = 80, num_samples_control= 95, delta = e, method="ofaat")
#'   }
#' }
#'
#' A single scenario, correlated data analyzed with lasso method
#' \dontrun{
#' if(interactive()){
#' set.alpha(0.05)
#' e <- effect_size(c(rep(0.3, 10), rep(0, 60), rep(0.1, 30), rep(0, 20), rep(0.7, 10)))
#' corr_structure <- create_cor_matrix(130, strong = 15, med = 10, weak = 90)
#' single_scenario_simulation(num_samples_trt = 30, num_samples_control= 30,
#' delta = e, cor_matrix = corr_structure, method="lasso")
#'   }
#' }
#'
#' Multiple scenarios, variable sample sizes
#'  \dontrun{
#' if(interactive()){
#' set.alpha(0.05)
#'
#' e <- effect_size(c(0, 0.7, 0.3))
#' corr_structure <- create_cor_matrix(3, med = 2)
#'
#' effect_list <- c(e)
#' sample_size_list <- c(20, 30, 40, 50, 60, 70)
#'
#' simulate_scenarios(effect_list, sample_size_list
#' cor_matrix = corr_structure, method="lasso")
#'   }
#' }
#' @seealso \code{\link{set.alpha}} for additional information
#' @seealso \code{\link{single_scenario_simulation}}
#' @seealso \code{\link{simulate_scenarios}}
#' @seealso \code{\link{analyze}} for details on method parameter
#' @seealso \code{\link{p.adjust}} for p value adjustment arguments

#' @docType package
#' @name mvsimstudy
NULL
