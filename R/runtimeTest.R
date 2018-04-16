#' @title For run-time testing
#' @description Sys.time() is used to track how long each piece of the
#' simulation needs to run. Please read ?simulate_single_scenario for detailed
#' walkthrough of parameters.
#' @param num_samples_trt Number of samples in group 1. Default: NULL
#' @param num_samples_control Number of samples in group 2. Default: NULL
#' @param num_genes Number of features being tested. Default: NULL
#' @param means A vector of the means for each feature. Default: NULL
#' @param delta An effect_size object containing a vector of effects.
#' Default: NULL
#' @param cor_matrix A matrix specifying the correlations. Default: NULL
#' @param num_reps Number of iterations of the simulation overall. Default: 1000
#' @param method Analysis method used. See ?analyze Default: NULL
#' @param p_adjust p value adjustment method used. See ?p.adjust Default: NULL
#' @return Summary of how long each piece of the simulation runs as well as
#' what kind of analysis method used.
#' @details For testing purposes only.
#' @rdname test_single_scenario_simulation
#' @export

test_single_scenario_simulation <- function(num_samples_trt = NULL,
    num_samples_control = NULL, num_genes = NULL, means = NULL, delta = NULL,
    cor_matrix = NULL, num_reps = 1000, method = NULL, p_adjust = NULL) {

    start_total_time <- Sys.time()
    start_dg <- Sys.time()

    generated_data <- generate_data(num_samples_trt, num_samples_control,
        means, delta@effects, cor_mat, num_reps = num_reps)

    data_generation <- Sys.time() - start_dg

    start_ad <- Sys.time()

    analyzed_data <- analyze(generated_data, method = method)
    analyzing_data <- Sys.time() - start_ad

    start_sd <- Sys.time()
    summarized_data <- summary(analyzed_data, delta, analysis_method = method)

    summarized_data$overall_summary_table <- cbind(summarized_data$overall_summary_table,
        Treatment_group = num_samples_trt, Control_group = num_samples_control,
        delta = delta@id)
    summarized_data$number_repetitions <- num_reps
    summarizing_data <- Sys.time() - start_sd

    total_time <- Sys.time() - start_total_time

    time_summary <- cbind(Num_samples = num_samples_trt + num_samples_control,
        Num_genes = length(means), num_iterations = num_reps, Data_generation = round(data_generation,
            2), Data_analysis = round(analyzing_data, 2), Analysis_method = method,
        Data_summarization = round(summarizing_data, 2), Total_time = round(total_time,
            2))
    return(time_summary)


}
