#' @title For run-time testing
#' @description FUNCTION_DESCRIPTION
#' @param num_samples_trt PARAM_DESCRIPTION, Default: NULL
#' @param num_samples_control PARAM_DESCRIPTION, Default: NULL
#' @param num_genes PARAM_DESCRIPTION, Default: NULL
#' @param means PARAM_DESCRIPTION, Default: NULL
#' @param delta PARAM_DESCRIPTION, Default: NULL
#' @param cor_matrix PARAM_DESCRIPTION, Default: NULL
#' @param num_reps PARAM_DESCRIPTION, Default: 1000
#' @param method PARAM_DESCRIPTION, Default: NULL
#' @param p_adjust PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
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
