#' @title Summary
#' @description Called after the analyze function within the simulation.
#' Though analysis is finished this function still requires the method of
#' analysis used in order to summarize the output correctly.
#' @param data_to_summarize Results of analyze function.
#' @param e The effect_size object used in data generation. The 'answer key'
#' to significant features found through analysis.
#' @param analysis_method Analysis method used, see analyze documentation for
#' more information. Default: NULL
#' @return A list of 2. Includes the overall summary dataframe and the
#' individual feature summary dataframe. Overall summary includes number of
#' true positives, false positives, true negatives, false negatives, FDR,
#' TPR/Sensitivity , and where applicable Model_p.value, AUC, and AUC_sd.
#' Individual feature summary includes the term, proportion of times it was
#' selected, and the user-specified effect size.
#' @details Best used from within the simulation functions to avoid any
#' confusion with arguments.
#' @examples
#' 
#' set.alpha(0.05)
#' e <- effect_size(c(0.3, 0, 0.7))
#' example <- generate_data(50, 50, c(1, 2, 3), c(0.3, 0, 0.7))
#' analyzed_example <- analyze(example, method= "ofaat")
#' summarized_example <- summary(analyzed_example, e, "ofaat")
#' 
#' 
#' @seealso
#'  \code{\link[dplyr]{select}}
#' @rdname summary
#' @export
#' @importFrom dplyr select
summary <- function(data_to_summarize, e, analysis_method = NULL) {
    if (is.null(analysis_method)) {
        stop("Please specify analysis method to be used.")
    }
    gene_names <- unique(as.vector(data_to_summarize$term))
    if (analysis_method == "ofaat") {

        summarized_by_gene <- data_to_summarize %>% dplyr::select(reps,
            term, selected) %>% group_by(term) %>% summarise(prop_p_under_alpha = sum(selected)/n()) %>%
            slice(match(gene_names, term))

        in_process <- data_to_summarize %>% dplyr::select(reps, term,
            selected) %>% spread(key = term, value = selected) %>%
            dplyr::select(gene_names)

        # NA for the purpose of avoiding errors downstream
        model_avg <- NA
        auc_avg <- NA
        auc_sd <- NA
    } else if (analysis_method == "mv_glm" | analysis_method == "lasso") {
        # in_process for the purpose of summarizing across rows
        in_process <- data_to_summarize %>% dplyr::select(reps, term,
            selected) %>% spread(key = term, value = selected) %>%
            dplyr::select(gene_names) %>% dplyr::select(-`(Intercept)`)

        # summarized_by_gene for the purpose of summarizing by columns
        summarized_by_gene <- data_to_summarize %>% dplyr::select(reps,
            term, selected) %>% group_by(term) %>% summarise(proportion_selected = sum(selected)/n()) %>%
            slice(match(gene_names, term)) %>% filter(term != "Model") %>%
            filter(term != "AUC") %>% filter(term != "(Intercept)")

        # store average Model p value (if applicable)
        model_avg <- mean(in_process$Model, na.rm = TRUE)

        # store average AUC and standard deviation
        auc_avg <- mean(in_process$AUC, na.rm = TRUE)
        auc_sd <- sd(in_process$AUC, na.rm = TRUE)

        # remove 'whole model' parameters before checking individual gene
        # results
        in_process <- dplyr::select(in_process, -Model, -AUC)
    }
    real_effects <- as.logical(e@effects)
    selection_summary <- t(t(in_process) == real_effects)
    selection_summary <- as.data.frame(selection_summary)
    colnames(selection_summary) <- 1:length(real_effects)
    summarized <- summarize_overall_power(selection_summary, real_effects)

    # summary of single scenario
    summarized <- cbind(summarized, Model_p.value = model_avg, AUC = auc_avg,
        AUC_sd = auc_sd)
    summarized_by_gene <- cbind(summarized_by_gene, specified_effect_size = e@effects)
    results <- list(overall_summary_table = summarized, individual_gene_summary = summarized_by_gene)


    return(results)
}

#' @title Summarize overall power
#' @description An internal function called by summary.
#' Processes each row to calculate the false discovery rate in an
#' accurate way (rather than making the anticipated calculations with the
#' overall results, the false discoveries for each row are noted and averaged
#' to get the final result.)
#' @param x A dataframe which contains the summary of selected (1) and not
#' selected (0) features.
#' @param truth A logical vector that states at which positions the user
#' specified there would be an effect (TRUE), and at which the effect was to be
#' 0 (FALSE).
#' @return A dataframe containing the overall summary results. See ?summary for
#' more details.
#' @details Internal function.
#' @rdname summarize_overall_power
summarize_overall_power <- function(x, truth) {

    true_neg = 0
    false_neg = 0
    true_pos = 0
    false_pos = 0
    # summarize by column
    for (i in 1:length(truth)) {

        if (truth[[i]] == FALSE) {
            true_neg <- true_neg + sum(x[[i]])
            false_pos <- false_pos + sum(!x[[i]])
        }
        if (truth[[i]] == TRUE) {
            true_pos <- true_pos + sum(x[[i]])
            false_neg <- false_neg + sum(!x[[i]])
        }
    }
    # insufficiently accurate false positive rate fdr <- false_pos/
    # (false_pos+true_pos)
    tpr <- true_pos/(false_neg + true_pos)

    col_names <- c("True Positive", "False Positive", "True Negative",
        "False Negative", "FDR", "TPR/Sensitivity")

    values <- c(true_pos, false_pos, true_neg, false_neg, tpr, tpr)

    summarized <- as.data.frame(matrix(nrow = 0, ncol = length(values)))
    colnames(summarized) <- col_names
    summarized[1, ] <- values


        # Rob wants to fix ; not sure tally would work here as a
        # comparison must be made for each value (ex. in a given row you
        # will have 1s and 0s; if the corresponding column to a 1 is TRUE
        # then true_pos+1 but if it's FALSE then false_pos+1). does tally
        # have a way to support this logic?

        # if FDR= 1 - positive predictive value, then does 1-FDR = PPV ?
        mv_summary_table <- data.frame(matrix(ncol = 5, nrow = nrow(x)))
        cols <- c("reps", "true_pos", "false_pos", "true_neg", "false_neg")
        colnames(mv_summary_table) <- cols
        mv_summary_table$reps <- 1:nrow(x)
        mv_summary_table[is.na(mv_summary_table)] <- 0

        for (sample in 1:nrow(x)) {
            for (i in 1:length(truth)) {
                if (truth[[i]] == TRUE) {
                  if (x[[sample, i]] == TRUE) {
                    mv_summary_table$true_pos[[sample]] <- mv_summary_table$true_pos[[sample]] +
                      1
                  }
                  if (x[[sample, i]] == FALSE) {
                    mv_summary_table$false_neg[[sample]] <- mv_summary_table$false_neg[[sample]] +
                      1
                  }
                }
                if (truth[[i]] == FALSE) {
                  if (x[[sample, i]] == TRUE) {
                    mv_summary_table$true_neg[[sample]] <- mv_summary_table$true_pos[[sample]] +
                      1
                  }
                  if (x[[sample, i]] == FALSE) {
                    mv_summary_table$false_pos[[sample]] <- mv_summary_table$false_neg[[sample]] +
                      1
                  }
                }
            }


        mv_summary_table <- mv_summary_table %>% mutate(avgd_fdr = false_pos/(false_pos +
            true_pos))
        # proper calculation for FDR
        fdr_weighted <- mean(mv_summary_table$avgd_fdr, na.rm = TRUE)
        summarized$FDR <- fdr_weighted
    }

    return(summarized)
}
