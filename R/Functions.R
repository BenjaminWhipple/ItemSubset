#' Get best n items
#' Uses the maximum average error across quantiles as the performance metric. Employs exhaustive search across item combinations. Is predictably slow for large numbers of items.
#' @param n The size of the subset to choose.
#' @param items vector of column names (strings)
#' @param data data.frame from which to pull values
#' @param actual_scores actual scores to compare performance to
#' @param scale The size of the scale to be used. Defualts to 5.
#' @param resamples The number of resample replicates by which to compute errors. Defaults to 0.
#' @param quantile_groups The number of quantile groups across which to compute the average error. Defaults to 4.
#' @export
#' get_best_n_items()
get_best_n_items <- function(n, items, data, actual_scores, scale=5, quantile_groups=4, resamples = 0){
  #Given items, actual scores, and a size of subset, find the subset that best approximates the actual scores. Note that this function only really works for relatively small size of items. items choose n has brutal scaling.
  #1. Generate all possible subsets of size n
  #2. for each subset, compute error of choice of subset.
  #3. choose subset that has minimal n
  #4. optionally, serialize all results
  #5. return best subset and associated errors.

  #Note: This procedure uses a minimax method across percentile groups.

  item_combinations <- combn(x=items, m=n)
  if (resamples == 0){

    inputs <- list()
    error <- list()
    average_errors <- list()

    #We compute the percentiles here in order to better allow for bootstrapping.
    actual_percentiles <- ecdf(actual_scores)(actual_scores)
    quantile_average_errors <- matrix(0,dim(item_combinations)[2],quantile_groups)
    #Compute error
    for (i in 1:dim(item_combinations)[2]){
      quantile_avg_err <- list()

      new_scores <- rowSums(data[item_combinations[,i]])
      new_percentiles <- ecdf(new_scores)(new_scores)
      absolute_errors <- abs(new_percentiles-actual_percentiles)

      temp_data <- data.frame("abs_err" = absolute_errors, "percentile"=actual_percentiles)

      for (j in 1:quantile_groups){
        data_percentile_subset <- subset(temp_data, percentile<(1/quantile_groups)*j & percentile>((1/quantile_groups)*(j-1)))

        subset_average_error <- sum(subset(data_percentile_subset, percentile<(1/quantile_groups)*j & percentile>((1/quantile_groups)*(j-1)))$abs_err)/nrow(data_percentile_subset)

        quantile_avg_err <- append(quantile_avg_err, subset_average_error)

        quantile_average_errors[i,j] <- subset_average_error
      }


      max_avg_err <- quantile_avg_err[which.max(quantile_avg_err)]

      inputs <- append(inputs, list(item_combinations[,i]))

      average_errors <- append(average_errors, max_avg_err)
    }

    unlisted_average_err <- unlist(average_errors)

    mins <- which.min(unlisted_average_err)

    quantile_average_error <- quantile_average_errors[mins,]

    chosen_inputs <- item_combinations[,mins]

    avg_error_mins <- unlist(unlisted_average_err[mins])

    result <- list(chosen_inputs,c(quantile_average_error),avg_error_mins)

    return(result)
  }

  if (resamples > 0){
    subset_counts <- data.frame('combination'=1:dim(item_combinations)[2],'count'=numeric(dim(item_combinations)[2]))

    complete_data <- data.frame(data[items],'score'=actual_scores)

    #Store averaged quantile errors
    subset_averaged_quantile_errors <- matrix(0,dim(item_combinations)[2],quantile_groups)

    for (k in 1:(resamples)){
      print(c(n,k))
      #resample
      boot_indices <- sample(1:nrow(complete_data), replace = TRUE)
      boot_data <- complete_data[boot_indices,]

      inputs <- list()
      error <- list()
      average_errors <- list()

      #We compute the percentiles here in order to better allow for bootstrapping.
      actual_percentiles <- ecdf(boot_data$score)(boot_data$score)

      quantile_average_errors <- matrix(0,dim(item_combinations)[2],quantile_groups)
      #Compute error
      for (i in 1:dim(item_combinations)[2]){
        quantile_avg_err <- list()

        new_scores <- rowSums(boot_data[item_combinations[,i]])
        new_percentiles <- ecdf(new_scores)(new_scores)
        absolute_errors <- abs(new_percentiles-actual_percentiles)

        temp_data <- data.frame("abs_err" = absolute_errors, "percentile"=actual_percentiles)

        for (j in 1:quantile_groups){
          data_percentile_subset <- subset(temp_data, percentile<(1/quantile_groups)*j & percentile>((1/quantile_groups)*(j-1)))

          subset_average_error <- sum(subset(data_percentile_subset, percentile<(1/quantile_groups)*j & percentile>((1/quantile_groups)*(j-1)))$abs_err)/nrow(data_percentile_subset)

          quantile_avg_err <- append(quantile_avg_err, subset_average_error)

          quantile_average_errors[i,j] <- subset_average_error
        }


        max_avg_err <- quantile_avg_err[which.max(quantile_avg_err)]

        inputs <- append(inputs, list(item_combinations[,i]))

        average_errors <- append(average_errors, max_avg_err)
      }


      unlisted_average_err <- unlist(average_errors)

      mins <- which.min(unlisted_average_err)

      quantile_average_error <- quantile_average_errors[mins,]

      arg_mins <- inputs[mins]

      avg_error_mins <- unlist(unlisted_average_err[mins])


      subset_averaged_quantile_errors = ((k-1)/k)*subset_averaged_quantile_errors + (1/k)*quantile_average_errors

    }

    maxes <- apply(subset_averaged_quantile_errors,1,max)
    mins <- which.min(maxes)

    chosen_inputs <- item_combinations[,mins]
    chosen_quantiles_error <- subset_averaged_quantile_errors[mins,]

    print(chosen_quantiles_error)

    max_err <- max(chosen_quantiles_error)


    result <- list(chosen_inputs, chosen_quantiles_error, max_err)

    return(result)
  }

}

#' Reorganize Output
#' This function is for internal use within get_best_item_subsets. Perhaps it should not be a function.
#' @param best_item_subsets_result Again, this is for internal use.
#' reorganize_output()
reorganize_output <- function(best_item_subsets_result){

  subsets <- list()
  n <- 1:length(best_item_subsets_result)
  absolute_errors <- list()
  average_errors <- list()
  quantile_errors <- list()

  for (i in 1:length(best_item_subsets_result)){
    subsets <- append(subsets,list(best_item_subsets_result[[i]][[1]]))
    quantile_errors <-append(quantile_errors, list(best_item_subsets_result[[i]][[2]]))
    average_errors <- append(average_errors,best_item_subsets_result[[i]][[3]])
  }

  result <- list('Subsets'=subsets,'QuantileErrors'=quantile_errors,'MaxAverageErrors'=unlist(average_errors))

  return(result)
}

#' Get Best Item Subsets
#' This is the main function of the program
#' @param items vector of column names (strings)
#' @param data data.frame from which to pull values
#' @param actual_scores actual scores to compare performance to
#' @param scale The size of the scale to be used. Defualts to 5.
#' @param resamples The number of resample replicates by which to compute errors
#' @param quantile_groups The number of quantile groups across which to compute the average error. Defaults to 4.
#' @export
#' get_best_item_subsets()
get_best_item_subsets <- function(items,data,actual_scores,scale=5,resamples=0,quantile_groups=4){
  #Expect n<= items
  #Default to 5 point scale.

  res <- list()

  for (i in 1:(length(items)-1)){
    print(i)

    if (resamples == 0){
      current_res <- get_best_n_items(i,items,data,actual_scores,scale=scale)
    }

    if (resamples > 0){
      print('Entering computation')
      current_res <- get_best_n_items(i,items,data,actual_scores,scale=scale,resamples = resamples,quantile_groups = quantile_groups)
    }

    res <- append(res, list(current_res))
  }

  print(res)

  return(reorganize_output(res))
}
