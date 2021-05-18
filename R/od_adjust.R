


ODlimits <- function(.data, data_type = "SR", sr_method = "SHMI"
                     , outlier_method="truncate", outlier_prop=0.1, ...){

  .data %>%
    transformed_zscore(.data, data_type, sr_method, ...) %>%
    outlier_trim(.data, outlier_method, outlier_prop, ...)

  n %>%


   phi_func()

   tau_func()


}


#' Transformation function for z-scoring
#'
#' @description Internal function to perform the transformations for data types.
#'
#' @param .data a tibble or data frame with columns 'numerator' and 'denominator', such as those produced by `aggregate_func`.
#' @param data_type Type of data for adjustment and plotting: Indirectly Standardised ratio (\"SR\"), proportion (\"PR\"), or ratio of counts (\"RC\").
#' @param sr_method Adjustment method, can take the value \"SHMI\" or \"CQC\". \"SHMI\" is default.
#'
#' @return A data.frame of original, aggregated data plus transformed z-score (unadjusted for overdispersion)
#' @keywords internal
#'

transformed_zscore<-function(.data, data_type = "SR", sr_method = "SHMI", ...){

  if(data_type == "SR"){

    # log-transformed SHMI version
    if(sr_method == "SHMI"){
      .data$target_transformed<- 0
      .data$Y <- log(.data$numerator / .data$denominator)
      .data$s <- 1 / (sqrt(.data$denominator))

      # SQRT-transformed CQC version
    } else if(sr_method == "CQC"){
      .data$target_transformed<- 1
      .data$Y <- sqrt(.data$numerator / .data$denominator)
      .data$s  <- 1 / (2 * sqrt(.data$denominator))

    }
  }

  if(data_type == "PR"){

    # use average proportion as target_transformed
    .data$target_transformed<- asin(sqrt(sum(.data$numerator)/ sum(.data$denominator)))

    .data$Y <- asin(sqrt(.data$numerator / .data$denominator))
    .data$s  <- 1 / (2 * sqrt(.data$denominator))


  }

  if(data_type=="RC"){

    # use average proportion as target_transformed
    .data$target_transformed<- log(sum(.data$numerator)/ sum(.data$denominator))

    .data$Y <- log((.data$numerator +0.5) / (.data$denominator +0.5))
    .data$s  <-
      sqrt(
        (.data$numerator/((.data$numerator + 0.5 )^2))
        +
          (.data$denominator/((.data$denominator +0.5)^2))
      )

    #.data$rrS2 = .data$s^2


  }


  if(data_type == "SR" & sr_method=="SHMI"){

    .data$Uzscore <-  sqrt(.data$denominator) * log(.data$numerator / .data$denominator)

  } else {

    .data$Uzscore <- (.data$Y - .data$target_transformed) / .data$s
  }


  return(.data)

}


outlier_trim <- function(.data, outlier_method="truncate", outlier_prop=0.1, ...){

  if(outlier_method=="truncate"){
    truncation(.data, outlier_prop)
  }

  if(outlier_method=="winsorize"){method <- "w"}

  if(outlier_method=="winsorise"){
    winsorisation(.data, outlier_prop)
  } else {
    stop("Unknown truncation method. Choose from 'truncate' or 'winsorise'")
  }

  return(.data)
}


#' Winsorisation function
#'
#' @description Internal function to perform the Winsorisation.
#'
#' @param .data Aggregated model input data
#' @param trim_by The amount to Winsorise the distribution by, prior to transformation. 0.1 means 10\% (at each end).
#'
#' @return A data.frame with winsorised z-scores returned added
#' @keywords internal
#'
#'
winsorisation <- function(.data, trim_by = 0.1){

  lz <- quantile(x = .data$Uzscore, trim_by, na.rm = TRUE)
  uz <- quantile(x = .data$Uzscore, (1 - trim_by), na.rm = TRUE)
  .data$winsorised = ifelse(.data$Uzscore > lz & .data$Uzscore < uz, 0, 1)

  .data$Wuzscore = ifelse(.data$Uzscore < lz
                                 , lz
                                 , ifelse(.data$Uzscore > uz
                                          , uz
                                          , .data$Uzscore))

  return(.data)
}


#' Truncation function for NHSD method
#'
#' @description Internal function to perform the truncation.
#'
#' @param .data Aggregated model input data
#' @param trim_by The amount to truncate the distribution by, prior to transformation. 0.1 means 10\% (at each end).
#'
#' @return A data.frame with truncated z-scores added
#' @keywords internal
#'
#'
truncation <- function(.data = .data, trim_by = 0.1){

  # How many groups for truncation
  k <- 1/trim_by
  maxk<- k-1
  mink<-min(k/k)-1

  .data$rk <- rank(.data$Uzscore, ties.method = "average")
  .data$sp <- floor(.data$rk  * k / (length(.data$rk) + 1))

  .data$truncated = ifelse(.data$sp > mink & .data$sp < maxk, 0, 1)

  .data$Wuzscore = ifelse(.data$truncated == 1
                                 , NA
                                 , .data$Uzscore)

  return(.data)
}


#' Calculate overdispersion ratio
#'
#' @description Internal function to perform the transformations for data types.
#'
#' @param n Single numeric value for the count of the number of groups (and therefore z-scores)
#' @param zscores Vector of z-scores z-scores to be used.  Commonly, this would be 'winsorised' first to remove impact of extreme outliers.  SHMI truncates instead, but this simply reduced the n as well as the z-score.
#'
#' @return A numeric phi value
#' @keywords internal
#'
#'
phi_func <- function(n, zscores){
  phi <- (1 / n) * sum(zscores^2)

  return(phi)
}



#' Calculate the between group standard error (tau2) using a dispersion factor
#'
#' @description Internal function to calculate the additional, between group, standard error (tau2) to add to S2.
#'
#' @param n The number of groups for data items, e.g. hospitals trusts that z-scores are calculated at.
#' @param phi The dispersion ratio, where > 1 means overdispersion
#' @param S Standard error (within cluster, calculated in z-score process)
#'
#' @return A numeric Tau2 value
#' @keywords internal
#'
#'
tau_func <- function(n,  phi, S){

  if(length(S) == 0){
    Tau2 <- 0
  } else {
    if((n*phi) < (n - 1)){
      Tau2 <- 0
    } else {

      Tau2 <- max(0, ((sum(n) * sum(phi)) - (sum(n) - 1)) /
                    (sum(1/(S^2)) - (sum((1/(S^2))^2) / sum(1/(S^2)))))
    }
  }

  return(Tau2)
}
