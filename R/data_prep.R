

target <- function(.data, data_type, ...){

  .data$target <- ifelse(data_type == "SR", 1, sum(.data$numerator)/ sum(.data$denominator))

  return(.data)

}

