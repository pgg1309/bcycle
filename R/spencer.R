#' Spencer filter
#'
#' Filter a series using a 15-term two-sided Spencer filter.
#'
#' This filter is used in the Bry-Boschan algorithm.
#'
#' @param x A time series (seasonally adjusted).
#' @param lambda Box-Cox transformation parameter. If lambda = "auto",
#'     then a transformation is automatically selected using [BoxCox.lambda()].
#' @return A time series after apllying the two-sided Spencer filter.
#'
spencer <- function(x, lambda = NULL) {
  s <- c(-3, -6, -5, 3, 21, 46, 67, 74, 67, 46, 21, 3, -5, -6, -3)
  s <- s/sum(s)

  # extends ts backward and fwd to allow two-sided filter
  fit <-
    forecast::auto.arima(
      x,
      stationary = FALSE,
      stepwise = FALSE,
      seasonal = FALSE,
      approximation = TRUE,
      trace = FALSE,
      lambda = lambda,
      biasadj = TRUE
    )
  xfc <- forecast::forecast(fit,7)$mean
  # reverses time series
  revx <- ts(rev(x), frequency = frequency(x))
  fitr <-
    forecast::auto.arima(
      revx,
      stationary = FALSE,
      stepwise = FALSE,
      seasonal = FALSE,
      approximation = TRUE,
      trace = FALSE,
      lambda = lambda,
      biasadj = TRUE
    )
  xbc <- forecast::forecast(fitr,7)$mean
  # Extends ts with forecast and backcast
  xpad <-
    ts(c(rev(xbc), x, xfc),
       start = tsp(x)[1] - (7 / frequency(x)),    # 7 datapoints to add
       frequency = frequency(x))

  # calculates the filtered ts
  xsp <- filter(xpad, s, sides = 2)
  xsp <- ts(na.omit(xsp), start = tsp(x)[1], frequency = frequency(x))
  return(xsp)
}

