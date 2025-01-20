int_info.default <- function(
    x,
    min_time = 0L, ...) {

  x <- round(x, 6)
  # check inputs
  assert_numeric(x, lower = 0, any.missing = FALSE)
  assert_numeric(min_time, lower  = 0L)

  # sort x and add origin if necessary
  if (is.unsorted(x)) {
    x <- sort(x)
  }
  if (min(x != 0)) {
    x <- c(0, x)
  }

  intlen <- diff(x)
  intlen <- round(intlen, 6)
  tstart <- x[-length(x)]
  tstart <- round(tstart, 6)
  tend   <- tstart + intlen
  tend <- round(tend, 6)

  tdf <- data.frame(
    tstart = tstart,
    tend   = tend,
    intlen = intlen) %>%
    mutate(
      intmid = tstart + intlen / 2,
      interval = paste0("(", tstart, ",", tend, "]"),
      interval = factor(.data$interval, levels = unique(.data$interval))
    )
  filter(tdf, tstart >= min_time)
}