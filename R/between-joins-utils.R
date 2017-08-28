# between-joins-utils.R

map_common <- function(x, y, x_name, y_name) {
  if (is(x, "GRanges") & is(y, "DataFrame")) {
    intersect(mcols(x)[[x_name]], y[[y_name]])
  } else {
    intersect(mcols(x)[[x_name]], mcols(y)[[y_name]])
  }
}
index_common <- function(x, y, x_name, y_name) {
  mcols(x)[[x_name]] %in% map_common(x, y, x_name, y_name)
}

rows_common <- function(by, x, y) {
  lapply(seq_along(by),
         function(i) {
           x_name <- names(by)[[i]]
           y_name <- by[[i]]
           index_common(x, y, x_name, y_name)
         })
}

negate_rows_common <- function(by, x, y) {
  lapply(seq_along(by),
         function(i) {
           x_name <- names(by)[[i]]
           y_name <- by[[i]]
           !index_common(x, y, x_name, y_name)
         })
}

filter_common <- function(by, x, y) {
  Reduce(all, rows_common(by, x, y))
}
