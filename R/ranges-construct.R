# constructor-ranges.R
#' Contstruct a I/GRanges object from a tibble or data.frame
#'
#' @param .data the object to construct a Ranges object from
#' @param ... optional named arguments specifying which the columns in .data
#' containing the core components a Ranges object.
#' @param keep_mcols place the remaining columns into the metadata columns slot
#' (default=TRUE)
#'
#' @description The Ranges function looks for column names in .data called start,
#' end, width, seqnames and strand in order to construct an IRanges or GRanges
#' object. By default other columns in .data are placed into the mcols (
#' metadata columns) slot of the returned object.
#'
#' @return a \link[IRanges]{Ranges} or a \link[GenomicRanges]{GRanges} object.
#' @seealso \link[IRanges]{IRanges-class} \link[GenomicRanges]{GRanges-class}
#'
#' @importFrom rlang quos eval_tidy
#' @importFrom S4Vectors mcols metadata mcols<- metadata<-
#' @importFrom BiocGenerics start end width strand score start<- end<- width<- score<-
#' @importFrom GenomeInfoDb seqnames seqnames<- seqinfo<-
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#'
#' @examples
#' df <- data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0))
#' Ranges(df)
#'
#' df <- data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0), strand = "+")
#' # will return an IRanges object
#' Ranges(df)
#'
#' df <- data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0),
#' strand = "+", seqnames = "chr1")
#' Ranges(df)
#'
#' @export
Ranges <- function(.data, ..., keep_mcols) UseMethod("Ranges")

#' @export
Ranges.default <- function(.data, ..., keep_mcols = TRUE) {
  Ranges.data.frame(as.data.frame(.data), ...)
}

#' @export
Ranges.data.frame <- function(.data, ..., keep_mcols = TRUE) {

  dots <- quos(...)
  col_names <- names(.data)

  if (!is_empty_quos(dots)) {
    valid_args <- names(dots) %in%
      c("start", "end", "width", "seqnames", "strand")
    if (any(!valid_args)) {
      stop("Named arguments must be start, end, width, seqnames, strand",
           .call = FALSE)
    }
    rd <- eval_tidy(dots, .data)
  } else {
    rd <- NULL
  }

  core_i <- quos(start = .data$start, end = .data$end, width = .data$width)
  match_cols_i <- names(core_i) %in% col_names
  match_dots_i <- names(core_i) %in% names(rd)

  if (sum(c(match_cols_i, match_dots_i)) < 2) {
    stop("Unable to construct IRanges from .data must have at least two of
         start, end or width columns.",
         call. = FALSE)
  } else {
    remain_cols <- match_cols_i & !match_dots_i
    remain_core <- core_i[remain_cols]
    if (length(remain_core) > 0) {
      ir <- c(eval_tidy(core_i[match_cols_i], .data),
              rd[names(rd) %in% names(core_i)])
    } else {
      ir <- rd[names(rd) %in% names(core_i)]
    }

    ir <- IRanges(start = ir[["start"]],
                  end = ir[["end"]],
                  width = ir[["width"]])
  }

  # Creating a GRanges object requires a seqnames column
  core_g <- quos(seqnames = .data$seqnames, strand = .data$strand)
  match_core_g <- names(core_g) %in% col_names
  match_dots_g <- names(core_g) %in% names(rd)

  if (all(match_dots_g)) {
    ir <- GRanges(seqnames = rd[["seqnames"]],
                  ranges = ir,
                  strand = rd[["strand"]])
  } else if (all(match_core_g)) {
    gr <- eval_tidy(core_g, .data)
    ir <- GRanges(seqnames = gr[["seqnames"]],
                  ranges = ir,
                  strand = gr[["strand"]])
  } else if (match_dots_g[1] & match_core_g[2]) {
    ir <- GRanges(seqnames = rd[["seqnames"]],
                  ranges = ir,
                  strand = unlist(eval_tidy(core_g[2], .data)))
  } else if (match_core_g[1] & match_dots_g[2]) {
    ir <- GRanges(seqnames = unlist(eval_tidy(core_g[1], .data)),
                  ranges = ir,
                  strand = rd[["strand"]])
  } else if (any(match_core_g)) {
    if (match_core_g[1]) {
      ir <- GRanges(seqnames = unlist(eval_tidy(core_g[1], .data)),
                    ranges = ir)
    }

    else if (match_core_g[2]) {
      mcols(ir) <- unlist(eval_tidy(core_g[2], .data))
      names(mcols(ir)) <- "strand"
    }

  } else if (any(match_dots_g)) {
    if (match_dots_g[1]) {
      ir <- GRanges(seqnames = rd[["seqnames"]],
                    ranges = ir)
    }

    else if (match_dots_g[2]) {
      mcols(ir) <- rd[["strand"]]
      names(mcols(ir)) <- "strand"
    }
  }

  if (keep_mcols) {
    old_cols <- unlist(lapply(dots, quo_name))
    remain_cols <- !(col_names %in%
                       c(old_cols, names(core_i), names(core_g)))

    if (length(names(mcols(ir))) == 0) {
      mcols(ir) <- .data[, remain_cols]
      names(mcols(ir)) <- col_names[remain_cols]
    } else {
      mcols(ir) <- cbind(mcols(ir), .data[, remain_cols])
      names(mcols(ir)) <- c(names(mcols(ir)), col_names[remain_cols])
    }
  }
  ir
}
