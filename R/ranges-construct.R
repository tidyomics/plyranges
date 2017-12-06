#' Construct a I/GRanges object from a tibble or data.frame
#'
#' @param .data the object to construct a Ranges object from
#' @param ... optional named arguments specifying which the columns in .data
#' containing the core components a Ranges object.
#' @param keep_mcols place the remaining columns into the metadata columns slot
#' (default=TRUE)
#'
#' @description The as_i(g)ranges function looks for column names in .data called start,
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
#' @rdname ranges-construct
#' @examples
#' df <- data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0))
#' as_iranges(df)
#'
#' df <- data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0), strand = "+")
#' # will return an IRanges object
#' as_iranges(df)
#'
#' df <- data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0),
#' strand = "+", seqnames = "chr1")
#' as_granges(df)
#'
#' @export
as_iranges <- function(.data, ..., keep_mcols = TRUE) UseMethod("as_iranges")

#' @export
as_iranges.default <- function(.data, ..., keep_mcols = TRUE) {
  as_iranges.data.frame(as.data.frame(.data), ...)
}

#' @export
as_iranges.data.frame <- function(.data, ..., keep_mcols = TRUE) {
  dots <- quos(...)

  col_names <- names(.data)

  check_names(dots, c("start", "end", "width"))

  if (length(dots) > 0) {
    rd <- lapply(dots, eval_tidy,  data = .data)
  } else {
    rd <- NULL
  }
  # IRanges constructor generate quos for core parts of class
  core_ir <- quos(start = .data$start, end = .data$end, width = .data$width)

  ir <- irng_construct(.data, rd, col_names, core_ir)

  if (keep_mcols) {
    return(make_mcols(.data, ir, col_names, dots, core_ir))
  }
  ir
}

#' @rdname ranges-construct
#' @export
as_granges <- function(.data, ..., keep_mcols = TRUE) UseMethod("as_granges")

#' @export
as_granges.default <- function(.data, ..., keep_mcols = TRUE) {
  as_granges.data.frame(as.data.frame(.data), ...)
}

#' @export
as_granges.data.frame <- function(.data, ..., keep_mcols = TRUE) {
  dots <- quos(...)

  col_names <- names(.data)

  valid_names <- c("start", "end", "width", "seqnames", "strand")
  check_names(dots, valid_names)

  if (length(dots) > 0) {
    rd <- lapply(dots, eval_tidy,  data = .data)
  } else {
    rd <- NULL
  }

  if (!(any(names(rd) %in% "seqnames") | any(names(.data) %in% "seqnames"))) {
    stop("seqnames column is required for GRanges.", call. = FALSE)
  }
  # IRanges constructor generate quos for core parts of class
  core_ir <- quos(start = .data$start, end = .data$end, width = .data$width)
  ir <- irng_construct(.data, rd, col_names, core_ir)

  # GRanges constructor generate quos for core parts of class
  core_gr <- quos(seqnames = .data$seqnames, strand = .data$strand)
  ir <- grng_construct(.data, rd, ir, col_names, core_gr)

  if (keep_mcols) {
    return(make_mcols(.data, ir, col_names, dots, c(core_ir, core_gr)))
  }

  ir
}

check_names <- function(dots, valid_names) {
  if (length(dots) > 0) {
    valid_args <- names(dots) %in% valid_names
    if (any(!valid_args)) {
      stop(paste("Named arguments must be",
                 paste(valid_names, collapse = ","), "."),
           .call = FALSE)
    }
  }
}

make_mcols <- function(.data, ir, col_names, dots, core) {
  # remaining columns in data
  old_cols <- unlist(lapply(dots, quo_name))
  remain_cols <- !(col_names %in% c(old_cols, names(core)))
  if (any(remain_cols)) {
    mcols(ir) <- .data[, remain_cols]
    names(mcols(ir)) <- col_names[remain_cols]
  }
  ir
}

irng_construct <- function(.data, rd, col_names, core_ir) {

  match_cols_i <- names(core_ir) %in% col_names
  match_dots_i <- names(core_ir) %in% names(rd)

  if (sum(c(match_cols_i, match_dots_i)) < 2) {
    stop("Unable to construct IRanges from .data must have at least two of
         start, end or width columns.",
         call. = FALSE)
  } else {
    remain_cols <- match_cols_i & !match_dots_i
    remain_core <- core_ir[remain_cols]
    if (length(remain_core) > 0) {
      ir <- lapply(core_ir[match_cols_i], eval_tidy, data = .data)

      ir <- c(ir, rd[names(rd) %in% names(core_ir)])
    } else {
      ir <- rd[names(rd) %in% names(core_ir)]
    }

    ir <- IRanges(start = ir[["start"]],
                  end = ir[["end"]],
                  width = ir[["width"]])
  }

  return(ir)

}

grng_construct <- function(.data, rd, ir, col_names, core_gr) {

  if (length(rd[["seqnames"]]) > 0) {
    ir <- GRanges(seqnames = rd[["seqnames"]],
                  ranges = ir)
  } else {
    ir <- GRanges(seqnames = unlist(eval_tidy(core_gr[[1]], .data)),
                  ranges = ir)
  }

  if (length(rd[["strand"]]) > 0) {
    strand(ir) <- rd[["strand"]]
  } else {
    if (any(col_names %in% "strand")) {
      strand(ir) <- unlist(eval_tidy(core_gr[[2]], .data))
    }
  }
  return(ir)

}
