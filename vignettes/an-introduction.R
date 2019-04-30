## --------------------------------------------------------------------------
suppressPackageStartupMessages(library(plyranges))
set.seed(100)
df <- data.frame(start=c(2:-1, 13:15), 
                 width=c(0:3, 2:0))

# produces IRanges
rng <- df %>% as_iranges()
rng


## --------------------------------------------------------------------------
# seqname is required for GRanges, metadata is automatically kept
grng <- df %>% 
  transform(seqnames = sample(c("chr1", "chr2"), 7, replace = TRUE),
         strand = sample(c("+", "-"), 7, replace = TRUE),
         gc = runif(7)) %>% 
  as_granges()

grng

## ---- echo = FALSE, out.width="400px", fig.align="center"------------------
knitr::include_graphics("anchors.png", dpi = 150)

