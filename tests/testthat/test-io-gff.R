context("test GFF files")

createCorrectGFF <- function(version) {
  space <- c(rep("chr10", 15), rep("chr12", 16))
  start <- c(rep(92828, 4), 92997, rep(94555, 4), rep(94744, 3), rep(95122, 2),
             95348, rep(87984, 4),
             rep(c(88257, 88570, 88860, 89675, 90587, 90796), each = 2))
  end <- c(95504, 95178, 95504, rep(94054, 2), rep(94665, 2), 94615, 94665,
           rep(94852, 3), rep(95178, 2), 95504,
           rep(c(91263, 88017, 88392, 88771, 89018, 89827, 90655, 91263),
               each = 2))
  type <- c("gene", "mRNA", "mRNA", "exon", "CDS", "exon", "exon",
            "CDS", "CDS", "exon", "exon", "CDS", "exon", "CDS", "exon",
            "gene", "mRNA", "exon", "CDS", "exon", "CDS", "exon", "CDS",
            "exon", "CDS", "exon", "CDS", "exon", "CDS", "exon", "CDS")
  type <- factor(type, unique(type))
  source <- factor("rtracklayer")
  phase <- NA_integer_
  score <- c(5, rep(NA, length(type) - 1L))
  strand <- strand(c(rep("-", 14), "*", rep("+", 15), "*"))
  Alias <- CharacterList(c(list(c("FLJ40100", "TUBB8")),
                           rep(list(character()), 14), "LOC100288778",
                           rep(list(character()), 15)))
  ID <- c("GeneID:347688", "873", "872", rep(NA, 12), "GeneID:100288778",
          "4644", rep(NA, 14))
  Name <- c(rep("TUBB8", 3), rep(NA, 12), rep("LOC100288778", 2),
            rep(NA, 14))
  Parent <- CharacterList(c(list(character()), rep("GeneID:347688", 2),
                            rep(list(c("872", "873")), 2),
                            rep(c("872", "873"), 3), rep("873", 3), "872",
                            list(character()), "GeneID:100288778",
                            rep("4644", 14)))
  geneName <- c("tubulin, beta 8", rep(NA, 14),
                "WAS protein family homolog 1; pseudogene", rep(NA, 15))
  genome <- c("hg19", rep(NA, length(geneName) - 1))

  correct_gff3 <- GRanges(space, IRanges(start, end), strand,
                          source, type, score, phase,
                          ID, Name, geneName, Alias, genome, Parent)
  seqinfo(correct_gff3) <- Seqinfo(c("chr10", "chr12"))

  if (version == 3) {
    return(correct_gff3)
  } else if (version == 2) {
    correct_gff2 <- correct_gff3
    toCSV <- function(x) {
      csv <- sapply(x, paste, collapse = ",")
      csv[nchar(csv) == 0] <- NA
      csv
    }
    correct_gff2$Alias <- toCSV(correct_gff2$Alias)
    correct_gff2$Parent <- toCSV(correct_gff2$Parent)
    return(correct_gff2)
  } else {
    correct_gff1 <- correct_gff3[,c("source", "type", "score", "phase")]
    correct_gff1$group <- as.factor(seqnames(correct_gff3))
    return(correct_gff1)
  }
}

test_path <- system.file("tests", package = "rtracklayer")

test_that("reading GFF files returns correct GRanges", {

  correct_gff3 <- createCorrectGFF(3)

  # basic import
  test_gff3 <- file.path(test_path, "genes.gff3")
  test_gr <- read_gff3(test_gff3)
  expect_identical(correct_gff3, test_gr)

  # as GFF(3)File
  test_gff_file <- rtracklayer::GFF3File(test_gff3)
  test_gr <- read_gff3(test_gff_file)
  expect_identical(correct_gff3, test_gr)
  test_gff_file <- rtracklayer::GFFFile(test_gff3, version = "3")
  test_gr <- read_gff3(test_gff_file)
  expect_identical(correct_gff3, test_gr)
  test_gff_file <- rtracklayer::GFF2File(test_gff3)
  expect_error(read_gff3(test_gff_file))
  expect_warning(read_gff2(test_gff3))


  # genome with character arg
  si_hg19 <- rtracklayer::SeqinfoForBSGenome("hg19")
  correct_hg19 <- correct_gff3
  seqlevels(correct_hg19) <- seqlevels(si_hg19)
  seqinfo(correct_hg19) <- si_hg19

  test_gr <- read_gff(test_gff3, genome_info = "hg19")
  expect_identical(test_gr, correct_hg19)

  # genome as genome_info
  # si_hg19_gr <- get_genome_info(si_hg19)
  # test_gr <- read_gff(test_gff3, genome_info = si_hg19_gr)
  # expect_identical(test_gr, correct_hg19)

  # overlaps
  which <- GRanges("chr10:90000-93000")
  which_target <- filter_by_overlaps(correct_gff3, which)
  test_gr <- read_gff3(test_gff3, overlap_ranges = which)
  expect_identical(which_target, test_gr)

  # colnames := "geneName"
  test_gr <- read_gff(test_gff3, col_names = "geneName")
  target <- select(correct_gff3, geneName)
  expect_identical(target, test_gr)

  # exporting gff3
  test_gff3_out <- file.path(tempdir(), "genes.gff3")
  on.exit(unlink(test_gff3_out))
  correct_genome_hg19 <- correct_gff3
  genome(correct_genome_hg19) <- "hg19"
  write_gff3(correct_genome_hg19, test_gff3_out)
  test_gr <- read_gff3(test_gff3_out)
  expect_identical(test_gr, correct_hg19)


  # indexing a gff file
  write_gff3(correct_gff3, test_gff3_out, index = TRUE)
  test_gff_bgz <- paste(test_gff3_out, ".bgz", sep = "")
  on.exit(unlink(test_gff_bgz))
  on.exit(unlink(paste(test_gff_bgz, ".tbi", sep = "")))
  test_gr <- read_gff3(test_gff_bgz, overlap_ranges = which)
  expect_identical(which_target, test_gr)

})

test_that("writing/reading other GFF files", {
  correct_gff3 <- createCorrectGFF(3)
  # 'gff' extension
  test_gff_out <- file.path(tempdir(), "genes.gff")
  on.exit(unlink(test_gff_out))

  write_gff1(correct_gff3, test_gff_out)
  test_gr <- read_gff1(test_gff_out)
  correct_gff1 <- createCorrectGFF(1)
  expect_identical(test_gr, correct_gff1)

  write_gff2(correct_gff3, test_gff_out)
  test_gr <- read_gff2(test_gff_out)
  correct_gff2 <- createCorrectGFF(2)
  expect_identical(test_gr, correct_gff2)

  write_gff3(correct_gff3, test_gff_out)
  test_gr <- read_gff3(test_gff_out)
  expect_identical(test_gr, correct_gff3)

  # 'gff2' extension
  test_gff2_out <- file.path(tempdir(), "genes.gff2")
  write_gff2(correct_gff3, test_gff2_out)
  test_gr <- read_gff2(test_gff2_out)
  expect_identical(test_gr, correct_gff2)

  #  'gff1' extension
  test_gff1_out <- file.path(tempdir(), "genes.gff1")
  write_gff1(correct_gff3, test_gff1_out)
  test_gr <- read_gff1(test_gff1_out)
  expect_identical(test_gr, correct_gff1)
})
