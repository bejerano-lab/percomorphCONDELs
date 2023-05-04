library(TFBSTools)
library(JASPAR2022)
opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["all_versions"]] <- FALSE
PFMatrixList22 <- getMatrixSet(JASPAR2022, opts)
pwmList22 <- toPWM(PFMatrixList22)

pCONDEL1189 <- "TCCAGTGACCACTGGCAGATGGGCTGATGATCTGCAGCTCCAGCTCTGCTAAATTAGGAAGAGAGAGTCTGTTCTCTAAATGACCTCCTCATCAACAGGATGCTCCCTGCAGCCACTCGTCTGTTTGTGAAAACA"
pCONDEL329 <- "ATCCTTTCACTGTGTAATCAGCTCAAACAAAGATAAAGAGTTGTTACGAC"

pCONDEL1189resultList <- searchSeq(pwmList22, pCONDEL1189, seqname="pCONDEL.1189", strand="*", min.score="80%")
pCONDEL1189Results <- as(pCONDEL1189resultList, "data.frame")
pCONDEL1189Results$start <- pCONDEL1189Results$start + 31129643
pCONDEL1189Results$end <- pCONDEL1189Results$end + 31129643
pCONDEL1189Results$chrom <- "chr4"
pCONDEL1189Results$relScore <-pCONDEL1189Results$relScore * 1000
pCONDEL1189ResultsSubset <- pCONDEL1189Results[c("chrom", "start", "end", "TF", "relScore", "strand", "ID", "class")]
pCONDEL1189ResultsSubset <-  pCONDEL1189ResultsSubset[pCONDEL1189ResultsSubset$relScore >= 925,]
write.table(pCONDEL1189ResultsSubset, file="pCONDEL1189tfbsResults925.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

pCONDEL329resultList <- searchSeq(pwmList22, pCONDEL329, seqname="pCONDEL.329", strand="*", min.score="80%")
pCONDEL329results <- as(pCONDEL329resultList, "data.frame")
pCONDEL329results$start <- pCONDEL329results$start + 3548841
pCONDEL329results$end <- pCONDEL329results$end + 3548841
pCONDEL329results$chrom <- "chr14"
pCONDEL329results$relScore <- pCONDEL329results$relScore * 1000
pCONDEL329resultsSubset <- pCONDEL329results[c("chrom", "start", "end", "TF", "relScore", "strand", "ID", "class")]
pCONDEL329resultsSubset <- pCONDEL329resultsSubset[pCONDEL329resultsSubset$relScore >= 925,]
write.table(pCONDEL329resultsSubset, file="pCONDEL329tfbsResults925.tsv", quote = FALSE, sep = "\t", row.names = FALSE)