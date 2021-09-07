library(stringr)
library(tidyverse)
iterations <- c(1,2,3,4,5,6,7,8,9,10,25,50,100,150,200,250,300,400)
runs <- 5

x <- list.files("/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/")
wbpath <- "/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/"
wbnum <- sum(str_detect(x, "blood"))
wbfiles <- vector(mode = "character", length = wbnum)
bc <- 1
for (a in 1:length(x)){
  if (str_detect(x[a], "blood")){
    wbfiles[bc] <- stringr::str_c(wbpath,x[a])
    bc <- bc + 1
  }
}
f_output <- vector(mode = "list", length = 5)

set.seed(123458)

gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
gtf_path <- ODER:::.file_cache(gtf_url)
rc = 1
for (b in 1:runs){
  wbsample <- sample(wbfiles,iterations[1])
  
  opt_ers <- ODER(
    bw_paths = wbsample, auc_raw = 40e6 * 100,  #gtex_metadata[["auc"]][1], # auc_example,
    auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
    genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
    gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
    exons_no_overlap = NULL, bw_chr = "chr"
  )
  
  opt_mcc <-  as.integer(stringr::str_remove(opt_ers[["opt_mcc_mrg"]][[1]],"mcc_"))
  opt_mrg <-  as.integer(stringr::str_remove(opt_ers[["opt_mcc_mrg"]][[2]],"mrg_"))
  deltas <- opt_ers[["deltas"]]
  med <- deltas[["median"]][deltas[["mcc"]]==opt_mcc & deltas[["mrg"]] ==opt_mrg]
  row <- list(iterations[1],b,med) #rc rc+=1
  names(row) <- c("Iteration","Run","Optimum Median Exon delta")
  f_output[[b]] <- row #rc rc+=1
}


print(f_output)
#wbfiles

wb1 <- sample(wbfiles,1)
wb2 <- sample(wbfiles,1)



opt_ers <- ODER(
  bw_paths = wb1, auc_raw = 40e6 * 100,  #gtex_metadata[["auc"]][1], # auc_example,
  auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
  genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
  gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
  exons_no_overlap = NULL, bw_chr = "chr"
)

opt_mcc <-  as.integer(stringr::str_remove(opt_ers[["opt_mcc_mrg"]][[1]],"mcc_"))
opt_mrg <-  as.integer(stringr::str_remove(opt_ers[["opt_mcc_mrg"]][[2]],"mrg_"))
#mrg <- opt_ers[["opt_mcc_mrg"]][[2]]
deltas <- opt_ers[["deltas"]]
med <- deltas[["median"]][deltas[["mcc"]]==opt_mcc & deltas[["mrg"]] ==opt_mrg]
print(med)