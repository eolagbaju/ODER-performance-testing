library(stringr)
library(tidyverse)

#gtex_metadata <- recount::all_metadata("gtex")

start <- stringr::str_c(Sys.time(), " - Start time")
#iterations <- c(1,2,3,4,5,6,7,8,9,10,25,50,100,150,200,250,300,400)
iterations <- c(1,2,3,4,5) # takes an hour - "2021-08-11 15:07:44 - Start time" : "2021-08-11 16:05:35 - Finish time"
#also got some NAs for 5
#iterations <- c(1,2) 
runs <- 5

x <- list.files("/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/")
#match(x[1],gtex_metadata$bigwig_file)
wbpath <- "/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/"
wbnum <- sum(str_detect(x, "blood"))
wbfiles <- vector(mode = "list", length = wbnum)
bc <- 1
for (a in 1:length(x)){
  if (str_detect(x[a], "blood")){
    #print(x[a])
    gtex_index <- match(x[a],gtex_metadata$bigwig_file)
    #print(gtex_index)
    auc <- gtex_metadata$auc[gtex_index]
   # print(auc)
    wbfiles[[bc]] <- list(stringr::str_c(wbpath,x[a]),auc)
    bc <- bc + 1
  }
}
f_output <- vector(mode = "list", length = 5 * length(iterations))



gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
gtf_path <- ODER:::.file_cache(gtf_url)
rc = 1
set.seed(123458) #https://stackoverflow.com/questions/44013535/why-set-seed-affects-sample-in-r
for (c in 1:length(iterations)){
  for (b in 1:runs){
    wbsample <- sample(wbfiles,iterations[c])
    for (d in 1:iterations[c]){
      opt_ers <- ODER(
        bw_paths = wbsample[[d]][[1]], auc_raw = wbsample[[d]][[2]], #40e6 * 100,  #gtex_metadata[["auc"]][1], # auc_example,
        auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
        genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
        gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
        exons_no_overlap = NULL, bw_chr = "chr")
        
        opt_mcc <-  as.integer(stringr::str_remove(opt_ers[["opt_mcc_mrg"]][[1]],"mcc_"))
        opt_mrg <-  as.integer(stringr::str_remove(opt_ers[["opt_mcc_mrg"]][[2]],"mrg_"))
        deltas <- opt_ers[["deltas"]]
        med <- deltas[["median"]][deltas[["mcc"]]==opt_mcc & deltas[["mrg"]] ==opt_mrg]
        row <- list(iterations[c],b,med) #rc rc+=1
        names(row) <- c("Iteration","Run","Optimum Median Exon delta")
        f_output[[rc]] <- row #rc rc+=1
        rc <- rc + 1
      
    }
  }
}

df_out <- as.data.frame(do.call(rbind,f_output))

print(start)
print(stringr::str_c(Sys.time(), " - Finish time"))

# [[54]]
# [[54]][[1]]
# [1] "/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/SRR814844_SRS408481_SRX260936_female_whole.blood.bw"
# 
# > match("SRR814844_SRS408481_SRX260936_female_whole.blood.bw",gtex_metadata$bigwig_file)
# [1] 3924

# [[54]][[2]]
# [1] 7761956972

# [[2]]
# [[2]][[1]]
# [1] "/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/SRR660957_SRS389789_SRX222921_male_whole.blood.bw"
# 
# [[2]][[2]]
# [1] 6000185418

# [[2]]
# [[2]][[1]]
# [1] "/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/SRR601263_SRS333488_SRX199066_male_whole.blood.bw"
# 
# [[2]][[2]]
# [1] 5736344125
# 


test_opt_ers <- ODER(
  bw_paths = c("/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/SRR601263_SRS333488_SRX199066_male_whole.blood.bw", 
               "/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/SRR660957_SRS389789_SRX222921_male_whole.blood.bw",
               "/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/SRR814844_SRS408481_SRX260936_female_whole.blood.bw"),
  auc_raw = c(5736344125,6000185418,7761956972), #40e6 * 100,  #gtex_metadata[["auc"]][1], # auc_example,
  auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
  genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
  gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
  exons_no_overlap = NULL, bw_chr = "chr")


set.seed(123458)
fir1 <- sample(wbfiles,1)
s1 <- sample(wbfiles,1)
t1 <- sample(wbfiles,1)
fo1 <- sample(wbfiles,1)
fif1 <- sample(wbfiles,1)
fir2 <- sample(wbfiles,2)
s2 <- sample(wbfiles,2)
t2 <- sample(wbfiles,2)
fo2 <- sample(wbfiles,2)
fif2 <- sample(wbfiles,2)
fir3 <- sample(wbfiles,3)
s3 <- sample(wbfiles,3)
t3 <- sample(wbfiles,3)
fo3 <- sample(wbfiles,3)
fif3 <- sample(wbfiles,3)
fir4 <- sample(wbfiles,4)
s4 <- sample(wbfiles,4)
t4 <- sample(wbfiles,4)
fo4 <- sample(wbfiles,4)
fif4 <- sample(wbfiles,4)
f5 <- sample(wbfiles,5)