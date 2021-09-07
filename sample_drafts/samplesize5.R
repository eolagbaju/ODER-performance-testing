library(stringr)
library(tidyverse)

sstesting <- function(sample_sizes,iterations,seed=NULL,tissue = "blood") {}


gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
gtf_path <- ODER:::.file_cache(gtf_url)


gtex_metadata <- recount::all_metadata("gtex")

start <- stringr::str_c(Sys.time(), " - Start time")
#iterations <- c(1,2,3,4,5,6,7,8,9,10,25,50,100,150,200,250,300,400)
iterations <- c(12,14,16,18,20)
#iterations <- c(12)



#"2021-08-12 16:29:21 - Start time" : "2021-08-12 23:30:13 - Finish time" for chrs 21-22 sp roughly 7hrs

#run started Friday afternoonish - still going
#250 started at - "2021-08-15 20:56:20 - Obtaining mean coverage across 250 samples"
#finished at - "2021-08-17 03:52:56 - Obtaining optimal set of ERs..."  
#took roughly 31 hrs altogether
#finally finished - "2021-08-21 08:15:40 

#iterations <- c(1,2,3,4,5,6,7,8,9,10,25) #[1] "2021-08-12 13:21:07 - Start time" : "2021-08-12 14:13:31 - Finish time"
#iterations <- c(1,2,3,4,5) # ~21mins - "2021-08-12 12:42:49 - Start time": "2021-08-12 13:03:40 - Finish time"
#also got some NAs for 5
#iterations <- c(1,2) 
runs <- 5

x <- list.files("/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/")
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
f_output <- vector(mode = "list", length = runs * length(iterations))




rc = 1
set.seed(123458) #https://stackoverflow.com/questions/44013535/why-set-seed-affects-sample-in-r
for (c in 1:length(iterations)){
  for (b in 1:runs){
    wbsample <- sample(wbfiles,iterations[c])
    opt_ers <- ODER(
      bw_paths = unlist(purrr::map(wbsample,1)), auc_raw = unlist(purrr::map(wbsample,2)), #40e6 * 100,  #gtex_metadata[["auc"]][1], # auc_example,
      auc_target = 40e6 * 100, chrs = "",
      genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
      gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
      exons_no_overlap = NULL, bw_chr = "chr")
    
    opt_mcc <-  as.integer(stringr::str_remove(opt_ers[["opt_mcc_mrg"]][[1]],"mcc_"))
    opt_mrg <-  as.integer(stringr::str_remove(opt_ers[["opt_mcc_mrg"]][[2]],"mrg_"))
    deltas <- opt_ers[["deltas"]]
    med <- deltas[["median"]][deltas[["mcc"]]==opt_mcc & deltas[["mrg"]] ==opt_mrg]
    neq0 <- deltas[["n_eq_0"]][deltas[["mcc"]]==opt_mcc & deltas[["mrg"]] ==opt_mrg]
    row <- list(iterations[c],b,med,neq0) #rc rc+=1
    names(row) <- c("Iteration","Run","Optimum_Median_Exon_delta","Number_of_zero_deltas")
    f_output[[rc]] <- row #rc rc+=1
    rc <- rc + 1
    
    
  }
}

# df_out <- as.data.frame(do.call(rbind,f_output))
# df_fin <- as.data.frame(t(apply(df_out, 1, unlist)))
# 
# df_plot <- df_fin %>% group_by(Iteration) %>% 
#   summarise(mean_delta = mean(Optimum_Median_Exon_delta), sd = sd(Optimum_Median_Exon_delta))

df_out2 <- as.data.frame(do.call(rbind,f_output))
df_fin2 <- as.data.frame(t(apply(df_out2, 1, unlist)))

df_plot2 <- df_fin2 %>% group_by(Iteration) %>% 
  summarise(mean_delta = mean(Optimum_Median_Exon_delta), sd = sd(Optimum_Median_Exon_delta))

df_plot3 <- df_fin2 %>% group_by(Iteration) %>% 
  summarise(mean_delta = mean(Number_of_zero_deltas), sd = sd(Number_of_zero_deltas))

# combiplot <- rbind(df_plot,df_plot2)
#sortedcombiplot <- combiplot[order(combiplot$Iteration),]
# save(df_out, file = stringr::str_c("/home/eolagbaju/projects/data", "/1to400run.Rdata"))
# save(df_fin, file = stringr::str_c("/home/eolagbaju/projects/data", "/1to400runsum.Rdata"))
# save(df_plot, file = stringr::str_c("/home/eolagbaju/projects/data", "/1to400runplot.Rdata"))
# save(df_out2, file = stringr::str_c("/home/eolagbaju/projects/data", "/12to20run.Rdata"))
# save(df_fin2, file = stringr::str_c("/home/eolagbaju/projects/data", "/12to20runsum.Rdata"))
# save(df_plot2, file = stringr::str_c("/home/eolagbaju/projects/data", "/12to20plotmed.Rdata"))
# save(df_plot3, file = stringr::str_c("/home/eolagbaju/projects/data", "/12to20plotnum.Rdata"))

#save(sortedcombiplot, file = stringr::str_c("/home/eolagbaju/projects/data", "/tot1to400medplot.Rdata"))
load("/home/eolagbaju/projects/data/tot1to400medplot.Rdata")

# load("/home/eolagbaju/projects/data/1to400run.Rdata")
# load("/home/eolagbaju/projects/data/1to400runsum.Rdata")
# load("/home/eolagbaju/projects/data/1to400runplot.Rdata")

# 
# ggplot(
#   data = df_fin,
#   mapping = aes(x = Iteration , y = Optimum_Median_Exon_delta)
# ) + geom_line()
# 
# ggplot(data = df_fin) + 
#   geom_point(mapping = aes(x = Iteration, y = Optimum_Median_Exon_delta))

# ggplot( #to do: add tick marks for the points
#   data = df_plot,
#   mapping = aes(x = Iteration, y = mean_delta)
# ) + 
#   geom_line() + 
#   geom_point(shape=4, size=3) +
#   geom_errorbar(aes(ymin=mean_delta - sd, ymax=mean_delta + sd), width=0.05, color = "red") +
#   scale_x_continuous(trans='log10') #scale_x_continuous(limits = c(0,25))
# 
tick_marks <- c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,50,100,150,200,250,300,400)

oneto400plot <- ggplot( #to do: add tick marks for the points
  data = sortedcombiplot,
  mapping = aes(x = Iteration, y = mean_delta)
) + 
  geom_line() + 
  geom_point(shape=4, size=3) +
  geom_errorbar(aes(ymin=mean_delta - sd, ymax=mean_delta + sd), width=0.05, color = "red") +
  scale_x_continuous(trans='log10') #, labels = tick_marks
#geom_errorbar(aes(ymin=mean.temp-sd, ymax=mean.temp+sd), width=.2)
#write.csv(sortedcombiplot, file = "/home/eolagbaju/projects/data/1to400completeplot.csv", row.names = FALSE)
ggsave(filename = "/home/eolagbaju/projects/data/oneto400plot.png" , plot = oneto400plot)

print(start)
print(stringr::str_c(Sys.time(), " - Finish time"))



# 
# test_opt_ers <- ODER(
#   bw_paths = unlist(purrr::map(f5,1)),
#   auc_raw = unlist(purrr::map(f5,2)), #40e6 * 100,  #gtex_metadata[["auc"]][1], # auc_example,
#   auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
#   genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
#   gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
#   exons_no_overlap = NULL, bw_chr = "chr")
