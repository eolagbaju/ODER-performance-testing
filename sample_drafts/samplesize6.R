library(stringr)
library(tidyverse)

sstesting <- function(sample_sizes,iterations,seed=NULL,tissue = "blood") {
  gtex_metadata <- recount::all_metadata("gtex")
  
  x <- list.files("/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/")
  tpath <- "/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/"
  tnum <- sum(str_detect(x, tissue))
  tfiles <- vector(mode = "list", length = tnum)
  tc <- 1 #tissue count
  for (a in 1:length(x)){
    if (str_detect(x[a], tissue)){
      gtex_index <- match(x[a],gtex_metadata$bigwig_file)
      auc <- gtex_metadata$auc[gtex_index]
      tfiles[[tc]] <- list(stringr::str_c(tpath,x[a]),auc)
      tc <- tc + 1
    }
  }
  f_output <- vector(mode = "list", length = iterations * length(sample_sizes))
  
  rc = 1
  if (!is.null(seed)) {set.seed(seed) }
  
    for (c in 1:length(sample_sizes)){
    for (b in 1:iterations){
      tsample <- sample(tfiles,sample_sizes[c])
      opt_ers <- ODER(
        bw_paths = unlist(purrr::map(tsample,1)), auc_raw = unlist(purrr::map(tsample,2)), #40e6 * 100,  #gtex_metadata[["auc"]][1], # auc_example,
        auc_target = 40e6 * 100, chrs = "",
        genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
        gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
        exons_no_overlap = NULL, bw_chr = "chr")
      
      opt_mcc <-  as.integer(stringr::str_remove(opt_ers[["opt_mcc_mrg"]][[1]],"mcc_"))
      opt_mrg <-  as.integer(stringr::str_remove(opt_ers[["opt_mcc_mrg"]][[2]],"mrg_"))
      deltas <- opt_ers[["deltas"]]
      med <- deltas[["median"]][deltas[["mcc"]]==opt_mcc & deltas[["mrg"]] ==opt_mrg]
      neq0 <- deltas[["n_eq_0"]][deltas[["mcc"]]==opt_mcc & deltas[["mrg"]] ==opt_mrg]
      row <- list(sample_sizes[c],b,med,neq0) #rc rc+=1
      names(row) <- c("Iteration","Run","Optimum_Median_Exon_delta","Number_of_zero_deltas")
      f_output[[rc]] <- row #rc rc+=1
      rc <- rc + 1
      
      
    }
    }
  
  return(f_output)
}



gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
gtf_path <- ODER:::.file_cache(gtf_url)
ns <- c(1,2)  #c(1,2,3,4,5,6,7,8,9,10,25)
reps <- 2#10
f_output <- sstesting(sample_sizes=ns,iterations=reps,seed=123458,tissue = "blood")

#############################################################

df_out <- as.data.frame(do.call(rbind,f_output))
df_fin <- as.data.frame(t(apply(df_out, 1, unlist)))

df_plot <- df_fin %>% group_by(Iteration) %>%
  summarise(mean_delta = mean(Optimum_Median_Exon_delta), sd = sd(Optimum_Median_Exon_delta))

# combiplot <- rbind(df_plot,df_plot2)
#sortedcombiplot <- combiplot[order(combiplot$Iteration),]
# save(df_out, file = stringr::str_c("/home/eolagbaju/projects/data", "/1to400run.Rdata"))
# save(df_fin, file = stringr::str_c("/home/eolagbaju/projects/data", "/1to400runsum.Rdata"))
# save(df_plot, file = stringr::str_c("/home/eolagbaju/projects/data", "/1to400runplot.Rdata"))
# save(df_out2, file = stringr::str_c("/home/eolagbaju/projects/data", "/12to20run.Rdata"))
# save(df_fin2, file = stringr::str_c("/home/eolagbaju/projects/data", "/12to20runsum.Rdata"))
# save(df_plot2, file = stringr::str_c("/home/eolagbaju/projects/data", "/12to20plotmed.Rdata"))
# save(df_plot3, file = stringr::str_c("/home/eolagbaju/projects/data", "/12to20plotnum.Rdata"))
save(sortedcombiplot, file = stringr::str_c("/home/eolagbaju/projects/data", "/tot1to400medplot.Rdata"))
load("/home/eolagbaju/projects/data/1to400run.Rdata")
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

ggplot( #to do: add tick marks for the points
  data = df_plot,
  mapping = aes(x = Iteration, y = mean_delta)
) +
  geom_line() +
  geom_point(shape=4, size=3) +
  geom_errorbar(aes(ymin=mean_delta - sd, ymax=mean_delta + sd), width=0.05, color = "red") +
  scale_x_continuous(trans='log10') #scale_x_continuous(limits = c(0,25))

#tick_marks <- c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,50,100,150,200,250,300,400)

oneto400plot <- ggplot( #to do: add tick marks for the points
  data = sortedcombiplot,
  mapping = aes(x = Iteration, y = mean_delta)
) + 
  geom_line() + 
  geom_point(shape=4, size=3) +
  geom_errorbar(aes(ymin=mean_delta - sd, ymax=mean_delta + sd), width=0.05, color = "red") +
  scale_x_continuous(trans='log10') #, labels = tick_marks
#geom_errorbar(aes(ymin=mean.temp-sd, ymax=mean.temp+sd), width=.2)
write.csv(sortedcombiplot, file = "/home/eolagbaju/projects/data/1to400completeplot.csv", row.names = FALSE)
ggsave(filename = "/home/eolagbaju/projects/data/oneto400plot.png" , plot = oneto400plot)

print(start)
print(stringr::str_c(Sys.time(), " - Finish time"))
