library(tidyverse)

unilanes <- function(lane){
  uniq_lanes <- 
    lane %>% 
    unique() %>% 
    str_c(collapse = ",") %>% 
    str_split(",") %>% 
    unlist() %>% 
    unique() %>% 
    as.integer() %>% 
    sort()
  return(uniq_lanes)
}

download.file(url = "https://github.com/LieberInstitute/BrainSeq_Phase1/blob/master/phenotype_annotated_szControlEqtl_DLPFC.rda?raw=true",  
              destfile = file.path(tempdir(), "data.rda"))

load(file = file.path(tempdir(), "data.rda"))

pd_tidy <- pd %>% 
  mutate(Lane_1 = Lane_1 %>% stringr::str_remove("L")) %>% 
  mutate(Lane_2 = Lane_2 %>% stringr::str_remove("L")) %>% 
  mutate(Lane_3 = Lane_3 %>% stringr::str_remove("L"))

uniq_lanes <- 
  pd_tidy[["Lane_1"]] %>% 
  unique() %>% 
  str_c(collapse = ",") %>% 
  str_split(",") %>% 
  unlist() %>% 
  unique() %>% 
  as.integer() %>% 
  sort()

uniq_lanes2 <- unilanes(pd_tidy[["Lane_2"]])
uniq_lanes3 <- unilanes(pd_tidy[["Lane_3"]])

lanes <- vector(mode = "list", length(uniq_lanes))
lanes2 <- vector(mode = "list", length(uniq_lanes2))
lanes3 <- vector(mode = "list", length(uniq_lanes3))

for(i in seq_along(uniq_lanes)){
  
  lanes[[i]] <- 
  pd_tidy %>% 
    filter(str_detect(Lane_1, as.character(uniq_lanes[i]))) %>% 
    mutate(Lane_1 = uniq_lanes[i])
  
}

for(i in seq_along(pd_tidy[["Lane_2"]])){
  
  lanes2[[i]] <- 
    pd_tidy %>% 
    filter(str_detect(Lane_2, as.character(uniq_lanes2[i]))) %>% 
    mutate(Lane_2 = uniq_lanes2[i])
  
}

for(i in seq_along(pd_tidy[["Lane_3"]])){
  
  lanes3[[i]] <- 
    pd_tidy %>% 
    filter(str_detect(Lane_3, as.character(uniq_lanes3[i]))) %>% 
    mutate(Lane_3 = uniq_lanes3[i])
  
}

pd_tidy_l1 <- do.call(bind_rows, lanes)
pd_tidy_l2 <- do.call(bind_rows, lanes2)
pd_tidy_l3 <- do.call(bind_rows, lanes3)

# str_c(pd_tidy_l1[["RNum"]][[1]], "_",pd_tidy_l1[["Flowcell_1"]][[1]],
#       "_",pd_tidy_l1[["Index_1"]][[1]],"_L00", pd_tidy_l1[["Lane_1"]][[1]] )

pd_tidy_l1$bigwigname <- str_c(pd_tidy_l1[["RNum"]], "_",pd_tidy_l1[["Flowcell_1"]],"_",pd_tidy_l1[["Index_1"]],
                               "_L00", pd_tidy_l1[["Lane_1"]], "!lieber_phase1_dlpfc!hg38!local.all.bw" )

pd_tidy_l2$bigwigname <- str_c(pd_tidy_l2[["RNum"]], "_",pd_tidy_l2[["Flowcell_2"]],"_",pd_tidy_l2[["Index_2"]],
                               "_L00", pd_tidy_l2[["Lane_2"]], "!lieber_phase1_dlpfc!hg38!local.all.bw" )

pd_tidy_l3$bigwigname <- str_c(pd_tidy_l3[["RNum"]], "_",pd_tidy_l3[["Flowcell_3"]],"_",pd_tidy_l3[["Index_3"]],
                               "_L00", pd_tidy_l3[["Lane_3"]], "!lieber_phase1_dlpfc!hg38!local.all.bw" )

tot_tidy_pd <- rbind(pd_tidy_l1,pd_tidy_l2,pd_tidy_l3)

n_occur <- data.frame(table(tot_tidy_pd$bigwigname))
#n_occur[n_occur$Freq > 1,]

# 
# 
# pd_l2 <- pd_tidy_l1 %>% 
#   filter(!str_detect(Lane_2, "NUL"))
# 
# count <- 0
# 
# for (r in 1:nrow(pd_l2)){
# 
#   if (str_detect(pd_l2[["Lane_1"]][[r]] , pd_l2[["Lane_2"]][[r]]  ))  {
#     count <- count + 1
#     print("repeat")
#     print(pd_l2[[1]][[r]])
#     print(pd_l2[["Lane_1"]][[r]])
#     print(pd_l2[["Lane_2"]][[r]])
#   }
# }
# 
# for (ro in 1:nrow(pd)){
#   
# }
# 
# 
# 
