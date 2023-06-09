## R script running Susie + Annotate on SparsePro example data
## S. Ning
## 20230602
## SuSie with Annotated Prior
## Using GLM to model the prior prob given annotation covariates

#library(mrpipeline)
library(susieR)
library(stringr)

## data dir
path_data <- c("./dat/")

## read in data
# gwas/z_score
if(!file.exists("SparsePro_test1.Rdata")){
  gwas_id <- "C1"
  gwas_file <- paste0(gwas_id, ".txt")
  dat_gwas <- read.table(file.path(path_data, gwas_file), header = F, row.names = 1)
  colnames(dat_gwas) <- "z.score"
  list_snp <- rownames(dat_gwas)
  
  # ld matrix
  dat_ld <- read.table(file.path(path_data, "ld.txt"))
  dat_ld <- as.matrix(dat_ld)
  colnames(dat_ld) <- list_snp
  rownames(dat_ld) <- list_snp
  
  # annotation
  dat_annot <- read.table(file.path(path_data, "anno.txt"), header = T, row.names = 1)
  
}

#save(gwas_id, dat_gwas, dat_ld, dat_annot, list_snp, file = "SparsePro_test1.Rdata")
load("SparsePro_test1.Rdata")

list_harm <- list(dat_gwas = dat_gwas, tab_lookup=NULL)

## running susie by batch
## batch by snp_size
batch_snp_size <- 2500
total_snp <- nrow(list_harm$dat_gwas)

#save time for each batch prep
time_log_susie_prep <- numeric(ceiling(total_snp/batch_snp_size))

#prep data format for each batch
list_dat_gwas <- list()
list_snp_gwas <- list()

for(i in 1:length(time_log_susie_prep)){
    time_cur <- Sys.time()
    idx_start <- batch_snp_size * (i-1) + 1
    idx_stop <- min(batch_snp_size * (i), length(list_snp))
    dat_maf <- list(z.score = dat_gwas$z.score[idx_start:idx_stop], ld.cor = dat_ld[idx_start:idx_stop, idx_start:idx_stop])
    idx_include <- rownames(dat_maf$ld.cor)%in%list_snp
    dat_maf$z.score <- dat_maf$z.score[idx_include]
    dat_maf$ld.cor <- dat_maf$ld.cor[idx_include, idx_include]
    list_dat_gwas[[i]] <- dat_maf
    list_snp_gwas[[i]] <- rownames(dat_maf$ld.cor)
    time_prev <- time_cur
    time_cur <- Sys.time()
    time_log_susie_prep[i] <-  difftime(time_cur, time_prev, units = "min")
    if(i%%10 == 0){
      print(i)
    }
  }

n_batch <- length(list_dat_gwas)
list_snp_susie <- (unlist(list_snp_gwas))

## re-order dat_annot to match with susie data
dat_annot_susie <- dat_annot[match(list_snp_susie, rownames(dat_annot)),]
rownames(dat_annot_susie) <- list_snp_susie

## saving GWAS data + annot
file_name <- paste0(paste("SparsePro", gwas_id, "test", sep = "_"), 
                    ".Rdata")
save(batch_snp_size, list_dat_gwas, gwas_id, list_snp_gwas, time_log_susie_prep, dat_annot_susie, list_harm, file=file_name)

source("annotate_susie_functions.R")
#### running susie
###initialize
n_iter <- 100
seed <- 1001
options(warn=0)

test <- annotate_susie_beta(list_dat_gwas=list_dat_gwas, list_snp_gwas=list_snp_gwas, dat_annot_susie=dat_annot_susie, n_iter = 20)
