# ---------------------------------------------------------------------------------------------------------------
# Module 1. Identifying discrete and continuous patterns of expression heterogeneity within cell lines and tumors
# ---------------------------------------------------------------------------------------------------------------

###########################################################################
###########################################################################
#     NOTE: The following code will take several hours to run!!
###########################################################################
###########################################################################


# **************************************************************************
# Basic setup

# load necessary R packages and functions
source("dbscan_programs.R")
source("dbscan_tsne_plot.R")
source("nmf_programs.R")


### New code added by DRT 2021-12-02
# This is to enable loading the data from the Tempo server's data directory
# and to define a new path for the output (outside of the Git repo)

OUTPUT_PATH <- "~/Kinker_output"
DATAPATH <- "/data/QSBSC/Kinker_et_al"

# read scRNA-seq data from cell lines and tumors
expr_ccle <- readRDS(file.path(DATAPATH,"CCLE_heterogeneity_Rfiles/CCLE_scRNAseq_CPM.RDS")) # CCLE cell lines
expr_tumor <- readRDS(file.path(DATAPATH,"CCLE_heterogeneity_Rfiles/tumors_scRNAseq_logTPM.RDS")) # human tumors

# ****************************************************************************************** 
# Identifying heterogeneity programs using tSNE and density-based clustering - cell lines
# (i.e. discrete heterogeneity) 

# get dbscan clusters
discr_clusters_minpt5_1.8_ccle <- lapply(expr_ccle, function(x) dbscan_programs(cpm =  x, eps = 1.8, minpts=5, min_size= 0.02))
discr_clusters_minpt5_1.8_ccle <- discr_clusters_minpt5_1.8_ccle[sapply(discr_clusters_minpt5_1.8_ccle, is.list)]

# save output
if(!file.exists(file.path(OUTPUT_PATH,"module1/discr_clusters_minpt5_eps1.8_ccle.RDS")))
  saveRDS(discr_clusters_minpt5_1.8_ccle, 
          file.path(OUTPUT_PATH,"module1/discr_clusters_minpt5_eps1.8_ccle.RDS"))


## Added by DRT 2021-12-03
## Needed to create empty directories for where output files will be saved
## Checks first to see whether they already exist
invisible(lapply(names(discr_clusters_minpt5_1.8_ccle), function(n) 
{
  if(!dir.exists(file.path(OUTPUT_PATH,"module1",n)))
    dir.create(file.path(OUTPUT_PATH,"module1",n))
}))

# tSNE plot showing dbscan clusters
lapply(names(discr_clusters_minpt5_1.8_ccle), function(x) 
  dbscan_tsne_plot(cpm =  expr_ccle[[x]], eps = 1.8, minpts=5, min_size= 0.02, save.plot = T, 
                   dir=file.path(OUTPUT_PATH,"module1", plot.title = x)))


# ******************************************************************************************** 
# Identifying heterogeneity programs using nonnegative matrix factorization - cell lines
# (i.e. continuous + discrete heterogeneity) 

# perform NMF with ranks ranging from 6 to 9       
w_basis_ccle <- list() # nmf gene scores
h_coef_ccle <- list() # nmf cell scores
       
for(i in names(expr_ccle)) {
    w <- NULL
    h <- NULL
    for(j in 6:9) {
        n <- nmf_programs(expr_ccle[[i]], rank=j)
        colnames(n$w_basis) <- paste0(i, "_", j, ".", 1:j)
        colnames(n$h_coef) <- paste0(i, "_", j, ".", 1:j)
        w <- cbind(w, n$w_basis)
        h <- cbind(h, n$h_coef)
    }
    w_basis_ccle[[i]] <- w
    h_coef_ccle[[i]] <- h
}

# save output
if(!file.exists(file.path(OUTPUT_PATH,"module1/nmf_w_basis_ccle.RDS")))
{
    saveRDS(w_basis_ccle, file.path(OUTPUT_PATH,"module1/nmf_w_basis_ccle.RDS"))
    saveRDS(h_coef_ccle, file.path(OUTPUT_PATH,"module1/nmf_h_coef_ccle.RDS"))
}
       
      
# ******************************************************************************************** 
# Identifying heterogeneity programs using nonnegative matrix factorization - human tumors
# (i.e. continuous + discrete heterogeneity) 

# perform NMF with ranks ranging from 6 to 9       
w_basis_tumor <- list() # nmf gene scores
h_coef_tumor <- list() # nmf cell scores
       
for(i in names(expr_tumor)) {
    w <- NULL
    h <- NULL
    for(j in 6:9) {
        n <- nmf_programs(expr_tumor[[i]], is.log=T, rank=j)
        colnames(n$w_basis) <- paste0(i, "_", j, ".", 1:j)
        colnames(n$h_coef) <- paste0(i, "_", j, ".", 1:j)
        w <- cbind(w, n$w_basis)
        h <- cbind(h, n$h_coef)
    }
    w_basis_tumor[[i]] <- w
    h_coef_tumor[[i]] <- h
}

# save output
if(!file.exists(file.path(OUTPUT_PATH,"module1/nmf_w_basis_tumor.RDS")))
{
    saveRDS(w_basis_tumor, file.path(OUTPUT_PATH,"module1/nmf_w_basis_tumor.RDS"))
    saveRDS(h_coef_tumor, file.path(OUTPUT_PATH,"module1/nmf_h_coef_tumor.RDS"))    
}

