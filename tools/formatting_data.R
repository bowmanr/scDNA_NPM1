library(Seurat)
library(dplyr)
scale_data<-readRDS("~/Projects/R_packages/BowmanLab_Dashboard/data/scale_treatment.rds")

count_RNA<-FetchData(object = scale_data,layer = "count",vars=rownames(scale_data))
genes_of_interest<- colnames(count_RNA)[colSums(count_RNA)<5]
norm_RNA<-FetchData(object = scale_data,layer = "data",vars=genes_of_interest)

final<- norm_RNA%>%
          tibble::rownames_to_column(var = "Cell")%>%
          tidyr::pivot_longer(cols = !Cell,names_to = "Gene",values_to = "Expression")%>%
          inner_join(scale_data@meta.data%>%
                      dplyr::select(Group,Cluster=cluster_annotations,Broad_cluster=broad_cell)%>%
                      tibble::rownames_to_column(var = "Cell"),
                      by="Cell")%>%
          dplyr::select(-Cell)
saveRDS(final,"~/Projects/R_packages/BowmanLab_Dashboard/data/scale_treatment.rds")
