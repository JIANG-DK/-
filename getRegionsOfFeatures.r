#'@param file a bed file get from ehh i forgot about it,which contains chrnames, gene start sites,gene end sites and geneid.
# the function of the function is to get the positional infomation of promoter,5'UTR,genebody and 3'UTR. be careful! however, some definitions may be incorrect cause i'm not sure about my needs
getRegionsOfFeatures<-function(file){
  a<-matrix(nrow = length(gene.bed[[1]]),ncol = 8)
colnames(a)<-c("promoter_start","promoter_end","five_gene_start","five_gene_end","tes_prox_strat","tes_prox_end","genebody_start","genebody_end")
 
for (i in 1:length(gene.bed[[1]])) {
    promoter_start<-max(1,gene.bed[i,2]-1000)
    promoter_end  <-gene.bed[i,2]
    five_gene_start  <-gene.bed[i,2]
    if (gene.bed[i,7]<1000) {
      five_gene_end  <-gene.bed[i,3]
      tes_prox_strat <-gene.bed[i,3]
      tes_prox_end   <-gene.bed[i,3]+1000
      genebody_start<-NA
      genebody_end  <-NA
    }
    if (gene.bed[i,7]>1000|gene.bed[i,7]<2000) {
      five_gene_end  <-gene.bed[i,2]+1000
      tes_prox_strat<-gene.bed[i,2]+1000
      tes_prox_end  <-gene.bed[i,3]+1000
      genebody_start<-NA
      genebody_end  <-NA
    }
    if (gene.bed[i,7]>2000) {
      five_gene_end  <-gene.bed[i,2]+1000
      tes_prox_strat<-gene.bed[i,3]-1000
      tes_prox_end  <-gene.bed[i,3]+1000
      genebody_start<-gene.bed[i,2]+1000
      genebody_end  <-gene.bed[i,3]-1000
    }
    tes_prox_strat<-max(gene.bed[i,2],gene.bed[i,3]-1000)
    tes_prox_end  <-gene.bed[i,3]+1000
    a[i,]<-c(promoter_start,promoter_end,five_gene_start,five_gene_end,tes_prox_strat,tes_prox_end,genebody_start,genebody_end)
 
  }
return(a)
}
