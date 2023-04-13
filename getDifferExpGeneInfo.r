getDifferExpGeneInfo<-function(reads_count="~/jdk/rna-seq/matrix/stem.bed",highvalue=5,mod=0){
  #'@ function: return a list contains high/mod/low expression geneid
  #'@ highvalue: a value about tpm, which was defined as a high_expression gene
  #'@ mod:a value about tpm, which was defined as a mod expression gene
  
  counts<-read.table(file = reads_count,header = TRUE,sep = "\t")         
  counts<-counts[,c(1,6,7)]
  
  
  counts[,4]<-(counts[,3]*10^9)/(as.numeric(sum(counts[,3]))*counts[,2])
  counts[,5]<-counts[,4]/sum(counts[,4])*10^6
  counts[,6]<-log(counts[,5],2)
  colnames(counts)[c(3,4,5,6)]<- c("count","fpkm","tpm","log2(tpm)")
  counts[,1]<-gsub(pattern = "\\.Araport.*447$",replacement = "",counts[,1])
 
  alist<-list()
  
  a<-which(counts[,6]<mod) #low gene ==0
  b<- which(counts[,6]>highvalue) #high expression  , customization
  c<-which(counts[,6]>=mod & counts[,6]<=highvalue)#mod exp
  high<- counts[b,]
  low<-counts[a,]
  mod<-counts[c,]
  alist[["high"]]<-high
  alist[["mod"]]<-mod
  alist[["low"]]<-low
  return(alist)
}
