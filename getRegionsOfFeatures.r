#'@param file a bed file get from ehh i forgot about it,which contains chrnames, gene start sites,gene end sites and geneid.
# the function of the function is to get the positional infomation of promoter,5'UTR,genebody and 3'UTR. be careful! however, some definitions may be incorrect cause i'm not sure about my needs
library(TxDb.Athaliana.BioMart.plantsmart28)
getRegionsOfFeatures<-function(upstream = 1000,downstream = 0,txdb=TxDb.Athaliana.BioMart.plantsmart28){
  #'@description a function for getting location info of promoter,5'UTR and so on, return a dataframe lastly.

  promoter<- GenomicFeatures::promoters(txdb,upstream = upstream,downstream = downstream,use.names = T)#get the promoter regions
  proinfo<-as.data.frame(promoter@ranges)
  a<-c(rep(c("+","-","+","-","+","-","+","-","+","-","+","-","+","-"),
           c("5409", "5237" ,"3497" ,"3358" ,"4084" ,"4135" ,"3191", "3193", "4675" ,"4613","59"  , "87"  , "55"  , "78")))#不同染色体正负链信息，还没改 流程化
  proinfo<-cbind(proinfo,a)
   proinfo <- proinfo[!grepl("ATCG|ATMG", proinfo[,4]),]#去除线粒体叶绿体基因
   proinfo<-proinfo[grepl(".*\\.1",proinfo[,4]),]#保留唯一转录本
  n1<-length(which(grepl(pattern = "^..{1}1",proinfo[,4])))  #加chrname
  n2<-length(which(grepl(pattern = "^..{1}2",proinfo[,4])))
  n3<-length(which(grepl(pattern = "^..{1}3",proinfo[,4])))
  n4<-length(which(grepl(pattern = "^..{1}4",proinfo[,4])))
  n5<-length(which(grepl(pattern = "^..{1}5",proinfo[,4])))
  a<-c(rep(c("Chr1","Chr2","Chr3","Chr4","Chr5"),c(n1,n2,n3,n4,n5)))
  proinfo<-cbind(proinfo,a)
  
  #获取基因信息
  gene<-genes(TxDb.Athaliana.BioMart.plantsmart28)
  genesinfo<-as.data.frame(gene@ranges)
  genesinfo<-genesinfo[!grepl("ATCG|ATMG", genesinfo[,4]),]
  
  proinfo[,4]<-gsub(pattern = "\\.1",replacement = "",x = proinfo$names)
  genesinfo<-genesinfo[(genesinfo[,4] %in% proinfo[,4]),] #基因信息
  proinfo  <-proinfo[proinfo[,4]%in% genesinfo[,4],]
  
  proinfo<-proinfo[order(proinfo[,4]),]
  features<-cbind(proinfo,genesinfo)#排序合并
  features[,4]<-NULL
  names(features)<-c("start","end","promoter_length","strand","chrname","gene_start","gene_end","gene_length","gene_name")
  
  
  p<-matrix(nrow = length(features[,1]),ncol = 6)
  for (i in 1:length(features[,1])) {
    five_start   <-  features[i,6]
    if (features[i,8]<2000) {
      
      five_end   <-  features[i,6]+0.2*features[i,8]
      genebody_s <-  features[i,6]+0.2*features[i,8]
      genebody_e <-  features[i,7]-0.2*features[i,8]
      three_start<-  features[i,7]-0.2*features[i,8]
    }else{
      five_end    <- features[i,6]+1000  
      genebody_s  <- features[i,6]+1000 
      genebody_e  <- features[i,7]-1000
      three_start <- features[i,7]-1000
    }
    three_end<-features[i,7]+1000
    p[i,]<-c(five_start,five_end,genebody_s,genebody_e,three_start,three_end)
    
  }
  p<-as.data.frame(p)
  names(p)<-c("five_start","five_end","genebody_s","genebody_e","three_start","three_end")
  features<-cbind(features,p)
  return(features)
  
}

