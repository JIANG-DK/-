system("Macs3 callpeak -f BAM -t xxx.bam --outdir x/x/x/ -n name --nomodel --shift 37 --extsize 73",intern = T)
system("macs3 callpeak -f BAM -t /wrk/wenchenjin/work/Chenjin_merge_1__MNase/pre_processed/bamfiles/A6_6_MN_Ara_stem_128x8ul__merge_data0202.sorted_rmdup.bam --outdir /wrk/jiangdingkun/jdk/MNase-seq1/MNase-seq1/callpeak/ -n Ara_stem_8ul --nomodel --shift 37 --extsize 73",intern = T)
system("macs3 callpeak -f BAM -t /wrk/wenchenjin/work/Chenjin_merge_1__MNase/pre_processed/bamfiles/A7_7_MN_Ara_stem_128x10ul__merge_data0202.sorted_rmdup.bam --outdir /wrk/jiangdingkun/jdk/MNase-seq1/MNase-seq1/callpeak/ -n Ara_stem_10ul --nomodel --shift 37 --extsize 73",intern = T)

peak_bed_filelist <-list.files(path = "~/jdk/MNase-seq1/MNase-seq1/callpeak/",pattern = "narrowPeak",full.names = T,)#获取路径下的包含pattern的文件名列表，是macs3输出的narrowpeak文件
peak_Anno_df<-function(peak_bed_file,tssregion=1000){
  library(clusterProfiler)
  library(ChIPseeker)
  library(GenomicFeatures)
  library(TxDb.Athaliana.BioMart.plantsmart28)
  
  txdb <- TxDb.Athaliana.BioMart.plantsmart28#构建txdb文件
  
  peak <- ChIPseeker::readPeakFile(peak_bed_file)#注释peak
  
  seqlevels(peak)[1:5] <- sub("Chr", "", seqlevels(peak)[1:5])
  seqlevels(peak)[seqlevels(peak) == "mitochondria"] <- "Mt"
  seqlevels(peak)[seqlevels(peak) == "chloroplast"] <- "Pt"#将染色体名称对应上（线粒体叶绿体或应该删除?）
  
  peakAnno <- ChIPseeker::annotatePeak(peak = peak,TxDb = txdb,tssRegion = c(-tssregion,tssregion))
  peakAnno_df<-as.data.frame(peakAnno)
  return(peakAnno_df)
} 
tmp <- lapply(peak_bed_filelist,peak_Anno_df)#批量对peak文件进行注释，tmp中含有多个peak注释结果。
df<-do.call(rbind,lapply(tmp,function(x){ 
  
  num1 = length(grep('Promoter', x$annotation))
  num2 = length(grep("5' UTR", x$annotation))
  num3 = length(grep('Exon', x$annotation))
  num4 = length(grep('Intron', x$annotation))
  num5 = length(grep("3' UTR", x$annotation))
  num6 = length(grep('Intergenic', x$annotation))
  return(c(num1,num2,num3,num4,num5,num6 ))
})) #提取不同peak的结构注释信息构建矩阵
colnames(df) <- c('Promoter',"5' UTR",'Exon','Intron',"3' UTR",'Intergenic')#赋值列名
rownames(df) <- unlist(lapply(peak_bed_filelist,function(y){basename(strsplit(y,'\\.')[[1]][1])}))#basename获取路径中的文件名
df1 <- apply(df, 1, function(x) x/sum(x))#转百分数
# 将数据转换成 ggpurb 所需要的格式
library(reshape2)
df2 <- reshape2::melt(df1)
colnames(df2) <- c('dis','sample','fraction')
df2[,2]<- sub("Ara_stem_10ul_peaks","stem_10ul",df2[,2])
df2[,2]<- sub("Ara_stem_8ul_peaks","stem_8ul",df2[,2])
df2[,2]<- sub("Ara_stem_4ul_peaks","stem_4ul",df2[,2])
df2[,2]<- sub("Ara_stem_1ul_peaks","stem_1ul",df2[,2])
library(ggpubr)#多个样本作图
my_color = c("#DD2A18", # 红色
             "#41678B", # 蓝色
             "#42A441", # 绿色
             "#8C4498", # 紫色
             "#FE7400", # 橙色
             "#FEFE2D") # 黄色
ggpubr::ggbarplot(df2, "sample", "fraction",
                  fill = "dis", color = "dis",palette=my_color)+
  guides(fill = guide_legend(nrow = 1))
