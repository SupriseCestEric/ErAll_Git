
library("AnnotationHub")
library("ggbio")
library("rtracklayer")
ahub = AnnotationHub()
ahub = subset(ahub, species == "Homo Sapiens")

qry = query(ahub, c("NFKB", "MEC1"))

peaks = qry["index"]


Rgs = ranges(peaks)
plotGrandLinear(peaks, aes(y = peaks$signalValue), geom = "point")


extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

list.filenames<-list.files()
#problematic samples
list.filenames1 <-list.filenames[-c(7,8)]
list.data<-list()
for (i in 1:length(list.filenames1))
{
  list.data[[i]]<- import.bed(list.filenames1[[i]], extraCols = extraCols_narrowPeak)
}
names(list.data)<-list.filenames1
peaks1<-list()

#Set target region
for (i in 1:length(list.data)){
  list.data[[i]]<-subset(list.data[[i]], seqnames == "4")
  peaks1[[i]]<-list.data[[i]][start(list.data[[i]]) > 69402903 & end(list.data[[i]]) <69446848]
}
names(peaks1)<-list.filenames1

#Get a list of files from all samples
summed_Gr<-unlist(peaks1)

#set colours
colours1<-summed_Gr$name
tmp_frame = data.frame("data"=colours1,"data2"=1:63, stringsAsFactors = FALSE)
tmp_frame$data = substr(tmp_frame$data,1,nchar(tmp_frame$data)-5)
summed_Gr$name<-tmp_frame$data

Graph_Frame = data.frame("data1"=start(summed_Gr),"data2"=summed_Gr$singnalValue, "data3" = colours1, stringsAsFactors = FALSE)

#plot signal intensity
gr<-ggplot(Graph_Frame, aes(x=data1, y=data2, fill = data3))+ geom_boxplot(aes(group = cut_width(data1, 1000))) +facet_wrap(~data3, scale = "free") +xlim(69435000,69442500) + ylim(2.5,9)+ylab("")+xlab("")+theme_bw()
peaks4 = subset(peaks, seqnames == "chr4")


#Example
#Find overlaps between region of interest
UGT_peaks = GRanges(seqnames = Rle("chr4"),ranges = IRanges(start=69402903, end = 69446848), strand = "*")
findOverlaps(UGT_peaks, peaks4)
over = peaks4[start(peaks4)> 69402903 & end(peaks4) <69446848]

scatter.smooth(start(over), over$signalValue, ylab = "H3K4me3 peak intensity",xlab = "chromosome 4 coordinates", xlim = c(69402900,69445000))


#to build GRanges with two isoforms
V1 = GRanges(seqnames = Rle("chr4"),ranges = IRanges(start=c(69403343,69416395,69417542,69426255,69431290,69433479), end = c(669403622,69416614,69417629,69426386,69431438,69434202), names = rep("V1",6)),strand = Rle(strand(c("-","-","-","-","-","-"))))
n2 = GRanges(seqnames = Rle("chr4"),ranges = IRanges(start=c(69403343,69416395,69417542,69426255,69431290,69433479,69441669), end = c(669403622,69416614,69417629,69426386,69431438,69434202,69442018),names = rep("n2",7)),strand = Rle(strand(c("-","-","-","-","-","-","-"))))
m<-GRangesList(V1,n2)
autoplot(m)




