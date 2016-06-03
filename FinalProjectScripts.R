library(dplyr)
library(reshape2)
library(ggplot2)
dta<-read.csv("Plate diversity_final.csv")
plateloc<-read.csv("Plate Locations.csv")
status<-read.csv("Status.csv")
point<-read.csv("PointCount.csv")
names(dta)<-c("CollectDate","PlateNo","VialID","Box","Row","Column","Taxa","Species","Color","VialType","Fixative","Notes") #rename data frame columns
names(plateloc)<-c("PlateNo","PlateLoc","PlateTrt","Depth","Lat","Long","Dist")
fulldta<-inner_join(dta,plateloc,by="PlateNo")
fulldta<-full_join(fulldta,status,by="Species")
fulldta<-fulldta[,c(1,2,13:18,7:9,19,3:6,10:12)]
count<-dcast(fulldta,PlateLoc+PlateTrt+PlateNo+Dist~PlateTrt,fun.aggregate=length)
names(count)<-c("PlateLoc","PlateTrt","PlateNo","Dist","Fixed","Floating")
count$count<-count$Fixed+count$Floating
count2<-data.frame(count[,-5:-6])

count3<-dcast(fulldta,PlateLoc+PlateTrt+PlateNo+Dist+Status~PlateTrt,fun.aggregate=length)
names(count3)<-c("PlateLoc","PlateTrt","PlateNo","Dist","Status","Fixed","Floating")
count3$count<-count3$Fixed+count3$Floating
count3<-data.frame(count3[,-6:-7])

statuscount<-dcast(fulldta,PlateLoc+PlateTrt+PlateNo+Dist~Status,fun.aggregate=length)
statuscount$total<-statuscount$Invasive+statuscount$Native+statuscount$Unknown

##to find exaxt distance of outside plates
library(sp)
plateloc<-plateloc[,-2:-4]
plateloc<-plateloc[,-4]
pts<-plateloc[,2:3]
pts<-pts[,c(2,1)]
pts2<-data.matrix(pts)

km<-spDistsN1(pts2,pts2[1,],longlat=TRUE)
m<-km*1000
platedist<-data.frame(plateloc$PlateNo,m)
names(platedist)<-c("PlateNo","Dist_m")

distdata<-inner_join(count2,platedist,by="PlateNo")
distdata<-distdata[,c(1:3,6,5)]

distdata2<-inner_join(count3,platedist, by="PlateNo")
distdata2<-distdata2[,c(1:3, 7,5:6)]

##shannon index: relative abundance of each sp (number of points for each sp divided by points of all spp on plate) multiplied by natural log of relative abundance and summed across all spp on a plate, then multiplied by -1

uniqueplates <- unique(point$PlateNo)
rich<-0
div<-0

Richness<-function(point,plate){
  platedat<-subset(point, PlateNo==plate)
  platesum<-data.frame(table(platedat$Primary))
  platesum<-subset(platesum, Freq > 0)
  platesum<-subset(platesum, Var1 != "")
  totalsp <- nrow(platesum)
  return(totalsp)
}

Shannon<-function(point,plate){
  platedat<-subset(point, PlateNo==plate)
  platesum<-data.frame(table(platedat$Primary))
  platesum<-subset(platesum, Freq > 0)
  platesum<-subset(platesum, Var1 != "")
  totalsp <- nrow(platesum)
  totalcount <- sum(platesum$Freq)
  platesum$relFreq <- platesum$Freq/totalcount
  platesum$X <- platesum$relFreq*log(platesum$relFreq)
  div<-sum(platesum$X)*-1
  return(div)
}

for (i in 1:length(uniqueplates)){
  rich[i]<-Richness(point, uniqueplates[i])
}
for (i in 1:length(uniqueplates)){
  div[i]<-Shannon(point, uniqueplates[i])
}

platerich<-data.frame(uniqueplates, rich)

platediv<-data.frame(uniqueplates, div)
names(platediv)<-c("PlateNo","Diversity")
divdata<-inner_join(count2,platediv, by="PlateNo")
divdata2<-inner_join(count3,platediv, by="PlateNo")

p<-ggplot(divdata,aes(PlateLoc,count))
p+geom_boxplot()+theme_bw()+labs(x="Plate Location", y="Species Richness/Plate")
p2<-ggplot(divdata,aes(PlateLoc,Diversity))
p2+geom_boxplot()+theme_bw()+labs(x="Plate Location", y="Shannon Diversity Index/Plate")
p3<-ggplot(divdata,aes(PlateTrt,count))
p3+geom_boxplot()+theme_bw()+labs(x="Plate Treatment", y="Species Richness/Plate")
p4<-ggplot(divdata,aes(PlateTrt,Diversity))
p4+geom_boxplot()+theme_bw()+labs(x="Plate Treatment", y="Shannon Diversity Index/Plate")
p5<-ggplot(data=divdata,aes(x=Dist,y=count))+geom_boxplot(aes(group=Dist))
p5+facet_grid(PlateTrt~Dist,scales="free",as.table=FALSE)+labs(y="Species Richness/Plate",x="Distance from Marina (m)")+ggtitle("Distance from Marina (m)")+theme_bw()+theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_text(size=10),axis.title.x=element_blank())+theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10))
p6<-ggplot(data=divdata,aes(x=Dist,y=Diversity))+geom_boxplot(aes(group=Dist))
p6+facet_grid(PlateTrt~Dist,scales="free",as.table=FALSE)+labs(y="Shannon Diversity Index/Plate",x="Distance from Marina (m)")+ggtitle("Distance from Marina (m)")+theme_bw()+theme(plot.title=element_text(size=10), axis.text.y=element_text(size=10),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_text(size=10),axis.title.x=element_blank())+theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10))
p7<-ggplot(data=divdata2,aes(x=Dist,y=count))+geom_boxplot(aes(group=Dist))
p7+facet_grid(PlateTrt~Status,scales="free",as.table=FALSE)+labs(y="Species Richness/Plate",x="Distance from Marina (m)")+theme_bw()+theme(axis.text=element_text(size=10),axis.title=element_text(size=10))+theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10))
p8<-ggplot(data=divdata2,aes(x=Dist,y=Diversity))+geom_boxplot(aes(group=Dist))
p8+facet_grid(PlateTrt~Status,scales="free",as.table=FALSE)+labs(y="Shannon Diversity Index/Plate",x="Distance from Marina (m)")+theme_bw()+theme(axis.text=element_text(size=10),axis.title=element_text(size=10))+theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=10))

