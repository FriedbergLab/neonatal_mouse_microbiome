#pca analysis for Meghan
## setwd('/home/nzhou/git/meghan')
setwd('/home/idoerg/soft/meghan')
library(scatterplot3d)
colon_temp<-read.csv('colon_new_noneo.csv')
  macro_temp<-read.csv('macro_new_noneo.csv')
histo_temp<-read.csv('histo_new_noneo.csv')
#concatenate
colon=c()
class=c()
macro=c()
histo=c()


# for (i in 1:8){
for (i in 1:6){
  colon = c(colon,colon_temp[,i])
  macro = c(macro,macro_temp[,i])
  temp_vec = c(histo_temp[,i],rep(NA,(32-length(histo_temp[,i]))))
  histo = c(histo,temp_vec)
  class = c(class,rep(colnames(colon_temp)[i],32))
}
new = data.frame(cbind(as.numeric(as.character(colon)),as.numeric(as.character(macro)),as.numeric(as.character(histo))))
data = cbind(new,class)
#data = data[complete.cases(data),]
#Do not omit NA in data frame
##dev.off()
#par(xpd=T)
#par(mfrow=c(1,2))
plot(data[data$class=="LF82A.DSS",c(2,1)],xlab = "Macroscopic Score", ylab="Colon Length",pch=5,col=4,ylim=c(5,10.5))
# points(data[data$class=="G0Neo.LF82..DSS",c(2,1)],col=2,pch=3)
points(data[data$class=="LF82N.DSS",c(2,1)],col=3,pch=4)
points(data[data$class=="DSS",c(2,1)],col=1,pch=20)
# legend("topright",c("G0.LF82.DSS","G0Neo.LF82..DSS","G1.LF82..DSS","Control..DSS"),pch=c(5,3,4,20),col=c(4,2,3,1))
legend("topright",legend=c(expression('LF82'[A]*'+DSS'),expression("LF82"[N]*"+DSS"),"Control+DSS"),pch=c(5,4,20),col=c(4,3,1))

temp = new[,1:2]
new1 = temp[complete.cases(temp),]
a=prcomp(new1)
new_pc=as.matrix(new1)%*%a$rotation
new_class=class[complete.cases(new[,1:2])]
plot(new_pc[new_class=="LF82A.DSS",c(2,1)],pch=5,col=4,xlab="PC2",ylab="PC1")
# points(new_pc[new_class=="G0Neo.LF82..DSS",c(2,1)],col=2,pch=3)
points(new_pc[new_class=="LF82N.DSS",c(2,1)],col=3,pch=4)
points(new_pc[new_class=="DSS",c(2,1)],col=1,pch=20)
legend("bottomleft",c(expression("LF82"[A]*"+DSS"),expression("LF82"[N]*"+DSS"),"Control+DSS"),pch=c(5,4,20),col=c(4,3,1))


#20180123
newdata=aggregate(data[,1:3],by=list(class),FUN=mean,na.rm=T)[1:4]
par(mar=c(5.1,4.1,4.1,9.1),xpd=T) # margins
#colon length vs macroscopic score
plot(newdata[newdata$Group.1=="Control",c(3,2)],xlab="Macroscopic Score",ylab="Colon Length",
     pch=20,col=1,xlim=c(0,5),ylim=c(6,10))
groups = newdata$Group.1
#for (i in 2:8){
for (i in 2:6){
  if (i %% 2==0){
    color=2
  }
  else { 
    color=1
  }
  if (i<3) ch=20
  else if (i<5) ch=2
  else if (i<7) ch=3
  # else ch=4
  points(newdata[newdata$Group.1==groups[i],c(3,2)],pch=ch,col=color) 
}
# legend('topright',inset=c(-0.5,0),groups,pch=rep(c(20,2,3,4),each=2),col=rep(c(1,2),4))
#legend('topright',inset=c(-0.5,0),groups,pch=rep(c(20,2,3),each=2),col=rep(c(1,2),4))
groups2 = c("Control","DSS",expression("LF82"[A]), expression("LF82"[A]*"+DSS"),expression("LF82"[N]), expression("LF82"[N]*"+DSS"))
legend('topright',inset=c(-0.35,0),groups2,pch=rep(c(20,2,3),each=2),col=rep(c(1,2),4))

#colon length vs histological score
plot(newdata[newdata$Group.1=="Control",c(4,2)],xlab="Histology Score",ylab="Colon Length",
     pch=20,col=1,xlim=c(0,11),ylim=c(6,10))
groups = newdata$Group.1
for (i in 2:8){
# for (i in 2:6){
  if (i %% 2==0){
    color=2
  }
  else { 
    color=1
  }
  if (i<3) ch=20
  else if (i<5) ch=2
  else if (i<7) ch=3
  else ch=4
  points(newdata[newdata$Group.1==groups[i],c(4,2)],pch=ch,col=color) 
}
# legend('topright',inset=c(-0.5,0),groups,pch=rep(c(20,2,3,4),each=2),col=rep(c(1,2),4))
legend('topright',inset=c(-0.35,0),groups2,pch=rep(c(20,2,3,4),each=2),col=rep(c(1,2),4))

#PCA
b = prcomp(newdata[,2:4],scale. = T)
new_new_pc=b$x

plot(x=new_new_pc[groups=="Control",1],y=new_new_pc[groups=="Control",2],xlab="PC1",ylab="PC2",xlim=c(-2,3),ylim=c(-1,1),pch=20,col=1)
for (i in 2:8){
  if (i %% 2==0){
    color=2
  }
  else { 
    color=1
  }
  if (i<3) ch=20
  else if (i<5) ch=2
  else if (i<7) ch=3
  else ch=4
  points(x=new_new_pc[groups==groups[i],1],y=new_new_pc[groups==groups[i],2],pch=ch,col=color) 
}
# legend('topright',inset=c(-0.5,0),groups,pch=rep(c(20,2,3),each=2),col=rep(c(1,2),4))
legend('topright',inset=c(-0.35,0),groups2,pch=rep(c(20,2,3),each=2),col=rep(c(1,2),4))


#import cytokines values
il17_temp<-read.csv('Il17_new_noneo.csv')
ifng_temp<-read.csv('Ifng_new_noneo.csv')
il1b_temp<-read.csv('il1b_new_noneo.csv')

il17=c()
ifng=c()
il1b=c()


for (i in 1:6){
  il17 = c(il17,il17_temp[,i])
  ifng = c(ifng,ifng_temp[,i])
  temp_vec = c(il1b_temp[,i],rep(NA,(24-length(il1b_temp[,i]))))
  il1b = c(il1b,temp_vec)
  class = rep(colnames(il17_temp),each=24)
}
cytodata = cbind(il17,ifng,il1b) 
temp = data.frame(cytodata)
cytodata<-cbind(temp,class)
cyto_new = aggregate(cytodata[,1:3],by=list(class),FUN=mean,na.rm=T)

# alldata = cbind(colnames(il17_temp),newdata[c(1,2,3,4,7,8),2:4],cyto_new[,2:4])
alldata = cbind(colnames(il17_temp),newdata[c(1,2,3,4,5,6),2:4],cyto_new[,2:4])

colnames(alldata)[1:4]<-c('class','colonLength','macro','histo')
#pca with all data and 6 groups
c = prcomp(alldata[,2:7],scale. = T)
# c = prcomp(alldata[,2:4],scale. = T)
all_new_pc=c$x

plot(x=all_new_pc[1,1],y=all_new_pc[1,2],xlab="PC1",ylab="PC2",xlim=c(-2,4),ylim=c(-1,1),pch=20,col=1)
for (i in 2:6){
  if (i %% 2==0){
    color=2
  }
  else { 
    color=1
  }
  if (i<3) ch=20
  else if (i<5) ch=2
  else  ch=3

  points(x=all_new_pc[i,1],y=all_new_pc[i,2],pch=ch,col=color) 
}
# legend('topright',inset=c(-0.5,0),colnames(il17_temp),pch=rep(c(20,2,4),each=2),col=rep(c(1,2),4))
# legend('topright',inset=c(-0.52,0),colnames(il17_temp),pch=rep(c(20,2,3),each=2),col=rep(c(1,2),4))
legend('topright',inset=c(-0.35,0),groups2,pch=rep(c(20,2,3),each=2),col=rep(c(1,2),4))

