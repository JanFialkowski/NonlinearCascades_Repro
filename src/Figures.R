source("./src/Ecuador_Pipeline.R")
library(Matrix)
library(latex2exp)

Marketshares <- function(Graph,IDData){
  Inds <- IDData[.(as.numeric(V(Graph)$name)),ciiu4n4]
  Industryoutputs <- rep(0,length(unique(Inds)))
  names(Industryoutputs) <- unique(Inds) 
  s_out <- strength(Graph,mode="out")
  for(Ind in names(Industryoutputs)){
    Industryoutputs[Ind] <- sum(s_out[Inds==Ind])
  }
  marketshares <- rep(0,vcount(Graph))  
  
  for(node in 1:vcount(Graph)){
    marketshares[node] <- s_out[node]/Industryoutputs[Inds[node]]
  }
  return(marketshares)
}

PairDataPipeline <- function(Inputfiles,Outputfiles,OG_ESRI){
  output <- createPairDataTableFromFileList(Inputfiles,Outputfiles)
  output <- unique(output,by=c("i","j"))
  output[,SumESRI := OG_ESRI$ESRI[i,1]+OG_ESRI$ESRI[j,1]]
  output[,SumESRI_Down := OG_ESRI$ESRI[i,2]+OG_ESRI$ESRI[j,2]]
  output[,MaxESRI := pmax(OG_ESRI$ESRI[i,1],OG_ESRI$ESRI[j,1])]
  output[,MaxESRI_down := pmax(OG_ESRI$ESRI[i,2],OG_ESRI$ESRI[j,2])]
  output[,Amp:=ESRI/SumESRI]
  return(output)
}
createPairDataTableFromFileList <- function(Inputfiles,Outputfiles){
  output <- createPairDataTable(Inputfiles[1],Outputfiles[1])
  for(i in seq(length(Inputfiles))[-1]){
    output <- rbind(output,createPairDataTable(Inputfiles[i],Outputfiles[i]))
  }
  return(output)
}
createPairDataTable <- function(Inputfile, Outputfile){
  inp <- readRDS(Inputfile)
  res <- readRDS(Outputfile)
  scenarios <- which(colSums(inp$psi_mat)==2)
  FirstFirm <- numeric(0)
  SecondFirm <- numeric(0)
  PairESRIs <- numeric(0)
  PairESRIs_Down <- numeric(0)
  if(length(scenarios)==dim(inp$psi_mat)[2]){
    FirstFirm <- inp$psi_mat@i[seq(1,2*length(scenarios),2)]+1
    SecondFirm <- inp$psi_mat@i[seq(2,2*length(scenarios),2)]+1
    PairESRIs <- res$ESRI[,1]
    PairESRIs_Down <- res$ESRI[,2]
  } else {
    for(scen in scenarios){
      Firms <- inp$psi_mat@i[(inp$psi_mat@p[scen]+1):inp$psi_mat@p[scen+1]]+1
      FirstFirm <- append(FirstFirm,min(Firms))
      SecondFirm <- append(SecondFirm,max(Firms))
      PairESRIs <- append(PairESRIs,res$ESRI[scen,1])
      PairESRIs_Down <- append(PairESRIs_Down,res$ESRI[scen,2])
    }
  }
  output <- data.table(i=FirstFirm,j=SecondFirm,ESRI=PairESRIs,ESRI_Down=PairESRIs_Down)
  return(output)
}

FirmcountPipe <- function(Inputfiles,Outputfiles,OG_ESRI){
  SinglePipe <- function(ifile,ofile,OG_ESRI){
    inp <- readRDS(ifile)
    out <- readRDS(ofile)
    nfirms <- colSums(inp$psi_mat)
    PairESRIs <- out$ESRI[,1]
    SumESRIs <- rep(0,length(PairESRIs))
    for(i in 1:length(PairESRIs)){
      firms <- which(inp$psi_mat[,i]==1)
      SumESRIs[i] <- sum(OG_ESRI$ESRI[firms,1])
    }
    output <- data.table(nfirms=nfirms,SumESRI=SumESRIs,PairESRIs=PairESRIs) 
    return(output)
  }
  
  output = data.table()
  nfiles <- length(Inputfiles)
  for(f in 1:nfiles){
    output <- rbind(output,SinglePipe(Inputfiles[f],Outputfiles[f],OG_ESRI))
  }
  return(output)
}

Graph <- ReadRawGraph("./data/EL.csv")
Graph <- ProcessRawGraph(Graph)

OG_ESRI <- readRDS("./OG_ESRI.rds")

# Read in the random data and deduplicate it ####
RandomInputs <- c("./data/RandomPairs_ESRIInput.rds")
RandomResults <- c("./data/RandomResults.rds")
FullRandomData <- PairDataPipeline(RandomInputs,RandomResults,OG_ESRI)

# Read in the Subsetdata, select only the pairs and dedupe them ####
ConfirmInputs <- grepv("ESRI",list.files("./Confirm",full.names = T))
ConfirmResults <- grepv("result",list.files("./Confirm",full.names = T))
SubsetInputs <- grepv("thresh",list.files("./Subsets",full.names = T))
SubsetResults <- grepv("result",list.files("./Subsets",full.names=T))

Inputs <- c(ConfirmInputs,SubsetInputs)
Results <- c(ConfirmResults,SubsetResults)
Pair1002Data <- PairDataPipeline(grepv("1002",Inputs),grepv("1002",Results),OG_ESRI) 
Pair1005Data <- PairDataPipeline(grepv("1005",Inputs),grepv("1005",Results),OG_ESRI) 
Pair5002Data <- PairDataPipeline(grepv("5002",Inputs),grepv("5002",Results),OG_ESRI) 
Pair5005Data <- PairDataPipeline(grepv("5005",Inputs),grepv("5005",Results),OG_ESRI) 

PairDataFusion <- rbind(Pair1002Data,Pair1005Data,Pair5002Data,Pair5005Data)
PairDataFusion <- unique(PairDataFusion,by=c("i","j"))
FullPairData <- rbind(PairDataFusion,FullRandomData)
FullPairData <- unique(FullPairData, by = c("i","j"))

FullPairData[,Amp:=(ESRI)/SumESRI]
FullRandomData[,Amp:=(ESRI)/SumESRI]
PairDataFusion[,Amp:=(ESRI)/SumESRI]

# Figure 2: Scatter Plot + subnetworks ####
png("./Figure2.png",width=8.7,height=8.7*1.5,res=600,units="cm")

# Figure 2c: Big scatterplots
bla100 <- FirmcountPipe(grepv("1002_thresh01",ConfirmInputs),grepv("1002_thresh01",ConfirmResults),OG_ESRI) 
bla500 <- FirmcountPipe(grepv("5002_thresh01",ConfirmInputs),grepv("5002_thresh01",ConfirmResults),OG_ESRI) 
Higherorderdata <- FirmcountPipe(ConfirmInputs,ConfirmResults,OG_ESRI)

cols = c("black","blue","green","purple","red")
symbs = c(20,20,17,18,8)
par(mar=c(4,3,0.5,1)+0.1)
plot(FullRandomData[runif(nrow(FullRandomData))<=pointfrac,.(SumESRI,ESRI)],xlim=c(0,0.2),ylim=c(0,1.0), xlab="", ylab="",main="",pch=symbs[1],cex=pointsize)
abline(0,1,lty="dashed")
text(0.02,0.03,adj=c(0,0),TeX("ESRI$(\\Sigma_i \\psi_i)=\\Sigma_i$ESRI$(\\psi_i$)"),srt=atan(93/480)/(2*pi)*360)
points(Higherorderdata[nfirms==2,SumESRI],Higherorderdata[nfirms==2,PairESRIs],col=cols[2],pch=symbs[2],cex=pointsize)
points(Higherorderdata[nfirms==3,SumESRI],Higherorderdata[nfirms==3,PairESRIs],col=cols[3],pch=symbs[3],cex=pointsize)
points(Higherorderdata[nfirms==4,SumESRI],Higherorderdata[nfirms==4,PairESRIs],col=cols[4],pch=symbs[4],cex=pointsize)
points(Higherorderdata[nfirms==5,SumESRI],Higherorderdata[nfirms==5,PairESRIs],col=cols[5],pch=symbs[5],cex=pointsize)
grid()
title(ylab=TeX("ESRI($\\Sigma_i \\psi_i$)"),line=2)
title(xlab=TeX("$\\Sigma_i$ESRI($\\psi_i$)"),line=2)
text(0,1.0,labels="C",adj=c(0.,1))
legend(x=0.2,y=max(Higherorderdata$PairESRIs),legend=c("Random pairs","Extracted pairs","Extracted triples","Extracted quadruples","Extracted quintuples"),pch=symbs,col=cols,bty="n",ncol=1, x.intersp=0.75,xjust=1,yjust=1,pt.cex=pointsize)
par(new=T)
SubPlotPosition <- c(0.56,0.56)
SubPlotSize <- c(0.3,0.3/1.425)*1.1
par(mar=c(3,3,0,0)+0.1)
layout(matrix(c(0,0,0,0,1,0,0,0,0),nrow=3),widths = c(SubPlotPosition[1],SubPlotSize[1],1-SubPlotPosition[1]-SubPlotSize[1]),heights = c(SubPlotPosition[2],SubPlotSize[2],1-SubPlotPosition[2]-SubPlotSize[2]))
plot(0,0,xlim=c(0,0.0125),ylim=c(0,0.05), xlab="", ylab="",main="",pch=20,bg="white",col="white",cex=99999999,yaxt="n",xaxt="n")
axis(2,mgp=c(3,0.5,0),tcl=-0.3)
axis(1,mgp=c(3,0.5,0),tcl=-0.3)
points(0,0,pch=20,cex=9999999,col="white")
grid()
points(FullRandomData[runif(nrow(FullRandomData))<=pointfrac,.(SumESRI,ESRI)],pch=symbs[1],cex=pointsize)
points(Higherorderdata[nfirms==2,SumESRI],Higherorderdata[nfirms==2,PairESRIs],col=cols[2],pch=symbs[2],cex=pointsize)
points(Higherorderdata[nfirms==3,SumESRI],Higherorderdata[nfirms==3,PairESRIs],col=cols[3],pch=symbs[3],cex=pointsize)
points(Higherorderdata[nfirms==4,SumESRI],Higherorderdata[nfirms==4,PairESRIs],col=cols[4],pch=symbs[4],cex=pointsize)
points(Higherorderdata[nfirms==5,SumESRI],Higherorderdata[nfirms==5,PairESRIs],col=cols[5],pch=symbs[5],cex=pointsize)
title(xlab=TeX("$\\Sigma_i$ESRI($\\psi_i$)"),line=1.5)
title(ylab=TeX("ESRI($\\Sigma_i\\psi_i$)"),line=1.5)
abline(0,1,lty="dashed",untf=T)
abline(0,3,lty="dashed",col="grey",untf=T)
dev.off()

#Figure 3: Survival curve ####
pdf("./Figure3.pdf",width=3.43,height=3.43,pointsize=9)
par(mar=c(4,3,1,1)+0.1,lwd=1.5)
plot(function(v){(1-ecdf(FullRandomData$Amp)(v))*nrow(FullRandomData)},log="xy",ylim=c(1,10^6),from=10^-0,ylab="",xlab="",to=260,add=T,col="darkgrey",lwd=1.5)
grid()
plot(function(v){(1-ecdf(FullPairData$Amp)(v))*nrow(FullPairData)},log="xy",ylim=c(1,10^6),from=10^-0,ylab="",xlab="",to=260,add=T,lwd=1.5)
title(ylab=TeX("Number of pairs with amplification $>\\alpha$"),line=2)
title(xlab=TeX("Amplification factor $\\alpha$"),line=2.5)
axis(2,at=10^seq(0,6),labels=TeX(paste0("$10^",seq(0,6),"$")))
fit <- fit_power_law(FullRandomData[Amp>3,Amp])
plot(function(x){(10^5)*x^-(fit$alpha-1)},from=4,to=250,add=T,col="black",lty="dashed")
#text(40,10^2.25,labels = TeX(paste0("# pairs $\\propto\\alpha^{-",format(fit$alpha-1,digits=3),"}$")),adj=c(0,0),cex=1.25)
text(15,10^3,labels = TeX(paste0("# pairs$\\propto\\alpha^{-",format(fit$alpha-1,digits=3),"}$")),adj=c(0,0),cex=1.25)
legend("topright",bty="n",legend=c("crustacean SCN","softdrink SCN","full SCN, random","full SCN, random + extracted"),col=c("blue","green","darkgrey","black"),pch=20)
dev.off()
# Figure Sensitivity_SI ####

Baseline <- FirmcountPipe(grepv("5002_thresh01",ConfirmInputs),grepv("5002_thresh01",ConfirmResults),OG_ESRI) 
png("./pics/Paperfigs/Sensitivity.png",height=8.7*3*0.75,width=8.7*0.9,units="cm",res=600)
layout(matrix(c(1,2,3),nrow=3),heights  = c(1,1,1.1))

# Number of Firms
Thresh2 <- unique(FirmcountPipe(grepv("1002_thresh01",ConfirmInputs),grepv("1002_thresh01",ConfirmResults),OG_ESRI))
par(mar=c(0.5,4,1,1)+0.1)
plot(Thresh2[,.(SumESRI,PairESRIs)],pch=20,col=rgb(0,1,0,1),cex=2,
     xlab="", ylab="",main="",xlim=c(0,0.15),ylim=c(0,1.2),xaxt="n")
axis(1,labels=F)
title(ylab=TeX("ESRI($\\Sigma_i\\psi_i$)"),line=2.5)
grid()
points(unique(Baseline)[,.(SumESRI,PairESRIs)],pch=20,col=rgb(0,0,0,1),cex=1)
legend("topright",legend=c(TeX(c("$N=500$","$N=100$"))),pch=20,col=c("black","green"),bty="n")
abline(0,1,lty="dashed")
text(0,1.2,labels="A")

# First threshold
Thresh1 <- unique(FirmcountPipe(grepv("5005_thresh01",ConfirmInputs),grepv("5005_thresh01",ConfirmResults),OG_ESRI))
par(mar=c(0.5,4,1,1)+0.1)
plot(Thresh1[,.(SumESRI,PairESRIs)],pch=20,col=rgb(0,1,0,1),cex=2,
     xlab="", ylab="",main="",xlim=c(0,0.15),ylim=c(0,1.2),xaxt="n")
title(ylab=TeX("ESRI($\\Sigma_i\\psi_i$)"),line=2.5)
axis(1,labels=F)
grid()
points(unique(Baseline)[,.(SumESRI,PairESRIs)],pch=20,col=rgb(0,0,0,1),cex=1)
legend("topright",legend=c(TeX(c("$\\Theta_1=3$","$\\Theta_1=6$"))),pch=20,col=c("black","green"),bty="n")
abline(0,1,lty="dashed")
text(0,1.2,labels="B")

# Second threshold
Thresh2_1 <- unique(FirmcountPipe(grepv("5002_thresh025",ConfirmInputs),grepv("5002_thresh025",ConfirmResults),OG_ESRI))
Thresh2_2 <- unique(FirmcountPipe(grepv("5002_thresh05",ConfirmInputs),grepv("5002_thresh05",ConfirmResults),OG_ESRI))
par(mar=c(4,4,1,1)+0.1)
plot(Thresh2_2[,.(SumESRI,PairESRIs)],pch=20,col=rgb(0,0,1,1),cex=2,xlab="", ylab="",main="",xlim=c(0,0.15),ylim=c(0,1.2))#,yaxt="n")
title(xlab=TeX("$\\Sigma_i$ ESRI($\\psi_i$)"),line=2.5)
title(ylab=TeX("ESRI($\\Sigma_i\\psi_i$)"),line=2.5)
grid()
points(Thresh2_1[,.(SumESRI,PairESRIs)],pch=20,col=rgb(0,1,0,1),cex=1.5)
points(unique(Baseline)[,.(SumESRI,PairESRIs)],pch=20,col=rgb(0,0,0,1),cex=1)
legend("topright",legend=c(TeX(c("$\\Theta_2=0.9$","$\\Theta_2=0.75$","$\\Theta_2=0.5$"))),pch=20,col=c("black","green","blue"),bty="n")
abline(0,1,lty="dashed")
text(0,1.2,labels="C")
dev.off()