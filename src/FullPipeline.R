# install packages from source
# setwd("~/NonlinearCascades_Repro/") # Change this to match your folder structure if needed

install.packages("./GLcascade_0.9.3.1.zip", 
                 repos = NULL, type = "win.binary")
install.packages("./fastcascade_0.9.3.1.zip", 
                 repos = NULL, type = "win.binary")

source("./src/Ecuador_Pipeline.R")
library(Matrix)
library(GLcascade)
# This is file takes a graph and performs all supersampling steps on that graph

PrepareScenarios <- function(input,OG_ESRI,nfirms,nscenarios,nruns,Name,Folder,Thresh=0.1){
  fcount <- length(input$p)
  ScenariosPerRun <- nscenarios
  Runs <- nruns
  psi_mat <- Matrix(0,nrow=fcount,ncol=ScenariosPerRun)
  Bad_Firms <- which((OG_ESRI$ESRI[,1]>Thresh) | OG_ESRI$ESRI[,2]==0)
  FirmsPerRun <- nfirms
  input$track_h <- F
  for(Run in 1:Runs){
    print(paste0("Run: ",Run))
    psi_mat <- Matrix(0,nrow=fcount,ncol=ScenariosPerRun)
    for(i in 1:ScenariosPerRun){
      Firms <- sample(seq(fcount)[-Bad_Firms],FirmsPerRun)
      psi_mat[Firms,i] <- 1
    }
    print("Finished psi_mat")
    input$psi_mat <- psi_mat
    input$run_id <- paste0(Name,Run)
    saveRDS(input,paste0(Folder,input$run_id,"_ESRIInput.rds"))
  }
}
CalcDiffsFromFiles <- function(Inputfiles,Resultfiles,OG_ESRI){
  print("This function assumes the same number of runs for each input file")
  nfiles <- length(Inputfiles)
  n <- dim(readRDS(Resultfiles[1])$ESRI)[1]
  RelDiffs <- rep(0,nfiles*n)
  RelDiffs_Down <- rep(0,nfiles*n)
  ESRIs <- rep(0,nfiles*n)
  for(i in seq(nfiles)){
    print(paste("Progress:",i/nfiles))
    result <- readRDS(Resultfiles[i])
    psi_mat <- readRDS(Inputfiles[i])$psi_mat
    for(j in 1:dim(psi_mat)[2]){
      Firms <- psi_mat@i[(psi_mat@p[j]+1):psi_mat@p[j+1]]+1
      TotalDamage <- result$ESRI[j,1]
      TD_Down <- result$ESRI[j,2]
      SumDamage <- sum(OG_ESRI$ESRI[Firms,1])
      SD_Down <- sum(OG_ESRI$ESRI[Firms,2])
      ESRIs[(i-1)*n+j] <- TotalDamage
      RelDiffs[(i-1)*n+j] <- (TotalDamage-SumDamage)/SumDamage
      RelDiffs_Down[(i-1)*n+j] <- (TD_Down-SD_Down)/SD_Down
    }
  }
  return(list(ESRIs = ESRIs,RelDiffs=RelDiffs,RelDiffs_Down=RelDiffs_Down))
}
DoubleCheckScenarios <- function(Inputfiles,RelDiffs,Thresh,Name,Folder){
  Interestingscenarios <- which(RelDiffs > Thresh)
  nscen <- dim(readRDS(Inputfiles[1])$psi_mat)[2]
  totalfirms <- dim(readRDS(Inputfiles[1])$psi_mat)[1]
  Firms <- list()
  for(scenario in Interestingscenarios){
    fileid <- ceiling(scenario/nscen)
    entryid <- scenario %% nscen
    input <- readRDS(Inputfiles[fileid])
    psi_mat <- input$psi_mat
    fs <- which(psi_mat[,entryid]==1)
    Firms <- append(Firms,list(fs))
    newpsi_mat <- Matrix(0,nrow=totalfirms,ncol=length(fs))
    for(i in 1:length(fs)){
      newpsi_mat[fs,i] <- 1
      newpsi_mat[fs[i],i] <- 0
    }
    input$psi_mat <- newpsi_mat
    input$run_id <- paste0(Name,scenario)
    saveRDS(input,paste0(Folder,input$run_id,"_ESRIInput.rds"))
    saveRDS(Firms,paste0(Folder,Name,"Firms.rds"))
  }
}
ConfirmPairs <- function(Outputfiles,Inputfiles,Folder,Name,thresh=0.1){
  Pairs <- list()
  for(i in 1:length(Outputfiles)){
    NewCheck <- readRDS(Inputfiles[i])
    test <- readRDS(Outputfiles[i])
    firmids <- which(((max(test$ESRI[,1]) - test$ESRI[,1])/max(test$ESRI[,1]))>thresh)
    firms <- which(rowSums(as.matrix(NewCheck$psi_mat))!=0)[firmids]
    if(length(firms)>1){Pairs <- append(Pairs,list(firms))}
  }
  NewCheck$psi_mat <- Matrix(0,nrow=dim(NewCheck$psi_mat)[1],ncol=length(Pairs))
  for(i in 1:length(Pairs)){
    NewCheck$psi_mat[Pairs[[i]],i] <- 1
  }
  NewCheck$run_id <- paste0(Name)
  saveRDS(NewCheck,paste0(Folder,Name,"_ESRIInput.rds"))
}
RunESRIs <- function(Folder){
  files <- list.files(Folder,full.names=T)
  for(file in files){
    if(grepl("ESRIInput.rds",file)){
      Input <- readRDS(file)
      if(!paste0(Input$run_id,"_results.rds") %in% files){
        ESRI <- do.call("GL_cascade",args=Input)
        saveRDS(ESRI,paste0(Folder,"/",Input$run_id,"_results.rds"))
      }
    }
  }
}

IDData <- ReadIDData("./data/IDData.csv")
Graph <- ReadRawGraph("./data/EL.csv")
Graph <- ProcessRawGraph(Graph)
ESRIInp <- DataPrepESRI(Graph,IDData,"Single Firm ESRI on random network",Targetfile="./data/OG_Input.rds")
OG_ESRI <- do.call(GLcascade,args=ESRIInp)
saveRDS(OG_ESRI,"./OG_ESRI.rds")
OG_ESRI <- readRDS("./OG_ESRI.rds")
input <- readRDS("./data/OG_Input.rds")
input$ncores <- 50

# Here we create the sets of 500 and 100 firms, and calculate their ESRIs.
set.seed(314159265)
dir.create("Supersamples")
PrepareScenarios(input,OG_ESRI,500,5000,20,"500_Firms_Run_","./Supersamples/")
PrepareScenarios(input,OG_ESRI,100,5000,20,"100_Firms_Run_","./Supersamples/")
RunESRIs("Supersamples")

# This part extracts the small subsets from the large 500 firm sets with a signal
Inputs500 <- grepv("500",list.files("./Supersamples",full.names=T))
Inputs100 <- grepv("100",list.files("./Supersamples",full.names=T))
Results500 <- grepv("500Firms_Run/w*results.rds",list.files("./Supersamples",full.names=T))
Results100 <- grepv("100Firms_Run/w*results.rds",list.files("./Supersamples",full.names=T))

RelDiffs100 <- CalcDiffsFromFiles(Inputfiles100,Resultfiles100,OG_ESRI)
RelDiffs500 <- CalcDiffsFromFiles(Inputfiles500,Resultfiles500,OG_ESRI)

dir.create("Doublecheck")
DoubleCheckScenarios(Inputfiles100,RelDiffs100$RelDiffs,5,"DoubleCheck100_thresh5_","./Doublecheck/")
DoubleCheckScenarios(Inputfiles100,RelDiffs100$RelDiffs,2,"DoubleCheck100_thresh2_","./Doublecheck/")
DoubleCheckScenarios(Inputfiles500,RelDiffs500$RelDiffs,5,"DoubleCheck500_thresh5_","./Doublecheck/")
DoubleCheckScenarios(Inputfiles500,RelDiffs500$RelDiffs,2,"DoubleCheck500_thresh2_","./Doublecheck/")
RunESRIs("Doublecheck")

DoublecheckInputs1002 <- grepv("DoubleCheck10_thresh2\\w*Input.rds",list.files("./Doublecheck/",full.names=T))
DoublecheckInputs1005 <- grepv("DoubleCheck10_thresh5\\w*Input.rds",list.files("./Doublecheck/",full.names=T))
DoublecheckInputs5002 <- grepv("DoubleCheck5_thresh2\\w*Input.rds",list.files("./Doublecheck/",full.names=T))
DoublecheckInputs5005 <- grepv("DoubleCheck5_thresh5\\w*Input.rds",list.files("./Doublecheck/",full.names=T))

Doublecheckresults1002 <- grepv("DoubleCheck100_thresh2\\w*results.rds",list.files("./Doublecheck/",full.names=T))
Doublecheckresults1005 <- grepv("DoubleCheck100_thresh5\\w*results.rds",list.files("./Doublecheck/",full.names=T))
Doublecheckresults5002 <- grepv("DoubleCheck500_thresh2\\w*results.rds",list.files("./Doublecheck/",full.names=T))
Doublecheckresults5005 <- grepv("DoubleCheck500_thresh5\\w*results.rds",list.files("./Doublecheck/",full.names=T))

dir.create("Confirm")
ConfirmPairs(Doublecheckresults1002,DoublecheckInputs1002,"./Confirm/","Pairs1002_thresh01",thresh=0.1)
ConfirmPairs(Doublecheckresults1005,DoublecheckInputs1005,"./Confirm/","Pairs1005_thresh01",thresh=0.1)
ConfirmPairs(Doublecheckresults5002,DoublecheckInputs5002,"./Confirm/","Pairs5002_thresh01",thresh=0.1)
ConfirmPairs(Doublecheckresults5005,DoublecheckInputs5005,"./Confirm/","Pairs5005_thresh01",thresh=0.1)

ConfirmPairs(Doublecheckresults1002,DoublecheckInputs1002,"./Confirm/","Pairs1002_thresh025",thresh=0.25)
ConfirmPairs(Doublecheckresults1005,DoublecheckInputs1005,"./Confirm/","Pairs1005_thresh025",thresh=0.25)
ConfirmPairs(Doublecheckresults5002,DoublecheckInputs5002,"./Confirm/","Pairs5002_thresh025",thresh=0.25)
ConfirmPairs(Doublecheckresults5005,DoublecheckInputs5005,"./Confirm/","Pairs5005_thresh025",thresh=0.25)

ConfirmPairs(Doublecheckresults1002,DoublecheckInputs1002,"./Confirm/","Pairs1002_thresh05",thresh=0.5)
ConfirmPairs(Doublecheckresults1005,DoublecheckInputs1005,"./Confirm/","Pairs1005_thresh05",thresh=0.5)
ConfirmPairs(Doublecheckresults5002,DoublecheckInputs5002,"./Confirm/","Pairs5002_thresh05",thresh=0.5)
ConfirmPairs(Doublecheckresults5005,DoublecheckInputs5005,"./Confirm/","Pairs5005_thresh05",thresh=0.5)
RunESRIs("Confirm")

ConfirmInputs <- grepv("ESRI",list.files("./Confirm",full.names = T))
ConfirmResults <- grepv("result",list.files("./Confirm",full.names = T))

dir.create("Subsets") # Here we check for all sets of n>2 firms, all the possible subsets of n'<n firms.
for(i in 1:length(ConfirmInputs)){
  Input <- readRDS(ConfirmInputs[i])
  RelDiffs <- CalcDiffsFromFiles(ConfirmInputs[i],ConfirmResults[i],OG_ESRI)
  file <- basename(ConfirmInputs[i])
  file <- substr(file,1,nchar(file)-4) #Remove file extension
  base <- dirname(ConfirmResults[i])
  saveRDS(RelDiffs,paste0(base,"/",file,"_RelDiffs.rds"))
  Input$run_id <- paste0(Input$run_id,"_Subsets")
  scenarios <- 0
  for(j in 1:dim(Input$psi_mat)[2]){
    nfirms <- colSums(Input$psi_mat)[j]
    scenarios <- scenarios + dim(hcube(rep(2,nfirms)))[1]-2-nfirms
  }
  newpsi <- Matrix(0,nrow=dim(Input$psi_mat)[1],ncol=scenarios)
  k <- 1
  for(j in 1:dim(Input$psi_mat)[2]){
    firms <- which(Input$psi_mat[,j]!=0)
    if(length(firms)>2){
      for(n in 2:(length(firms)-1)){
        chub <- combn(firms,n)
        for(column in 1:dim(chub)[2]){
          newpsi[chub[,column],k] <- 1
          k <- k+1
        }
      }
    }
  }
  Input$psi_mat <- newpsi
  if(scenarios > 0){
    saveRDS(Input,paste0("./Subsets/",Input$run_id,"_ESRIInput.rds"))
  }
}
RunESRIs("Subsets")

source("./src/RandomPairs.R")
