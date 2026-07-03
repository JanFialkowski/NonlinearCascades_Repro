input <- readRDS("./data/OG_Input.rds")
input$psi_mat <- NULL
input$track_h <- NULL
input$run_id <- "RandomPairs"
OG_ESRI <- readRDS("./OG_ESRI.rds")
Good_Firms <- which(!((OG_ESRI$ESRI[,1]>0.1) | OG_ESRI$ESRI[,2]==0))
n <- 2000000
j <- rep(1:n,each=2)
i <- sample(Good_Firms,2*n,replace=T)
while(any(i[2:(2*n)]-i[1:(2*n-1)]==0)){ # A crime is when a firm is shocked twice in a single scenario
  crimes <- which(i[2:(2*n)]-i[1:(2*n-1)]==0)
  i[crimes] <- sample(Good_Firms,length(crimes),replace=T)
}
input$psi_mat <- sparseMatrix(i=i,j=j,x=1,dims=c(dim(OG_ESRI$ESRI)[1],n))
saveRDS(input,"./data/RandomPairs_ESRIInput.rds")

Randomout <- do.call(GLcascade,args=input)
saveRDS(Randomout,"./data/RandomResults.rds")