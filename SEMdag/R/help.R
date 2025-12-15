# Helper functions for SEMdag paper

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

complete_matrix <- function(mat, ref) {
  dif <- setdiff(rownames(ref), rownames(mat))
  mat <- rbind(mat, matrix(0, length(dif), ncol(mat), dimnames = list(dif, NULL)))
  mat <- cbind(mat, matrix(0, nrow(mat), length(dif), dimnames = list(NULL, dif)))
  return(mat)
}

hammingDist <- function(G1,G2, allMistakesOne = TRUE)
  # Computes Hamming Distance between DAGs G1 and G2 with SHD(->,<-) = 1 if allMistakesOne == TRUE
  # INPUT:  G1, G2  adjacency graph containing only zeros and ones: (i,j)=1 means edge from X_i to X_j.
  # OUTPUT: hammingDis Hamming Distance between G1 and G2
  # Copyright (c) 2012-2013  Jonas Peters [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms.
{
  if(allMistakesOne) {
    Gtmp <- (G1+G2)%%2
    Gtmp <- Gtmp + t(Gtmp)
    nrReversals <- sum(Gtmp == 2)/2
    nrInclDel <- sum(Gtmp == 1)/2
    hammingDis <- nrReversals + nrInclDel
  }else{        
    hammingDis <- sum(abs(G1 - G2))
    # correction: dist(-,.) = 1, not 2
    hammingDis <- hammingDis - 0.5*sum(G1 * t(G1) * (1-G2) * t(1-G2) + G2 * t(G2) * (1-G1) * t(1-G1))
  }    
  return(hammingDis)
}

# gg=dag_ig; df_train=df_train_all; group_train=group_train
# df_test=df_test_all; group_test=group_test

run_fda <- function(gg, df_train, group_train, df_test, group_test, ...)
{
	all_res <- c()

	 for (g in 1:length(gg)){ #g=1
	  print(names(gg[g]))
	  try(res<- quiet(predictY(graph=gg[[g]], data=df_train, group=group_train, 
								newdata=df_test, newgroup=group_test, 
								source = FALSE, verbose=FALSE)))
	  yobs<- res$yobs
	  ypred<- res$ypred
	  CT <- table(yobs, ypred)

	  dag<- gg[[g]]
	  din<- igraph::degree(dag, mode = "in")
	  root<- V(dag)$name[din == 0]
	  
	  res1 <- quiet(performance(CT=CT))
	  res2<- data.frame(
		vcount= vcount(dag),
		ecount= ecount(dag),
		root= length(root),
		rec= res1[1],
		pre= res1[2],
		f1= res1[3],
		mcc= res1[4]
	  )
	  all_res <- rbind(all_res, res2)
	  rownames(all_res)[g] <-names(gg)[g]
	 }

	return(all_res)
}

run_rf <- function(df_train, group_train, df_test, group_test, ...)
{
	all_res <- c()
	df_all <- rbind(df_train, df_test)
	group_all <- c(group_train, group_test)
	train <- 1:nrow(df_train)
	 
	res<- CMA::rfCMA(X=df_all, y=group_all, learnind=train, models = TRUE)
	yobs<- res@y
	ypred<- res@yhat
	ypred<- factor(ypred, levels=c(0,1))
	yobs<- factor(yobs, levels=c(0,1))
	CT <- table(yobs, ypred) + 0.001
	res1<- quiet(performance(CT=CT))
	res2<- data.frame(
	   vcount= ncol(df_all),
	   ecount= 0,
	   root= 0,
	   rec= res1[1],
	   pre= res1[2],
	   f1= res1[3],
	   mcc= res1[4]
	 )
	all_res <- rbind(all_res, res2)

	return(all_res)
}

quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    tmpf <- tempfile()
    sink(tmpf)
    on.exit({sink(); file.remove(tmpf)})
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}
