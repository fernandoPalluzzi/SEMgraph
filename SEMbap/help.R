### SIMULATION DESIGN (Helper functions)

#Sest=solve(S); n=nrow(X); p=ncol(X)

generate.data <- function(Sest, n, p, ...)
{
	#Sest = the covariance matrix = solve(lrpsadmm$S)
	fake.data <- mvtnorm::rmvnorm(n, sigma = Sest)
	e <- eigen(Sest)
	sqrt.true.cov.mat <- e$vectors%*%sqrt(diag(e$values))
	
	samp.cov.mat <- cov(fake.data)
	if(!corpcor::is.positive.definite(samp.cov.mat)){
	 samp.cov.mat<- corpcor::cov.shrink(samp.cov.mat, verbose = TRUE)
	}
	
	e <- eigen(samp.cov.mat)
	sqrt.samp.cov.mat <- e$vectors%*%sqrt(diag(e$values))
	fake.data <- t(sqrt.true.cov.mat%*%solve(sqrt.samp.cov.mat,t(fake.data)))
	fake.data <- as.data.frame(fake.data)
	
	return(fake.data)
}

#graph = g0 #image(Beta); V(g0)$name

generateBeta <- function(graph, ...) 
{
	ftm<- as_data_frame(graph)
	
	num_edges <- ecount(graph)
	Betas <- matrix(runif(num_edges,0.1,1),nrow=num_edges)
	Sign2 <- matrix(rbinom(num_edges,1,0.5),nrow=num_edges)
	Betas[Sign2==1] <- -Betas[Sign2==1]
	
	ftm$weight <- Betas
	ig <- graph_from_data_frame(ftm)
	Beta <- as_adj(ig, attr = "weight", sparse = FALSE)
	
	return(Beta)
}

#dag=g0; size = 80; p = 32 #image(Omega)

generateBAP <- function(dag, size, p, ...)
{
	A1 <- diag(p)
	ri <- sample(seq(from = 1, to = p),
		   size = size,
		   replace = TRUE)
	index <- abs(seq(
		     from = -size,
		     to = -1,
		     size = size
	))
	rj <- ri[order(index)]
	for (k in 1:size) {
		i <- ri[k]
		j <- rj[k]
		A1[i, j] <- 1
		A1[j, i] <- 1
	}
	A2<- as_adj(as.undirected(dag), sparse=FALSE) 
	Omega<- ifelse(A2 == 1, 0, ifelse (A1 == 1, 1, 0))
	
	return(Omega)
}

#dag = g0; p = vcount(dag); nei = 5 #image(Omega)

generateSW <- function(dag, p=vcount(dag), nei, ...)
{
	guu<- igraph::sample_smallworld(1, size=p, nei=nei, p=0.9)
	A1<- as_adj(simplify(guu), sparse=FALSE) # image (A1)
	A2<- as_adj(as.undirected(dag), sparse=FALSE) # image (A2)
	Omega<- ifelse(A2 == 1, 0, ifelse (A1 == 1, 1, 0)) #image(Omega)
	diag(Omega) <- 1
	
	return(Omega)
}

#Omega=adj; O=c(0.2,0.7) #image(Omega)

generateOmega <- function(Omega, O=c(0,1), ...)
{
	# Set covariances to uniform random numbers
	ind <- lower.tri(Omega) & (Omega==1)
	Omega[ind] <- runif(length(which(ind)),O[1],O[2])
	Omega[upper.tri(Omega)] <- t(Omega)[upper.tri(Omega)]
	
	# Set variances to rowsum of abs values plus U(0,1)
	#diag(Omega) <- rowMeans(abs(Omega)) + runif(nrow(Omega),0.1,0.9)
	Ojj <- diag(runif(nrow(Omega), 0.1, 0.9))
	Omega <- Omega - Ojj #diag(Omega)

	#d <- svd(Omega, nu=0, nv=0)$d
	#Omega <- Omega + (0.1 - min(d))*diag(nrow(Omega))
	#Omega<- corpcor::cor.shrink(Omega, verbose=TRUE)
	#print(corpcor::is.positive.definite(Omega))
	
	return(Omega)
}

# W=pre; W=cov

compute.SqrtW <- function(W, pre = TRUE, ...)
{
	# Eigenvalues-eigenvectors of W
	if (!corpcor::is.positive.definite(W)) {
	 W <- corpcor::cov.shrink(W, verbose = FALSE)
	}
	E<- eigen(W)
	if (pre) {
	 R<- E$vectors%*%diag(sqrt(E$values))%*%t(E$vectors)
	 #sum(w - R %*% R)
	}else{
	 R<- E$vectors%*%diag(1/sqrt(E$values))%*%t(E$vectors)
	 #sum(solve(w) - R %*% R)
	}
	return(R)
}	     

#true.mat = A0; est.mat = A1

compute.metrics <- function(true.mat, est.mat, ...)
{
	prec <- sum(true.mat * est.mat) / sum(est.mat)
	rec <- sum(true.mat * est.mat) / sum(true.mat)
	f1 <- 2*(prec*rec)/(prec+rec)
	tpr <- rec
	fpr <- 1 - (sum((true.mat ==0) * (est.mat ==0)) / sum(true.mat == 0))
	return(list(prec = prec, rec = rec, f1 = f1, tpr = tpr, fpr = fpr))
}

estimate_latent_count <- function(X, method="auto", seed = 1, ...)
{
	n <- nrow(X)
	p <- ncol(X)
	r <- min(n-1,p)
	svdX <- svd(scale(X), nu=0, nv=r)
		
	if (method == "auto") {
	# Parallel analys by N permutation (Dobriban, 2020)
	 permutation_thresholding <- function(X, quantile = 0.99, seed = 1)
	 {
		N <- 50
		X <- scale(X)
		X <- X[, !is.na(apply(X, 2, sum))]
		X <- X[, sort(apply(X, 2, sd),
					index.return = TRUE,
					decreasing = TRUE)$ix[seq_len(min(1000, ncol(X)))]]
		r <- min(dim(X))
		
		evals <- matrix(0, nrow = N, ncol = r)
		set.seed(seed)
		for (i in seq_len(N)) {
		 X_perm <- apply(X, 2, function(xx) sample(xx))
		 evals[i, ] <- svd(X_perm, nu = 0, nv = 0)$d[seq_len(r)]
		}
		thresholds <- apply(
		 evals, 2,
		 function(xx) quantile(xx, probs = quantile)
		)

		# limit number of confounders to at most 10 LVs
		limit <- ifelse(r > 10, 10, r)
		# last which crosses the threshold
		crosses <- (svd(X, nu = 0, nv = 0)$d > thresholds)#[1:limit]
		return(max(which(c(TRUE, crosses)) - 1))
	 }
	 q <- permutation_thresholding(X, seed = seed)
	 idx <- ceiling(q)
	}

	return(idx)
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
