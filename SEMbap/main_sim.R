### SIMULATION DESIGN (SEMbap)

#load libraries
rm(list = ls(all = TRUE)); graphics.off(); cat("\014")
#setwd("~/Desktop/SEMgraph/SEMbap")
library(lrpsadmm)
library(mvtnorm)
library(SEMgraph)
#library(SEMdata)

source("help.R")
options(warn = -1)

getwd()


### Generate results' table

simulation <- "sparse"
#simulation <- "dense"
if (simulation == "sparse") {
  lv_vec <- c("HDLVS_sporadic", "HDLVS_interconnected", "DAG")
} else{
  lv_vec <- c("1LV_all", "3LVS_cluster", "3LVS_over")
}

graph <- c("small", "large")
#seed_vec <- c(1:100)
seed_vec <- 1

methods <- c(
  "simulated",
  "dsep_cggm",
  "dsep_glpc",
  "glasso",
  "d_pcss",
  "d_trim",
  "d_lava",
  "d_rsvp",
  "dce",
  "efa",
  "lrpsadmm"
)


nset <- length(graph)*length(lv_vec) # design gxv
nmet <- length(methods) # n. of methods
k <- length(graph) #n. of graphs
l <- length(seed_vec)
 
results <- data.frame(
  iter = rep(seed_vec, each = nset * nmet),
  design = rep(rep(lv_vec, each = nmet * k), l),
  graph = rep(rep(graph, each = nmet), l),
  method = rep(methods, nset * l),
  K = rep(NA, nset * nmet * l),
  srmr = rep(NA, nset * nmet * l),
  node_act = rep(NA, nset * nmet * l),
  node_inh = rep(NA, nset * nmet * l),
  nlog10P = rep(NA, nset * nmet * l),
  vcountP = rep(NA, nset * nmet * l),
  prec = rep(NA, nset * nmet * l),
  rec = rep(NA, nset * nmet * l),
  f1 = rep(NA, nset * nmet * l),
  tpr = rep(NA, nset * nmet * l),
  fpr = rep(NA, nset * nmet * l)
  #PDO = rep(NA, nset * nmet * l),
  #PDR = rep(NA, nset * nmet * l)
)
dim(results); colnames(results)


# Initialise graph (dense or sparse)
small <- upgrade_graph(alsData$graph); is.dag(small) # TRUE
j <- which(names(kegg.pathways) == "Amyotrophic lateral sclerosis");j
ig <- upgrade_graph(kegg.pathways[[j]])
ig <- properties(ig)[[1]]
large <- ig - E(ig)[which_mutual(ig)]; is.dag(large) #TRUE

# threshold factor loadings (Widaman, SEM 2018; 25: 829-847)
# High l=0.8; l^2=0.64
# Mod  l=0.6; l^2=0.36
# Low  l=0.4: l^2=0.16

# threshold Beta* (???)
# High abs(B)=0.8
# Mod  abs(B)=0.5
# Low  abs(B)=0.2


# Run simulations

start <- Sys.time()
results <-
  run_sim(
    simulation,
    seed_vec,
    vec = lv_vec,
    graph = graph,
	methods = methods,
	n_size = 100
  )
end <- Sys.time()
print(end - start)

results[1:33,]

write.csv(results, file="xxx.csv")


run_sim <- function(simulation, seed_vec, vec, graph, methods, n_size, ...)
{
  # simulation: type of simulation (dense/sparse)
  # seed_vec: sequence of numbers of iterations
  # vec: type of simulation design (subcategory of dense/sparse)
  # graph: initial KEGG graphs (small/large)
  # methods: running methods
  # n_size: number of samples

  idx <- 1

  for (j in seed_vec) {
    #j=1
	
	# set seed and result objects = NULL
	set.seed(j)
	X=Y1=Y2=Y3=Y4=Y5=Y6=Y7=Y8=Y9=Y10=NULL
	A0=A1=A2=A3=A4=A5=A6=A7=A8=A9=A10=NULL
	g0=gH2=gH8=gH9=NULL
	
    for (q in vec) {
      #q="HDLVS_sporadic"
 
      for (g in graph) {
        #g="small"
        
        cat("\n", j, "/", q, "/", g)

        if (g == "small") g0 <- small else g0 <- large
        
        # Initialise matrix of coefficients (B)
        #B <- 0.5 * as_adj(g0, attr = "weight", sparse = FALSE)# row:X, col:Y
        B <- generateBeta(g0)
		B <- B[V(g0)$name, V(g0)$name]
        p <- nrow(B)
        n0 <- 0.25 * n_size
        n1 <- n_size - n0
        md <- round(p / 6)
        mu <- runif(md, 0.05, 0.75)
        group <- c(rep(0, n0), rep(1, n1))

        if (g == "small") size <- 80 else size <- 800
        if (simulation == "sparse"){
          if (g == "small") alpha <- 5E-2 else alpha <- 5E-3#sparse
        } else{
          if (g == "small") alpha <- 5E-3 else alpha <- 5E-4#dense
        }
        if (g == "small") thrL <- 0.16 else thrL <- 0.16 #sqrt(0.16)=0.4
        if (g == "small") thrW <- c(0.1, 0.05) else thrW <- c(0.1, 0.05)


        # Generate Omega
        
        if (simulation == "sparse") {

            if (q == "HDLVS_sporadic") {
              adj <- generateBAP(dag = g0, size = size, p = p)
              Omega0 <- generateOmega(adj, O=c(0,1))
			  colnames(Omega0) <- rownames(Omega0) <- colnames(B)
            } 

            else if (q == "HDLVS_interconnected") {
              adj <- generateSW(dag = g0, p = p, nei = 5)
              Omega0 <- generateOmega(adj, O=c(0,1))
			  colnames(Omega0) <- rownames(Omega0) <- colnames(B)
            }

            else if (q == "DAG") {
              adj <- diag(p)
			  Omega0 <- diag(runif(p, 0.1, 0.9))
              rownames(Omega0) <- colnames(Omega0) <- colnames(B)
            }

        } else if (simulation == "dense") {

            if (q == "1LV_all") {
              Lambda <- rep(1, p)
              adj <- Lambda %*% t(Lambda)
              Omega0 <- generateOmega(adj, O = c(0.64, 0.81))
			  colnames(Omega0) <- rownames(Omega0) <- colnames(B)
            }

            else if (q == "3LVS_cluster") {
              #a<- 10; b<- 6; c<- 11; d <- 21
              a <- round(p / 3)
              b <- 0
              c <- a + 1
              d <- 2 * a + 1
              l1 <- c(rep(1, a), rep(0, p - a))
              l2 <- c(rep(0, a), rep(1, c), rep(0, p - (a + c)))
              l3 <- c(rep(0, d), rep(1, p - d))
              Lambda <- cbind(l1, l2, l3)
              adj <- Lambda %*% t(Lambda)
              Omega0 <- generateOmega(adj, O = c(0.2, 0.7))
			  colnames(Omega0) <- rownames(Omega0) <- colnames(B)
            }

            else if (q == "3LVS_over") {
              #a<- 15; b<- 8; c<- 15; d <- 13
              a <- round(p / 3)
              b <- round(a / 2)
              c <- a + 1
              d <- 2 * a + 1
              l1 <- c(rep(1, a + b), rep(0, p - (a + b)))
              l2 <- c(rep(0, a), rep(1, c + b), rep(0, b))
              l3 <- c(rep(0, d), rep(1, p - d))
              Lambda <- cbind(l1, l2, l3)[1:p, ]
              adj <- Lambda %*% t(Lambda)
			  adj<- ifelse(adj == 2, 1, adj)
              Omega0 <- generateOmega(adj, O = c(0.2, 0.7))
			  colnames(Omega0) <- rownames(Omega0) <- colnames(B)
            }
        }


        # Generate simulated data X = E*(I-B)^-1
        # mu1 <- c(rep(mu, md), rep(-mu, md), rep(0, p - 2 * md))
        mu1 <- c(mu, -mu, rep(0, p - 2 * md))
		din <- igraph::degree(g0, mode = "in") #din
        diag(Omega0)[which(din == 0)] <- 1
		
		set.seed(j)
        
		E1 <- mvtnorm::rmvnorm(n = n0,
                               mean = rep(0, p),
							   sigma = Omega0) #dim(E1)
		E2 <- mvtnorm::rmvnorm(n = n1,
                               mean = mu1,
                               sigma = Omega0) #dim(E2)
        
		X1 <- E1 %*% solve(diag(p) - B) #dim(X1)
        X2 <- E2 %*% solve(diag(p) - B) #dim(X2)
        X <- rbind(X1, X2) #dim(X)
		#image(cor(X)); image(Omega0); image(adj)
            
		V <- colnames(Omega0)


        ### Sparse
        cat("\n...Sparse")  
				
        #...with CGGM
		if ("dsep_cggm" %in% methods) {
		 if (g == "small") {
          res<- tryCatch(quiet(SEMbap(g0, X, group = NULL, dalgo = "cggm",
                        method="BH", alpha=alpha, cmax=Inf)),
                        error = function(err) NA)			
		 }else{	
		  res<- tryCatch(quiet(SEMbap(g0, X, group = NULL, dalgo = "cggm",
                        method="BH", alpha=alpha, cmax=5)),
                        error = function(err) NA)
         }
		 if (is.list(res)){
            Y1 <- res$data
            g1 <- res$dag
            A1 <- res$adj[V, V]
            # image(A1)
           }else{
            Y1 <- X
            g1 <- g0
            A1 <- matrix(0, p, p)
           }
		}
		
		#...with gLPCA
		if ("dsep_glpc" %in% methods) {
         if (g == "small") {
          res<- tryCatch(quiet(SEMbap(g0, X, group = NULL, dalgo = "glpc",
                        method="BH", alpha=alpha, cmax=Inf)),
                        error = function(err) NA) #str(res, max.level=1)
		 }else{	
		  res<- tryCatch(quiet(SEMbap(g0, X, group = NULL, dalgo = "glpc",
                        method="BH", alpha=alpha, cmax=5)),
                        error = function(err) NA) #str(res, max.level=1)
         }
         if (is.list(res)){
            Y2 <- res$data
			gH2 <- res$dag
            A2 <- res$adj[V, V]
            # image(A2)
           }else{
            Y2 <- X
            gH2 <- g0
            A2 <- matrix(0, p, p)
           }
		}

		#... with glasso
		if ("glasso" %in% methods) { 
		  res<- tryCatch(quiet(SEMbap(g0, X, group = NULL, limit=1)),
                        error = function(err) NA) #str(res, max.level=1)
		 if (is.list(res)){
            Y3 <- res$data
            g3 <- res$dag
            A3 <- res$adj[V, V]
            # image(A3)
           }else{
            Y3 <- X
            g3 <- g0
            A3 <- matrix(0, p, p)
           }
		}
		
        ### SVD
        cat("...SVD")

		if ("d_pcss" %in% methods) { 
         #...with svd (pcss)
         Y4 <- quiet(SEMbap(g0, X, group = NULL, dalgo = "pcss")$data)
         A4 <- NA
        }
		if ("d_trim" %in% methods) { 
         #...with svd (trim)
         Y5 <- SEMbap(g0, X, group = NULL, dalgo = "trim", hcount=0)$data
         A5 <- NA
		}
		if ("d_lava" %in% methods) { 
         #...with svd (lava)
         Y6 <- SEMbap(g0, X, group = NULL, dalgo = "lava", hcount=0)$data
         A6 <- NA
		}
		if ("d_rsvp" %in% methods) { 
         #...with svd (rsvp)
         Y7 <- SEMbap(g0, X, group = NULL, dalgo = "rsvp", hcount=0)$data
         A7 <- NA
		}

        ### Dense
        cat("...Dense")

        # ...with estimated hidden variables (dce=pc)
        if ("dce" %in% methods) { 
		 res <- quiet(SEMbap(g0, X, group = NULL, dalgo = "pc"))
         Y8 <- res$data
         if (!identical(colnames(Y8)[grepl("LV", colnames(Y8))], character(0))) {
          gH8 <- res$dag
          fit <- quiet(SEMrun(gH8, Y8, algo = "ricf", n_rep = 0))
          K <- length(colnames(Y8)[grepl("LV", colnames(Y8))])
          L <- as.matrix(fit$fit$Beta[c(1:K), -c(1:K)])
          LL <- t(L) %*% L
          A8 <- ifelse(abs(LL) > thrL, 1, 0)
          A8 <- A8[V, V]
          diag(A8) <- 0
          # image(A8)
         }else{
          gH8 <- g0
          A8 <- matrix(0, p, p)
         }
		}
		
        # ...with estimated hidden variables (efa=fa)
		if ("efa" %in% methods) { 
         res <- quiet(SEMbap(g0, X, group = NULL, dalgo = "fa"))
         Y9 <- res$data
         if (!identical(colnames(Y9)[grepl("LV", colnames(Y9))], character(0))) {
          gH9 <- res$dag
          fit <- quiet(SEMrun(gH9, Y9, algo = "ricf", n_rep = 0))
          K <- length(colnames(Y9)[grepl("LV", colnames(Y9))])
          L <- as.matrix(fit$fit$Beta[c(1:K), -c(1:K)])
          LL <- t(L) %*% L
          A9 <- ifelse(abs(LL) > thrL, 1, 0)
          A9 <- A9[V, V]
          diag(A9) <- 0
          # image(A9)
         }else{
          gH9 <- g0
          A9 <- matrix(0, p, p)
         }
		}
		
        #...with lrpsadmm
		cat("...lrpsadmm")
		
		if ("lrpsadmm" %in% methods) { 
		 lambda <- sqrt(log(ncol(X)) / nrow(X))
         gamma <- 0.1 # A small value of gamma favours an L with a small rank
         l1 <- lambda * gamma
         l2 <- lambda * (1 - gamma)
         
		 fit <- lrpsadmm::lrpsadmm(
          Sigma = cor(X),
          Lambda1 = l1,
          Lambda2 = l2,
          maxiter = 1000,
          print_progress = FALSE,
          zeros = NULL,
          backend = "RcppEigen"
         )

         S <- fit$S # The estimated Shat (sparse precision matrix)
         L <- fit$L # The estimated Lhat (dense low-rank matrix)
         K <- Matrix::rankMatrix(L)[1] # K
         A10 <- ifelse(abs(L) > 0.05, 1, 0) #fixed thrW=0.05
         colnames(A10) <- rownames(A10) <- colnames(X)
         A10 <- A10[V, V]
         diag(A10) <- 0
         # image(A10)
         Y10 <- quiet(generate.data(
          Sest = solve(S),
          n = nrow(X),
          p = ncol(X)
         ))
         colnames(Y10) <- colnames(X)
		 Y10<- as.matrix(Y10)
		}
				
        # SEM fitting + Classification metrics
        
        Y <- list(
          list(g0, X),
          list(g0, Y1),
          list(gH2, Y2),
          list(g0, Y3),
          list(g0, Y4),
          list(g0, Y5),
          list(g0, Y6),
          list(g0, Y7),
          list(gH8, Y8),
          list(gH9, Y9),
          list(g0, Y10)
		)

        names(Y) <- c(
          "simulated",
          "dsep_cggm",
          "dsep_glpc",
          "glasso",
          "d_pcss",
          "d_trim",
          "d_lava",
          "d_rsvp",
          "dce",
          "efa",
          "lrpsadmm"
        )

		
		#A0 <- ifelse(Omega0 == 0, 0, 1)
		A0 <- adj; diag(A0) <- 0
        A <- list(A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10)
        
	    I <- which(names(Y) %in% methods)
		Y <- Y[I]
		A <- A[I]
		        
        for (i in 1:length(I)) {
          #cat("\n", q, names(Y)[i], "\n")
          if(!is.matrix(Y[[i]][[2]])){ 
           idx <- idx + 1
           next
          }

          fit <- quiet(SEMrun(
              graph = Y[[i]][[1]],
              data = Y[[i]][[2]],
              group = group,
              algo = "ricf"
            ))
          fitIdx <- fit$fit$fitIdx

		  results$srmr[idx] <- fitIdx[3]
          results$node_act[idx] <- fit$pval[1]
          results$node_inh[idx] <- fit$pval[2]
		  results$nlog10P[idx] <- -log10(2*min(fit$pval))
          if (i == 3 | i == 9 | i == 10) {
           fit$graph <- delete_vertices(fit$graph,
           V(fit$graph)$name[grep("LV", V(fit$graph)$name)])
          }
          results$vcountP[idx] <- sum(V(fit$graph)$color != "white")

          metrics <- compute.metrics(true.mat = A0, est.mat = A[[i]])
          results$prec[idx] <- metrics$prec
          results$rec[idx] <- metrics$rec
          results$f1[idx] <- metrics$f1
          results$tpr[idx] <- metrics$tpr
          results$fpr[idx] <- metrics$fpr
          results$K[idx] <- tryCatch(estimate_latent_count(
                                   X %*% (diag(p) - B),
								   method = "auto", seed = j),
								   error = function(err) NA)					   
		  
		  #results$PDO[idx] <- corpcor::is.positive.definite(Omega0)
		  #results$PDR[idx] <- corpcor::is.positive.definite(cor(Y[[i]][[2]]))
		  
		  idx <- idx + 1  #results[idx,]
        }
      }
    }
  }

  return(results)
}

