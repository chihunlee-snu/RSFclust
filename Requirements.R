### requirements


### Install packages
# install.packages("fda")
# install.packages("lqmm")
# install.packages("mclust")
# install.packages("glmnet")
# install.packages("gam")
# install.packages("quantreg")
# install.packages("pdfCluster")
# install.packages("arrangements")
# install.packages("expectreg")
# install.packages("LaplacesDemon")
# install.packages("multisensi")
# install.packages("scatterplot3d")
# install.packages("gridExtra")
# install.packages("quadprog")
# install.packages("Matrix")
# install.packages("sparcl")
# install.packages("RSKC")
# install.packages("https://cran.r-project.org/src/contrib/Archive/Funclustering/Funclustering_1.0.2.tar.gz", repos = NULL, type="source")
# install.packages("lpSolve")
# install.packages("readxl")
# install.packages("maps")

library(maps)
library(readxl)
library(lpSolve)
library(Funclustering)
library(RSKC)
library(sparcl)
library(Matrix)
library(quadprog)
library(gridExtra)
library(scatterplot3d)
library(multisensi)
library(LaplacesDemon)
library(expectreg)
library(arrangements)
library(pdfCluster)
library(quantreg)
library(fda)
library(lqmm)
library(mclust)
library(glmnet)
library(gam)
library(lubridate)

### Some functions
mquant <- function(t, tau, M = 2){
  if(t>=0) return(tau*abs(t)^M)
  else return((1-tau)*abs(t)^M)
}

permute_c_true <- function(c_true){
  n_cluster = length(unique(c_true))
  per = permutations(n_cluster, n_cluster)
  ret = matrix(nrow = nrow(per), ncol = length(c_true))
  for(i in 1:nrow(per)){
    for(j in 1:length(c_true)){
      ret[i, j] = per[i, c_true[j]]
    }
  }
  return(ret)
}

### Performance measure
perform_score <- function(c_result, c_true, measure){
  
  if(length(c_result) != length(c_true)){
    print("clustering error(type : length) ")
    return(NA)
  }
  n = length(c_result)  # data 개수
  K = length(unique(c_true))   # cluster 개수
  K_data = length(unique(c_result))
  
  permute_matrix = permute_c_true(c_true)
  
  TP = matrix(nrow = nrow(permute_matrix), ncol = K)
  FP = matrix(nrow = nrow(permute_matrix), ncol = K)
  FN = matrix(nrow = nrow(permute_matrix), ncol = K)
  
  for(i in 1:nrow(permute_matrix)){
    ## calculate TP_k, FP_k, FN_k
    for(k in 1:K){
      true = (permute_matrix[i, ] == rep(k, n))  # type : boolean
      false = (permute_matrix[i, ] != rep(k, n))
      pos = as.numeric(c_result == rep(k, n))
      neg = as.numeric(c_result != rep(k, n))
      
      TP[i, k] = sum(pos[true])
      FP[i, k] = sum(pos[false])
      FN[i, k] = sum(neg[true])
    }
  }
  
  if(measure %in% c("CSR", "csr")){
    return(max(rowSums(TP)/n))
  }
  else if(measure %in% c("ARI", "ari")){
    return(adj.rand.index(c_result, c_true))
  }
  else if(measure %in% c("MAP", "map")){
    return(max(rowMeans(TP/(TP+FP), na.rm = TRUE)) * K_data / K)
  }
  else if(measure %in% c("MAF", "maf")){
    return(max(rowMeans(2*TP/(2*TP+FP+FN))))
  }
  else{
    print("error(type : non existing performance measure) ")
    return(NA)
  }
  
}


num.knots = 10
n_breaks = 100
b.sp0 = bs(seq(0, 1, by = 1/100), knots = seq(0, 1, length.out = 11))

### calculate values at break points
basis_value_point <- function(b.weight, x, b.sp = b.sp0){
  ret = 0
  nb = nrow(b.sp)-1
  for(j in 1:num.knots+3) ret = ret + b.weight[j]*b.sp[round(nb*x)+1, j][-(num.knots+1)]
  return(ret)
}

basis_value_vector <- function(b.weight, b.sp = b.sp0){
  nb = nrow(b.sp)-1
  ret = rep(0, nb)
  for(j in 1:(num.knots+3)) ret = ret + b.weight[j]*b.sp[, j][-(nb+1)]
  return(ret)
}


### Quantile curve
coeff_quantile_vector <- function(curve, tau, n_breaks = 100, lambda = 0.1, num.knots = 10){
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
  index <- c(1:(num.knots+1))
  D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1
  y = curve
  A1 = matrix(0, nrow = (num.knots+3), ncol = (length(y) + num.knots + 3))
  A1[1:(num.knots+3), 1:(num.knots+3)] = diag(num.knots+3)
  A2 = matrix(0, nrow = length(y), ncol = length(y) + num.knots + 3)
  A2[1:length(y), (num.knots+4):(num.knots+3+length(y))] = diag(length(y))
  dvec = (tau-0.5)*t(A1)%*%t(b.basis)%*%rep(1, length(y)) - 0.5*t(A2)%*%rep(1, length(y))
  Dmat = 2*lambda*t(A1)%*%t(D)%*%D%*%A1
  Dmat = as.matrix(nearPD(Dmat)$mat)
  Amat = rbind(A2+b.basis%*%A1, A2-b.basis%*%A1)
  bvec = c(y, -y)
  
  Q = solve.QP(Dmat, dvec, t(Amat), bvec)
  coeff = Q$solution[1:(num.knots+3)]
  y.fit = b.basis%*%coeff
  weight = tau*(y>=y.fit)+(1-tau)*(y<y.fit)
  coeff <- list(coeff=coeff, weight=weight)
  return(coeff)
}
quant.aocv <- function(curve, tau, n_breaks = 100, lambda = 0.1, num.knots = 10){
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
  index <- c(1:(num.knots+1))
  D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1
  aocv = 0
  for(j in 1:length(y)){
    y = curve[-j]
    b.basis.j = b.basis[-j, ]
    A1 = matrix(0, nrow = (num.knots+3), ncol = (length(y) + num.knots + 3))
    A1[1:(num.knots+3), 1:(num.knots+3)] = diag(num.knots+3)
    A2 = matrix(0, nrow = length(y), ncol = length(y) + num.knots + 3)
    A2[1:length(y), (num.knots+4):(num.knots+3+length(y))] = diag(length(y))
    dvec = (tau-0.5)*t(A1)%*%t(b.basis.j)%*%rep(1, length(y)) - 0.5*t(A2)%*%rep(1, length(y))
    Dmat = 2*lambda*t(A1)%*%t(D)%*%D%*%A1
    Dmat = as.matrix(nearPD(Dmat)$mat)
    Amat = rbind(A2+b.basis.j%*%A1, A2-b.basis.j%*%A1)
    bvec = c(y, -y)
    
    Q = solve.QP(Dmat, dvec, t(Amat), bvec)
    coeff = Q$solution[1:(num.knots+3)]
    y.fit = b.basis%*%coeff
    weight = tau*(curve>=y.fit)+(1-tau)*(curve<y.fit)
    aocv = aocv + abs(curve[j]-y.fit[j])*weight[j]
  }
  return(aocv/length(y))
}
quant.aocv.modified <- function(curve, tau, n_breaks = 100, lambda = 0.1, num.knots = 10){
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
  index <- c(1:(num.knots+1))
  D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1
  
  miss_index = seq(10, 90, 10)
  aocv = 0
  y = curve[-miss_index]
  b.basis.j = b.basis[-miss_index, ]
  A1 = matrix(0, nrow = (num.knots+3), ncol = (length(y) + num.knots + 3))
  A1[1:(num.knots+3), 1:(num.knots+3)] = diag(num.knots+3)
  A2 = matrix(0, nrow = length(y), ncol = length(y) + num.knots + 3)
  A2[1:length(y), (num.knots+4):(num.knots+3+length(y))] = diag(length(y))
  dvec = (tau-0.5)*t(A1)%*%t(b.basis.j)%*%rep(1, length(y)) - 0.5*t(A2)%*%rep(1, length(y))
  Dmat = 2*lambda*t(A1)%*%t(D)%*%D%*%A1
  Dmat = as.matrix(nearPD(Dmat)$mat)
  Amat = rbind(A2+b.basis.j%*%A1, A2-b.basis.j%*%A1)
  bvec = c(y, -y)
  
  Q = solve.QP(Dmat, dvec, t(Amat), bvec)
  coeff = Q$solution[1:(num.knots+3)]
  y.fit = b.basis%*%coeff
  weight = tau*(curve>=y.fit)+(1-tau)*(curve<y.fit)
  for(j in miss_index) aocv = aocv + abs(curve[j]-y.fit[j])*weight[j]
  return(aocv/length(miss_index))
}
quant.gcv.select <- function(y, tau, lambda=exp(seq(-2, 4, length=7)), n_breaks=100, num.knots = 10){
  gcv.result <- vector(length=length(lambda))
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
  index <- c(1:(num.knots+1))
  D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1
  for(i in 1:length(lambda)){
    lamb <- lambda[i]
    q.fit <- coeff_quantile_vector(y, tau, n_breaks, lambda=lamb, num.knots = num.knots)
    y.fit <- b.basis %*% q.fit$coeff
    W <- diag(as.vector(q.fit$weight))
    S.mat <- b.basis %*% solve(t(b.basis) %*% W %*% b.basis+lamb*t(D)%*%D) %*% t(b.basis) %*% W 
    gcv.result[i] <- sum(q.fit$weight*(y-y.fit)^2)/(length(y)-sum(diag(S.mat)))^2
  }
  #		print(gcv.result)
  return(lambda[which.min(gcv.result)])
}
quant.aocv.select <- function(y, tau, lambda=exp(seq(-2, 4, length=7)), n_breaks=100, num.knots = 10){
  gcv.result <- vector(length=length(lambda))
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
  index <- c(1:(num.knots+1))
  D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1
  for(i in 1:length(lambda)){
    lamb <- lambda[i]
    gcv.result[i] <- quant.aocv.modified(y, tau, n_breaks, lamb, num.knots)
  }
  #		print(gcv.result)
  return(lambda[which.min(gcv.result)])
}
quant.curve.pspline <- function(y, tau, n_breaks = 100, num.knots = 10, forcompute = F, gcv = F){
  if(!forcompute & !gcv) lambda <- quant.aocv.select(y, tau, n_breaks = n_breaks, num.knots = num.knots)
  else if(!forcompute) lambda <- quant.gcv.select(y, tau, n_breaks = n_breaks, num.knots = num.knots)
  else lambda <- 0.1
  coeff <- coeff_quantile_vector(y, tau, n_breaks, lambda=lambda, num.knots = num.knots)$coeff
  return(coeff)
}


### Expectile curve
coeff_vector <- function(curve, tau, n_breaks = 100, lambda = 0.1, num.knots = 10){
  D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
  index <- c(1:(num.knots+1))
  D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1
  
  y = curve
  weight.new <- rep(tau, n_breaks); weight <- rep(0, n_breaks)
  iter <- 1
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  
  while(sum((weight-weight.new)^2)>1e-6){
    weight <- weight.new
    # print(paste(iter, "th iteration"))
    W <- diag(weight)
    coeff <- solve(t(b.basis) %*% W %*% b.basis+lambda*t(D)%*%D, t(b.basis) %*% W %*% y)
    iter <- iter + 1 
    y.fit <- b.basis %*% coeff
    weight.new <- as.vector(tau * (y > y.fit) + (1-tau) * (y <= y.fit))
    if(iter == 1000){
      print(paste("procedure does not converge"))
      break
    }
  }
  coeff <- list(coeff=coeff, weight=weight)
  return(coeff)
}
gcv.select <- function(y, tau, lambda=exp(seq(-2, 4, length=7)), n_breaks=100, num.knots = 10){
  gcv.result <- vector(length=length(lambda))
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
  index <- c(1:(num.knots+1))
  D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1
  for(i in 1:length(lambda)){
    lamb <- lambda[i]
    exp.fit <- coeff_vector(y, tau, n_breaks, lambda=lamb, num.knots = num.knots)
    y.fit <- b.basis %*% exp.fit$coeff
    W <- diag(exp.fit$weight)
    S.mat <- b.basis %*% solve(t(b.basis) %*% W %*% b.basis+lamb*t(D)%*%D) %*% t(b.basis) %*% W 
    gcv.result[i] <- sum(exp.fit$weight*(y-y.fit)^2)/(length(y)-sum(diag(S.mat)))^2
  }
  #		print(gcv.result)
  return(lambda[which.min(gcv.result)])
}
expec.curve.pspline <- function(y, tau, n_breaks = 100, num.knots = 10, forcompute = F){
  if(!forcompute) lambda <- gcv.select(y, tau, n_breaks = n_breaks, num.knots = num.knots)
  else lambda <- 0.1
  coeff <- coeff_vector(y, tau, n_breaks = n_breaks, lambda=lambda, num.knots = num.knots)$coeff
  return(coeff)
}


### l2-distance
l2_distance <- function(coeff){
  funvec = rep(0, n_breaks)
  t = seq(0, 1, 1/n_breaks)[-(n_breaks+1)]
  for(i in 1:num.knots){
    funvec = funvec + coeff[i]*b.sp[, i][-(n_breaks+1)]
  }
  return(sqrt(sum(funvec**2)))
}

### Clustering with l2-distance
distance_clustering <- function(n_cluster, data, print = FALSE){
  # data : (data dimension) * (amount of data)
  n = ncol(data)
  a = rep(1/n_cluster, n_cluster-1)
  cluster = rep(0, n)
  new.cluster = sample(1:n_cluster, size=n, rep=T, prob=c(a, 1-sum(a)))
  cluster_num = 5
  iter = 0
  final.cluster = new.cluster
  obj.fun = 100000
  
  for(z in 1:cluster_num){
    new.cluster = sample(1:n_cluster, size=n, rep=T, prob=c(a, 1-sum(a)))
    error_count = 0  # cluster 개수가 n_cluster보다 작아지면 new.cluster를 다시 sampling
    error_bound = 10
    
    while(sum((cluster-new.cluster)^2) != 0){
      # print(new.cluster)
      iter = iter + 1
      cluster = new.cluster
      mu = matrix(nrow = nrow(data), ncol = n_cluster)
      
      cluster_set <- 1:n_cluster
      
      if(length(unique(cluster)) < n_cluster & error_count < error_bound){
        cluster = sample(1:n_cluster, size=n, rep=T, prob=c(a, 1-sum(a)))
        error_count <- error_count + 1
      }
      
      if(length(unique(cluster)) < n_cluster & error_count == error_bound){
        cluster_set <- unique(cluster)
      }
      
      
      for(k in cluster_set){
        if(sum(cluster==k) > 1) mu[, k] = rowMeans(data[, cluster == k])
        else mu[, k] = data[, cluster == k]
      }
      
      new.cluster = vector(length = n)
      for(i in 1:n){
        dist = vector(length = length(cluster_set))
        for(j in 1:length(cluster_set)){
          dist[j] = l2_distance(data[, i] - mu[, cluster_set[j]])
        }
        new.cluster[i] = which.min(dist)
      }
    }
    
    obj = 0
    for(i in 1:n) obj <- obj + l2_distance(data[, i] - mu[, cluster[i]])^2
    if(obj < obj.fun){
      obj.fun <- obj
      final.cluster <- cluster
    }
  }
  
  print(paste(iter, " iterations for l2-distance clustering"))
  if(!print) return(final.cluster)
  else return(list(final.cluster, iter))
}


y = curvelist_normal(N = 1)[[1]]
quant.curve.pspline(y, 0.5)
expec.curve.pspline(y, 0.5)
