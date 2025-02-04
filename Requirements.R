### requirements

### Install packages
install.packages("fda")
install.packages("lqmm")
install.packages("mclust")
install.packages("glmnet")
install.packages("gam")
install.packages("quantreg")
install.packages("pdfCluster")
install.packages("arrangements")
install.packages("expectreg")
install.packages("LaplacesDemon")
install.packages("multisensi")
install.packages("scatterplot3d")
install.packages("gridExtra")
install.packages("quadprog")
install.packages("Matrix")
install.packages("sparcl")
install.packages("RSKC")
install.packages("lpSolve")
install.packages("readxl")
install.packages("maps")
install.packages("sn")
install.packages("clustEff")

library(cluster)
library(stringr)
library(maps)
library(readxl)
library(lpSolve)
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
library(sn)
library(clustEff)

### Some functions
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

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

### Some functions

f1 <- function(x) return(3*sin(0.5*pi*x))
f2 <- function(x) return(3*x)
f3 <- function(x) return(10*x*(1-x))


alpha_f <- function(t){
  return(9^(2*t-1)/(1+9^(2*t-1)))
}

alpha_minus_f <- function(t){
  return(1-alpha_f(t))
}

minus_alpha_f <- function(t){
  return(-alpha_f(t))
}

### Generating functional data

curvelist_normal <- function(type, n_breaks = 100, N = 100, sd = 0.1){
  data = list()
  t = seq(0, 1, length.out = n_breaks+1)[1:n_breaks]
  for(i in 1:N){
    if(type == 0) data = c(data, list(rnorm(length(t), 0, sd)))
    else if(type == 1) data = c(data, list(f1(t) + rnorm(length(t), 0, sd)))
    else if(type == 2) data = c(data, list(f2(t) + rnorm(length(t), 0, sd)))
    else if(type == 3) data = c(data, list(f3(t) + rnorm(length(t), 0, sd)))
  }
  return(data)
}

curvelist_alaplace <- function(type, n_breaks = 100, N = 100, sigma = 0.1, alpha = 0.3, hetero = F){
  data = list()
  t = seq(0, 1, length.out = n_breaks+1)[1:n_breaks]
  if(!hetero){ # which means alpha is scalar
    for(i in 1:N){
      lam = sqrt(2)*sigma/sqrt(alpha*(1-alpha))
      kappa = sqrt(alpha/(1-alpha))
      mu = -lam*(1/kappa-kappa)/sqrt(2)
      if(type == 0) newdata = ralaplace(length(t), mu, lam, kappa)
      else if(type == 1) newdata = f1(t) + ralaplace(length(t), mu, lam, kappa)
      else if(type == 2) newdata = f2(t) + ralaplace(length(t), mu, lam, kappa)
      else if(type == 3) newdata = f3(t) + ralaplace(length(t), mu, lam, kappa)
      data = c(data, list(newdata))
    }
  }
  else{
    for(i in 1:N){
      newdata = c()
      for(j in 1:length(t)){
        k = (j-1)/length(t)
        lam = sqrt(2)*sigma/sqrt(alpha(k)*(1-alpha(k)))
        kappa = sqrt(alpha(k)/(1-alpha(k)))
        mu = -lam*(1/kappa-kappa)/sqrt(2)
        if(type == 0) newdata = c(newdata, ralaplace(1, mu, lam, kappa))
        else if(type == 1) newdata = c(newdata, f1(k) + ralaplace(1, mu, lam, kappa))
        else if(type == 2) newdata = c(newdata, f2(k) + ralaplace(1, mu, lam, kappa))
        else if(type == 3) newdata = c(newdata, f3(k) + ralaplace(1, mu, lam, kappa))
      }
      data = c(data, list(newdata))
    }
  }
  return(data)
}

curvelist_skewt <- function(type, n_breaks = 100, N = 100, omega = 0.1, alpha = pi, nu = 3, const = 5, hetero = F){
  # alpha > 0 : right skewed, alpha < 0 : left skewed
  data = list()
  t = seq(0, 1, length.out = n_breaks+1)[1:n_breaks]
  b.nu = sqrt(nu)*gamma((nu-1)/2)/(sqrt(pi)*gamma(nu/2))
  
  if(!hetero){ # which means alpha is scalar
    delta = alpha/sqrt(1+alpha**2)
    xi = -omega*b.nu*delta
    for(i in 1:N){
      if(type == 1) newdata = f1(t) + rst(length(t), dp = c(xi, omega, alpha, nu))
      else if(type == 2) newdata = f2(t) + rst(length(t), dp = c(xi, omega, alpha, nu))
      else if(type == 3) newdata = f3(t) + rst(length(t), dp = c(xi, omega, alpha, nu))
      data = c(data, list(newdata))
    }
  }
  else{
    for(i in 1:N){
      newdata = c()
      for(j in 1:length(t)){
        k = (j-1)/length(t)
        delta = const*alpha(k)/sqrt(1+(const*alpha(k))**2)
        xi = -omega*b.nu*delta
        if(type == 1) newdata = c(newdata, f1(k) + rst(1, dp = c(xi, omega, const*alpha(k), nu)))
        else if(type == 2) newdata = c(newdata, f2(k) + rst(1, dp = c(xi, omega, const*alpha(k), nu)))
        else if(type == 3) newdata = c(newdata, f3(k) + rst(1, dp = c(xi, omega, const*alpha(k), nu)))
      }
      data = c(data, list(newdata))
    }
  }
  return(data)
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

### calculate values at break points
basis_value_point <- function(b.weight, x, num.knots = 10, n_breaks = 100){
  ret = 0
  b.sp = bs(seq(0, 1, by = 1/n_breaks), knots = seq(0, 1, length.out = num.knots+1))
  nb = nrow(b.sp)-1
  for(j in 1:num.knots+3) ret = ret + b.weight[j]*b.sp[round(nb*x)+1, j][-(num.knots+1)]
  return(ret)
}

basis_value_vector <- function(b.weight, num.knots = 10, n_breaks = 100){
  b.sp = bs(seq(0, 1, by = 1/n_breaks), knots = seq(0, 1, length.out = num.knots+1))
  nb = nrow(b.sp)-1
  ret = rep(0, nb)
  for(j in 1:(num.knots+3)) ret = ret + b.weight[j]*b.sp[, j][-(nb+1)]
  return(ret)
}


### Quantile curve
coeff_quantile_vector <- function(curve, tau, num.knots = 10, n_breaks = 100, lambda = 0.1){
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
quant.aocv <- function(curve, tau, num.knots = 10, n_breaks = 100, lambda = 0.1){
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
quant.akcv <- function(curve, tau, num.knots = 10, n_breaks = 100, lambda = 0.1, k = 10){
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
  index <- c(1:(num.knots+1))
  D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1
  
  rand_index = sample(1:n_breaks)
  rand_chunk = chunk(rand_index, k)
  akcv = 0
  
  for(j in 1:k){
    miss_index = rand_chunk[[j]]
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
    for(m in miss_index) akcv = akcv + abs(curve[m]-y.fit[m])*weight[m]
  }
  return(akcv/length(n_breaks))
}

quant.gcv.select <- function(y, tau, num.knots = 10, n_breaks=100, lambda=exp(seq(-2, 4, length=7))){
  gcv.result <- vector(length=length(lambda))
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
  index <- c(1:(num.knots+1))
  D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1
  for(i in 1:length(lambda)){
    lamb <- lambda[i]
    q.fit <- coeff_quantile_vector(y, tau, num.knots = num.knots, n_breaks = n_breaks, lambda=lamb)
    y.fit <- b.basis %*% q.fit$coeff
    W <- diag(as.vector(q.fit$weight))
    S.mat <- b.basis %*% solve(t(b.basis) %*% W %*% b.basis+lamb*t(D)%*%D) %*% t(b.basis) %*% W 
    gcv.result[i] <- sum(q.fit$weight*(y-y.fit)^2)/(length(y)-sum(diag(S.mat)))^2
  }
  #		print(gcv.result)
  return(lambda[which.min(gcv.result)])
}
quant.akcv.select <- function(y, tau, num.knots = 10, n_breaks = 100, lambda = exp(seq(-2, 4, length=7)), k = 10){
  akcv.result <- rep(0, length=length(lambda))
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
  index <- c(1:(num.knots+1))
  D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1
  for(i in 1:length(lambda)){
    lamb <- lambda[i]
    akcv.result[i] <- quant.akcv(y, tau, num.knots, n_breaks, lamb, k)
  }
  return(lambda[which.min(akcv.result)])
}
quant.curve.pspline <- function(y, tau, num.knots = 10, n_breaks = 100, lambda = 0, k = 10){
  if(lambda == 0) lambda <- quant.akcv.select(y, tau, num.knots = num.knots, n_breaks = n_breaks, k = 10)
  coeff <- coeff_quantile_vector(y, tau, num.knots = num.knots, n_breaks = n_breaks, lambda=lambda)$coeff
  return(coeff)
}


### Expectile curve
coeff_vector <- function(curve, tau, num.knots = 10, n_breaks = 100, lambda = 0.1){
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
gcv.select <- function(y, tau, num.knots = 10, n_breaks = 100, lambda=exp(seq(-2, 4, length=7))){
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
expec.curve.pspline <- function(y, tau, num.knots = 10, n_breaks = 100, lambda = 0){
  if(lambda == 0) lambda <- gcv.select(y, tau, num.knots = num.knots, n_breaks = n_breaks)
  coeff <- coeff_vector(y, tau, num.knots = num.knots, n_breaks = n_breaks, lambda=lambda)$coeff
  return(coeff)
}


### l2-distance
l2_distance <- function(coeff, num.knots = 10, n_breaks = 100){
  b.sp = bs(seq(0, 1, by = 1/n_breaks), knots = seq(0, 1, length.out = num.knots+1))
  funvec = rep(0, n_breaks)
  t = seq(0, 1, 1/n_breaks)[-(n_breaks+1)]
  for(i in 1:(num.knots+3)){
    funvec = funvec + coeff[i]*b.sp[, i][-(n_breaks+1)]
  }
  return(sqrt(sum(funvec**2)))
}

### Clustering with l2-distance
distance_clustering <- function(n_cluster, data, num.knots = 10, n_breaks = 100, print = FALSE){
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
          dist[j] = l2_distance(data[, i] - mu[, cluster_set[j]],
                                num.knots = num.knots, n_breaks = n_breaks)
        }
        new.cluster[i] = which.min(dist)
      }
    }
    
    obj = 0
    for(i in 1:n) obj <- obj + l2_distance(data[, i] - mu[, cluster[i]],
                                           num.knots = num.knots, n_breaks = n_breaks)^2
    if(obj < obj.fun){
      obj.fun <- obj
      final.cluster <- cluster
    }
  }
  
  print(paste(iter, " iterations for l2-distance clustering"))
  if(!print) return(final.cluster)
  else return(list(final.cluster, iter))
}



