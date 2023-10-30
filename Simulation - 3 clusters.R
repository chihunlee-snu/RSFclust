### 3 cluster

### Some functions

func1 <- function(x, scale = 10, p = 0.5, ver = 0, leftp = TRUE){
  if(leftp) return(scale*x^(1+p)*(1-x) + ver)
  else return(scale*x*(1-x)^(1+p) + ver)
}

alpha_f <- function(t){
  return(9^(2*t-1)/(1+9^(2*t-1)))
}

alpha_minus_f <- function(t){
  return(1-alpha_f(t))
}

max_f <- function(pdiff, a){
  xmax = 0.5-pdiff/2
  deg = 2*pdiff/(1-pdiff)+1
  return(a*(1-xmax)^deg*xmax)
}


### Generating functional data

curvelist_normal <- function(param = c(10, 0.5, 0, 0, 2.5, TRUE), N = 100, sd = 0.1, type = 1){
  scale = param[1]
  p = param[2]
  ver = param[3]
  hor = param[4]
  peak = param[5]
  leftp = param[6]
  
  data = list()
  t = seq(0, 0.99, 0.01)
  set.seed(123)
  for(i in 1:N){
    if(type == 1) data = c(data, list(func1(t, scale, p, ver, leftp) + rnorm(length(t), 0, sd)))
    else data = c(data, list(rnorm(length(t), 0, sd)))
  }
  
  return(data)
}

curvelist_alaplace <- function(param = c(10, 0.5, 0, 0, 2.5, TRUE), N = 100, sd = 0.1, alpha = 0.3, hetero = F, type = 1){
  scale = param[1]
  p = param[2]
  ver = param[3]
  hor = param[4]
  peak = param[5]
  leftp = param[6]
  sigma = sd
  alpha1 = alpha
  if(!hetero){ # which means alpha is scalar
    for(i in 1:N){
      lam = sqrt(2)*sigma/sqrt(alpha*(1-alpha))
      kappa = sqrt(alpha/(1-alpha))
      mu = -lam*(1/kappa-kappa)/sqrt(2)
      if(type == 1) newdata = func1((0:99)*0.01, scale, p, ver, leftp) + ralaplace(100, mu, lam, kappa)
      else newdata = ralaplace(100, mu, lam, kappa)
      
      if(i == 1) data = list(newdata)
      else data = c(data, list(newdata))
    }
  }
  else{
    for(i in 1:N){
      newdata = c()
      for(j in 1:100){
        k = (j-1)/100
        lam = sqrt(2)*sigma/sqrt(alpha(k)*(1-alpha(k)))
        kappa = sqrt(alpha(k)/(1-alpha(k)))
        mu = -lam*(1/kappa-kappa)/sqrt(2)
        if(type == 1) newdata = c(newdata, func1(k, scale, p, ver, leftp) + ralaplace(1, mu, lam, kappa))
        else newdata = c(newdata, ralaplace(1, mu, lam, kappa))
      }
      
      if(i == 1) data = list(newdata)
      else data = c(data, list(newdata))
    }
  }
  return(data)
}


### Simulation for 3 clusters

final3 <- function(num, sc, param = c(10, 0.5, 0, 0, 2.5, TRUE), curvenum = 50, simulnum = 50, pdiff = 0.5, sig = 0.2){
  csr1 = rep(0, simulnum)
  ari1 = rep(0, simulnum)
  maf1 = rep(0, simulnum)
  csr2 = rep(0, simulnum)
  ari2 = rep(0, simulnum)
  maf2 = rep(0, simulnum)
  
  csr5 = rep(0, simulnum)
  ari5 = rep(0, simulnum)
  maf5 = rep(0, simulnum)
  csr6 = rep(0, simulnum)
  ari6 = rep(0, simulnum)
  maf6 = rep(0, simulnum)
  csr7 = rep(0, simulnum)
  ari7 = rep(0, simulnum)
  maf7 = rep(0, simulnum)
  csr8 = rep(0, simulnum)
  ari8 = rep(0, simulnum)
  maf8 = rep(0, simulnum)
  
  num.knots = 10
  n_breaks = 100
  N = 3*curvenum
  bsp = create.bspline.basis(breaks=seq(0, 1, by=1/num.knots)) 
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  t = seq(0, 0.99, 0.01)
  c_true.matrix = matrix(nrow = simulnum, ncol = N)
  
  for(i in 1:simulnum){
    set.seed(123)
    per = sample(1:N, replace=FALSE)
    c_true = c(rep(1, curvenum), rep(2, curvenum), rep(3, curvenum))[per]
    c_true.matrix[i, ] = c_true
    print(paste("num : ", num, "scenario : ", sc, "sigma : ", sig))
    print(paste("clustering for ", i, "th simulation"))
    if(num == 1){
      param2 = param
      param2[2] = 0
      param2[1] = 4*max_f(pdiff, param2[1])
      param[2] = 2*pdiff/(1-pdiff) 
      param3 = param
      param3[6] = FALSE
      if(sc == 1) curve_list = c(curvelist_alaplace(param, sd = sig, alpha = 0.3, N = curvenum), curvelist_alaplace(param2, sd = sig, alpha = 0.3, N = curvenum), 
                                 curvelist_alaplace(param3, sd = sig, alpha = 0.3, N = curvenum))
      if(sc == 2) curve_list = c(curvelist_alaplace(param, sd = sig, alpha = 0.1, N = curvenum), curvelist_alaplace(param2, sd = sig, alpha = 0.1, N = curvenum), 
                                 curvelist_alaplace(param3, sd = sig, alpha = 0.9, N = curvenum))
      if(sc == 3) curve_list = c(curvelist_alaplace(param, sd = sig, alpha = 0.1, N = curvenum), curvelist_alaplace(param2, sd = sig, alpha = 0.5, N = curvenum), 
                                 curvelist_alaplace(param3, sd = sig, alpha = 0.9, N = curvenum))
    }
    if(num == 2){
      param[2] = 2*pdiff/(1-pdiff)
      param_ = param
      param_[6] = FALSE
      if(sc == 1) curve_list = c(curvelist_alaplace(param, sd = sig, alpha = 0.1, N = curvenum), curvelist_alaplace(param, sd = sig, alpha = 0.9, N = curvenum), 
                                 curvelist_alaplace(param_, sd = sig, alpha = 0.9, N = curvenum))
      if(sc == 2) curve_list = c(curvelist_alaplace(param, sd = sig, alpha = 0.3, N = curvenum), curvelist_alaplace(param, sd = sig, alpha = 0.7, N = curvenum), 
                                 curvelist_alaplace(param_, sd = sig, alpha = 0.7, N = curvenum))
      if(sc == 3) curve_list = c(curvelist_alaplace(param, sd = sig, alpha = alpha_f, N = curvenum, hetero = T), curvelist_alaplace(param, sd = sig, alpha = alpha_minus_f, N = curvenum, hetero = T), 
                                 curvelist_alaplace(param_, sd = sig, alpha = alpha_minus_f, N = curvenum, hetero = T))
    }
    if(num == 3){
      if(sc == 1) curve_list = c(curvelist_normal(param, sd = 0.5, N = curvenum), curvelist_alaplace(param, sd = sig, alpha = 0.1, N = curvenum), 
                                 curvelist_alaplace(param, sd = sig, alpha = 0.9, N = curvenum))
      if(sc == 2) curve_list = c(curvelist_normal(param, sd = 0.1, N = curvenum), curvelist_normal(param, sd = 0.5, N = curvenum), 
                                 curvelist_alaplace(param, sd = sig, alpha = 0.1, N = curvenum))
      if(sc == 3) curve_list = c(curvelist_normal(param, sd = 0.5, N = curvenum), curvelist_alaplace(param, sd = sig, alpha = alpha_f, N = curvenum, hetero = T), 
                                 curvelist_alaplace(param, sd = sig, alpha = alpha_minus_f, N = curvenum, hetero = T))
    }
    
    curve_list0 = list()
    for(sk in 1:N){
      curve_list0[[sk]] = curve_list[[per[sk]]]
    }
    curve_list = curve_list0
    b.sp = bs(seq(0, 1, by = 1/n_breaks), knots = seq(0, 1, length.out = num.knots+1))
    
    bmatrix = matrix(nrow = length(curve_list), ncol = (num.knots+3))
    coef.mean = matrix(nrow = (num.knots+3), ncol = length(curve_list))
    
    coef.qmedian = coef.mean
    coef.q0.3 = coef.mean
    coef.q0.7 = coef.mean
    
    vec.mean = matrix(nrow = length(t), ncol = length(curve_list))
    vec.q0.3 = vec.mean
    vec.qmedian = vec.mean
    vec.q0.7 = vec.mean
    
    time0 = Sys.time()
    
    for(j in 1:length(curve_list)){
      y = curve_list[[j]]
      coef.mean[,j] = expec.curve.pspline(y, 0.5, n_breaks = n_breaks, num.knots = num.knots)
      
      coef.qmedian[,j] = quant.curve.pspline(y, 0.5, n_breaks = n_breaks, num.knots = num.knots)
      coef.q0.3[,j] = quant.curve.pspline(y, 0.25, n_breaks = n_breaks, num.knots = num.knots)
      coef.q0.7[,j] = quant.curve.pspline(y, 0.75, n_breaks = n_breaks, num.knots = num.knots)
      
      vec.mean[,j] = basis_value_vector(coef.mean[,j], b.sp = b.sp)
      
      vec.qmedian[,j] = basis_value_vector(coef.qmedian[,j], b.sp = b.sp)
      vec.q0.3[,j] = basis_value_vector(coef.q0.3[,j], b.sp = b.sp)
      vec.q0.7[,j] = basis_value_vector(coef.q0.7[,j], b.sp = b.sp)
      
      bbb = as.matrix(nearPD(t(b.basis)%*%b.basis)$mat)
      bmatrix[j,] = solve(bbb, t(b.basis)%*%y)
    }
    if(i == 1) bmatrix.matrix = bmatrix
    else bmatrix.matrix = rbind(bmatrix.matrix, bmatrix)
    
    basis = create.bspline.basis(rangeval = c(0, 1), num.knots+3)
    argvals = matrix(rep(t, length(curve_list)), nrow = length(t), ncol = length(curve_list))
    
    mean.obj = Data2fd(argvals = argvals, y = vec.mean, basisobj = basis, lambda = 0.1)
    
    qmedian.obj = Data2fd(argvals = argvals, y = vec.qmedian, basisobj = basis, lambda = 0.1)
    qobj.0.3 = Data2fd(argvals = argvals, y = vec.q0.3, basisobj = basis, lambda = 0.1)
    qobj.0.7 = Data2fd(argvals = argvals, y = vec.q0.7, basisobj = basis, lambda = 0.1)
    
    mean.pca = pca.fd(mean.obj, nharm = 3)
    
    qmedian.pca = pca.fd(qmedian.obj, nharm = 3)
    qpca.0.3 = pca.fd(qobj.0.3, nharm = 3)
    qpca.0.7 = pca.fd(qobj.0.7, nharm = 3)
    
    
    mean.score = mean.pca$scores[,1]
    qmedian.score = qmedian.pca$scores[,1]
    qscore.0.3 = qpca.0.3$scores[,1]
    qscore.0.7 = qpca.0.7$scores[,1]
    
    mean.score2 = mean.pca$scores[,2]
    qmedian.score2 = qmedian.pca$scores[,2]
    qscore.0.32 = qpca.0.3$scores[,2]
    qscore.0.72 = qpca.0.7$scores[,2]
    
    qsparse = cbind(qscore.0.3, qscore.0.32, qmedian.score, qmedian.score2, mean.score, mean.score2, qscore.0.7, qscore.0.72)
    
    RSQfunclust1 = RSKC(qsparse, 3, alpha = 0.1, L1 = 3)
    RSQfunclust2 = RSKC.multidim(qsparse, ncl = 3, l1bound = 3)
    
    print("RSQfunclust1 weight")
    print(RSQfunclust1$weights)
    print("RSQfunclust2 weight")
    print(RSQfunclust2$weights)
    
    c1 <- RSQfunclust1$labels
    c2 <- RSQfunclust2$labels
    
    ### Exfunclust
    domain = seq(from=0, to=1, by=0.001)
    obs.domain = seq(from=0, to=1, by=0.01)[-101]
    knots <- seq(from=0, to=1, by=1/num.knots)
    b.expand <- bsplineS(domain, knots)
    simul.norm <- matrix(nrow=length(obs.domain), ncol=length(curve_list))
    for(ss in 1:length(curve_list)) simul.norm[,ss] = curve_list[[ss]]
    
    coeff.tau0.5 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.5)
    coeff.tau0.01 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.01)
    coeff.tau0.99 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.99)
    fit.expec0.5 <- b.expand%*%coeff.tau0.5
    fit.expec0.01 <- b.expand%*%coeff.tau0.01
    fit.expec0.99 <- b.expand%*%coeff.tau0.99
    
    w1 <- w2 <- matrix(nrow=N, ncol=num.knots)
    weighted.expec <- matrix(nrow=length(domain), ncol=N)
    for(i2 in 1:N){
      for(is in 1:num.knots){
        w1[i2,is] <- abs(0.001*sum(fit.expec0.99[c((1000/num.knots*(is-1)+1):(1000/num.knots*is+1)), i2]-fit.expec0.5[c((1000/num.knots*(is-1)+1):(1000/num.knots*is+1)), i2]))
        w2[i2,is] <- abs(0.001*sum(fit.expec0.5[c((1000/num.knots*(is-1)+1):(1000/num.knots*is+1)), i2]-fit.expec0.01[c((1000/num.knots*(is-1)+1):(1000/num.knots*is+1)), i2]))
      }
      #w1[i2,] <- 1/w1[i2,]; w2[i2,] <- 1/w2[i2,]
      w11 <- w1[i2,]/(w1[i2,]+w2[i2,]); w22 <- w2[i2,]/(w1[i2,]+w2[i2,]); w1[i2,] <- w11; w2[i2,] <- w22
      poo1 <- c(w1[i2,1], rep(w1[i2,], each=1000/num.knots)); poo2 <- c(w2[i2,1], rep(w2[i2,], each=1000/num.knots))
      weighted.expec[,i2] <-poo1*fit.expec0.99[,i2]+poo2*fit.expec0.01[,i2]
    }
    coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
    for(is in 1:N){
      coeff.tau.al[is,] <- solve(t(b.expand)%*%b.expand, t(b.expand)%*%weighted.expec[,is])
    }     
    
    Exfunclust <- kmeans(coeff.tau.al, centers=3, nstart = 10)
    c5 <- Exfunclust$cluster
    
    ### Baseclust
    c6 <- kmeans(bmatrix, centers = 3, nstart = 10)$cluster
    
    ### Funclust
    fd.obj <- fd(coef=coeff.tau0.5, basisobj=bsp)
    c7 <- funclust(fd.obj, K=3)$cls
    
    ### L2 clustering (mean)
    c8 <- distance_clustering(3, coef.mean)
    
    time1 = Sys.time()
    print(paste("1 simulation time : ", time1 - time0, "(min)"))
    
    
    csr1[i] = 100*perform_score(c1, c_true, "CSR")
    ari1[i] = perform_score(c1, c_true, "ARI")
    maf1[i] = 100*perform_score(c1, c_true, "MAF")
    print(paste("RSQfunclust1       : ", csr1[i], ari1[i], maf1[i]))
    
    csr2[i] = 100*perform_score(c2, c_true, "CSR")
    ari2[i] = perform_score(c2, c_true, "ARI")
    maf2[i] = 100*perform_score(c2, c_true, "MAF")
    print(paste("RSQfunclust2       : ", csr2[i], ari2[i], maf2[i]))
    
    csr5[i] = 100*perform_score(c5, c_true, "CSR")
    ari5[i] = perform_score(c5, c_true, "ARI")
    maf5[i] = 100*perform_score(c5, c_true, "MAF")
    print(paste("Exfunclust         : ", csr5[i], ari5[i], maf5[i]))
    
    csr6[i] = 100*perform_score(c6, c_true, "CSR")
    ari6[i] = perform_score(c6, c_true, "ARI")
    maf6[i] = 100*perform_score(c6, c_true, "MAF")
    print(paste("Baseclust          : ", csr6[i], ari6[i], maf6[i]))
    
    csr7[i] = 100*perform_score(c7, c_true, "CSR")
    ari7[i] = perform_score(c7, c_true, "ARI")
    maf7[i] = 100*perform_score(c7, c_true, "MAF")
    print(paste("Funclust           : ", csr7[i], ari7[i], maf7[i]))
    
    csr8[i] = 100*perform_score(c8, c_true, "CSR")
    ari8[i] = perform_score(c8, c_true, "ARI")
    maf8[i] = 100*perform_score(c8, c_true, "MAF")
    print(paste("Mean               : ", csr8[i], ari8[i], maf8[i]))
  }
  
  ret = c(mean(csr1), mean(ari1), mean(maf1), mean(csr2), mean(ari2), mean(maf2), 
          mean(csr5), mean(ari5), mean(maf5), mean(csr6), mean(ari6), mean(maf6),
          mean(csr7), mean(ari7), mean(maf7), mean(csr8), mean(ari8), mean(maf8))
  ret = c(ret, sd(csr1), sd(ari1), sd(maf1), sd(csr2), sd(ari2), sd(maf2), 
          sd(csr5), sd(ari5), sd(maf5), sd(csr6), sd(ari6), sd(maf6),
          sd(csr7), sd(ari7), sd(maf7), sd(csr8), sd(ari8), sd(maf8))
  
  return(list(ret, bmatrix.matrix, c_true.matrix)) # B-spline basis coefficients to use in K-expectile clustering
}

### For K-expectile clustering, since the source code is written in Python, run the code separately

k.expectile3 <- function(bmatrix.matrix, c_true.matrix, num, sc, sig) {
  setwd("C://Users//user//Desktop//4학년 2학기//인턴//Code//bmatrix")
  filename = paste("bmatrix3_", num, "_", sc, "_", 10*sig, ".csv", sep = '')
  write.csv(bmatrix.matrix, file=filename)
  filename = paste("c_true3_", num, "_", sc, "_", 10*sig, ".csv", sep = '')
  write.csv(c_true.matrix, file=filename)
}

a = final3(1, 1, simulnum = 2, sig = 0.1)
k.expectile3(a[[2]], a[[3]], 1, 1, 0.1) # save b-spline coefficient matrix data


### Example 1

final_score1 = matrix(nrow = 9, ncol = 36)

for(i in 1:3){
  for(j in 1:3){
    a = final3(1, i, simulnum = 50, sig = 0.1*j)
    k.expectile3(a[[2]], a[[3]], 1, i, 0.1*j)
    final_score1[(i-1)*3+j, ] = a[[1]]
  }
}


### Example 2

final_score2 = matrix(nrow = 9, ncol = 64)

final_score2[1, ] = final3(2, 1, simulnum = 50, sig = 0.1)
final_score2[2, ] = final3(2, 1, simulnum = 50, sig = 0.2)
final_score2[3, ] = final3(2, 1, simulnum = 50, sig = 0.3)

final_score2[4, ] = final3(2, 2, simulnum = 50, sig = 0.1)
final_score2[5, ] = final3(2, 2, simulnum = 50, sig = 0.2)
final_score2[6, ] = final3(2, 2, simulnum = 50, sig = 0.3)

final_score2[7, ] = final3(2, 3, simulnum = 50, sig = 0.1)
final_score2[8, ] = final3(2, 3, simulnum = 50, sig = 0.2)
final_score2[9, ] = final3(2, 3, simulnum = 50, sig = 0.3)


### Example 3

final_score3 = matrix(nrow = 9, ncol = 64)

final_score3[1, ] = final3(3, 1, simulnum = 50, sig = 0.1)
final_score3[2, ] = final3(3, 1, simulnum = 50, sig = 0.2)
final_score3[3, ] = final3(3, 1, simulnum = 50, sig = 0.3)

final_score3[4, ] = final3(3, 2, simulnum = 50, sig = 0.1)
final_score3[5, ] = final3(3, 2, simulnum = 50, sig = 0.2)
final_score3[6, ] = final3(3, 2, simulnum = 50, sig = 0.3)

final_score3[7, ] = final3(3, 3, simulnum = 50, sig = 0.1)
final_score3[8, ] = final3(3, 3, simulnum = 50, sig = 0.2)
final_score3[9, ] = final3(3, 3, simulnum = 50, sig = 0.3)




