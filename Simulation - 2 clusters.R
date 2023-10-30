### 2 cluster

final2 <- function(num, sc, errornum = 1, param = c(10, 0.5, 0, 0, 2.5, TRUE), curvenum = 50, simulnum = 50, pdiff = 0.2, sig = 0.2){
  csr1 = rep(0, simulnum)
  ari1 = rep(0, simulnum)
  maf1 = rep(0, simulnum)

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
  
  bsp = create.bspline.basis(breaks=seq(0, 1, by=1/num.knots))
  N = 2*curvenum
  t = seq(0, 0.99, 0.01)
  c_true.matrix = matrix(nrow = simulnum, ncol = N)
  
  for(i in 1:simulnum){
    set.seed(123)
    per = sample(1:N, replace=FALSE)
    c_true = c(rep(1, curvenum), rep(2, curvenum))[per]
    c_true.matrix[i, ] = c_true
    print(paste("num : ", num, "scenario : ", sc, "sigma : ", sig))
    print(paste("clustering for ", i, "th simulation"))
    if(num == 1 & sc %in% c(1, 2)){
      param_ = param
      param_[3] = 0.2
      if(sc == 1) curve_list = c(curvelist_normal(param, sd = sig, N = curvenum), curvelist_normal(param_, sd = sig, N = curvenum))
      if(sc == 2) curve_list = c(curvelist_alaplace(param, sd = sig, alpha = 0.3, N = curvenum), curvelist_alaplace(param_, sd = sig, alpha = 0.3, N = curvenum))
    }
    if(num == 1 & sc %in% c(3, 4)){
      param[2] = 2*pdiff/(1-pdiff)
      param_ = param
      param_[6] = FALSE
      if(sc == 3) curve_list = c(curvelist_normal(param, sd = sig, N = curvenum), curvelist_normal(param_, sd = sig, N = curvenum))
      if(sc == 4) curve_list = c(curvelist_alaplace(param, sd = sig, alpha = 0.3, N = curvenum), curvelist_alaplace(param_, sd = sig, alpha = 0.3, N = curvenum))
    }
    
    if(num == 2 & sc == 1) curve_list = c(curvelist_alaplace(param, sd = 0.1, alpha = 0.3, N = curvenum), curvelist_alaplace(param, sd = sig, alpha = 0.3, N = curvenum))
    if(num == 2 & sc == 2) curve_list = c(curvelist_normal(param, sd = 0.1, N = curvenum), curvelist_normal(param, sd = sig, N = curvenum))
    if(num == 2 & sc == 3){
      if(errornum == 1) curve_list = c(curvelist_alaplace(param, sd = sig, alpha = 0.1, N = curvenum), curvelist_alaplace(param, sd = sig, alpha = 0.9, N = curvenum))
      if(errornum == 2) curve_list = c(curvelist_alaplace(param, sd = sig, alpha = 0.3, N = curvenum), curvelist_alaplace(param, sd = sig, alpha = 0.7, N = curvenum))
      if(errornum == 3) curve_list = c(curvelist_alaplace(param, sd = sig, alpha = alpha_f, N = curvenum, hetero = T), 
                                       curvelist_alaplace(param, sd = sig, alpha = alpha_minus_f, N = curvenum, hetero = T))
    }
    if(num == 2 & sc == 4){
      if(errornum == 1) curve_list = c(curvelist_normal(param, sd = 0.5, N = curvenum), curvelist_alaplace(param, sd = 0.1, alpha = 0.3, N = curvenum))
      if(errornum == 2) curve_list = c(curvelist_normal(param, sd = 0.5, N = curvenum), curvelist_alaplace(param, sd = 0.3, alpha = 0.3, N = curvenum))
    }
    
    curve_list0 = list()
    for(sk in 1:N){
      curve_list0[[sk]] = curve_list[[per[sk]]]
    }
    curve_list = curve_list0
    
    bmatrix = matrix(nrow = length(curve_list), ncol = num.knots+3)
    coef.mean = matrix(nrow = num.knots+3, ncol = length(curve_list))
    
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
      coef.mean[,j] = expec.curve.pspline(y, 0.5)
      
      coef.qmedian[,j] = quant.curve.pspline(y, 0.5)
      coef.q0.3[,j] = quant.curve.pspline(y, 0.25)
      coef.q0.7[,j] = quant.curve.pspline(y, 0.75)
      
      vec.mean[,j] = basis_value_vector(coef.mean[,j])
      
      vec.qmedian[,j] = basis_value_vector(coef.qmedian[,j])
      vec.q0.3[,j] = basis_value_vector(coef.q0.3[,j])
      vec.q0.7[,j] = basis_value_vector(coef.q0.7[,j])
      
      bmatrix[j,] = solve(t(b.basis)%*%b.basis, t(b.basis)%*%y)
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
    
    qsparse = cbind(qscore.0.3, qmedian.score, mean.score, qscore.0.7)
    
    RSQfunclust2 = RSKC(qsparse, 2, alpha = 0.1, L1 = 3)
    
    print("RSQfunclust2 weight")
    print(RSQfunclust2$weights)
    
    c1 <- RSQfunclust2$labels
    
    ### Exfuclust
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
    
    efckim <- kmeans(coeff.tau.al, centers=2, nstart = 10)
    c5 <- efckim$cluster
    
    ### Baseclust
    c6 <- kmeans(bmatrix, centers = 2, nstart = 10)$cluster
    
    ### Funclust
    fd.obj <- fd(coef=coeff.tau0.5, basisobj=bsp)
    c7 <- funclust(fd.obj, K=2)$cls
    
    ### Mean
    c8 <- distance_clustering(2, coef.mean)
    
    time1 = Sys.time()
    print(paste("1 simulation time : ", time1 - time0, "(s)"))
    
    
    csr1[i] = 100*perform_score(c1, c_true, "CSR")
    ari1[i] = perform_score(c1, c_true, "ARI")
    maf1[i] = 100*perform_score(c1, c_true, "MAF")
    print(paste("RSQfunclust2       : ", csr1[i], ari1[i], maf1[i]))
    
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
  
  ret = c(mean(csr1), mean(ari1), mean(maf1), 
          mean(csr5), mean(ari5), mean(maf5), mean(csr6), mean(ari6), mean(maf6),
          mean(csr7), mean(ari7), mean(maf7), mean(csr8), mean(ari8), mean(maf8))
  ret = c(ret, sd(csr1), sd(ari1), sd(maf1),
          sd(csr5), sd(ari5), sd(maf5), sd(csr6), sd(ari6), sd(maf6),
          sd(csr7), sd(ari7), sd(maf7), sd(csr8), sd(ari8), sd(maf8))
  
  return(list(ret, bmatrix.matrix, c_true.matrix))
}

k.expectile2 <- function(bmatrix.matrix, c_true.matrix, num, sc, j) {
  setwd("C://Users//user//Desktop//4학년 2학기//인턴//Code//bmatrix")
  filename = paste("bmatrix2_", num, "_", sc, "_", j, ".csv", sep = '')
  write.csv(bmatrix.matrix, file=filename)
  filename = paste("c_true2_", num, "_", sc, "_", j, ".csv", sep = '')
  write.csv(c_true.matrix, file=filename)
}

a = final2(1, 1, simulnum = 2, sig = 0.1)
k.expectile2(a[[2]], a[[3]], 1, 1, 0.1) # save b-spline coefficient matrix data


### Example 1

score.2.1 = matrix(nrow = 8, ncol = 30)

for(i in 1:4){
  for(j in 1:2){
    if(i%%2 == 1 & j == 1) sig = 0.3
    if(i%%2 == 1 & j == 2) sig = 0.7 
    if(i%%2 == 0 & j == 1) sig = 0.1 
    if(i%%2 == 0 & j == 2) sig = 0.3 
    a = final2(1, i, simulnum = 50, sig = sig)
    k.expectile2(a[[2]], a[[3]], 1, i, j)
    score2.1[(i-1)*2+j, ] = a[[1]]
  }
}


### Example 2

score.2.2 = matrix(nrow = 12, ncol = 30)

for(i in 1:4){
 if(i != 3){ 
    for(j in 1:2){
      if(i == 1) sig = 0.3*(j==1) + 0.5*(j==2)
      if(i == 2) sig = 0.5*(j==1) + 1.0*(j==2)
      if(i == 4) sig = 0.1 # In this case, we don't use $sig$
      if(i < 3){
        a = final2(2, i, simulnum = 50, sig = sig)
        k.expectile2(a[[2]], a[[3]], 2, i, j)
        score2.2[(i-1)*2+j, ] = a[[1]]
      }
      else{
        a = final2(2, i, errornum = j, simulnum = 50)
        k.expectile2(a[[2]], a[[3]], 2, i, j)
        score2.2[10+j, ] = a[[1]]
      }
    }
 }
 else{
   for(j in 1:6){
     sig = 0.1 + 0.2*(j%%2 == 0)
     a = final2(2, 3, errornum = 1+(j-1)%/%2, simulnum = 50, sig = sig)
     k.expectile2(a[[2]], a[[3]], 2, i, j)
     score2.2[4+j, ] = a[[1]]
   }
 }
}


















