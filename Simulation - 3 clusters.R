### 3 cluster

sig.simul.set = rep(c(0.1, 0.2, 0.3), 6)
submodel.simul.set = rep(c(rep(1, 3), rep(2, 3)), 3)
simul1.set = cbind(rep(1, 18), c(rep(1, 6), rep(2, 6), rep(3, 6)), submodel.simul.set, sig.simul.set)
simul2.set = cbind(rep(2, 18), c(rep(1, 6), rep(2, 6), rep(3, 6)), submodel.simul.set, sig.simul.set)
simul3.set = cbind(rep(3, 18), c(rep(1, 6), rep(2, 6), rep(3, 6)), submodel.simul.set, sig.simul.set)
simul.set = rbind(simul1.set, simul2.set, simul3.set)

### Simulation function (K=3)
simulation3 <- function(sim, model, submodel, sig, curvenum = 100, simulnum = 200, num.knots = 10, n_breaks = 100, kexp = F){
  csr1.a = rep(0, simulnum)
  ari1.a = rep(0, simulnum)
  maf1.a = rep(0, simulnum)
  csr1.b = rep(0, simulnum)
  ari1.b = rep(0, simulnum)
  maf1.b = rep(0, simulnum)
  csr1.c = rep(0, simulnum)
  ari1.c = rep(0, simulnum)
  maf1.c = rep(0, simulnum)
  
  csr2.a = rep(0, simulnum)
  ari2.a = rep(0, simulnum)
  maf2.a = rep(0, simulnum)
  csr2.b = rep(0, simulnum)
  ari2.b = rep(0, simulnum)
  maf2.b = rep(0, simulnum)
  csr2.c = rep(0, simulnum)
  ari2.c = rep(0, simulnum)
  maf2.c = rep(0, simulnum)
  
  csr4 = rep(0, simulnum)
  ari4 = rep(0, simulnum)
  maf4 = rep(0, simulnum)
  csr5 = rep(0, simulnum)
  ari5 = rep(0, simulnum)
  maf5 = rep(0, simulnum)
  csr6 = rep(0, simulnum)
  ari6 = rep(0, simulnum)
  maf6 = rep(0, simulnum)
  csr7 = rep(0, simulnum)
  ari7 = rep(0, simulnum)
  maf7 = rep(0, simulnum)
  
  N = 3*curvenum
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  t = seq(0, 1, length.out = n_breaks+1)[1:n_breaks]
  c_true.matrix = matrix(nrow = simulnum, ncol = N)
  bmatrix.matrix = c()
  quant.levels = c(0.25, 0.5, 0.75)
  expec.levels = c(0.25, 0.5, 0.75)
  
  q.num = length(quant.levels)
  e.num = length(expec.levels)
  K = 3
  
  for(i in 1:simulnum){
    start = Sys.time()
    if(i != 197) set.seed(i*100)
    if(i == 197) set.seed(197)
    per = sample(1:N, replace=FALSE)
    c_true = c(rep(1, curvenum), rep(2, curvenum), rep(3, curvenum))[per]
    c_true.matrix[i,] = c_true
    
    dist = ifelse(submodel == 1, "asymmetric laplace", "Skew-t")
    print(paste("Simulation :", sim, ", Model :", model, ",  Sig :", sig, ",  Error :", dist))
    print(paste("clustering for", i, "th simulation set"))
    
    if(sim == 1){
      if(submodel == 1){
        if(model == 1) curve_list = c(curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.1),
                                      curvelist_alaplace(2, n_breaks, curvenum, sig, alpha = 0.1),
                                      curvelist_alaplace(3, n_breaks, curvenum, sig, alpha = 0.1))
        if(model == 2) curve_list = c(curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.1),
                                      curvelist_alaplace(2, n_breaks, curvenum, sig, alpha = 0.1),
                                      curvelist_alaplace(3, n_breaks, curvenum, sig, alpha = 0.9))
        if(model == 3) curve_list = c(curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.1),
                                      curvelist_alaplace(2, n_breaks, curvenum, sig, alpha = 0.5),
                                      curvelist_alaplace(3, n_breaks, curvenum, sig, alpha = 0.9))
      }
      else if(submodel == 2){
        if(model == 1) curve_list = c(curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = 5),
                                      curvelist_skewt(2, n_breaks, curvenum, omega = sig, alpha = 5),
                                      curvelist_skewt(3, n_breaks, curvenum, omega = sig, alpha = 5))
        if(model == 2) curve_list = c(curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = 5),
                                      curvelist_skewt(2, n_breaks, curvenum, omega = sig, alpha = 5),
                                      curvelist_skewt(3, n_breaks, curvenum, omega = sig, alpha = -5))
        if(model == 3) curve_list = c(curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = 5),
                                      curvelist_skewt(2, n_breaks, curvenum, omega = sig, alpha = 0),
                                      curvelist_skewt(3, n_breaks, curvenum, omega = sig, alpha = -5))
      }
    }
    if(sim == 2){
      if(submodel == 1){
        if(model == 1) curve_list = c(curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.1),
                                      curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.9),
                                      curvelist_alaplace(3, n_breaks, curvenum, sig, alpha = 0.9))
        if(model == 2) curve_list = c(curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.3),
                                      curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.7),
                                      curvelist_alaplace(3, n_breaks, curvenum, sig, alpha = 0.7))
        if(model == 3) curve_list = c(curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = alpha_f, hetero = T),
                                      curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = alpha_minus_f, hetero = T),
                                      curvelist_alaplace(3, n_breaks, curvenum, sig, alpha = alpha_minus_f, hetero = T))
      }
      else if(submodel == 2){
        if(model == 1) curve_list = c(curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = 5),
                                      curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = -5),
                                      curvelist_skewt(3, n_breaks, curvenum, omega = sig, alpha = -5))
        if(model == 2) curve_list = c(curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = 1),
                                      curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = -1),
                                      curvelist_skewt(3, n_breaks, curvenum, omega = sig, alpha = -1))
        if(model == 3) curve_list = c(curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = alpha_f, hetero = T),
                                      curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = minus_alpha_f, hetero = T),
                                      curvelist_skewt(3, n_breaks, curvenum, omega = sig, alpha = minus_alpha_f, hetero = T))
      }
    }
    if(sim == 3){
      if(submodel == 1){
        if(model == 1) curve_list = c(curvelist_normal(1, n_breaks, curvenum, 0.5),
                                      curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.1),
                                      curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.9))
        if(model == 2) curve_list = c(curvelist_normal(1, n_breaks, curvenum, 0.1),
                                      curvelist_normal(1, n_breaks, curvenum, 0.5),
                                      curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.1))
        if(model == 3) curve_list = c(curvelist_normal(1, n_breaks, curvenum, 0.5),
                                      curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = alpha_f, hetero = T),
                                      curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = alpha_minus_f, hetero = T))
      }
      else if(submodel == 2){
        if(model == 1) curve_list = c(curvelist_normal(1, n_breaks, curvenum, 0.5),
                                      curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = 5),
                                      curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = -5))
        if(model == 2) curve_list = c(curvelist_normal(1, n_breaks, curvenum, 0.1),
                                      curvelist_normal(1, n_breaks, curvenum, 0.5),
                                      curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = 5))
        if(model == 3) curve_list = c(curvelist_normal(1, n_breaks, curvenum, 0.5),
                                      curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = alpha_f, hetero = T),
                                      curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = minus_alpha_f, hetero = T))
      }
    }
    
    curve_list0 = list()
    for(s in 1:N) curve_list0[[s]] = curve_list[[per[s]]]
    curve_list = curve_list0
    
    bmatrix = matrix(nrow = length(curve_list), ncol = num.knots+3)
    coef.mean = matrix(nrow = (num.knots+3), ncol = length(curve_list))
    
    for(j in 1:length(curve_list)){
      y = curve_list[[j]]
      bmatrix[j,] = solve(t(b.basis)%*%b.basis, t(b.basis)%*%y)
    }
    
    bmatrix.matrix = rbind(bmatrix.matrix, bmatrix)
    if(kexp) next
    
    for(j in 1:length(curve_list)){
      y = curve_list[[j]]
      if(j == 1) lambda.mean = gcv.select(y, 0.5, num.knots = num.knots, n_breaks = n_breaks)
      coef.mean[,j] = expec.curve.pspline(y, 0.5, num.knots = num.knots, n_breaks = n_breaks, lambda = lambda.mean)
    }
    
    ### Extract FPCA input features
    score.matrix = FPC.score.matrix(curve_list, quant.levels, expec.levels, K, num.knots, n_breaks)
    
    score.matrix.a = score.matrix[, c(1, 2, 3, 4, 5, 6)]
    score.matrix.b = score.matrix[, c(1, 2, 3, 4, 5, 6, 9, 10)]
    score.matrix.c = score.matrix
    
    RSF1.a <- RSKC(score.matrix.a, 3, alpha = 0.1, L1 = 3)
    RSF1.b <- RSKC(score.matrix.b, 3, alpha = 0.1, L1 = 3)
    RSF1.c <- RSKC(score.matrix.c, 3, alpha = 0.1, L1 = 3)
    RSF2.a <- RSKC.multidim(score.matrix.a, 3)
    RSF2.b <- RSKC.multidim(score.matrix.b, 3)
    RSF2.c <- RSKC.multidim(score.matrix.c, 3)
    
    c1.a <- RSF1.a$labels
    c1.b <- RSF1.b$labels
    c1.c <- RSF1.c$labels
    c2.a <- RSF2.a$labels
    c2.b <- RSF2.b$labels
    c2.c <- RSF2.c$labels
    
    ### Baseclust
    c4 <- kmeans(bmatrix, centers = 3, nstart = 10)$cluster
    
    ### FPCAC
    X = matrix(nrow = length(curve_list), ncol = n_breaks)
    for(j in 1:length(curve_list)) X[j,] = curve_list[[j]]
    c5 <- fpcac(t(X), K = 3, alpha = 0.1)$clusters
    
    ### Mean
    c6 <- distance_clustering(3, coef.mean, num.knots, n_breaks)
    
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
      # w1[i2,] <- 1/w1[i2,]; w2[i2,] <- 1/w2[i2,]
      w11 <- w1[i2,]/(w1[i2,]+w2[i2,]); w22 <- w2[i2,]/(w1[i2,]+w2[i2,]); w1[i2,] <- w11; w2[i2,] <- w22
      poo1 <- c(w1[i2,1], rep(w1[i2,], each=1000/num.knots)); poo2 <- c(w2[i2,1], rep(w2[i2,], each=1000/num.knots))
      weighted.expec[,i2] <-poo1*fit.expec0.99[,i2]+poo2*fit.expec0.01[,i2]
    }
    coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
    for(is in 1:N){
      coeff.tau.al[is,] <- solve(t(b.expand)%*%b.expand, t(b.expand)%*%weighted.expec[,is])
    }     
    
    efckim <- kmeans(coeff.tau.al, centers=3, nstart = 10)
    c7 <- efckim$cluster
    
    end = Sys.time()
    print(difftime(end, start, units = "secs"))
    
    csr1.a[i] = 100*perform_score(c1.a, c_true, "CSR")
    ari1.a[i] = perform_score(c1.a, c_true, "ARI")
    maf1.a[i] = 100*perform_score(c1.a, c_true, "MAF")
    print(paste("RSFclust1 (a)     : ", csr1.a[i]))
    
    csr1.b[i] = 100*perform_score(c1.b, c_true, "CSR")
    ari1.b[i] = perform_score(c1.b, c_true, "ARI")
    maf1.b[i] = 100*perform_score(c1.b, c_true, "MAF")
    print(paste("RSFclust1 (b)     : ", csr1.b[i]))
    
    csr1.c[i] = 100*perform_score(c1.c, c_true, "CSR")
    ari1.c[i] = perform_score(c1.c, c_true, "ARI")
    maf1.c[i] = 100*perform_score(c1.c, c_true, "MAF")
    print(paste("RSFclust1 (c)     : ", csr1.c[i]))
    
    csr2.a[i] = 100*perform_score(c2.a, c_true, "CSR")
    ari2.a[i] = perform_score(c2.a, c_true, "ARI")
    maf2.a[i] = 100*perform_score(c2.a, c_true, "MAF")
    print(paste("RSFclust2 (a)     : ", csr2.a[i]))
    
    csr2.b[i] = 100*perform_score(c2.b, c_true, "CSR")
    ari2.b[i] = perform_score(c2.b, c_true, "ARI")
    maf2.b[i] = 100*perform_score(c2.b, c_true, "MAF")
    print(paste("RSFclust2 (b)     : ", csr2.b[i]))
    
    csr2.c[i] = 100*perform_score(c2.c, c_true, "CSR")
    ari2.c[i] = perform_score(c2.c, c_true, "ARI")
    maf2.c[i] = 100*perform_score(c2.c, c_true, "MAF")
    print(paste("RSFclust2 (c)     : ", csr2.c[i]))
    
    csr4[i] = 100*perform_score(c4, c_true, "CSR")
    ari4[i] = perform_score(c4, c_true, "ARI")
    maf4[i] = 100*perform_score(c4, c_true, "MAF")
    print(paste("Baseclust         : ", csr4[i]))
    
    csr5[i] = 100*perform_score(c5, c_true, "CSR")
    ari5[i] = perform_score(c5, c_true, "ARI")
    maf5[i] = 100*perform_score(c5, c_true, "MAF")
    print(paste("FPCAC             : ", csr5[i]))
    
    csr6[i] = 100*perform_score(c6, c_true, "CSR")
    ari6[i] = perform_score(c6, c_true, "ARI")
    maf6[i] = 100*perform_score(c6, c_true, "MAF")
    print(paste("L2-dist           : ", csr6[i]))
    
    csr7[i] = 100*perform_score(c7, c_true, "CSR")
    ari7[i] = perform_score(c7, c_true, "ARI")
    maf7[i] = 100*perform_score(c7, c_true, "MAF")
    print(paste("Exfunclust        : ", csr7[i]))
  }
  ret = c(mean(csr1.a), mean(ari1.a), mean(maf1.a), mean(csr1.b), mean(ari1.b), mean(maf1.b), 
          mean(csr1.c), mean(ari1.c), mean(maf1.c), mean(csr2.a), mean(ari2.a), mean(maf2.a), 
          mean(csr2.b), mean(ari2.b), mean(maf2.b), mean(csr2.c), mean(ari2.c), mean(maf2.c), 
          mean(csr4), mean(ari4), mean(maf4), mean(csr5), mean(ari5), mean(maf5), mean(csr6), 
          mean(ari6), mean(maf6), mean(csr7), mean(ari7), mean(maf7))
  ret = c(ret, sd(csr1.a), sd(ari1.a), sd(maf1.a), sd(csr1.b), sd(ari1.b), sd(maf1.b),
          sd(csr1.c), sd(ari1.c), sd(maf1.c), sd(csr2.a), sd(ari2.a), sd(maf2.a), sd(csr2.b), sd(ari2.b), sd(maf2.b),
          sd(csr2.c), sd(ari2.c), sd(maf2.c), sd(csr4), sd(ari4), sd(maf4),
          sd(csr5), sd(ari5), sd(maf5), sd(csr6), sd(ari6), sd(maf6),
          sd(csr7), sd(ari7), sd(maf7))
  
  if(!kexp) return(list(ret, bmatrix.matrix, c_true.matrix)) # B-spline basis coefficients to use in K-expectile clustering
  if(kexp) return(list(bmatrix.matrix, c_true.matrix))
}

### For K-expectile clustering, since the source code is written in Python, run the code separately
### Refer to https://github.com/QuantLet/KEC/tree/master/

k.expectile3 <- function(bmatrix.matrix, c_true.matrix, sim, model, submodel, sig) {
  filename = paste("bmatrix3_", sim, "_", model, "_", submodel, "_", sig, ".csv", sep = '')
  write.csv(bmatrix.matrix, file=filename)
  filename = paste("c_true3_", sim, "_", model, "_", submodel, "_", sig, ".csv", sep = '')
  write.csv(c_true.matrix, file=filename)
}

a = simulation3(2, 2, 2, 0.1, simulnum = 1) # 33 seconds per simulation
k.expectile3(a[[2]], a[[3]], 1, 1, 0.1) # save b-spline coefficient data, and upload it to the python code of K-expectie clustering


