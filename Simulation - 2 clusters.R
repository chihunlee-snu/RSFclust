### 2 cluster

### Simulation function (K=2)
simulation2 <- function(sim, model, submodel, sig, curvenum = 100, simulnum = 200, num.knots = 10, n_breaks = 100){
  # sim : simulation number (e.g. (S1) -> 1, (S2) -> 2)
  # model : model number (e.g. (M1) -> 1, (M2) -> 2, (M3-3) -> 3)
  # submodel : model (sub)number (e.g. (M1) -> 1, (M3-1) -> 1, (M3-2)-> 2, (M3-3) -> 3)
  csr1 = rep(0, simulnum)
  ari1 = rep(0, simulnum)
  maf1 = rep(0, simulnum)

  csr2 = rep(0, simulnum)
  ari2 = rep(0, simulnum)
  maf2 = rep(0, simulnum)
  csr3 = rep(0, simulnum)
  ari3 = rep(0, simulnum)
  maf3 = rep(0, simulnum)
  csr4 = rep(0, simulnum)
  ari4 = rep(0, simulnum)
  maf4 = rep(0, simulnum)
  csr5 = rep(0, simulnum)
  ari5 = rep(0, simulnum)
  maf5 = rep(0, simulnum)
  
  b.basis = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
  N = 2*curvenum
  K = 2
  t = seq(0, 1, length.out = n_breaks+1)[1:n_breaks]
  c_true.matrix = matrix(nrow = simulnum, ncol = N)
  
  quant.levels = c(0.25, 0.5, 0.75)
  expec.levels = c(0.5)
  
  q.num = length(quant.levels)
  e.num = length(expec.levels)
  
  for(i in 1:simulnum){
    start = Sys.time()
    set.seed(100*i)
    per = sample(1:N, replace=FALSE)
    c_true = c(rep(1, curvenum), rep(2, curvenum))[per]
    c_true.matrix[i, ] = c_true
    print(paste("Simulation :", sim, "model :", model, "submodel :", submodel, "sig :", sig))
    print(paste("clustering for", i, "th simulation"))
    
    if(sim == 1){
      if(model == 1) curve_list = c(curvelist_normal(1, n_breaks, curvenum, sig), curvelist_normal(2, n_breaks, curvenum, sig))
      if(model == 2) curve_list = c(curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.1), curvelist_alaplace(2, n_breaks, curvenum, sig, alpha = 0.1))
      if(model == 3) curve_list = c(curvelist_normal(1, n_breaks, curvenum, sig), curvelist_normal(3, n_breaks, curvenum, sig))
      if(model == 4) curve_list = c(curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.1), curvelist_alaplace(3, n_breaks, curvenum, sig, alpha = 0.1))
    }
    if(sim == 2){
      if(model == 1) curve_list = c(curvelist_normal(1, n_breaks, curvenum, 0.1), curvelist_normal(1, n_breaks, curvenum, sig))
      if(model == 2) curve_list = c(curvelist_alaplace(1, n_breaks, curvenum, sigma = 0.05, alpha = 0.1), curvelist_alaplace(1, n_breaks, curvenum, sigma = sig, alpha = 0.1))
      if(model == 3){
        if(submodel == 1) curve_list = c(curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.1), curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.9))
        if(submodel == 2) curve_list = c(curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.3), curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = 0.7))
        if(submodel == 3) curve_list = c(curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = alpha_f, hetero = T), curvelist_alaplace(1, n_breaks, curvenum, sig, alpha = alpha_minus_f, hetero = T))
      }
      if(model == 4){
        if(submodel == 1) curve_list = c(curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = 5), curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = -5))
        if(submodel == 2) curve_list = c(curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = 1), curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = -1))
        if(submodel == 3) curve_list = c(curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = alpha_f, hetero = T), curvelist_skewt(1, n_breaks, curvenum, omega = sig, alpha = minus_alpha_f, hetero = T))
      }
    }
    
    curve_list0 = list()
    for(sk in 1:N) curve_list0[[sk]] = curve_list[[per[sk]]]
    curve_list = curve_list0
    
    bmatrix = matrix(nrow = length(curve_list), ncol = num.knots+3)
    coef.mean = matrix(nrow = (num.knots+3), ncol = length(curve_list))

    for(j in 1:length(curve_list)){
      y = curve_list[[j]]
      if(j == 1) lambda.mean = gcv.select(y, 0.5, num.knots = num.knots, n_breaks = n_breaks)
      coef.mean[,j] = expec.curve.pspline(y, 0.5, num.knots = num.knots, n_breaks = n_breaks, lambda = lambda.mean)
      bmatrix[j,] = solve(t(b.basis)%*%b.basis, t(b.basis)%*%y)
    }
    
    if(i == 1) bmatrix.matrix = bmatrix
    else bmatrix.matrix = rbind(bmatrix.matrix, bmatrix)
    
    ### Extract FPCA input features
    score.matrix = FPC.score.matrix(curve_list, quant.levels, expec.levels, K, num.knots, n_breaks)
    
    RSF <- RSKC(score.matrix, 2, alpha = 0.1, L1 = 3)
    c1 <- RSF$labels
    
    ### Baseclust
    c2 <- kmeans(bmatrix, centers = 2, nstart = 10)$cluster
    
    ### FPCAC
    X = matrix(nrow = length(curve_list), ncol = n_breaks)
    for(j in 1:length(curve_list)) X[j,] = curve_list[[j]]
    c3 <- fpcac(t(X), K = 2, alpha = 0.1)$clusters
    
    ### Mean
    c4 <- distance_clustering(2, coef.mean, num.knots, n_breaks)
        
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
    
    efckim <- kmeans(coeff.tau.al, centers=2, nstart = 10)
    c5 <- efckim$cluster
    
    end = Sys.time()
    print(paste("1 simulation time : ", end - start, "(s)"))
    
    csr1[i] = 100*perform_score(c1, c_true, "CSR")
    ari1[i] = perform_score(c1, c_true, "ARI")
    maf1[i] = 100*perform_score(c1, c_true, "MAF")
    print(paste("RSFclust2         : ", csr1[i], ari1[i], maf1[i]))
    
    csr2[i] = 100*perform_score(c2, c_true, "CSR")
    ari2[i] = perform_score(c2, c_true, "ARI")
    maf2[i] = 100*perform_score(c2, c_true, "MAF")
    print(paste("Baseclust         : ", csr2[i], ari2[i], maf2[i]))
    
    csr3[i] = 100*perform_score(c3, c_true, "CSR")
    ari3[i] = perform_score(c3, c_true, "ARI")
    maf3[i] = 100*perform_score(c3, c_true, "MAF")
    print(paste("FPCAC             : ", csr3[i], ari3[i], maf3[i]))
    
    csr4[i] = 100*perform_score(c4, c_true, "CSR")
    ari4[i] = perform_score(c4, c_true, "ARI")
    maf4[i] = 100*perform_score(c4, c_true, "MAF")
    print(paste("L2-dist           : ", csr4[i], ari4[i], maf4[i]))
    
    csr5[i] = 100*perform_score(c5, c_true, "CSR")
    ari5[i] = perform_score(c5, c_true, "ARI")
    maf5[i] = 100*perform_score(c5, c_true, "MAF")
    print(paste("Exfunclust        : ", csr5[i], ari5[i], maf5[i]))
  }
  ret = c(mean(csr1), mean(ari1), mean(maf1), mean(csr2), mean(ari2), mean(maf2),
          mean(csr3), mean(ari3), mean(maf3), mean(csr4), mean(ari4), mean(maf4),
          mean(csr5), mean(ari5), mean(maf5))
  ret = c(ret, sd(csr1), sd(ari1), sd(maf1), sd(csr2), sd(ari2), sd(maf2),
          sd(csr3), sd(ari3), sd(maf3), sd(csr4), sd(ari4), sd(maf4),
          sd(csr5), sd(ari5), sd(maf5))
  
  return(list(ret, bmatrix.matrix, c_true.matrix))
}

### Run K-expectile clustering separately, since the source code of K-expectile clustering is written in Python code
### Refer to https://github.com/QuantLet/KEC/tree/master/

k.expectile2 <- function(bmatrix.matrix, c_true.matrix, sim, model, submodel, sig) {
  filename = paste("bmatrix2_", sim, "_", model, "_", submodel, "_", sig, ".csv", sep = '')
  write.csv(bmatrix.matrix, file=filename)
  filename = paste("c_true2_", sim, "_", model, "_", submodel, "_", sig, ".csv", sep = '')
  write.csv(c_true.matrix, file=filename)
}

a = simulation2(2, 1, 2, 0.1, simulnum = 1) # 16 seconds per simulation
k.expectile2(a[[2]], a[[3]], 2, 3, 2, 0.1) # save b-spline coefficient data, and upload it to the python code of K-expectie clustering












