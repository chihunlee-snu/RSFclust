### 3 cluster

sig.simul.set = rep(c(0.1, 0.2, 0.3), 6)
submodel.simul.set = c(rep(1, 9), rep(2, 9)) # 1 for asymmteric laplace, 2 for skew-t
simul1.set = cbind(rep(3, 18), rep(c(rep(1, 3), rep(2, 3), rep(3, 3)), 2), submodel.simul.set, sig.simul.set)
simul2.set = cbind(rep(4, 18), rep(c(rep(1, 3), rep(2, 3), rep(3, 3)), 2), submodel.simul.set, sig.simul.set)
simul3.set = cbind(rep(5, 18), rep(c(rep(1, 3), rep(2, 3), rep(3, 3)), 2), submodel.simul.set, sig.simul.set)
simul.set = rbind(simul1.set, simul2.set, simul3.set)

# (e.g. (S3)-(A1) : sim = 3, model = 1, submodel = 1)
# (e.g. (S3)-(T1) : sim = 3, model = 1, submodel = 2)
# (e.g. (S4)-(T2) : sim = 4, model = 2, submodel = 2)

### Simulation function (K=3)
simulation3 <- function(sim, model, submodel, sig, curvenum = 100, simulnum = 200, num.knots = 10, n_breaks = 100, kexp = F){
  csr1 = rep(0, simulnum)
  ari1 = rep(0, simulnum)
  maf1 = rep(0, simulnum)
  csr2 = rep(0, simulnum)
  ari2 = rep(0, simulnum)
  maf2 = rep(0, simulnum)
  
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
  expec.levels = c(0.5)
  
  q.num = length(quant.levels)
  e.num = length(expec.levels)
  K = 3
  
  for(i in 1:simulnum){
    start = Sys.time()
    set.seed(i)
    per = sample(1:N, replace=FALSE)
    c_true = c(rep(1, curvenum), rep(2, curvenum), rep(3, curvenum))[per]
    c_true.matrix[i,] = c_true
    
    error_dist = ifelse(submodel == 1, "A", "T")
    print(paste("(S", sim, ") - (", error_dist, model, "), sig : ", sig, sep = ""))
    print(paste("clustering for", i, "th simulation set"))
    
    if(sim == 3){
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
    if(sim == 4){
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
    if(sim == 5){
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
    
    RSF1 <- RSKC(score.matrix, 3, alpha = 0.1, L1 = 3)
    RSF2 <- RSKC.multidim(score.matrix, 3)
    
    c1 <- RSF1$labels
    c2 <- RSF2$labels
    
    ### Baseclust
    c7 <- kmeans(bmatrix, centers = 3, nstart = 10)$cluster
    
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
    c4 <- efckim$cluster
    
    end = Sys.time()
    print(paste("1 simulation time : ", end - start, "(s)"))
    print("---------- ARI (%) for each method ------------")
    
    csr1[i] = 100*perform_score(c1, c_true, "CSR")
    ari1[i] = 100*perform_score(c1, c_true, "ARI")
    maf1[i] = 100*perform_score(c1, c_true, "MAF")
    print(paste("RSFclust1          : ", round(ari1[i], 2)))
    
    csr2[i] = 100*perform_score(c2, c_true, "CSR")
    ari2[i] = 100*perform_score(c2, c_true, "ARI")
    maf2[i] = 100*perform_score(c2, c_true, "MAF")
    print(paste("RSFclust2          : ", round(ari2[i], 2)))
    
    csr4[i] = 100*perform_score(c4, c_true, "CSR")
    ari4[i] = 100*perform_score(c4, c_true, "ARI")
    maf4[i] = 100*perform_score(c4, c_true, "MAF")
    print(paste("Exfunclust         : ", round(ari4[i], 2)))
    
    csr5[i] = 100*perform_score(c5, c_true, "CSR")
    ari5[i] = 100*perform_score(c5, c_true, "ARI")
    maf5[i] = 100*perform_score(c5, c_true, "MAF")
    print(paste("FPCAC              : ", round(ari5[i], 2)))
    
    csr6[i] = 100*perform_score(c6, c_true, "CSR")
    ari6[i] = 100*perform_score(c6, c_true, "ARI")
    maf6[i] = 100*perform_score(c6, c_true, "MAF")
    print(paste("L2-dist            : ", round(ari6[i], 2)))
    
    csr7[i] = 100*perform_score(c7, c_true, "CSR")
    ari7[i] = 100*perform_score(c7, c_true, "ARI")
    maf7[i] = 100*perform_score(c7, c_true, "MAF")
    print(paste("Baseclust          : ", round(ari7[i], 2)))
    print("-----------------------------------------------")
  }
  ret = c(mean(csr1), mean(ari1), mean(maf1), mean(csr2), mean(ari2), mean(maf2), 
          mean(csr4), mean(ari4), mean(maf4), mean(csr5), mean(ari5), mean(maf5), mean(csr6), 
          mean(ari6), mean(maf6), mean(csr7), mean(ari7), mean(maf7))
  ret = c(ret, sd(csr1), sd(ari1), sd(maf1), sd(csr2), sd(ari2), sd(maf2), sd(csr4), sd(ari4), sd(maf4),
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

a = simulation3(5, 2, 1, 0.2, simulnum = 200)
k.expectile3(a[[2]], a[[3]], 1, 1, 0.1) # save b-spline coefficient data, and upload it to the python code of K-expectie clustering






