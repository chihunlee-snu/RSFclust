### RSFclust2 code

l2n <- function(vec){
  return(sqrt(sum(vec**2)))
}

Bin.search <- function (argu, sumabs){
  if (l2n(argu) == 0 || sum(abs(argu/l2n(argu))) <= sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu)) - 1e-05
  iter <- 1
  while (iter <= 15 && (lam2 - lam1) > (1e-04)) {
    su <- soft(argu, (lam1 + lam2)/2)
    if (sum(abs(su/l2n(su))) < sumabs) {
      lam2 <- (lam1 + lam2)/2
    }
    else {
      lam1 <- (lam1 + lam2)/2
    }
    iter <- iter + 1
  }
  return((lam1 + lam2)/2)
}

soft <- function (x, d){
  return(sign(x) * pmax(0, abs(x) - d))
}

bcss <- function(data, cluster){ # data : n*p
  k = length(unique(cluster))
  if(length(cluster) != nrow(data)) print("cluster, data dimension error")
  p = ncol(data)
  n = nrow(data)
  ret = rep(0, p)
  for(i in 1:p){
    pvec = data[,i]
    for(c in 1:k){
      pvec_c = pvec[cluster==c]
      nk = length(pvec_c)
      ck = 0
      for(j in pvec_c){
        for(s in pvec_c) ck = ck + (j-s)**2
      }
      ck = ck/nk
      ret[i] = ret[i] - ck
    }
    
    for(j in 1:n){
      for(s in 1:n){
        ret[i] = ret[i] + (pvec[j] - pvec[s])**2/n
      }
    }
    
  }
  return(ret/2)
}

trimmed.Kmeans <- function(x, w, ncl, alpha){
  # x : n*p(K-1), w : p(K-1) length vector
  z = sweep(x, 2, sqrt(w), "*")
  ifbreak = 1
  
  while(ifbreak == 1){
    centre = z[sample(1:nrow(z), ncl), ]
    old.centre = centre + 1
    ifbreak = 0
    iter = 0
    
    while(sum((old.centre - centre)**2) > 1e-6){
      ifstart = 1
      iter = iter + 1
      old.centre = centre
      kobj = kmeans(z, centers = old.centre)
      cl = kobj$cluster
      centre = kobj$centers
      mudist = vector(length = nrow(x))
      for(i in 1:nrow(x)) mudist[i] = sum((centre[cl[i], ] - z[i, ])**2)
      
      sort_mudist = sort(mudist, decreasing = TRUE)
      cut = sort_mudist[round(nrow(x)*alpha)+1]
      for(k in 1:ncl){
        if(sum(cl==k & mudist<=cut) > 1) centre[k, ] = colMeans(z[(cl==k & mudist<=cut),])
        else if(sum(cl==k & mudist<=cut) == 1){
          print("sum(cl==k & mudist<=cut) == 1")
          print(c(sum(cl==k), sum(cl==k & mudist<=cut)))
          centre[k, ] = z[(cl==k & mudist<=cut),]
          print(centre[k,])
        }
        else{
          print("Initial cluster center reset")
          ifbreak = 1
          break
        }
      }
      if(ifbreak == 1) break
      ow = (1:nrow(x))[mudist>cut]
      if(iter > 500) break
    }
  }
  unweighted.centre = sweep(centre, 2, sqrt(w), "/")
  
  mudist = vector(length = nrow(x))
  for(i in 1:nrow(x)) mudist[i] = sum((unweighted.centre[cl[i], ] - x[i, ])**2)
  
  sort_mudist = sort(mudist, decreasing = TRUE)
  cut = sort_mudist[round(nrow(x)*alpha)+1]
  
  oe = (1:nrow(x))[mudist>cut]
  o = union(ow, oe)
  o = o[!is.na(o)]
  
  return(list(weighted.centre = centre, unweighted.centre = unweighted.centre,
              ow = ow, oe = oe, o = o, cluster = cl, iter = iter))
}

weight.update <- function(x, o, cluster, unweighted.centre, l1bound){
  x = x[-o, ]
  K = length(unique(cluster))
  if(ncol(x)%%(K-1) != 0) print("input data dimension error")
  else p = ncol(x)/(K-1)
  cluster = cluster[-o]
  
  bc = bcss(x, cluster)
  bc = matrix(bc, nrow = K-1, ncol = p)
  bc_tilda = colSums(bc)
  
  l1bound = l1bound/sqrt(K-1)
  # argu = bc_tilda
  # sumabs = l1bound
  lam <- Bin.search(bc_tilda, l1bound)
  ws.unscaled <- soft(bc_tilda, lam)
  beta = ws.unscaled/sqrt(sum(ws.unscaled**2))
  ws = beta/sqrt(K-1)
  ret = c()
  for(i in 1:p) ret = c(ret, rep(ws[i], K-1))
  return(ret)
}

### Modification of RSKC, run a single clustering
RSKC.multidim.onerep <- function(x, ncl, l1bound, alpha = 0.1){
  ifstart = 0
  if(ncol(x)%%(ncl-1) != 0) print("input data dimension error")
  else p = ncol(x)/(ncl-1)
  w = rep(1/sqrt(ncol(x)), ncol(x))
  w.old = w
  start = Sys.time()
  while(ifstart == 0 | sum((w-w.old)**2) > 1e-6){
    ifstart = ifstart + 1
    end = Sys.time()
    w.old = w
    trim = trimmed.Kmeans(x = x, w = w, ncl = ncl, alpha = alpha)
    o = trim$o
    un.centre = trim$unweighted.centre
    cl = trim$cluster
    w = weight.update(x, o = o, cluster = cl, unweighted.centre = un.centre, l1bound = l1bound)
    if(ifstart > 1000 | end-start > 10) break
  }
  w = w[seq(from = ncl-1, to = ncol(x), by = ncl-1)]*sqrt(ncl-1)
  return(list(Cs = cl, w = w))
}

### Function to detect mode value
detect.mode <- function(modevec){
  uniq = sum(modevec>0)
  trial = sum(modevec)
  detect = modevec[modevec>trial*0.7]
  if(length(detect) == 0) return(F)
  else return(T)
}


### Modification of RSKC, run $nrep$ clustering and find mode value
RSKC.multidim <- function(x, ncl, l1bound = 3, alpha = 0.1, nrep = 20){
  modevec = rep(0, nrep)
  for(i in 1:nrep){
    b = RSKC.multidim.onerep(x, ncl, l1bound, alpha)
    cl = b$Cs
    w = b$w
    if(i == 1){
      cl_list = list(cl)
      w_list = list(w)
      modevec[1] = 1
    }
    else{
      whatif = 0
      for(j in 1:length(cl_list)){
        if(perform_score(cl_list[[j]], cl, "CSR") == 1){
          modevec[j] = modevec[j] + 1
          whatif = whatif + 1
        }
      }
      if(whatif == 0){
        modevec[length(cl_list)+1] = 1
        cl_list = c(cl_list, list(cl))
        w_list = c(w_list, list(w))
      }
    }
    if(i >= 3 & detect.mode(modevec)) break
  }
  mode = (1:nrep)[modevec == max(modevec)]
  cat("RSFclust2 mode vector :", modevec, "\n")
  if(length(mode)>1){
    mode = mode[1]
    print("mode unstable")
  }
  cl = cl_list[[mode]]
  w = w_list[[mode]]
  return(list(labels = cl, weights = w))
}


### Return (length(curve_list)) * (|Q|+|E|)(K-1) FPC score matrix
FPC.score.matrix <- function(curve_list, quant.levels, expec.levels, K, 
                             num.knots = 10, n_breaks = 100){
  t = seq(0, 1, length.out = n_breaks+1)[1:n_breaks]
  q.num = length(quant.levels)
  e.num = length(expec.levels)
  
  quant.coef = list()
  quant.vec = list()
  expec.coef = list()
  expec.vec = list()
  
  if(q.num > 0){
    for(i in 1:q.num){
      quant.coef[[i]] = matrix(nrow = (num.knots+3), ncol = length(curve_list))
      quant.vec[[i]] = matrix(nrow = length(t), ncol = length(curve_list))
    }
  }
  if(e.num > 0){
    for(j in 1:e.num){
      expec.coef[[j]] = matrix(nrow = (num.knots+3), ncol = length(curve_list))
      expec.vec[[j]] = matrix(nrow = length(t), ncol = length(curve_list))
    }
  }

  for(j in 1:length(curve_list)){
    y = curve_list[[j]]
    if(j == 1){
      quant.lambda = c()
      expec.lambda = c()
      if(q.num > 0){
        for(i in 1:q.num) quant.lambda = c(quant.lambda, quant.akcv.select(y, quant.levels[i], num.knots = num.knots, n_breaks = n_breaks))
      }
      if(e.num > 0){
        for(i in 1:e.num) expec.lambda = c(expec.lambda, gcv.select(y, expec.levels[i], num.knots = num.knots, n_breaks = n_breaks))
      }
    }
    if(q.num > 0){
      for(i in 1:q.num){
        quant.coef[[i]][,j] = quant.curve.pspline(y, quant.levels[i], num.knots = num.knots, n_breaks = n_breaks, lambda = quant.lambda[i])
        quant.vec[[i]][,j] = basis_value_vector(quant.coef[[i]][,j], num.knots = num.knots, n_breaks = n_breaks)
      }
    }
    if(e.num > 0){
      for(i in 1:e.num){
        expec.coef[[i]][,j] = expec.curve.pspline(y, expec.levels[i], num.knots = num.knots, n_breaks = n_breaks, lambda = expec.lambda[i])
        expec.vec[[i]][,j] = basis_value_vector(expec.coef[[i]][,j], num.knots = num.knots, n_breaks = n_breaks)
      }
    }
  }
  
  basis = create.bspline.basis(rangeval = c(0, 1), num.knots+3)
  argvals = matrix(rep(t, length(curve_list)), nrow = length(t), ncol = length(curve_list))
  quant.obj = list()
  quant.score = list()
  expec.obj = list()
  expec.score = list()
  score.matrix = c()
  
  if(q.num > 0){
    for(i in 1:q.num){
      quant.obj[[i]] = Data2fd(argvals = argvals, y = quant.vec[[i]], basisobj = basis, lambda = 0.1)
      pca = pca.fd(quant.obj[[i]], nharm = K)
      quant.score[[i]] = pca$scores[, 1:(K-1)]
      score.matrix = cbind(score.matrix, quant.score[[i]])
    }
  }
  if(e.num > 0){
    for(i in 1:e.num){
      expec.obj[[i]] = Data2fd(argvals = argvals, y = expec.vec[[i]], basisobj = basis, lambda = 0.1)
      pca = pca.fd(expec.obj[[i]], nharm = K)
      expec.score[[i]] = pca$scores[, 1:(K-1)]
      score.matrix = cbind(score.matrix, expec.score[[i]])
    }
  }
  return(score.matrix)
}



