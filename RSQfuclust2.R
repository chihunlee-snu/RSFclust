### RSQfunclust2 code

l2n <- function(vec){
  return(sqrt(sum(vec**2)))
}

Bin.search <- function (argu, sumabs){
  if (l2n(argu) == 0 || sum(abs(argu/l2n(argu))) <= sumabs) 
    return(0)
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
  old.centre = z[sample(1:nrow(z), ncl), ]
  centre = old.centre
  ifstart = 0
  iter = 0
  start = Sys.time()
  
  while(ifstart == 0 | sum((old.centre - centre)**2) > 1e-6){
    ifstart = 1
    end = Sys.time()
    if(end-start > 5) break
    iter = iter + 1
    old.centre = centre
    kobj = kmeans(z, centers = old.centre)
    cl = kobj$cluster
    centre = kobj$centers
    mudist = vector(length = nrow(x))
    for(i in 1:nrow(x)) mudist[i] = sum((centre[cl[i], ] - z[i, ])**2)
    
    sort_mudist = sort(mudist, decreasing = TRUE)
    cut = sort_mudist[round(nrow(x)*alpha)]
    for(k in 1:ncl){
      if(sum(cl==k & mudist<cut) > 1) centre[k, ] = colMeans(z[(cl==k & mudist<cut),])
      else centre[k, ] = z[(cl==k & mudist<cut),]
    }
    ow = (1:nrow(x))[mudist>=cut]
  }
  unweighted.centre = sweep(centre, 2, sqrt(w), "/")
  
  mudist = vector(length = nrow(x))
  for(i in 1:nrow(x)) mudist[i] = sum((unweighted.centre[cl[i], ] - x[i, ])**2)
  
  sort_mudist = sort(mudist, decreasing = TRUE)
  cut = sort_mudist[round(nrow(x)*alpha)]
  
  oe = (1:nrow(x))[mudist>=cut]
  o = union(ow, oe)
  
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
  lam <- Bin.search(bc_tilda, l1bound)
  ws.unscaled <- soft(bc_tilda, lam)
  beta = ws.unscaled/sqrt(sum(ws.unscaled)^2)
  ws = beta/sqrt(K-1)
  ret = c()
  for(i in 1:p) ret = c(ret, rep(ws[i], K-1))
  return(ret)
}

RSKC.multidim.onerep <- function(x, ncl, l1bound, alpha = 0.1){
  ifstart = 0
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
  w = w[c(1, ncl, 2*ncl-1, 3*ncl-2)]*sqrt(ncl-1)
  return(list(Cs = cl, w = w))
}

### Modification of RSKC
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
  }
  mode = (1:nrep)[modevec == max(modevec)]
  cat("RSQfunclust2 mode vector :", modevec, "\n")
  if(length(mode)>1){
    mode = mode[1]
    print("mode unstable")
  }
  cl = cl_list[[mode]]
  w = w_list[[mode]]
  return(list(labels = cl, weights = w))
}





