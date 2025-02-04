# Canadian weather data
library("maps")

# Settings
num.knots = 10
n_breaks = 12
b.sp.c = bs(seq(0, 1, by = 1/n_breaks), knots = seq(0, 1, length.out = num.knots+1))
b.basis.c = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))
tc = seq(from = 0, to = 1, length.out = 13)[-13]

# Creating canada.curve.list
canada <- CanadianWeather$dailyAv[,,1]
canada.month <- matrix(nrow=12, ncol=35) # row : month, col : place
canada.month[1,] <- apply(canada[c(1:31),], 2, mean)
canada.month[2,] <- apply(canada[seq(from=31+1, by=1, length=28),], 2, mean)
canada.month[3,] <- apply(canada[seq(from=31+28+1, by=1, length=31),], 2, mean)
canada.month[4,] <- apply(canada[seq(from=31+28+31+1, by=1, length=30),], 2, mean)
canada.month[5,] <- apply(canada[seq(from=31+28+31+30+1, by=1, length=31),], 2, mean)
canada.month[6,] <- apply(canada[seq(from=31+28+31+30+31+1, by=1, length=30),], 2, mean)
canada.month[7,] <- apply(canada[seq(from=31+28+31+30+31+30+1, by=1, length=31),], 2, mean)
canada.month[8,] <- apply(canada[seq(from=31+28+31+30+31+30+31+1, by=1, length=31),], 2, mean)
canada.month[9,] <- apply(canada[seq(from=31+28+31+30+31+30+31+31+1, by=1, length=30),], 2, mean)
canada.month[10,] <- apply(canada[seq(from=31+28+31+30+31+30+31+31+30+1, by=1, length=31),], 2, mean)
canada.month[11,] <- apply(canada[seq(from=31+28+31+30+31+30+31+31+30+31+1, by=1, length=30),], 2, mean)
canada.month[12,] <- apply(canada[seq(from=31+28+31+30+31+30+31+31+30+31+30+1, by=1, length=31),], 2, mean)

for(i in 1:35){
  if(i == 1) canada.curve.list = list(canada.month[, i])
  else canada.curve.list = c(canada.curve.list, list(canada.month[, i]))
}

############# K = 2, 3 (RSFclust, Funclust, L2-dist, Baseclust, Exfunclust, FPCAC) #############

bmatrix = matrix(nrow = length(canada.curve.list), ncol = num.knots+3)
vec.original = matrix(nrow = n_breaks, ncol = length(canada.curve.list))
coef.mean = t(bmatrix)
coef.q0.3 = coef.mean

for(j in 1:length(canada.curve.list)){
  y = canada.curve.list[[j]]
  vec.original[,j] = canada.curve.list[[j]]
  if(j == 1){
    lambda.q0.3 = quant.akcv.select(y, 0.25, num.knots = num.knots, n_breaks = n_breaks)
    lambda.mean = gcv.select(y, 0.5, num.knots = num.knots, n_breaks = n_breaks)
  }
  coef.mean[,j] = expec.curve.pspline(y, 0.5, num.knots = num.knots, n_breaks = n_breaks, lambda = lambda.mean)
  coef.q0.3[,j] = quant.curve.pspline(y, 0.25, num.knots = num.knots, n_breaks = n_breaks, lambda = lambda.q0.3)
  bmat = as.matrix(nearPD(t(b.basis.c)%*%b.basis.c)$mat)
  bmatrix[j,] = solve(bmat, t(b.basis.c)%*%y)
}

N = length(canada.curve.list)
domain = seq(from=0, to=1, by=0.001)
obs.domain = seq(from=0, to=1, by=1/n_breaks)[-(n_breaks+1)]
knots <- seq(from=0, to=1, by=1/num.knots)
b.expand <- bsplineS(domain, knots)
simul.norm <- matrix(nrow=length(obs.domain), ncol=length(canada.curve.list))
for(ss in 1:length(canada.curve.list)) simul.norm[,ss] = canada.curve.list[[ss]]
bsp = create.bspline.basis(breaks=seq(0, 1, by=1/num.knots))
coeff.tau0.5 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.5, n_breaks = n_breaks, num.knots = num.knots)
coeff.tau0.01 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.01, n_breaks = n_breaks, num.knots = num.knots)
coeff.tau0.99 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.99, n_breaks = n_breaks, num.knots = num.knots)
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
  w11 <- w1[i2,]/(w1[i2,]+w2[i2,]); w22 <- w2[i2,]/(w1[i2,]+w2[i2,]); w1[i2,] <- w11; w2[i2,] <- w22
  poo1 <- c(w1[i2,1], rep(w1[i2,], each=1000/num.knots)); poo2 <- c(w2[i2,1], rep(w2[i2,], each=1000/num.knots))
  weighted.expec[,i2] <-poo1*fit.expec0.99[,i2]+poo2*fit.expec0.01[,i2]
}
coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
for(is in 1:N){
  coeff.tau.al[is,] <- solve(t(b.expand)%*%b.expand, t(b.expand)%*%weighted.expec[,is])
}   

# RSFclust
score.matrix2 = FPC.score.matrix(canada.curve.list, c(0.25, 0.5, 0.75), c(0.5), K = 2, num.knots, n_breaks)
RSFclust.2 = RSKC(score.matrix2, 2, alpha = 0.1, L1 = 3) # 2 clusters
r2 = RSFclust.2$labels

score.matrix3 = FPC.score.matrix(canada.curve.list, c(0.25, 0.5, 0.75), c(0.5), K = 3, num.knots, n_breaks)
RSFclust1.3 = RSKC(score.matrix3, 3, alpha = 0.1, L1 = 3) # 3 clusters
RSFclust2.3 = RSKC.multidim(score.matrix3, ncl = 3, l1bound = 3, alpha = 0.1) # 3 clusters
r1.3 = RSFclust1.3$labels
r2.3 = RSFclust2.3$labels

# Exfunclust
efckim2 <- kmeans(coeff.tau.al, centers = 2, nstart = 10)
e2 <- efckim2$cluster
efckim3 <- kmeans(coeff.tau.al, centers = 3, nstart = 10)
e3 <- efckim3$cluster

# Funclust (Unavailable now)
fd.obj <- fd(coef=coeff.tau0.5, basisobj=bsp)
f2 <- funclust(fd.obj, K=2)$cls

# FPCAC, Baseclust, L2-dist
fpc2 = fpcac(vec.original, K = 2, alpha = 0)$clusters
b2 = kmeans(bmatrix, centers = 2, nstart = 10)$cluster
l2 = distance_clustering(2, coef.mean, num.knots, n_breaks)

fpc3 = fpcac(vec.original, K = 3, alpha = 0)$clusters
b3 = kmeans(bmatrix, centers = 3, nstart = 10)$cluster
l3 = distance_clustering(3, coef.mean, num.knots, n_breaks)

# Match colors (1 = blue, 2 = red)
match_col_2 <- function(x){
  ret = x
  if(sum(x == 1) < sum(x == 2)){
    ret[x == 1] = 2
    ret[x == 2] = 1
  }
  return(ret)
}

r2 = match_col_2(r2)
e2 = match_col_2(e2)
fpc2 = match_col_2(fpc2)
b2 = match_col_2(b2)
l2 = match_col_2(l2)

# Match colors (1 = blue, 2 = red, 3 = green)
match_col_3 <- function(x){
  ret = x
  mean1 = mean(location[,1][x == 1])
  mean2 = mean(location[,1][x == 2])
  mean3 = mean(location[,1][x == 3])
  ord = order(c(mean1, mean2, mean3)) # ord = c(greennum, bluenum, rednum)
  ret[x == ord[1]] = 3
  ret[x == ord[2]] = 1
  ret[x == ord[3]] = 2
  return(ret)
}

r1.3 = match_col_3(r1.3)
r2.3 = match_col_3(r2.3)
e3 = match_col_3(e3)
fpc3 = match_col_3(fpc3)
b3 = match_col_3(b3)
l3 = match_col_3(l3)

### Figure 8-1
location = CanadianWeather$coordinates
location[, 2] <- 360 - location[, 2]
poo1 = b2 == 1
poo2 = b2 == 2

maps::map("world2", ylim=c(40, 80), xlim=c(200, 320))
points(location[poo1,2], location[poo1,1], col="blue", pch = 4)
points(location[poo2,2], location[poo2,1], col="red")
title("Baseclust result")


### Figure 8-2
location = CanadianWeather$coordinates
location[, 2] <- 360 - location[, 2]
poo1 = f2 == 1
poo2 = f2 == 2

maps::map("world2", ylim=c(40, 80), xlim=c(200, 320))
points(location[poo1,2], location[poo1,1], col="blue", pch = 4)
points(location[poo2,2], location[poo2,1], col="red")
title("Funclust result")


# Figure 11
poo1 = r2.3 == 1 # blue
poo2 = r2.3 == 2 # red
poo3 = r2.3 == 3 # green

pdf("figure11.pdf")
maps::map("world2", ylim=c(40, 80), xlim=c(200, 320))
points(location[poo1,2], location[poo1,1], col="blue", pch = 4)
points(location[poo2,2], location[poo2,1], col="red")
points(location[poo3,2], location[poo3,1], col="green", pch = 2)
title("RSFclust2 result")
dev.off()


### Figure 9, 10, 12

# Locations where funclust != RSFclust : Dawson, Yellowknife, Uranium Cty 

c2tot = r2
c2tot[rownames(location) %in% c("Dawson", "Yellowknife", "Uranium Cty")] = 3
c2tot

ID2 = c()
for(l in 1:length(canada.curve.list)){
  ID2 = c(ID2, rep(l, 100))
}
ID2 <- as.factor(ID2)

vec0.3.2 = c()
for(l in 1:length(canada.curve.list)) vec0.3.2 = c(vec0.3.2, basis_value_vector(coef.q0.3[,l]))
vec.y = as.vector(vec.original)
c_input2 = as.vector(matrix(rep(c2tot, 100), byrow = T, nrow = 100))
c_input2 = as.factor(c_input2)
c2tot = as.factor(c2tot)
c3r <- as.factor(r2.3)
x_grid2 = rep((1:100)*12/100, length(canada.curve.list))

dataset1.2 = data.frame(ID2, x_grid2, c_input2, vec0.3.2) # To plot 0.25th quantile curves

ID = c()
for(l in 1:length(canada.curve.list)){
  ID = c(ID, rep(l, 12))
}
vec.q0.3 = basis_value_vector(coef.q0.3, num.knots = num.knots, n_breaks = n_breaks)
vec0.3 = as.vector(vec.q0.3)
vec.y = as.vector(vec.original)
c_input = as.vector(matrix(rep(c2tot, 12), byrow = T, nrow = 12))
ID <- as.factor(ID)
c_input <- as.factor(c_input)
x_grid = rep(1:12, length(canada.curve.list))

dataset1.1 = data.frame(ID, x_grid, c_input, vec.y) # To plot raw data curves
dataset2 = data.frame(c2tot, c3r, score.matrix3)
colnames(dataset2)[3:10] = c("qscore.0.3", "qscore.0.32", "qmedian.score", "qmedian.score2",
                             "qscore.0.7", "qscore.0.72", "mean.score", "mean.score2")
# To plot FPC1, FPC2 score scatterplots

### Figure 9-1
pdf("figure9-1.pdf", width = 5, height = 5)
ggplot(data = dataset1.1, mapping = aes(x = x_grid, y = vec.y, group = ID, colour = c_input)) + 
  geom_line(size = 0.5) +
  geom_point(size = 1, shape = 15) +
  scale_color_manual(values = c('1'='blue', '2'='red', '3'='black')) +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  labs(title = "Raw data", x = "Month", y = "Temperature") +
  theme_test() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
dev.off()

### Figure 9-2
pdf("figure9-2.pdf", width = 5, height = 5)
ggplot(data = dataset1.2, mapping = aes(x = x_grid2, y = vec0.3.2, group = ID2, colour = c_input2)) + 
  geom_line(size = 0.5) +
  scale_color_manual(values = c('1'='blue', '2'='red', '3'='black')) +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  labs(title = "0.25th quantile curves", x = "Month", y = "Temperature") +
  theme_test() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
dev.off()

### Figure 10
pdf("figure10.pdf", width = 4.5, height = 4.5)
ggplot(data = dataset2, mapping = aes(x = qscore.0.3, y = qscore.0.32, colour = c2tot)) +
  geom_point(size = 2, shape = 16) +
  scale_color_manual(values = c('1'='blue', '2'='red', '3'='black')) +
  labs(title = "FPC scores of 0.25th quantile curves", x = "FPC 1 score", y = "FPC 2 score") +
  theme_test() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-25, 10)) +
  scale_y_continuous(limits = c(-10, 10))
dev.off()

### Figure 12-1
pdf("figure12-1.pdf", width = 4.5, height = 4.5)
ggplot(data = dataset2, mapping = aes(x = qscore.0.3, y = qscore.0.32, colour = c3r)) +
  geom_point(size = 2, shape = 16) +
  scale_color_manual(values = c('1'='blue', '2'='red','3'='green')) +
  labs(title = "FPC scores of 0.25th quantile curves", x = "FPC 1 score", y = "FPC 2 score") +
  theme_test() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-25, 10)) +
  scale_y_continuous(limits = c(-10, 10))
dev.off()

### Figure 12-2
pdf("figure12-2.pdf", width = 4.5, height = 4.5)
ggplot(data = dataset2, mapping = aes(x = qscore.0.7, y = qscore.0.72, colour = c3r)) +
  geom_point(size = 2, shape = 16) +
  scale_color_manual(values = c('1'='blue', '2'='red','3'='green')) +
  labs(title = "FPC scores of 0.75th quantile curves", x = "FPC 1 score", y = "FPC 2 score") +
  theme_test() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-25, 10)) +
  scale_y_continuous(limits = c(-10, 10))
dev.off()



