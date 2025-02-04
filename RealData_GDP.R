### Economic growth rate


### Draw the entire world plot
world <- map_data("world")
worldplot <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  theme_test() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank())
worldplot


### Preprocessing of the gdp growth rate data
library(readxl)
library(stringr)
library(dplyr)
gdp <- read_excel("realgdprate.xls")
gdp <- gdp[, -c(45, 46, 47, 48, 49, 50)]
gdp <- gdp[-c(1, 230, 231), ]
gdp <- gdp[1:196, ]

unlist(lapply(sapply(gdp, FUN = "class"), FUN = "[", 1))
gdp[, 2:ncol(gdp)] = sapply(gdp[, 2:ncol(gdp)], FUN = "as.numeric")
unlist(lapply(sapply(gdp, FUN = "class"), FUN = "[", 1))
colnames(gdp)[1] <- "region"
head(gdp)
tail(gdp)


### Matching the country name in gdp data and world map data
colnames(gdp)
gdp <- gdp %>% mutate(region = recode(str_trim(region), "United States" = "USA",
                                      "United Kingdom" = "UK",
                                      "Korea, Republic of" = "South Korea",
                                      "Congo, Dem. Rep. of the" = "Democratic Republic of the Congo",
                                      "Congo, Republic of" = "Republic of Congo",
                                      "Kyrgyz Republic" = "Kyrgyzstan",
                                      "Russian Federation" = "Russia",
                                      "Bahamas, The" = "Bahamas",
                                      "Brunei Darussalam" = "Brunei",
                                      "Cabo Verde" = "Cape Verde",
                                      "China, People's Republic of" = "China",
                                      "Côte d'Ivoire" = "Ivory Coast",
                                      "Gambia, The" = "Gambia",
                                      "Lao P.D.R." = "Laos",
                                      "Micronesia, Fed. States of" = "Micronesia",
                                      "Saint Vincent and the Grenadines" = "Saint Vincent",
                                      "Slovak Republic" = "Slovakia",
                                      "São Tomé and Príncipe" = "Sao Tome and Principe",
                                      "Taiwan Province of China" = "Taiwan",
                                      "Trinidad and Tobago" = "Tobago"))


### Asian countries
asia.region <- unique(c('Azerbaijan', 'Bahrain', 'Bangladesh', 'Bhutan', 'Brunei', 'India', 'Japan', 'Kyrgyzstan', 'Mongolia', 'Myanmar',
                        'Tajikistan', 'Turkmenistan', 'Vietnam', 'Cambodia', 'Hong Kong SAR', 'Indonesia', 'Israel', 'Jordan', 'South Korea',
                        'Laos', 'Philippines', 'Qatar', 'Sri Lanka', 'United Arab Emirates', 'Uzbekistan', 'Armenia', 'Georgia',
                        'Iran', 'Kazakhstan', 'Malaysia', 'Nepal', 'Oman', 'Saudi Arabia', 'Thailand', 'Turkey', 'Yemen',
                        'China', 'Iraq', 'Maldives', 'Pakistan', 'Russia', 'Singapore', 'Taiwan', 'Afghanistan', 'North Korea', 'Syria', 'Lebanon', 'Cyprus'))

gdp$region[!gdp$region %in% world$region]
gdp$region[182] <- "Turkey"
gdp$region[!gdp$region %in% world$region] # Now Turkiye included
# 196 countries, 44 columns (1 + 43years) -> Use 2001~2022 (23:44 columns)

gdp <- as.data.frame(gdp)
gdp[197, ] <- c("North Korea", rep(NA, 43))
asia_index <- gdp$region %in% asia.region
asia.gdp <- gdp[asia_index, ]

sum(rowSums(is.na(asia.gdp))==0) # 29/48 countries
sum(rowSums(is.na(asia.gdp[,23:44]))==0) # 44/48 countries

gdp$region[!gdp$region %in% world$region] # Asian countries now all have same country names with world data


### European countries
europe.region <- unique(c('Belarus', 'Czech Republic', 'France', 'Germany', 'Iceland', 'Lithuania', 'Malta', 'Montenegro', 'North Macedonia',
                          'Serbia', 'Slovakia', 'Albania', 'Belgium', 'Ireland', 'Luxembourg', 'Norway', 'Slovenia', 'UK', 'Bosnia and Herzegovina',
                          'Bulgaria', 'Croatia', 'Denmark', 'Greece', 'Moldova', 'Poland', 'Romania', 'Spain', 'Sweden', 'Andorra', 'Armenia',
                          'Azerbaijan', 'Estonia', 'Finland', 'Georgia', 'Hungary', 'Italy', 'Latvia', 'Netherlands', 'Portugal', 'Switzerland',
                          'Ukraine', 'Austria'))

eurasia.region <- unique(c(asia.region, europe.region))
europe_index <- gdp$region %in% europe.region
europe.gdp <- gdp[europe_index, ]

sum(rowSums(is.na(europe.gdp))==0) # 23/42 countries
sum(rowSums(is.na(europe.gdp[, 23:44]))==0) # 42/42 countries
# Asian countries now all have same country names with world data


### Creating gdp growth rate data of Eurasia countries
eurasia_index <- gdp$region %in% eurasia.region
eurasia.gdp <- gdp[eurasia_index, ]
sum(rowSums(is.na(eurasia.gdp))==0) # 52/87 countries
sum(rowSums(is.na(eurasia.gdp[, 23:44]))==0) # 83/87 countries

g.data <- eurasia.gdp[,c(1, 23:44)]
g.data # 87 eurasia countries
unlist(lapply(sapply(g.data, FUN = "class"), FUN = "[", 1))
g.data[, 2:ncol(g.data)] = sapply(g.data[, 2:ncol(g.data)], FUN = "as.numeric")

g.data$region[rowSums(is.na(g.data))>0]
# Afghanistan, Lebanon, Syria, North Korea -> cluster number 0

g.data.valid <- g.data[rowSums(is.na(g.data))==0, ]
unlist(lapply(sapply(g.data.valid, FUN = "class"), FUN = "[", 1))
g.data.valid[, 2:ncol(g.data.valid)] = sapply(g.data.valid[, 2:ncol(g.data.valid)], FUN = "as.numeric")
valid.index = rowSums(is.na(g.data))==0
nrow(g.data.valid) # 83 eurasia countries where all data in 2001~2022 is valid


### Creating Eurasia map
asia_map <- map_data(map = "world", region = asia.region)
europe_map <- map_data(map = "world", region = europe.region)
eurasia_map <- map_data(map = "world", region = eurasia.region)

asiaplot <- ggplot() +
  geom_polygon(data = asia_map, aes(x=long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  theme_test() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank())
asiaplot

europeplot <- ggplot() +
  geom_polygon(data = europe_map, aes(x=long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  theme_test() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank())
europeplot


eurasiaplot <- ggplot() +
  geom_polygon(data = eurasia_map, aes(x=long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  theme_test() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank())
eurasiaplot

# Finally, we use eurasia data in 2001 ~ 2022 of 83 countries




### Run RSFclust1, RSFclust2, FPCAC, Exfunclust, Baseclust, L2-dist
num.knots = 10
n_breaks = 22
b.sp.w = bs(seq(0, 1, by = 1/n_breaks), knots = seq(0, 1, length.out = num.knots+1))
b.basis.w = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))

tw = seq(from = 0, to = 1, length.out = n_breaks+1)[-(1+n_breaks)]
gdp_matrix = as.matrix(g.data.valid[,-1])

bmatrix = matrix(nrow = nrow(gdp_matrix), ncol = num.knots+3)
coef.mean = t(bmatrix)
coef.q0.75 = t(bmatrix)
coef.q0.25 = t(bmatrix)
vec.original = matrix(nrow = length(tw), ncol = nrow(gdp_matrix))
vec.mean = vec.original
vec.q0.75 = matrix(nrow = 100, ncol = nrow(gdp_matrix))
vec.q0.25 = matrix(nrow = 100, ncol = nrow(gdp_matrix))
gdp.list = list()
for(i in 1:nrow(gdp_matrix)) gdp.list = c(gdp.list, list(as.vector(gdp_matrix[i,])))
score.matrix2 = FPC.score.matrix(gdp.list, c(0.25, 0.5, 0.75), c(0.5), K = 2, num.knots, n_breaks)
score.matrix3 = FPC.score.matrix(gdp.list, c(0.25, 0.5, 0.75), c(0.5), K = 3, num.knots, n_breaks)
score.matrix4 = FPC.score.matrix(gdp.list, c(0.25, 0.5, 0.75), c(0.5), K = 4, num.knots, n_breaks)

for(j in 1:nrow(gdp_matrix)){
  y = gdp_matrix[j,]
  vec.original[,j] = gdp_matrix[j,]
  if(j == 1){
    lambda.q0.75 = quant.akcv.select(y, 0.5, num.knots = num.knots, n_breaks = n_breaks)
    lambda.mean = gcv.select(y, 0.5, num.knots = num.knots, n_breaks = n_breaks)
  }
  coef.mean[,j] = expec.curve.pspline(y, 0.5, num.knots = num.knots, n_breaks = n_breaks, lambda = lambda.mean)
  vec.mean[, j] = basis_value_vector(coef.mean[,j], num.knots = num.knots, n_breaks = n_breaks)
  coef.q0.75[,j] = quant.curve.pspline(y, 0.75, num.knots = num.knots, n_breaks = n_breaks)
  vec.q0.75[,j] = basis_value_vector(coef.q0.75[,j], num.knots = num.knots, n_breaks = 100)
  coef.q0.25[,j] = quant.curve.pspline(y, 0.25, num.knots = num.knots, n_breaks = n_breaks)
  vec.q0.25[,j] = basis_value_vector(coef.q0.25[,j], num.knots = num.knots, n_breaks = 100)
  bmat = as.matrix(nearPD(t(b.basis.w)%*%b.basis.w)$mat)
  bmatrix[j,] = solve(bmat, t(b.basis.w)%*%y)
}

# RSFclust1, RSFclust2 (K = 2, 3, 4)
RSF_2 <- RSKC(score.matrix2, 2, alpha = 0.1, L1 = 2)
RSF2_3 <- RSKC.multidim(score.matrix3, ncl = 3, l1bound = 2)
RSF2_4 <- RSKC.multidim(score.matrix4, ncl = 4, l1bound = 2)
RSF1_3 <- RSKC(score.matrix3, ncl = 3, alpha = 0.1, L1 = 2)
RSF1_4 <- RSKC(score.matrix4, ncl = 4, alpha = 0.1, L1 = 2)

r.2.c <- RSF_2$labels
r.2.w <- RSF_2$weights
r2.3.c <- RSF2_3$labels
r2.3.w <- RSF2_3$weights
r2.4.c <- RSF2_4$labels
r2.4.w <- RSF2_4$weights
r1.3.c <- RSF1_3$labels
r1.3.w <- RSF1_3$weights
r1.4.c <- RSF1_4$labels
r1.4.w <- RSF1_4$weights

# Exfunclust (K = 2, 3, 4)
domain = seq(from=0, to=1, by=0.001)
obs.domain = seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)]
knots <- seq(from=0, to=1, by=1/num.knots)
b.expand <- bsplineS(domain, knots)
simul.norm <- matrix(nrow=length(obs.domain), ncol=nrow(gdp_matrix))
for(ss in 1:nrow(gdp_matrix)) simul.norm[,ss] = gdp_matrix[ss,]

coeff.tau0.5 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.5, n_breaks = n_breaks)
coeff.tau0.01 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.01, n_breaks = n_breaks)
coeff.tau0.99 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.99, n_breaks = n_breaks)
fit.expec0.5 <- b.expand%*%coeff.tau0.5
fit.expec0.01 <- b.expand%*%coeff.tau0.01
fit.expec0.99 <- b.expand%*%coeff.tau0.99

w1 <- w2 <- matrix(nrow=nrow(gdp_matrix), ncol=num.knots)
weighted.expec <- matrix(nrow=length(domain), ncol=nrow(gdp_matrix))
for(i2 in 1:nrow(gdp_matrix)){
  for(is in 1:num.knots){
    w1[i2,is] <- abs(0.001*sum(fit.expec0.99[c((1000/num.knots*(is-1)+1):(1000/num.knots*is+1)), i2]-fit.expec0.5[c((1000/num.knots*(is-1)+1):(1000/num.knots*is+1)), i2]))
    w2[i2,is] <- abs(0.001*sum(fit.expec0.5[c((1000/num.knots*(is-1)+1):(1000/num.knots*is+1)), i2]-fit.expec0.01[c((1000/num.knots*(is-1)+1):(1000/num.knots*is+1)), i2]))
  }
  # w1[i2,] <- 1/w1[i2,]; w2[i2,] <- 1/w2[i2,]
  w11 <- w1[i2,]/(w1[i2,]+w2[i2,]); w22 <- w2[i2,]/(w1[i2,]+w2[i2,]); w1[i2,] <- w11; w2[i2,] <- w22
  poo1 <- c(w1[i2,1], rep(w1[i2,], each=1000/num.knots)); poo2 <- c(w2[i2,1], rep(w2[i2,], each=1000/num.knots))
  weighted.expec[,i2] <-poo1*fit.expec0.99[,i2]+poo2*fit.expec0.01[,i2]
}
coeff.tau.al <- matrix(nrow=nrow(gdp_matrix), ncol=num.knots+3)
for(is in 1:nrow(gdp_matrix)){
  coeff.tau.al[is,] <- solve(t(b.expand)%*%b.expand, t(b.expand)%*%weighted.expec[,is])
}     
e.2.c <- kmeans(coeff.tau.al, centers=2, nstart = 10)$cluster
e.3.c <- kmeans(coeff.tau.al, centers=3, nstart = 10)$cluster
e.4.c <- kmeans(coeff.tau.al, centers=4, nstart = 10)$cluster


# FPCAC (K = 2, 3, 4)
fpc.2.c <- fpcac(t(gdp_matrix), K = 2, alpha = 0)$clusters
fpc.3.c <- fpcac(t(gdp_matrix), K = 3, alpha = 0)$clusters
fpc.4.c <- fpcac(t(gdp_matrix), K = 4, alpha = 0)$clusters

# Baseclust (K = 2, 3, 4)
b.2.c = kmeans(bmatrix, centers = 2, nstart = 10)$cluster
b.3.c = kmeans(bmatrix, centers = 3, nstart = 10)$cluster
b.4.c = kmeans(bmatrix, centers = 4, nstart = 10)$cluster

# L2-dist (K = 2, 3, 4)
l.2.c = distance_clustering(2, coef.mean, num.knots, n_breaks)
l.3.c = distance_clustering(3, coef.mean, num.knots, n_breaks)
l.4.c = distance_clustering(4, coef.mean, num.knots, n_breaks)

### Choosing the number of the clusters

# Silhouette coefficient (The bigger, the better)

mean(silhouette(r1.3.c, dist(sweep(score.matrix3, 2, sqrt(r1.3.w), "*")))[,3]) # 0.473
mean(silhouette(r1.4.c, dist(sweep(score.matrix4, 2, sqrt(r1.4.w), "*")))[,3]) # 0.388

mean(silhouette(r.2.c, dist(sweep(score.matrix2, 2, sqrt(r.2.w), "*")))[,3]) # 0.413
mean(silhouette(r2.3.c, dist(sweep(score.matrix3, 2, sqrt(r2.3.w), "*")))[,3]) # 0.372
mean(silhouette(r2.4.c, dist(sweep(score.matrix4, 2, sqrt(r2.4.w), "*")))[,3]) # 0.395

mean(silhouette(fpc.2.c, dist(t(vec.mean)))[,3]) # 0.429
mean(silhouette(fpc.3.c, dist(t(vec.mean)))[,3]) # 0.339
mean(silhouette(fpc.4.c, dist(t(vec.mean)))[,3]) # 0.357

mean(silhouette(b.2.c, dist(bmatrix))[,3]) # 0.873
mean(silhouette(b.3.c, dist(bmatrix))[,3]) # 0.484
mean(silhouette(b.4.c, dist(bmatrix))[,3]) # 0.500

mean(silhouette(e.2.c, dist(coeff.tau.al))[,3]) # 0.791
mean(silhouette(e.3.c, dist(coeff.tau.al))[,3]) # 0.278
mean(silhouette(e.4.c, dist(coeff.tau.al))[,3]) # 0.267

mean(silhouette(l.2.c, dist(t(vec.mean)))[,3]) # 0.486
mean(silhouette(l.3.c, dist(t(vec.mean)))[,3]) # 0.357
mean(silhouette(l.4.c, dist(t(vec.mean)))[,3]) # 0.349


### Matching the color of the map

convert_label_3 <- function(x){
  goldnum = x[g.data.valid$region == 'Russia'] # Group 1
  orangenum = x[g.data.valid$region == 'India'] # Group 2
  khakinum = x[g.data.valid$region == 'France'] # Group 3
  x1 = x
  x1[x == goldnum] = 1
  x1[x == orangenum] = 2
  x1[x == khakinum] = 3
  return(x1)
}

convert_label_4 <- function(x){
  rednum = x[g.data.valid$region == 'Myanmar'] # Group 1
  orangenum = x[g.data.valid$region == 'India'] # Group 2
  khakinum = x[g.data.valid$region == 'France'] # Group 3
  goldnum = x[g.data.valid$region == 'Russia'] # Group 4
  x1 = x
  x1[x == rednum] = 1
  x1[x == orangenum] = 2
  x1[x == khakinum] = 3
  x1[x == goldnum] = 4
  return(x1)
}

convert_label_base <- function(x){
  khakinum = x[g.data.valid$region == 'France'] # Group 1
  orangenum = which.min(table(x))[[1]] # Group 2
  rednum = x[g.data.valid$region == 'Myanmar'] # Group 3
  goldnum = x[g.data.valid$region == 'Russia'] # Group 4
  x1 = x
  x1[x == khakinum] = 1
  x1[x == orangenum] = 2
  x1[x == rednum] = 3
  x1[x == goldnum] = 4
  return(x1)
}

convert_label_exf <- function(x){
  goldnum = x[g.data.valid$region == 'France'] # Group 1
  orangenum = x[g.data.valid$region == 'China'] # Group 2
  rednum = x[g.data.valid$region == 'Iraq'] # Group 3
  x1 = x
  x1[x == goldnum] = 1
  x1[x == orangenum] = 2
  x1[x == rednum] = 3
  return(x1)
}

fpc.4.c = convert_label_4(fpc.4.c)
r2.4.c = convert_label_4(r2.4.c)
r1.3.c = convert_label_3(r1.3.c)
l.3.c = convert_label_3(l.3.c)
b.4.c = convert_label_base(b.4.c)
e.3.c = convert_label_exf(e.3.c)


### Adding the cluster assignment to g.data

g.data.valid[,"r1.3"] = as.factor(r1.3.c)
g.data.valid[,"r2.4"] = as.factor(r2.4.c)
g.data.valid[,"e3"] = as.factor(e.3.c)
g.data.valid[,"fpc4"] = as.factor(fpc.4.c)
g.data.valid[,"l3"] = as.factor(l.3.c)
g.data.valid[,"b4"] = as.factor(b.4.c)

eurasiasubset <- inner_join(world, g.data.valid, by = "region")
head(eurasiasubset)


### Figure theme
plain <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank(),
  panel.background = element_rect(fill = "white"),
  plot.title = element_text(hjust = 0.5),
  legend.position = "bottom",
  legend.text = element_text(size = 10),
  legend.key.spacing.x = unit(0.5, "cm")
)


### Figure 4-rsf1 (RSFclust1)
pdf(file = "figure4-rsf1.pdf", width = 5.5, height = 3.9)
eurasiagdp.r1.3 <- ggplot(data = eurasiasubset, mapping = aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3) +
  geom_polygon(aes(fill = r1.3)) + 
  scale_fill_manual(values = c('1'='gold', '2'='orange', '3'='khaki3'),
                    guide = guide_legend(title = NULL, byrow = TRUE)) +
  ggtitle("RSFclust1 result (K=3)") +
  plain
eurasiagdp.r1.3
dev.off()

### Figure 4-rsf2 (RSFclust2)
pdf(file = "figure4-rsf2.pdf", width = 5.5, height = 3.9)
eurasiagdp.r2.4 <- ggplot(data = eurasiasubset, mapping = aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3) +
  geom_polygon(aes(fill = r2.4)) + 
  scale_fill_manual(values = c('1'='indianred2', '2'='orange', '3'='khaki3', '4'='gold'),
                    guide = guide_legend(title = NULL, byrow = TRUE)) +
  ggtitle("RSFclust2 result (K=4)") +
  plain
eurasiagdp.r2.4
dev.off()

### Figure 4-exf (Exfunclust)
pdf(file = "figure4-exf.pdf", width = 5.5, height = 3.9)
eurasiagdp.e3 <- ggplot(data = eurasiasubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = e3)) +
  scale_fill_manual(values = c('1'='gold', '2'='orange', '3'='indianred2'),
                    guide = guide_legend(title = NULL, byrow = TRUE)) +
  ggtitle("Exfunclust result (K=3)") +
  plain
eurasiagdp.e3
dev.off()

### Figure 4-fpc (FPCAC)
pdf(file = "figure4-fpcac.pdf", width = 5.5, height = 3.9)
eurasiagdp.f4 <- ggplot(data = eurasiasubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = fpc4)) +
  scale_fill_manual(values = c('1'='indianred2', '2'='orange', '3'='khaki3', '4'='gold'),
                    guide = guide_legend(title = NULL, byrow = TRUE)) +
  ggtitle("FPCAC result (K=4)") +
  plain
eurasiagdp.f4
dev.off()

### Figure 4-base (Baseclust)
pdf(file = "figure4-base.pdf", width = 5.5, height = 3.9)
eurasiagdp.b4 <- ggplot(data = eurasiasubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = b4)) +
  scale_fill_manual(values = c('1'='khaki3', '2'='orange', '3'='indianred2', '4'='gold'),
                    guide = guide_legend(title = NULL, byrow = TRUE)) +
  ggtitle("Baseclust result (K=4)") +
  plain
eurasiagdp.b4
dev.off()

### Figure 4-l2 (L2-dist)
pdf(file = "figure4-l2.pdf", width = 5.5, height = 3.9)
eurasiagdp.l3 <- ggplot(data = eurasiasubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = l3)) +
  scale_fill_manual(values = c('1'='gold', '2'='orange', '3'='khaki3'),
                    guide = guide_legend(title = NULL, byrow = TRUE)) +
  ggtitle("L2-dist result (K=3)") +
  plain
eurasiagdp.l3
dev.off()


### 0.25th quantile curves of each group (RSFclust2)
c_result = r2.4.c
years <- seq(from = 2001, to = 2022, length.out = 100)
c4_1 = rowMeans(vec.q0.25[, c_result == 1])
c4_2 = rowMeans(vec.q0.25[, c_result == 2])
c4_3 = rowMeans(vec.q0.25[, c_result == 3])
c4_4 = rowMeans(vec.q0.25[, c_result == 4])
dataq0.25 = data.frame(years, c4_1, c4_2, c4_3, c4_4)

pdf("figure6-q0.25.pdf", width = 4.5, height = 4.5)
ggplot(data = dataq0.25) +
  geom_line(aes(x = years, y = c4_1), color = 'indianred2') +
  geom_line(aes(x = years, y = c4_2), color = 'darkorange') +
  geom_line(aes(x = years, y = c4_3), color = 'khaki3') +
  geom_line(aes(x = years, y = c4_4), color = 'gold') +
  labs(title = '(RSFclust2) 0.25th quantile curves', x = 'years', y = 'Real GDP growth rate (%)') +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


### 0.75th quantile curves of each group (RSFclust2)
c_result = r2.4.c
years <- seq(from = 2001, to = 2022, length.out = 100)
c4_1 = rowMeans(vec.q0.75[, c_result == 1])
c4_2 = rowMeans(vec.q0.75[, c_result == 2])
c4_3 = rowMeans(vec.q0.75[, c_result == 3])
c4_4 = rowMeans(vec.q0.75[, c_result == 4])
dataq0.75 = data.frame(years, c4_1, c4_2, c4_3, c4_4)

pdf("figure6-q0.75.pdf", width = 4.5, height = 4.5)
ggplot(data = dataq0.75) +
  geom_line(aes(x = years, y = c4_1), color = 'indianred2') +
  geom_line(aes(x = years, y = c4_2), color = 'darkorange') +
  geom_line(aes(x = years, y = c4_3), color = 'khaki3') +
  geom_line(aes(x = years, y = c4_4), color = 'gold') +
  labs(title = '(RSFclust2) 0.75th quantile curves', x = 'years', y = 'Real GDP growth rate (%)') +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



### Mean curves of each group (if K = 4)
colname = 'b4'

c4matrix = matrix(nrow = 4, ncol = 22)
c4matrix[1,] = as.vector(as.matrix(colMeans(g.data.valid[g.data.valid[, colname]== '1', 2:23], na.rm = TRUE)))
c4matrix[2,] = as.vector(as.matrix(colMeans(g.data.valid[g.data.valid[, colname]== '2', 2:23], na.rm = TRUE)))
c4matrix[3,] = as.vector(as.matrix(colMeans(g.data.valid[g.data.valid[, colname]== '3', 2:23], na.rm = TRUE)))
c4matrix[4,] = as.vector(as.matrix(colMeans(g.data.valid[g.data.valid[, colname]== '4', 2:23], na.rm = TRUE)))

years <- 2001:2022
c4_1 = c4matrix[1,]
c4_2 = c4matrix[2,]
c4_3 = c4matrix[3,]
c4_4 = c4matrix[4,]
meandata4 <- data.frame(years, c4_1, c4_2, c4_3, c4_4)

# '1'='red2', '2'='orange', '3'='khaki4','4'='yellow2'

### Figure 6-rsf2
pdf("figure6-rsf2.pdf", width = 4.5, height = 4.5)
ggplot(data = meandata4) +
  geom_line(aes(x = years, y = c4_1), color = 'indianred2') +
  geom_line(aes(x = years, y = c4_2), color = 'darkorange') +
  geom_line(aes(x = years, y = c4_3), color = 'khaki3') +
  geom_line(aes(x = years, y = c4_4), color = 'gold') +
  labs(title = '(RSFclust2) mean curves', x = 'years', y = 'Real GDP growth rate (%)') +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### Figure 6-fpcac
pdf("figure6-fpcac.pdf", width = 4.5, height = 4.5)
ggplot(data = meandata4) +
  geom_line(aes(x = years, y = c4_1), color = 'indianred2') +
  geom_line(aes(x = years, y = c4_2), color = 'darkorange') +
  geom_line(aes(x = years, y = c4_3), color = 'khaki3') +
  geom_line(aes(x = years, y = c4_4), color = 'gold') +
  labs(title = '(FPCAC) mean curves', x = 'years', y = 'Real GDP growth rate (%)') +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### Figure 6-base
pdf("figure6-base.pdf", width = 4.5, height = 4.5)
ggplot(data = meandata4) +
  geom_line(aes(x = years, y = c4_1), color = 'khaki3') +
  geom_line(aes(x = years, y = c4_2), color = 'darkorange') +
  geom_line(aes(x = years, y = c4_3), color = 'indianred2') +
  geom_line(aes(x = years, y = c4_4), color = 'gold') +
  labs(title = '(Baseclust) mean curves', x = 'years', y = 'Real GDP growth rate (%)') +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


### Mean curves of each group (if K = 3)
colname = 'l3'

c3matrix = matrix(nrow = 3, ncol = 22)
c3matrix[1,] = as.vector(as.matrix(colMeans(g.data.valid[g.data.valid[, colname]== '1', 2:23], na.rm = TRUE)))
c3matrix[2,] = as.vector(as.matrix(colMeans(g.data.valid[g.data.valid[, colname]== '2', 2:23], na.rm = TRUE)))
c3matrix[3,] = as.vector(as.matrix(colMeans(g.data.valid[g.data.valid[, colname]== '3', 2:23], na.rm = TRUE)))

years <- 2001:2022
c3_1 = c3matrix[1,]
c3_2 = c3matrix[2,]
c3_3 = c3matrix[3,]
meandata3 <- data.frame(years, c3_1, c3_2, c3_3)

### Figure 6-rsf1
pdf("figure6-rsf1.pdf", width = 4.5, height = 4.5)
ggplot(data = meandata3) +
  geom_line(aes(x = years, y = c3_1), color = 'gold') +
  geom_line(aes(x = years, y = c3_2), color = 'darkorange') +
  geom_line(aes(x = years, y = c3_3), color = 'khaki3') +
  labs(title = '(RSFclust1) mean curves', x = 'years', y = 'Real GDP growth rate (%)') +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### Figure 6-exf
pdf("figure6-exf.pdf", width = 4.5, height = 4.5)
ggplot(data = meandata3) +
  geom_line(aes(x = years, y = c3_1), color = 'gold') +
  geom_line(aes(x = years, y = c3_2), color = 'darkorange') +
  geom_line(aes(x = years, y = c3_3), color = 'indianred2') +
  labs(title = '(Exfunclust) mean curves', x = 'years', y = 'Real GDP growth rate (%)') +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### Figure 6-l2
pdf("figure6-l2.pdf", width = 4.5, height = 4.5)
ggplot(data = meandata3) +
  geom_line(aes(x = years, y = c3_1), color = 'gold') +
  geom_line(aes(x = years, y = c3_2), color = 'darkorange') +
  geom_line(aes(x = years, y = c3_3), color = 'khaki3') +
  labs(title = '(L2-dist) mean curves', x = 'years', y = 'Real GDP growth rate (%)') +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



### Analysis of comparing FPCAC and RSFclust2

# FPCAC (1 : red, 2 : orange)
fpcac.orange = g.data.valid$region[fpc.4.c == 2] # 19 countries
fpcac.red = g.data.valid$region[fpc.4.c == 1] # 5 countries

# RSFclust2 (2 : orange, 1 : red)
r2.orange = g.data.valid$region[r2.4.c == 2] # 16 countries
r2.red = g.data.valid$region[r2.4.c == 1] # 9 countries

intersect.countries = intersect(union(fpcac.orange, fpcac.red), union(r2.orange, r2.red)) # 21 countries

orange.countries = intersect(fpcac.orange, r2.orange) # 13 countries
red.countries = intersect(fpcac.red, r2.red) # 4 countries
diff.countries = setdiff(intersect.countries, union(orange.countries, red.countries)) # 4 countries

g.data.valid.fpc = g.data.valid
g.data.valid.fpc$FPC1 = score.matrix3[, 5]
g.data.valid.fpc$FPC2 = score.matrix3[, 6]
x = rep(0, nrow(g.data.valid))
x[g.data.valid$region %in% orange.countries] = 1
x[g.data.valid$region %in% red.countries] = 2
g.data.valid.fpc$cluster = as.factor(x)
sub.data.fpc = g.data.valid.fpc[g.data.valid$region %in% intersect.countries, ]

iraq.fpc1 = sub.data.fpc[sub.data.fpc$region == 'Iraq', 'FPC1']
iraq.fpc2 = sub.data.fpc[sub.data.fpc$region == 'Iraq', 'FPC2']

pdf(file = "figure7-supp.pdf", width = 5.5, height = 5.5)
ggplot(data = sub.data.fpc, mapping = aes(x = FPC1, y = FPC2, colour = cluster)) +
  geom_point(size = 2, shape = 16) +
  geom_text(aes(label = ifelse(FPC1 == iraq.fpc1 & FPC2 == iraq.fpc2, region, "")), vjust = -1, color = "black") +
  labs(title = "FPC scores of 0.75 quantile curves", x = "FPC 1 score", y = "FPC 2 score", colour = "group") +
  scale_color_manual(values = c('0' = 'black', '1'='orange', '2'='red')) +
  theme_test() +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-1, 9)) +
  scale_y_continuous(limits = c(-5, 5))
dev.off()


### Scatterplot of FPC1, 2 scores
rrr <- as.factor(r2.4.c)
qscore.0.3 = score.matrix3[,1]
qscore.0.32 = score.matrix3[,2]
qscore.0.7 = score.matrix3[,5]
qscore.0.72 = score.matrix3[,6]
dataset = data.frame(qscore.0.3, qscore.0.32, qscore.0.7, qscore.0.72, rrr)

### Fiure 7-1
pdf("figure7-1.pdf", width = 5.5, height = 5.5)
ggplot(data = dataset, mapping = aes(x = qscore.0.7, y = qscore.0.72, colour = rrr)) +
  geom_point(size = 2, shape = 16) +
  labs(title = "FPC scores of 0.75 quantile curves", x = "FPC 1 score", y = "FPC 2 score", colour = "group") +
  scale_color_manual(values = c('1'='red2', '2'='orange', '3'='khaki4','4'='yellow2')) +
  theme_test() +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-6, 7)) +
  scale_y_continuous(limits = c(-6, 7))
dev.off()

### Figure 7-2
pdf("figure7-2.pdf", width = 5.5, height = 5.5)
ggplot(data = dataset, mapping = aes(x = qscore.0.3, y = qscore.0.32, colour = rrr)) +
  geom_point(size = 2, shape = 16) +
  labs(title = "FPC scores of 0.25 quantile curves", x = "FPC 1 score", y = "FPC 2 score", colour = "group") +
  scale_color_manual(values = c('1'='red2', '2'='orange', '3'='khaki4','4'='yellow2')) +
  theme_test() +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-7, 6)) +
  scale_y_continuous(limits = c(-6, 7))
dev.off()

### Economic growth rate graph for some selected countries
selected.region <- c('South Korea', 'China', 'Iraq', 'Vietnam', 'France', 'India', 'Japan', 'Maldives', 'Turkey')

graph.list = list(rep(1, length(selected.region)))
ze <- rep(0, 22)

for(i in 1:length(selected.region)){
  co <- selected.region[i]
  co.vec <- gdp_matrix[g.data.valid$region == co, ]
  gdata <- data.frame(years, co.vec, ze)
  if(co != 'Iraq' & co != 'Maldives'){
    graph.list[[i]] <- ggplot(data = gdata) +
      geom_line(mapping = aes(x = years, y = co.vec)) +
      geom_line(mapping = aes(x = years, y = ze), linetype = "dotted", linewidth = 0.3) +
      scale_y_continuous(limits = c(-12, 15)) +
      labs(title = co, x = "", y = "") +
      theme_test() +
      theme(legend.position = "none") +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_blank()) +
      theme(plot.title = element_text(size = 10, hjust = 0.5))
  }
  else{
    graph.list[[i]] <- ggplot(data = gdata) +
      geom_line(mapping = aes(x = years, y = co.vec)) +
      geom_line(mapping = aes(x = years, y = ze), linetype = "dotted", linewidth = 0.3) +
      labs(title = co, x = "", y = "") +
      theme_test() +
      theme(legend.position = "none") +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_blank()) +
      theme(plot.title = element_text(size = 10, hjust = 0.5, colour = "red"))
  }
}

### Figure 3
pdf("figure3.pdf")
g.grid = arrangeGrob(graph.list[[1]], graph.list[[2]], graph.list[[3]], graph.list[[4]], graph.list[[5]],
                      graph.list[[6]], graph.list[[7]], graph.list[[8]], graph.list[[9]], ncol = 3)
ggsave("figure3.pdf", g.grid)





















