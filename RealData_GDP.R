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


### Run RSFclust2, Baseclust, Funclust
num.knots = 10
n_breaks = 22
b.sp.w = bs(seq(0, 1, by = 1/n_breaks), knots = seq(0, 1, length.out = num.knots+1))
b.basis.w = bsplineS(seq(from=0, to=1, length.out = n_breaks+1)[-(n_breaks+1)], seq(from=0, to=1, by=1/num.knots))

tw = seq(from = 0, to = 1, length.out = n_breaks+1)[-(1+n_breaks)]
gdp_matrix = as.matrix(g.data.valid[,-1])

bmatrix = matrix(nrow = nrow(gdp_matrix), ncol = num.knots+3)
vec.original = matrix(nrow = length(tw), ncol = nrow(gdp_matrix))
gdp.list = list()
for(i in 1:nrow(gdp_matrix)) gdp.list = c(gdp.list, list(as.vector(gdp_matrix[i,])))
score.matrix2 = FPC.score.matrix(gdp.list, c(0.25, 0.5, 0.75), c(0.5), K = 2, num.knots, n_breaks)
score.matrix3 = FPC.score.matrix(gdp.list, c(0.25, 0.5, 0.75), c(0.5), K = 3, num.knots, n_breaks)
score.matrix4 = FPC.score.matrix(gdp.list, c(0.25, 0.5, 0.75), c(0.5), K = 4, num.knots, n_breaks)

for(j in 1:nrow(gdp_matrix)){
  y = gdp_matrix[j,]
  vec.original[,j] = gdp_matrix[j,]
  
  bmat = as.matrix(nearPD(t(b.basis.w)%*%b.basis.w)$mat)
  bmatrix[j,] = solve(bmat, t(b.basis.w)%*%y)
}


# Funclust (K = 2, 3, 4)
bsp = create.bspline.basis(breaks=seq(0, 1, by=1/num.knots))
fd.obj <- fd(coef=coef.mean, basisobj=bsp)
f2.c <- funclust(fd.obj, K=2)$cls
f3.c <- funclust(fd.obj, K=3)$cls
f4.c <- funclust(fd.obj, K=4)$cls

# RSFclust2 (K = 2, 3, 4)
R2g <- RSKC(score.matrix2, 2, alpha = 0.1, L1 = 2)
R3g <- RSKC.multidim(score.matrix3, ncl = 3, l1bound = 2)
R4g <- RSKC.multidim(score.matrix4, ncl = 4, l1bound = 2)

r2.cg <- R2g$labels
r2.wg <- R2g$weights
r3.cg <- R3g$labels
r3.wg <- R3g$weights
r4.cg <- R4g$labels
r4.wg <- R4g$weights


# Baseclust (K = 2, 3, 4)
b2.c = kmeans(bmatrix, centers = 2, nstart = 10)$cluster
b3.c = kmeans(bmatrix, centers = 3, nstart = 10)$cluster
b4.c = kmeans(bmatrix, centers = 4, nstart = 10)$cluster


### Choosing the number of the clusters

# Silhouette coefficient (The bigger, the better, K=3)

mean(silhouette(r2.cg, dist(sweep(qsparse, 2, sqrt(r2.wg), "*")))[,3]) # 0.41
mean(silhouette(r3.cg, dist(sweep(qsparse2, 2, sqrt(r3.wg), "*")))[,3]) # 0.36
mean(silhouette(r4.cg, dist(sweep(qsparse3, 2, sqrt(r4.wg), "*")))[,3]) # 0.38

mean(silhouette(f2.c, dist(t(vec.mean)))[,3]) # 0.02
mean(silhouette(f3.c, dist(t(vec.mean)))[,3]) # -
mean(silhouette(f4.c, dist(t(vec.mean)))[,3]) # -

mean(silhouette(b2.c, dist(bmatrix))[,3]) # 0.87
mean(silhouette(b3.c, dist(bmatrix))[,3]) # 0.43
mean(silhouette(b4.c, dist(bmatrix))[,3]) # 0.50


### Adding the cluster assignment to g.data
x = rep(0, 87)
x[valid.index] = r2.cg
x<- as.character(x)
for(i in 1:87){
  x[i] <- paste(x[i], " ")
}
g.data[,"r2"] = x

x = rep(0, 87)
x[valid.index] = r3.cg
x<- as.character(x)
for(i in 1:87){
  x[i] <- paste(x[i], " ")
}
g.data[,"r3"] = x

x = rep(0, 87)
x[valid.index] = r4.cg
x<- as.character(x)
for(i in 1:87){
  x[i] <- paste(x[i], " ")
}
g.data[,"r4"] = x


x = rep(0, 87)
x[valid.index] = f2.c
x<- as.character(x)
for(i in 1:87){
  x[i] <- paste(x[i], " ")
}
g.data[,"f2"] = x

x = rep(0, 87)
x[valid.index] = b4.c
x<- as.character(x)
for(i in 1:87){
  x[i] <- paste(x[i], " ")
}
g.data[,"b4"] = x

eurasiasubset <- inner_join(world, g.data, by = "region")
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
  legend.text = element_text(size = 10)
)

### Figure 4-1 (Funclust)
eurasiagdp2 <- ggplot(data = eurasiasubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = f2)) +
  scale_fill_manual(values = c('0  '='seashell2', '2  '='orange', '1  '='gold'),
                    name = "") +
  ggtitle("Funclust result (K=2)") +
  plain
eurasiagdp2

### Figure 4-2 (Baseclust)
eurasiagdp_b4 <- ggplot(data = eurasiasubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = b4)) +
  scale_fill_manual(values = c('0  '='seashell2', '1  '='khaki3', '2  '='orange', '3  '='indianred2', '4  '='gold'),
                    name = "") +
  ggtitle("Baseclust result (K=4)") +
  plain
eurasiagdp_b4

### Figure 5 (RSFclust2)
eurasiagdp4 <- ggplot(data = eurasiasubset, mapping = aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3) +
  geom_polygon(aes(fill = r4)) + 
  scale_fill_manual(values = c('0  '='seashell2', '1  '='indianred2', '2  '='khaki3', '3  '='orange', '4  '='gold'),
                    name = "") +
  ggtitle("RSFclust2 result (K=4)") +
  plain
eurasiagdp4



### Mean curves of Baseclust
c4matrix = matrix(nrow = 4, ncol = 22)

c4matrix[1,] = as.vector(as.matrix(colMeans(g.data[g.data$b4== '1  ', 2:23], na.rm = TRUE)))
c4matrix[2,] = as.vector(as.matrix(colMeans(g.data[g.data$b4== '2  ', 2:23], na.rm = TRUE)))
c4matrix[3,] = as.vector(as.matrix(colMeans(g.data[g.data$b4== '3  ', 2:23], na.rm = TRUE)))
c4matrix[4,] = as.vector(as.matrix(colMeans(g.data[g.data$b4== '4  ', 2:23], na.rm = TRUE)))

years <- 2001:2022
c4_1 = c4matrix[1,]
c4_2 = c4matrix[2,]
c4_3 = c4matrix[3,]
c4_4 = c4matrix[4,]

meandata4 <- data.frame(years, c4_1, c4_2, c4_3, c4_4)

### Figure 6-1
ggplot(data = meandata4) +
  geom_line(aes(x = years, y = c4_1), color = 'khaki3') +
  geom_line(aes(x = years, y = c4_2), color = 'orange') +
  geom_line(aes(x = years, y = c4_3), color = 'indianred2') +
  geom_line(aes(x = years, y = c4_4), color = 'gold') +
  labs(title = '(Baseclust) mean curves', x = 'years', y = 'Real GDP growth rate') +
  theme_test() +
  theme(plot.title = element_text(size = 12, hjust = 0.5))



### Mean curves of RSFclust2
c4matrix = matrix(nrow = 4, ncol = 22)

c4matrix[1,] = as.vector(as.matrix(colMeans(g.data[g.data$r4== '1  ', 2:23], na.rm = TRUE)))
c4matrix[2,] = as.vector(as.matrix(colMeans(g.data[g.data$r4== '2  ', 2:23], na.rm = TRUE)))
c4matrix[3,] = as.vector(as.matrix(colMeans(g.data[g.data$r4== '3  ', 2:23], na.rm = TRUE)))
c4matrix[4,] = as.vector(as.matrix(colMeans(g.data[g.data$r4== '4  ', 2:23], na.rm = TRUE)))

years <- 2001:2022
c4_1 = c4matrix[1,]
c4_2 = c4matrix[2,]
c4_3 = c4matrix[3,]
c4_4 = c4matrix[4,]

meandata4 <- data.frame(years, c4_1, c4_2, c4_3, c4_4)

### Figure 6-2
ggplot(data = meandata4) +
  geom_line(aes(x = years, y = c4_1), color = 'indianred2') +
  geom_line(aes(x = years, y = c4_2), color = 'khaki3') +
  geom_line(aes(x = years, y = c4_3), color = 'orange') +
  geom_line(aes(x = years, y = c4_4), color = 'gold') +
  labs(title = '(RSFclust2) mean curves', x = 'years', y = 'Real GDP growth rate') +
  theme_test() +
  theme(plot.title = element_text(size = 12, hjust = 0.5))





### Scatterplot of FPC1, 2 scores
rrr <- as.factor(r4.cg)
dataset = data.frame(qscore.0.3, qscore.0.32, qscore.0.7, qscore.0.72, rrr)

### Fiure 7-1
ggplot(data = dataset, mapping = aes(x = qscore.0.7, y = qscore.0.72, colour = rrr)) +
  geom_point(size = 2, shape = 16) +
  labs(title = "FPC scores of 0.75th quantile curves", x = "FPC 1 score", y = "FPC 2 score", colour = "group") +
  scale_color_manual(values = c('1'='red2', '2'='khaki4', '3'='orange','4'='yellow2')) +
  theme_test() +
  theme(legend.position = "right") +
  theme(legend.text = element_text(size = 10)) +
  scale_x_continuous(limits = c(-6, 7)) +
  scale_y_continuous(limits = c(-6, 7))

### Figure 7-2
ggplot(data = dataset, mapping = aes(x = qscore.0.3, y = qscore.0.32, colour = rrr)) +
  geom_point(size = 2, shape = 16) +
  labs(title = "FPC scores of 0.25th quantile curves", x = "FPC 1 score", y = "FPC 2 score", colour = "group") +
  scale_color_manual(values = c('1'='red2', '2'='khaki4', '3'='orange','4'='yellow2')) +
  theme_test() +
  theme(legend.position = "right") +
  theme(legend.text = element_text(size = 10)) +
  scale_x_continuous(limits = c(-7, 6)) +
  scale_y_continuous(limits = c(-6, 7))


### Economic growth rate graph for some selected countries
selected.region <- c('South Korea', 'China', 'Iraq', 'Vietnam', 'France', 'India', 'Japan', 'Maldives', 'Turkey')

graph.list = list(rep(1, 9))
ze <- rep(0, 22)

for(i in 1:9){
  co <- selected.region[i]
  co.vec <- gdp_matrix[g.data.valid$region == co, ]
  gdata <- data.frame(years, co.vec, ze)
  if(co != 'Iraq' & co != 'Maldives'){
    graph.list[[i]] <- ggplot(data = gdata) +
      geom_line(mapping = aes(x = years, y = co.vec)) +
      geom_line(mapping = aes(x = years, y = ze), linetype = "dotted", size = 0.3) +
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
      geom_line(mapping = aes(x = years, y = ze), linetype = "dotted", size = 0.3) +
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





















