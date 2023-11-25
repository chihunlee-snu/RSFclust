# RSFclust

Requirements.R includes basic functions to run simulations and real data analysis, including the codes to obtain quantile and expectile curves using B-spline fits. 
RSFclust2.R includes codes to run RSFclust2 algorithm, which incorporates a modification of RSK-means(Robust and Sparse K-means) algorithm.

Simulation - 2 clsuters.R includes the codes to run simulations with 2 clusters. It requires R functions in Requirements.R and RSFclust2.R.

Simulation - 3 clsuters.R includes the codes to run simulations with 3 clusters. It requires R functions in Requirements.R and RSFclust2.R.

Real data - Economic growth rate.R incudes the codes to apply RSFclust1, RSFclust2 to the economic growth rate data. It requires dataset realgdprate.csv, and R functions in Requirements.R and RSFclust2.R.

Real data - Canadian weather.R includes the codes to apply RSFclust1, RSFclust2 to the Canadian weather data. It requires R functions in Requirements.R and RSFclust2.R.

The results and figures of the simulation code and real data analysis are displayed in the paper "The use of quantile curves to extend the distribution scope of functional data in functional data clustering" by Chihun Lee, Hee-Seok Oh and Joonpyo Kim.
