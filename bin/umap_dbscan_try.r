setwd(".") #we use the current folder of the script as working directory
options(stringsAsFactors = FALSE) # we set the input strings not to be considered 
set.seed(10) # we set a seed to be able to replicate our tests
options(repos = list(CRAN="http://cran.rstudio.com/")) # we set the URL where to download the packages

# install.packages("pacman", dependencies = TRUE)
library("pacman")

p_load("dlookr", "dplyr", "ggplot2",  "pastecs", "tableone", "umap", "textshape", "factoextra", "ggdendro", "fpc", "cluster", "ggdendro", "clusterSim", "parameters", "dbscan", "umap")

cat(" ~ ~ ~ UMAP ~ ~ ~\n")


custom.settings <- umap.defaults
custom.settings$"n_neighbors" <- 8
custom.settings$"min_dist" <- 0.1

iris.umap <- umap(iris[,1:4], config=custom.settings)

# display object summary
iris.umap

# display embedding coordinates
head(iris.umap$layout)

cat("\n ~ ~ ~ DBSCAN ~ ~ ~\n")

n <- 100
x <- cbind(
x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
y = runif(10, 0, 10) + rnorm(n, sd = 0.2)
)

head(x)
dim(x)

res <- dbscan::dbscan(iris.umap$layout, eps = (1-custom.settings$min_dist), minPts = custom.settings$n_neighbors)
res %>% print()
res %>% str()  

table(res$cluster) 

cat(" ~ ~ ~ clustering statistics ~ ~ ~\n")

stats <- fpc::cluster.stats(dist(as.data.frame(iris.umap$layout)), res$cluster, silhouette = TRUE) 
stats 
     
