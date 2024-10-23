setwd(".") #we use the current folder of the script as working directory
options(stringsAsFactors = FALSE) # we set the input strings not to be considered 
set.seed(10) # we set a seed to be able to replicate our tests
options(repos = list(CRAN="http://cran.rstudio.com/")) # we set the URL where to download the packages

# install.packages("pacman", dependencies = TRUE)
library("pacman")

p_load("dlookr", "dplyr", "ggplot2",  "pastecs", "tableone", "umap", "textshape", "factoextra", "ggdendro", "fpc", "cluster", "ggdendro", "clusterSim", "parameters", "dbscan", "umap", "ggplot2", "ggtext")

source("./utils.r")


# data reading

# fileName <- "/home/davide/my_projects/UMAP_assessment_metric/data/journal.pone.0175818_S1Dataset_Spain_cardiac_arrest_EDITED.csv"
fileName <- "/home/davide/my_projects/UMAP_assessment_metric/data/journal.pone.0148699_S1_Text_Sepsis_SIRS_EDITED.csv"

patients_data_original <- read.csv(fileName, header = TRUE, sep =",")
patients_data <- patients_data_original[sample(1:nrow(patients_data_original)), ] %>% na.omit()  # shuffled rows
cat("\nread file ", fileName, "\n", sep="")


# setting DBSCAN hyperparameters

this_eps <- 0.17
these_minimal_points <- 10

# apply DBSCAN to the original dataset

dbscan_initial_dataset_results <- dbscan::dbscan(patients_data, eps = this_eps, minPts = these_minimal_points)

dbscan_initial_dataset_results %>% print()

# calculate Silhouette, Davies-Bouldin, Calinski-Harabasz, Gap statistic, Dunn index (pre-UMAP)

       one_cluster_stop <- FALSE
        if(max(dbscan_initial_dataset_results$"cluster") == 1) {
                one_cluster_stop <- TRUE
                cat("There's only 1 cluster: we'll skip this analysis configuration\n")
            }


        if(max(dbscan_initial_dataset_results$"cluster") >= 2 & one_cluster_stop  == FALSE) {

        clusters_labels <- NULL
        clusters_labels <- dbscan_initial_dataset_results$"cluster"

        cluster_statistics <- fpc::cluster.stats(dist(as.data.frame(patients_data)), clusters_labels, silhouette = TRUE)

        cat("\n ~ :  ~ :  ~ :  ~ : ~ DBSCAN results on the original  data  ~ :  ~ :  ~ :  ~ : ~\n\n")


        cat("DBSCAN hyperparameters: epsilon = ", this_eps , " & minimal points = ", these_minimal_points, "\n\n", sep="")
        cat("average Silhouette score = ", cluster_statistics$"avg.silwidth" %>% dec_three(), " in the [-1,+1] interval \n", sep="")
        cat("Calinski-Harabasz index = ", cluster_statistics$"ch" %>% dec_three(), " in the [0,+∞) interval\n", sep="")
        cat("Dunn index = ", cluster_statistics$"dunn" %>% dec_three(), " in the [0,+∞) interval\n", sep="")

        davies_bouldin_index_results <- index.DB(x=as.data.frame(patients_data), cl=clusters_labels)
        davies_bouldin_index <- davies_bouldin_index_results$"DB"

        cat("reciprocal Davies-Bouldin index = ",  (1/davies_bouldin_index)  %>% dec_three(), " in the [0,+∞) interval\n", sep="")

        bootstrap_parameter <- 60
        this_gap <- clusGap(x = patients_data, FUNcluster = dbscan, K.max = max(clusters_labels), B = bootstrap_parameter, eps = this_eps, verbose = FALSE)
        average_gap <- as.data.frame(this_gap$"Tab")$"gap" %>% mean()

        cat("average Gap statistic = ", average_gap %>% dec_three(), " in the [0,+∞) interval\n", sep="")
        cat("(all these metrics: the higher, the better)\n")



        cat("\nDBSCAN number of clusters = ", dbscan_initial_dataset_results$"cluster" %>% max(), " on ", patients_data %>% nrow() ," patients \n", sep="")

        clusters_over_patients_ratio <- dbscan_initial_dataset_results$"cluster" %>% max()  / patients_data %>% nrow()
        cat("clusters over patients ratio = ", clusters_over_patients_ratio %>% dec_three(), "\n", sep="")

}

# apply UMAP to the original dataset

# apply DBSCAN to UMAP's results

# calculate Silhouette, Davies-Bouldin, Calinski-Harabasz, Gap statistic, Dunn index (post-UMAP)

# compare the results of pre-UMAP and post-UMAP
