setwd(".") #we use the current folder of the script as working directory
options(stringsAsFactors = FALSE) # we set the input strings not to be considered 
set.seed(10) # we set a seed to be able to replicate our tests
options(repos = list(CRAN="http://cran.rstudio.com/")) # we set the URL where to download the packages

# install.packages("pacman", dependencies = TRUE)
library("pacman")

p_load("dlookr", "dplyr", "ggplot2",  "pastecs", "tableone", "umap", "textshape", "factoextra", "ggdendro", "fpc", "cluster", "ggdendro", "clusterSim", "parameters", "dbscan", "umap", "ggplot2")

source("./utils.r")

VERBOSE <- FALSE

fileName <- "../data/Takashi2019_diabetes_type1_dataset_preprocessed.csv"
patients_data <- read.csv(fileName, header = TRUE, sep =",")

if(VERBOSE) cat("read dataset from file ", fileName, "\n", sep="")

if(VERBOSE) patients_data %>% dim()
if(VERBOSE) patients_data %>% summary()
if(VERBOSE) patients_data %>% str()
if(VERBOSE) patients_data %>% colnames() %>% sort()

list_of_umap_results <- c()
list_of_dbscan_results <- c()

hyperparameters_and_metrics_names <- c("DBSCAN k", "UMAP neighbors", "UMAP min dist", "DBSCAN eps",  "DBSCAN min points", "Silhouette", "Calinski-Harabasz", "recipr Davies-Bouldin", "Dunn", "Gap")

NUM_VALUES <- hyperparameters_and_metrics_names %>% length()
resultsDataFrame <- as.data.frame(matrix(ncol=NUM_VALUES, nrow=1))
colnames(resultsDataFrame) <- hyperparameters_and_metrics_names

if(VERBOSE) cat(" ~ ~ ~ UMAP ~ ~ ~\n")

min_neighbors <- 2
max_neighbors <- 4
number_of_neighbors <- c(min_neighbors:max_neighbors)

min_min_distance <- 0.0001
max_min_distance <- 0.002
step <- 0.0001
min_distances <- seq(from=min_min_distance, to=max_min_distance, by=step)

# min_dist 0.02 and n_neighbors 2 produces Silhouette 0.97

cat(" :: :: DBSCAN applied to the UMAP reduced dataset of ", fileName, "  :: :: \n", sep="")

Sys.time() %>% print()

cat("UMAP miminal distances to try:")
min_distances %>% print()

cat("UMAP number of n_neighbors to try':")
number_of_neighbors %>% print()

exe_i <- 1

for(this_neigh_num in number_of_neighbors) {
    for(this_min_dist in min_distances) {

        custom.settings <- umap.defaults
        custom.settings$"n_neighbors" <- this_neigh_num
        custom.settings$"min_dist" <- this_min_dist

        patients_data_umap <- umap(patients_data, config=custom.settings)

        # display object summary
        if(VERBOSE) patients_data_umap %>% print()

        # display embedding coordinates
        if(VERBOSE) head(patients_data_umap$"layout")

        if(VERBOSE) cat("\n ~ ~ ~ DBSCAN ~ ~ ~\n")

        method <- "DBSCAN"

        this_eps <- (1-custom.settings$"min_dist")
        these_minimal_points <- custom.settings$"n_neighbors"
        
        dbscan_result <- dbscan::dbscan(patients_data_umap$"layout", eps = this_eps, minPts = these_minimal_points)
        
        if(VERBOSE) dbscan_result %>% print()
        if(VERBOSE) dbscan_result %>% str()

        if(VERBOSE) table(dbscan_result$"cluster") %>% print()
        
        noise_points <- dbscan_result$cluster[which(dbscan_result$cluster==0)] %>% length()
        cat("number of noise points in the 0 cluster = ", noise_points, "\n", sep="")

        perc_noise_cluster <- noise_points * 100 / (patients_data %>% nrow())
        cat("percentage of data in the 0 noise cluster = ", perc_noise_cluster %>% dec_two(), "%\n", sep="")
        
        noise_threshold <- 0
        noise_threshold_exceeded <- FALSE
        if(perc_noise_cluster > noise_threshold) { 
            noise_threshold_exceeded <- TRUE
            cat("The 0 noise cluster contains too much data: we'll skip this analysis configuration\n")
        }

        
        # an integer vector of length of the number of cases, which indicates a clustering. The clusters have to be numbered from 1 to the number of clusters.

        one_cluster_stop <- FALSE
        if(max(dbscan_result$"cluster") == 1) {
                one_cluster_stop <- TRUE 
                cat("There's only 1 cluster: we'll skip this analysis configuration\n")
            }
        
        
        if(max(dbscan_result$"cluster") >= 2 & noise_threshold_exceeded  == FALSE) {
        
        clusters_labels <- NULL

        if(min(dbscan_result$"cluster") == 0) {
            clusters_labels <- dbscan_result$"cluster" + 1
        } else { 
            clusters_labels <- dbscan_result$"cluster"
        }


        cluster_statistics <- fpc::cluster.stats(dist(as.data.frame(patients_data_umap$"layout")), clusters_labels, silhouette = TRUE) 


        cat(" ~ :  ~ :  ~ :  ~ : ~ DBSCAN results on the UMAP projection of the EHRs data  ~ :  ~ :  ~ :  ~ : ~\n\n")


        cat("UMAP hyperparameters: number of neighbors = ", custom.settings$n_neighbors , " &  min_dist = ", custom.settings$min_dist , "\n", sep="")
        cat("DBSCAN hyperparameters: epsilon = ", this_eps , " & minimal points = ", these_minimal_points, "\n\n", sep="")


        cat("average Silhouette score = ", cluster_statistics$"avg.silwidth" %>% dec_three(), " in the [-1,+1] interval \n", sep="")

        cat("Calinski-Harabasz index = ", cluster_statistics$"ch" %>% dec_three(), " in the [0,+∞) interval\n", sep="")

        cat("Dunn index = ", cluster_statistics$"dunn" %>% dec_three(), " in the [0,+∞) interval\n", sep="")

        davies_bouldin_index_results <- index.DB(x=as.data.frame(patients_data_umap$"layout"), cl=clusters_labels) 
        davies_bouldin_index <- davies_bouldin_index_results$"DB"

        cat("reciprocal Davies-Bouldin index = ",  (1/davies_bouldin_index)  %>% dec_three(), " in the [0,+∞) interval\n", sep="")

        bootstrap_parameter <- 60
        this_gap <- clusGap(x = patients_data_umap$layout, FUNcluster = dbscan, K.max = max(dbscan_result$cluster), B = bootstrap_parameter, eps = custom.settings$min_dist, verbose = FALSE)
        average_gap <- as.data.frame(this_gap$"Tab")$"gap" %>% mean()

        cat("average Gap statistic = ", average_gap %>% dec_three(), " in the [0,+∞) interval\n", sep="")
        cat("(all these metrics: the higher, the better)\n")
        
        # c("UMAP neighbors", "UMAP min dist", "DBSCAN eps",  "DBSCAN min points", 
        # "Silhouette", "Calinski-Harabasz", "Davies-Bouldin", "Dunn", "Gap")

        cat("DBSCAN number of clusters = ", dbscan_result$"cluster" %>% max(), "\n", sep="")
        
        outputDataframe <- as.data.frame(matrix(ncol=NUM_VALUES, nrow=1))
        colnames(outputDataframe) <- hyperparameters_and_metrics_names
        outputDataframe[,1] <- dbscan_result$"cluster" %>% max()
        outputDataframe[,2] <- custom.settings$"n_neighbors"
        outputDataframe[,3] <- custom.settings$"min_dist"
        outputDataframe[,4] <- this_eps
        outputDataframe[,5] <- these_minimal_points
        
        outputDataframe[,6] <- cluster_statistics$"avg.silwidth"
        outputDataframe[,7] <- cluster_statistics$"ch"
        outputDataframe[,8] <- (1/davies_bouldin_index)
        outputDataframe[,9] <- cluster_statistics$"dunn"
        outputDataframe[,10] <- average_gap

        if (exe_i == 1)  resultsDataFrame <-  outputDataframe
        else  resultsDataFrame <- rbind(resultsDataFrame, outputDataframe)
        
        list_of_umap_results[exe_i] <- patients_data_umap
        list_of_dbscan_results[exe_i] <- dbscan_result


        exe_i <- exe_i + 1
        cat(" ~ :  ~ :  ~ :  ~ : ~ The end  ~ :  ~ :  ~ :  ~ : ~\n")

        }

        else cat("Only 1 cluster found by DBSCAN or too much noise: we do not compute the clustering statistics in this case\n")
    }
}


resultsDataFrame %>% print()

cat("\n : : results with maximum Silhouette score: \n")

best_outcome_by_Silhouette <- resultsDataFrame[resultsDataFrame$"Silhouette" ==  max(resultsDataFrame$"Silhouette"),]
best_outcome_by_Silhouette[,c("DBSCAN k", "Silhouette")] %>% print()

cat("\n : : results with maximum Calinski-Harabasz index: \n")
resultsDataFrame[resultsDataFrame$"Calinski-Harabasz" ==  max(resultsDataFrame$"Calinski-Harabasz"), c("DBSCAN k", "Calinski-Harabasz")] %>% print()

cat("\n : : results with maximum reciprocal Davies-Bouldin index: \n")
resultsDataFrame[resultsDataFrame$"recipr Davies-Bouldin" ==  max(resultsDataFrame$"recipr Davies-Bouldin"), c("DBSCAN k", "recipr Davies-Bouldin")] %>% print()

cat("\n : : results with maximum Dunn index: \n")
resultsDataFrame[resultsDataFrame$"Dunn" ==  max(resultsDataFrame$"Dunn"), c("DBSCAN k", "Dunn")] %>% print()

cat("\n : : results with maximum Gap statistic: \n")
resultsDataFrame[resultsDataFrame$"Gap" ==  max(resultsDataFrame$"Gap"), c("DBSCAN k", "Gap")] %>% print()

SAVE_PLOT_BEST_OUTCOME <- TRUE
if(SAVE_PLOT_BEST_OUTCOME) {

    cat("\nWe save the UMAP and DBSCAN plot with the best hyperparameters according to the Silhouette\n")

    index_best_result_by_Silhouette <- best_outcome_by_Silhouette %>% rownames() %>% as.numeric()

    best_outcome_by_Silhouette %>% print()

    patients_data_umap_best <- list_of_umap_results[[index_best_result_by_Silhouette]]
    best_dbscan_result <- list_of_dbscan_results[[index_best_result_by_Silhouette]]

    thisK <- best_dbscan_result %>% max()

    data_for_ggplot2 <- as.data.frame(cbind(as.data.frame(patients_data_umap_best), as.data.frame(best_dbscan_result)))
    colnames(data_for_ggplot2) <- c("UMAP1", "UMAP2", "cluster")
    data_for_ggplot2$cluster <- as.factor(data_for_ggplot2$"cluster")
    this_plot <- ggplot(data_for_ggplot2, aes(x = UMAP1, y = UMAP2)) + geom_jitter(aes(color = `cluster`, size = 3), height = 0.9, width = 0.9) + theme(text = element_text(size = 20))

    pdfFileName <- paste0("../results/umap_dbscan_", toString(thisK), "clusters_", this_moment_formatted(), ".pdf")
    ggsave(this_plot, file=pdfFileName, width=20, height=10)
    cat("Saved file ", pdfFileName, "\n", sep="")

}



cat(" // The end of the script // \n")

     
