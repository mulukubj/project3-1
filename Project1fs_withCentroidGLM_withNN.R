# important -- this makes sure our runs are consistent and reproducible
set.seed(0)

library(neuralnet)

file <- "TCGA_breast_cancer_ERstatus_allGenes.txt"
nfold <- 5
sd_threshold <- 1
#sd_threshold <- 4
#top_num <- 10
#top_num <- 50
#top_num <- 100

header <- scan(file, nlines = 1, sep="\t", what = character())
data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names=FALSE)
names(data) <- header

header2 <- scan(file, skip = 1, nlines = 1, sep="\t", what = character())

# cleanup - remove genes with sd < 1
# compute sd for each gene
data_sd<-sapply(seq(nrow(data)), function(x) { as.numeric(sd(data[x,-1])) })

# add gene names to the sd list
data_sd_names<-cbind(data.frame(data_sd),data[,1])

# create an "include" list of all those genes where sd > threashold
include_list <- data_sd_names[data_sd_names[,1]>sd_threshold,2]

Positive <- data[data$id %in% include_list,header2=='Positive']
Negative <- data[data$id %in% include_list,header2=='Negative']

# define function cross_valid so we can rerun the cross validataion with various parameters
cross_validation <- function (nfold, top_num, alg) {
  
  Positive_groups <- split(sample(colnames(Positive)), 1+(seq_along(colnames(Positive)) %% nfold))
  Negative_groups <- split(sample(colnames(Negative)), 1+(seq_along(colnames(Negative)) %% nfold))
  
  result <- array()
  
  for (test_group in 1:nfold) {
    
    testA <- Positive[,colnames(Positive) %in% unlist(Positive_groups[test_group])]
    testB <- Negative[,colnames(Negative) %in% unlist(Negative_groups[test_group])]
    
    trainingA <- Positive[,!(colnames(Positive) %in% unlist(Positive_groups[test_group]))]
    trainingB <- Negative[,!(colnames(Negative) %in% unlist(Negative_groups[test_group]))]
    
    # Feature selection -- 
    
    # compute t-statistic for each row
    training_t_stat<-data.frame(sapply(seq(nrow(trainingA)), function(x) { abs(as.numeric(t.test(trainingA[x,], trainingB[x,])$statistic)) }))
    
    # add gene id column
    training_t_stat_geneid<-cbind(training_t_stat,rownames(trainingA))
    colnames(training_t_stat_geneid) <- c('t','id')
    
    # pick top 50 based on t-statistic
    selected_genes <- head(training_t_stat_geneid[order(-training_t_stat_geneid$t),],n=top_num)[,2]
    
    # narrow down the list of genes based on t-statistic
    testA <- testA[rownames(testA) %in% selected_genes,]
    testB <- testB[rownames(testB) %in% selected_genes,]
    trainingA <- trainingA[rownames(trainingA) %in% selected_genes,]
    trainingB <- trainingB[rownames(trainingB) %in% selected_genes,]
    
    if(alg == "centroid") {
      
      centroidA <- rowMeans(trainingA)
      centroidB <- rowMeans(trainingB)
      
      misclassifiedA <- sum(sapply(testA, function(x) { sqrt(sum((x-centroidA)^2))-sqrt(sum((x-centroidB)^2))>0 }))
      misclassifiedB <- sum(sapply(testB, function(x) { sqrt(sum((x-centroidA)^2))-sqrt(sum((x-centroidB)^2))<0 }))
    }
    
    if(alg == "glm") {
      trainingCombined <- rbind(cbind(data.frame(t(trainingA)),cancer=0),cbind(data.frame(t(trainingB)),cancer=1))
      testA0 <- data.frame(t(testA))
      testB0 <- data.frame(t(testB))
      
      model <- glm(cancer ~ ., data=trainingCombined, family=binomial, control = list(maxit=50))
      pA <- predict(model, newdata= testA0, type="response")
      pB <- predict(model, newdata= testB0, type="response")
      
      misclassifiedA <- sum(ifelse(pA<0.5,0,1))
      misclassifiedB <- sum(ifelse(pB>0.5,0,1))
    }
    
    if(alg == "nn") {
      trainingCombined <- rbind(cbind(data.frame(t(trainingA)),cancer=0),cbind(data.frame(t(trainingB)),cancer=1))
      testA0 <- data.frame(t(testA))
      testB0 <- data.frame(t(testB))
      
      ## 1-layer network with 2 neurons, startweights - random, linear.output = False for classification
      model <- neuralnet(cancer ~., data=trainingCombined, hidden=c(2), startweights = NULL, linear.output =F)
      pA <- predict(model, testA0, rep = 1, all.units = FALSE)
      pB <- predict(model, testB0, rep = 1, all.units = FALSE)
      
      misclassifiedA <- sum(ifelse(pA<0.5,0,1))
      misclassifiedB <- sum(ifelse(pB>0.5,0,1))
    }
    ## Plot the Neural Net model - Black lines = connnections between each layer and weights, Blue lines = bias
    plot(model)
    result[test_group] <- (misclassifiedA+misclassifiedB)/(ncol(testA)+ncol(testB))
  }
  
  paste0(round(mean(result),4)," sd=(",round(sd(result),4),")")
}

# centroid_50_all <- cross_validation(nfold=5, top_num = 50, alg= "centroid")
# centroid_100_all <- cross_validation(nfold=5, top_num = 100, alg= "centroid")

# glm_50_all <- cross_validation(nfold=5, top_num = 50, alg = "glm")
# glm_100_all <- cross_validation(nfold=5, top_num = 100, alg = "glm")

nn_50_all <- cross_validation(nfold=5, top_num = 50, alg = "nn")
#nn_100_all <- cross_validation(nfold=5, top_num = 100, alg = "nn")

# x<-data.frame("Centroid"=c(centroid_50_all),"GLM"=c(glm_50_all))
# rownames(x) <- c("50 genes")
# kable(x)  # Requires knitr package

#x<-data.frame("Centroid"=c(centroid_50_all,centroid_100_all),"GLM"=c(glm_50_all,glm_100_all),"NeuralNet"=c(nn_50_all,nn_100_all))
#rownames(x) <- c("50 genes","100 genes")


#x<-data.frame("Centroid"=c(centroid_50_all),"GLM"=c(glm_50_all), "NeuralNet"=c(nn_50_all))
#rownames(x) <- c("50 genes")
#kable(x)
