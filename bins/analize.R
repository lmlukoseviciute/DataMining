# This script is the analyses of mice protein expression data

# Not all libraries are needed I just don't know which does what so I 
# left them all
library(pacman)
p_load("dplyr", "ggplot2","foreach", "ggplot2", "datasets", "MASS", 
       "e1071", "class", "caret") 

# Loading data
df = read.csv("../inputs/Data_Cortex_Nuclear.csv")
df_copy = df
# Looking to the structure of the data 
summary (df)
# Figuring out how many NAs are there
sum(is.na(df))
# Applying NA filling alorithm
# Missing Value Imputation by Weighted Moving Average
library(imputeTS)
for (i in 1:length(colnames(df[2:78]))){
    df[i] = na_ma(df[i], k = 2, weighting = "exponential")
}
# just checking if everything was replaced
sum(is.na(df))

# I decided to work with 2 groups: control mice and trisomic mice 
# The variable for calsification is Ganotype

# Normalizing the data
library(BBmisc)
for (i in 2:78){
    df[i] = normalize(df[i], method = "range", range = c(0, 1),
                      margin = 1L, on.constant = "quiet")
}

summary(df)

# Preparing data for crossvalidation
# Sorry, I borrowed the code

# Let's try knn
library(class)
df_knn= df[,2:79]

cvKNN <- function(df_knn, kcv = 5, kNN = 2, label="Genotype") {
    # creating index to separate train and validation sets for k-fold validation 
    set.seed(1) 
    idx <- sample(rep_len(1:kcv, nrow(df_knn)))
    # loop to fit model for each fold
    res <- foreach (i = 1:kcv, .combine = rbind) %do% {
        # spliting datsets
        dTrain <- df_knn[idx != i, which(names(df_knn) != label)]
        dTest <- df_knn[idx == i, which(names(df_knn) != label)]
        vTrain <- df_knn[idx != i, which(names(df_knn) == label)]
        vTest <- df_knn[idx == i, which(names(df_knn) == label)]
        predictionKNN <- knn(train = dTrain, test = dTest, cl = vTrain, k = kNN)
        # it's also good to look at confusion table
        print(caret::confusionMatrix(predictionKNN, vTest)$overall[1])
        data.frame(
            foldID = i, 
            kNN = kNN, 
            validationAccuracy = caret::confusionMatrix(predictionKNN, vTest)$overall[1])
    }
    return(res)
}

# loop to find best k-neighbor
resultCVKNN <- foreach(kNN = 1:10, .combine = rbind) %do% {
    print(kNN)
    cvKNN(df_knn, kcv = 5, kNN = kNN)
}
resultCVKNN %>%
    group_by(kNN) %>% 
    summarise(meanAcc = mean(validationAccuracy)) %>%
    ggplot(., aes(as.factor(kNN), meanAcc, group = 1)) +
    geom_point() +
    geom_line() +
    theme_bw()

resultCVKNNb <- cvKNN(df_knn, kcv = 5, kNN = 1)

# Another method is naive bayes
# Sorry again, I borrowed the code :D

cvNBC <- function(df_knn, ModelFormula, kcv = 5, label="Genotype") {
    # creating index to separate train and validation sets for k-fold validation 
    set.seed(1) 
    idx <- sample(rep_len(1:kcv, nrow(df_knn)))
    # loop to fit model for each fold
    res <- foreach (i = 1:kcv, .combine = rbind) %do% {
        # spliting datsets
        dTrain <- df_knn[idx != i, ]
        dTest <- df_knn[idx == i, ]
        modelNBClass <- naiveBayes(ModelFormula, data = dTrain)
        predictionMBClass <- predict(modelNBClass, dTest)
        # it's also good to look at confusion table
        print(paste0("foldID ", i))
        print(caret::confusionMatrix(predictionMBClass, dTest$Genotype))
        data.frame(
            foldID = i, 
            validationAccuracy = caret::confusionMatrix(
                predictionMBClass, dTest$Genotype)$overall[1])
    }
    return(res)
}
resultCVNBCb <- cvNBC(df_knn, ModelFormula = formula(Genotype ~ .), kcv = 5)

# PCA for redusing dimentions 
library(stats)
library(pca3d)
dfPCA = prcomp(df[, 2:78])
summary(dfPCA)
plot(dfPCA)
biplot(dfPCA)
str(dfPCA)
pca3d(dfPCA, group= factor(df_copy$Genotype))
pca2d(dfPCA, group= factor(df_copy$Genotype))

# kmeans clustering
res = kmeans(df[, 2:78], centers = 2, iter.max = 100, nstart = 10)
table(res$cluster, df$Genotype)

# Let's try kmeans with pca object
dfkmeanspca = kmeans(dfPCA$x[, 1:3], centers = 2, iter.max = 1000,  nstart = 10)
table(dfkmeanspca$cluster, df$Genotype)
