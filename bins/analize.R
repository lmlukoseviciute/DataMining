p_load("dplyr", "ggplot2","foreach", "ggplot2", "datasets", "MASS", 
       "e1071", "class", "C50", "caret") 

# This script is the analyses of mice protein expression data
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
# Since the groups are in order let's mix the rows
df = df[sample(nrow(df)),]
write.csv(df, "../outputs/df.csv")

# And devide the data into 5 data sets
df_1 = df[1:216,]
df_2 = df[217:432,]
df_3 = df[433:648,]
df_4 = df[649:864,]
df_5 = df[865:1080,]
df_all = list(df_1, df_2, df_3, df_4, df_5)
df_index = list(1:216, 217:432, 433:648, 649:864, 865:1080)
df_index_ind = 1:5
write.csv(df_all, "../outputs/df_all.csv")

df_all_index = c(1:5)
kazkas1 = c('a','b','c','d','e')
kazkas2 = c(1:5)
kazkas = list(kazkas1,kazkas2)
kazkas_ind = (1:5)

# Let's try knn
library(class)
df_knn= df[,2:79]

for (i in 1:5){
    train_data = df[df_index != df_index[i] ,]
    write.csv(train_data, "../outputs/train_data.csv")
    test_data = df_all[df_all_index == i]
    knn_model <- knn(train = train_data[,2:78], test = test_data[],
                     cl = train_data$Genotype, k = 2)
}

knn_model <- knn(train = df[1:864,2:78], test = df[865:1080,2:78], cl = df[1:864,79], k = 1)
summary(knn_model)

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
        # print(caret::confusionMatrix(predictionKNN, vTest)$overall[1])
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
