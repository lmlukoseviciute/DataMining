# This script is the analyses of mice protein expression data
# Loading data
df = read.csv("../inputs/Data_Cortex_Nuclear.csv")
df_copy = df
# Looking to the structure of the data 
summary (df)
# Figuring out how much NAs are there
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
# So code below will add an additional column deviding data into 2 new classes

df$my_class [1:570] = 'control'
df$my_class [571:1080] = 'trisomic'
















