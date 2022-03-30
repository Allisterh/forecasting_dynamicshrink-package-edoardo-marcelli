# Mine Google Trends Data
#-------------------------
# In this script I download data for inflation and unemployment from Google Trend 
#
#
# Set working directory 
#-------------------------------

setwd("C:/Users/edoar/Dropbox/elementi salvati R")

# Load Package
#----------------

library(readr) 
library(gtrendsR) 
library(purrr) 
library(tidyverse)

# Select Trends
#---------------

wordlist<-c("what is inflation", 
            "what is deflation", 
            "us inflation",
            "us inflation rate", 
            "us inflation index", 
            "us deflation", 
            "united states inflation", 
            "real inflation",
            "rate of inflation",
            "price inflation",
            "price index",
            "national inflation",
            "inflation in us",
            "inflation graph", 
            "inflation forecast", 
            "inflation deflation", 
            "inflation definition",
            "inflation data", 
            "inflation chart", 
            "inflation calculator",
            "inflation and deflation",
            "historical inflation",
            "high inflation", 
            "deflation",
            "deflation rate", 
            "deflation interest rates",
            "deflation in us", 
            "deflation gold",
            "deflation economy", 
            "definition inflation",
            "definition deflation", 
            "define inflation", 
            "debt deflation", 
            "current inflation", 
            "current inflation rate"
)

country <- c('US') 
time <- ("all") 
channel <- 'web' 

trends <- gtrends("inflation", 
                  gprop = channel,
                  geo = country,
                  time = time ) 

out <- trends$interest_over_time 

GT_matrix = matrix(NA,nrow=dim(out)[1],ncol=length(wordlist)+1)
GT_matrix[,1] = out$hits

j=2
for(i in wordlist){
  
  country <- c('US') 
  time <- ("all") 
  channel <- 'web' 
  
  trends <- gtrends(i, 
                    gprop = channel,
                    geo = country,
                    time = time ) 
  
  results <- trends$interest_over_time
  
  GT_matrix[,j] = results$hits
  j=j+1
}

GT_matrix[GT_matrix<1] = 1
Inf_GT = data.frame(out$date,GT_matrix)

columns_names = c("Date","inflation",wordlist)
colnames(Inf_GT) = columns_names

write.csv(Inf_GT, "Inflation_GoogleTrends.csv")


# Download the dataframe "out" as a .csv file 
#---------------------------------------------
# Note: re-writing the dataset, it will include also the most recent searches
# therefore it would not be the same used in the in-depth study

write.csv(GT_trans, "GT_trans.csv")

################################################################################
################################################################################
################################################################################


# Now we repeat the analysis with Unemployment
#---------------------------------------------

# Select Trends
#---------------

wordlist<-c(
  "Careerbuilder",
  "Dice",
  "Unemployment compensation",
  "Glassdoor",
  "Unemployment agency",
  "Indeed.com",
  "LinkedIn",
  "Monster.com",
  "salary",
  "online jobs",
  "washington unemployment",
  "us unemployment rate",
  "unemployment statistics",
  "unemployment rate",
  "unemployment pa",
  "unemployment office",
  "unemployment michigan",
  "unemployment insurance",
  "unemployment great depression",
  "unemployment extension",
  "unemployment depression",
  "unemployment checks",
  "unemployment check",
  "unemployment benefits",
  "texas unemployment",
  "subsidies",
  "state compensation fund",
  "oregon unemployment",
  "ohio unemployment",
  "ny unemployment",
  "nj unemployment",
  "new york unemployment",
  "michigan works",
  "michigan works unemployment",
  "michigan state unemployment",
  "marvin unemployment",
  "marvin michigan unemployment",
  "job growth",
  "florida unemployment",
  "federal unemployment",
  "employee benefits",
  "depression",
  "compensation packages",
  "compensation package",
  "california unemployment" )

country <- c('US') 
time <- ("all") 
channel <- 'web' 

trends <- gtrends("unemployment", 
                  gprop = channel,
                  geo = country,
                  time = time ) 

out <- trends$interest_over_time 

GT_matrix = matrix(NA,nrow=dim(out)[1],ncol=length(wordlist)+1)
GT_matrix[,1] = out$hits

j=2
for(i in wordlist){
  
  country <- c('US') 
  time <- ("all") 
  channel <- 'web' 
  
  trends <- gtrends(i, 
                    gprop = channel,
                    geo = country,
                    time = time ) 
  
  results <- trends$interest_over_time
  
  GT_matrix[,j] = results$hits
  j=j+1
}

GT_matrix[GT_matrix<1] = 1
class(GT_matrix) <- "numeric"

columns_names = c("unemployment",wordlist)
colnames(GT_matrix) = columns_names

# Monthly GoogleTrend Dataset
#-----------------------------

Monthly_GoogleTrends = data.frame(out$date,GT_matrix)

write.csv(Monthly_GoogleTrends, "Monthly_GoogleTrends.csv")




# Quarterly transformation
#-------------------------



class(GT_matrix) <- "numeric"

GT_trans = matrix(NA, nrow=72,ncol=dim(GT_matrix)[2])
j=1
for(i in seq(1,214,by=3)){
  
  GT_trans[j,] = colMeans(GT_matrix[i:(i+2),])
  j=j+1
}

columns_names = c("Date","unemployment",wordlist)
colnames(GT_trans) = columns_names

GT_df = data.frame(GT_trans)
GT_df[,1] =  seq(as.Date("2004-01-01"), by = "quarter", length.out = 72)

# Download the dataframe "out" as a .csv file 
#---------------------------------------------
# Note: re-writing the dataset, it will include also the most recent searches
# therefore it would not be the same used in the in-depth study

write.csv(GT_trans, "GT_trans.csv")
