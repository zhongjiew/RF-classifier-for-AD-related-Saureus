#install packages
install.packages('randomForest')
install.packages('rfUtilities')
install.packages('caret')
install.packages('e1071')
install.packages('doParallel')


#library required packages
library(foreach)
library(doParallel)
library("randomForest")
library("plyr") # for the "arrange" function
library("rfUtilities") # to test model significance
library(caret) # to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function 
library(pROC)
library(dplyr)


########## read training data ##########

#set work directory
setwd("D:/helmholtz/staphylococcus/genomes_ncbi/public_strains/All_assemblies_good/10.augsburg48_public300_binary_splits/")


#read the training data and metada into R
table_staph300 <- read.table("training_data.txt", sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
#transform  the gene abundance table into gene presence-absence table
table_staph300[table_staph300>0] =1


########## Pre-processing training data ##########
#Removing Rare Features
nonzero_counts <- apply(table_staph300, 1, function(y) sum(length(which(y > 0))))
hist(nonzero_counts, breaks=1000, col="grey", main="", ylab="Number of ", xlab="Number of Non-Zero Values")

#This R function removes features present in less than a specified proportion
remove_rare <- function( table , cutoff_pro ) {
    row2keep <- c()
    cutoff <- ceiling( cutoff_pro * ncol(table) )  
    for ( i in 1:nrow(table) ) {
        row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
        if ( row_nonzero > cutoff ) {
            row2keep <- c( row2keep , i)
        }
    }
    return( table [ row2keep , , drop=F ])
    }

#remove gene families present in <= 10% of samples:
table_staph300_rare_removed <- remove_rare(table=table_staph300, cutoff_pro=0.1)

#read metadata
group <- read.table("training_metadata.txt", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

dim(table_staph300)
# [1] 5281  300
dim(table_staph300_rare_removed)
# [1] 2982  300
dim(group)
# [1] 300   1

#select AD and HE group, repectively.
group_AD <- group[group$HealthStatus=='AD', ,drop=F]
group_HE <- group[group$HealthStatus=='HE', ,drop=F]

dim(group_AD)
# [1] 150   1
dim(group_HE)
# [1] 150   1



#set the seed to make the workflow reproducible
set.seed(111)


###########################################################
##########              Part I                   ##########
##########                                       ##########
##########    Partition ratio determination      ##########
###########################################################


# split the dataset into different Train-Test proportions:9-1,...,5-5
# 10 repetitions for each splitting
for (p in 1:5){
    for (i in 1:10){
        # different splitting proportions
        test_AD_strains = dim(group_AD)[1]*p*0.1
        test_HE_strains = dim(group_HE)[1]*p*0.1

        #training dataset (270) and test dataset(30 subsampled genomes)
        #subsample 15 AD genomes and 15 HE genomes from table_staph300 as test dataset
        group_AD_subsample <- group_AD[sample(nrow(group_AD),test_AD_strains, replace=F), ,drop=F]
        group_HE_subsample <- group_HE[sample(nrow(group_HE),test_HE_strains, replace=F), ,drop=F]

        #get the metadata of training dataset and test dataset
        group_test <- rbind(group_AD_subsample, group_HE_subsample)
        group_training <- group[!rownames(group) %in% c(rownames(group_test)), ,drop=F]

        #get the tables of training and test datasets 
        table_staph_training <- select(table_staph300_rare_removed, -c(rownames(group_test)))
        table_staph_test <- select(table_staph300_rare_removed, c(rownames(group_test)))

        
        ###create RF model
        #to prep input tables
        table_training_group <- data.frame(t(table_staph_training))
        table_training_group$group_training <- group_training[rownames(table_training_group), "HealthStatus"]

        #determination of the number of tree
        #rf_ntree <- randomForest(group_training~.,data=table_training_group, ntree=5001)
        #plot(rf_ntree)


        ####################################
        #####    For classification    #####
        ####################################

        ########## create RF model ##########
        RF_classifier <- randomForest(x=table_training_group[,1:(ncol(table_training_group)-1)] , y=table_training_group[ ,ncol(table_training_group)], ntree=4000, mtry=500, importance=TRUE, proximities=TRUE )
        RF_classifier


        ########## Test data ###########
        ###Run Model with test data
        table_test_group <- data.frame(t(table_staph_test))  
        table_test_group$group_test <- group_test[rownames(table_test_group), "HealthStatus"]

        #predict test dataset with RF model
        prob<-predict(RF_classifier,table_test_group, type="prob")
        predict<-as.data.frame(predict(RF_classifier,table_test_group))


        # create and plot the ROC curve and calculate the AUC
        rf_pred <- predict(RF_classifier, table_test_group, type="prob")[,2]
        roc_obj <- roc(table_test_group$group_test, rf_pred, levels = rev(levels(table_test_group$group_test)))
        auc_area <- auc(roc_obj)

        plot(roc_obj, main="ROC curve")


        # write output results
        write.table(prob, file=paste0('split',10-p,'-',p,'_prob_repetition',i,'.txt'),sep = '\t')
        write.table(predict, file=paste0('split',10-p,'-',p,'_predict_repetition',i,'.txt'),sep = '\t')
        write.table(auc_area, file=paste0('split',10-p,'-',p,'_auc_area_repetition',i,'.txt'),sep = '\t')
    }
}


