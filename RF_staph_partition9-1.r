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
set.seed(12345) 


###########################################################
##########              Part II                  ##########
##########                                       ##########
##########    Partition ratio determined 9:1     ##########
###########################################################

# note: running the part of code above the repitition first, before running the code below. 

########## Train-Test partition 9:1 was selected ##########

test_AD_strains = dim(group_AD)[1]*0.1
test_HE_strains = dim(group_HE)[1]*0.1

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


# prepare input tables
table_training_group <- data.frame(t(table_staph_training))
table_training_group$group_training <- group_training[rownames(table_training_group), "HealthStatus"]

#determination of the number of tree
rf_ntree <- randomForest(group_training~.,data=table_training_group, ntree=5001)
plot(rf_ntree)
# best ntree = 1000


####################################
#####    For classification    #####
####################################

########## create RF model ##########
RF_classifier <- randomForest(x=table_training_group[,1:(ncol(table_training_group)-1)] , y=table_training_group[ ,ncol(table_training_group)], ntree=1000, mtry=500, importance=TRUE, proximities=TRUE )
RF_classifier

'''
Call:
 randomForest(x = table_training_group[, 1:(ncol(table_training_group) -      1)], y = table_training_group[, ncol(table_training_group)],      ntree = 1000, mtry = 500, importance = TRUE, proximities = TRUE) 
               Type of random forest: classification
                     Number of trees: 1000
No. of variables tried at each split: 500

        OOB estimate of  error rate: 10.74%
Confusion matrix:
    AD  HE class.error
AD 121  14   0.1037037
HE  15 120   0.1111111
'''

###Run Model with the test data
table_test_group <- data.frame(t(table_staph_test ))  
table_test_group$group_test <- group_test[rownames(table_test_group), "HealthStatus"]

#predict test dataset with RF model
prob_test_data <- predict(RF_classifier,table_test_group, type="prob")
predict_test_data <- as.data.frame(predict(RF_classifier,table_test_group))


# create and plot the ROC curve and calculate the AUC
rf_pred <- predict(RF_classifier, table_test_group, type="prob")[,2]
roc_obj <- roc(table_test_group$group_test, rf_pred, levels = rev(levels(table_test_group$group_test)))
auc_area <- auc(roc_obj)
auc_area
# Area under the curve: 0.9467
plot(roc_obj, main="ROC curve")


write.table(prob_test_data, file=paste0('prob_test_data_RF_classifier_split91.txt'), sep = '\t')
write.table(predict_test_data, file=paste0('predict_test_data_RF_classifier_split91.txt'),sep = '\t')


########## verify the model reliability based on the real-world test dataset ##########
# read real-world test dataset
table_test_staph48 <- read.table("test_data.txt", sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
table_test_staph48[table_test_staph48 > 0] = 1
group_staph48 <- read.table("test_metadata.txt", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

table_group_staph48 <- data.frame(t(table_test_staph48))  
table_group_staph48$group <- group_staph48[rownames(table_group_staph48), "HealthStatus"]

#predict the real-world dataset with RF model
prob_staph48 <- predict(RF_classifier,table_group_staph48, type="prob")
predict_staph48 <- as.data.frame(predict(RF_classifier,table_group_staph48))


# create and plot the ROC curve and calculate the AUC
rf_pred_staph48 <- predict(RF_classifier, table_group_staph48, type="prob")[,2]
roc_obj_staph48 <- roc(table_group_staph48$group, rf_pred_staph48, levels = rev(levels(table_group_staph48$group)))
auc_area_staph48 <- auc(roc_obj_staph48)
auc_area_staph48
# Area under the curve: 0.7636
plot(roc_obj_staph48, main="ROC curve")


write.table(prob_staph48, file=paste0('prob_staph48_RF_classifier_split91.txt'), sep = '\t')
write.table(predict_staph48, file=paste0('predict_staph48_RF_classifier_split91.txt'),sep = '\t')



####################################
#####    Model optimization    #####
####################################

########## Optimize RF model with the top important features ##########

# determine the number of top features by checking error rate
RF_classify_err_rate <- as.data.frame(RF_classifier$err.rate)
RF_classify_err_rate$features <- rownames( RF_classify_err_rate )
RF_classify_err_rate$features <- as.numeric(RF_classify_err_rate$features)
ggplot(RF_classify_err_rate[1:100,], aes(y=OOB, x=features)) + geom_point() + geom_smooth(color = "red", span=0.1, se=F) + scale_x_continuous(breaks = seq(min(RF_classify_err_rate$features), max(RF_classify_err_rate$features), by = 2))

'''
The lowest predicted error rate is achieved at the top 50 gene features
'''


# identifying the important Features
RF_classify_imp <- as.data.frame( RF_classifier$importance )
RF_classify_imp$features <- rownames( RF_classify_imp )
RF_classify_imp_sorted <- arrange( RF_classify_imp, desc(`MeanDecreaseAccuracy`)  )
barplot(RF_classify_imp_sorted$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")
barplot(RF_classify_imp_sorted[1:50,"MeanDecreaseAccuracy"], names.arg=RF_classify_imp_sorted[1:50,"features"] , ylab="Mean Decrease in Accuracy (Variable Importance)", las=2, ylim=c(0,0.08), main="RF Classification")  


#Select the top 50 features:
RF_table_top <- table_training_group[ , RF_classify_imp_sorted[1:50,"features"] ]
RF_table_top$group <- group_training[rownames(RF_table_top), "HealthStatus"]

# output thee MeanDecreaseAccuracy of all gene features
RF_classify_imp_sorted <- arrange(RF_classify_imp, desc(`MeanDecreaseAccuracy`)  )
write.table(RF_classify_imp_sorted, file='RF_classify_imp_sorted.txt',sep = '\t')


rf_ntree <- randomForest(group~.,data=RF_table_top, ntree=5001)
rf_ntree
# best ntree = 500

#Run the optimized RF model with test dataset
RF_classifier_top <- randomForest(x=RF_table_top[,1:(ncol(RF_table_top)-1)], y=RF_table_top[ , ncol(RF_table_top)] , ntree=500, mtry= 50, importance = TRUE )  
RF_classifier_top

'''
Call:
 randomForest(x = RF_table_top[, 1:(ncol(RF_table_top) - 1)],      y = RF_table_top[, ncol(RF_table_top)], ntree = 500, mtry = 50,      importance = TRUE) 
               Type of random forest: classification
                     Number of trees: 500
No. of variables tried at each split: 50

        OOB estimate of  error rate: 7.41%
Confusion matrix:
    AD  HE class.error
AD 127   8  0.05925926
HE  12 123  0.08888889
'''


#predict test dataset with the optimized RF model
prob_test_data_top <- predict(RF_classifier_top,table_test_group, type="prob")
predict_test_data_top <- as.data.frame(predict(RF_classifier_top,table_test_group))

# create and plot the ROC curve and calculate the AUC
rf_pred_top_test <- predict(RF_classifier_top, RF_table_top, type="prob")[,2]
roc_obj_top_test <- roc(RF_table_top$group, rf_pred_top_test, levels = rev(levels(RF_table_top$group)))
auc_area_top_test <- auc(roc_obj_top_test)
auc_area_top_test
# Area under the curve: 1
plot(roc_obj_top_test, main="ROC curve (Optimized RF)")


write.table(prob_test_data_top, file=paste0('prob_test_data_top_RF_classifier_top_split91.txt'), sep = '\t')
write.table(predict_test_data_top, file=paste0('predict_test_data_RF_classifier_top_split91.txt'),sep = '\t')


#predict the real-world dataset with RF model
prob_staph48_top <- predict(RF_classifier_top,table_group_staph48, type="prob")
predict_staph48_top <- as.data.frame(predict(RF_classifier_top,table_group_staph48))

# create and plot the ROC curve and calculate the AUC
rf_pred_top_staph48 <- predict(RF_classifier_top, table_group_staph48, type="prob")[,2]
roc_obj_top_staph48 <- roc(table_group_staph48$group, rf_pred_top_staph48, levels = rev(levels(table_group_staph48$group)))
auc_area_top_staph48 <- auc(roc_obj_top_staph48)
auc_area_top_staph48
# Area under the curve: 0.8182
plot(roc_obj_top_staph48, main="ROC curve (Optimized RF)")


# write output results
write.table(prob_staph48_top, file=paste0('prob_staph48_top_RF_classifier_top_split91.txt'), sep = '\t')
write.table(predict_staph48_top, file=paste0('predict_staph48_top_RF_classifier_top_split91.txt'),sep = '\t')



# generate the ROC curves in one figure
plot(roc_obj, col="skyblue", main="ROC curves")
lines(roc_obj_top_test, col="blue")
lines(roc_obj_staph48, col="red")
lines(roc_obj_top_staph48, col="darkred")
legend("bottomright",
legend = c("RF model (test data)",
    "Optimized RF model (test data)",
    "RF model (real-world data)",
    "Optimized RF model (real-world data)"),
    col = c("skyblue", "blue", "red", "darkred"),
    lty = 1, cex=0.5, lwd = 3)

