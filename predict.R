#.libPaths("caret/")
#library("caret")
library("randomForest")

data=read.table("descriptors.csv",header=T)
print(colnames(data))
load("finmodel_9.rda")
write.table(file="results.csv",as.data.frame(predict.randomForest(final_model2,data,type="prob")$X1),row.names=F)

load("fff.rda")
write.table(file="results2.csv",as.data.frame(predict.randomForest(final_model3,data,type="prob")$I),row.names=F)
write.table(file="results3.csv",as.data.frame(predict.randomForest(final_model3,data,type="prob")$II),row.names=F)
write.table(file="results4.csv",as.data.frame(predict.randomForest(final_model3,data,type="prob")$III),row.names=F)
write.table(file="results5.csv",as.data.frame(predict.randomForest(final_model3,data,type="prob")$IV),row.names=F)
