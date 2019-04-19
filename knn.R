# library
library(ggplot2)
library(tidyr)
library(dplyr)
library(parallel)

# functions
set.seed(1)
multi_roc <- function(roc_list) {
    for (ele in roc_list) {
        ele_temp <- eval(parse(text = ele))
        if (which(roc_list == ele) == 1) {
            df_temp <- data.frame(k_value = rep(ele, length(ele_temp@x.values[[1]])),
                                  x = ele_temp@x.values[[1]],
                                  y = ele_temp@y.values[[1]])
            x_name <- ele_temp@x.name[[1]]
            y_name <- ele_temp@y.name[[1]]
        } else {
            df_temp <- rbind(df_temp, data.frame(k_value = rep(ele, length(ele_temp@x.values[[1]])),
                                                 x = ele_temp@x.values[[1]],
                                                 y = ele_temp@y.values[[1]]))
        }
    }
    df_temp <- separate(df_temp, "k_value", c("temp", "k_value"), sep = "_")
    ggplot(df_temp) + geom_line(aes(x, y, col = k_value)) + labs(x = x_name, y = y_name)
}

k_cal <- function(k) {
    pred <- class::knn(train = train_in, test = test_in, cl = train_out, k = k)
    conf <- table(pred, test_out)
    acc <- sum(diag(conf)) / sum(conf)
    td_temp <- data.frame(k = k, acc = acc)
}

# info
path <- "/mnt/c/machinelearning/"
input_file <- "titanic.csv"
output_name <- "knn"
k_area <- 0.8

# load file
cat("------------loading file------------", "\n")
titanic <- readr::read_csv(file.path(path, input_file))

# dataset arrange
titanic_fil <- titanic[sample(nrow(titanic)),]
titanic_fil$Sex <- ifelse(titanic_fil$Sex == "male", 1, 0)
titanic$Age <- (titanic$Age - min(titanic$Age)) / (max(titanic$Age) - min(titanic$Age))
titanic_fil <- titanic_fil[,c(1,2,4,5,8)]

train <- titanic_fil[1:round(nrow(titanic_fil)*0.7),]
train_in <- train[,2:ncol(train)]
train_out <- train$Survived

test <- titanic_fil[(round(nrow(titanic_fil)*0.7)+1):nrow(titanic_fil),]
test_in <- test[,2:ncol(test)]
test_out <- test$Survived

# k value vs accuracy
cat("------------calculating k------------", "\n")
k_list <- 1:round(k_area * nrow(train))

# 2 threads process
cl <- makeCluster(2)
clusterExport(cl, c("train_in", "train_out", "test_in", "test_out"))
df_list <- parLapply(cl, k_list, k_cal)
acc_df <- Reduce("rbind", df_list)
stopCluster(cl)

ggplot(acc_df) + geom_point(aes(k, acc)) 
ggsave(paste0(path, "/", output_name, "_k.pdf"), dpi = "print")
# q()

# ROC curve
cat("-----------calculating ROC-----------", "\n")
acc_df <- acc_df[order(acc_df$acc, decreasing = TRUE),]
acc_df <- acc_df[!duplicated(acc_df$acc),][1:10,]
roc_list <- vector()
for (i in 1:nrow(acc_df)) {
    pred_temp <- class::knn(train = train_in, test = test_in, cl = train_out, k = acc_df$k[i])
    pred_temp2 <- ROCR::prediction(as.numeric(pred_temp), as.numeric(test_out))
    roc_name <- paste0("roc_", acc_df$k[i])
    assign(roc_name, ROCR::performance(pred_temp2, "tpr", "fpr"))
    roc_list <- append(roc_list, roc_name)
}
multi_roc(roc_list)

ggsave(paste0(path, "/", output_name, "_ROC.pdf"), dpi = "print")