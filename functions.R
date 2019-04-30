# ROC
multi_roc <- function(roc_list) {
    for (ele in roc_list) {
        ele_temp <- eval(parse(text = ele))
        if (which(roc_list == ele) == 1) {
            df_temp <- data.frame(item = rep(ele, length(ele_temp@x.values[[1]])),
                                  x = ele_temp@x.values[[1]],
                                  y = ele_temp@y.values[[1]])
            x_name <- ele_temp@x.name[[1]]
            y_name <- ele_temp@y.name[[1]]
        } else {
            df_temp <- rbind(df_temp, data.frame(item = rep(ele, length(ele_temp@x.values[[1]])),
                                                 x = ele_temp@x.values[[1]],
                                                 y = ele_temp@y.values[[1]]))
        }
    }
    ggplot(df_temp) + geom_line(aes(x, y, col = item)) + labs(x = x_name, y = y_name)
}

#single kNN Regression
kNN_Reg <- function(train_x, train_y, test_x, k) {
    pred_y <- rep(0, length(test_x))
    for (x in test_x) {
        dis_list <- abs(train_x - x)
        select_y <- train_y[order(dis_list)][1:k]
        pred_y <- append(pred_y, mean(select_y))
    }
    return(pred_y)
}

# multi kNN Regression
kNN_mulReg <- function(train_x, train_y, test_x, k) {
    pred_y <- rep(0, length(test_x))
    for (i in 1:nrow(test_x)) {
        # euclidean distanct
        temp_df <- data.frame(V1 = rep(0,nrow(train_x)))
        for (d in 1:ncol(train_x)) {
            temp_df[,d] <- (train_x[,d] - test_x[i,d])^2
        }
        temp_df$dE <- sqrt(do.call("+", temp[,1:ncol(temp_df)]))
        select_y <- train_y[order(temp_df$dE)][1:k]
        #nonweighted
        pred_y <- append(pred_y, mean(train_y))
    }
}