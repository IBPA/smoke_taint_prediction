#Functions for predicting smoke taint index
run_lasso = function(selected_sensory,selected_features){
    pcc_lasso = rep(0, ncol(selected_sensory))
    mse_lasso = rep(0, ncol(selected_sensory))
    lambda_lasso = rep(0, ncol(selected_sensory))
    lambdas = 10^seq(3,-3,-0.01)
    coef_lasso = NULL
    
    pred_Y = NULL
    truth_Y = NULL
    fold_id = NULL
    
    #fold_id = sample(1:5,nrow(selected_features),replace = T)
    for (i in 1 : ncol(selected_sensory)){
        scores = selected_sensory[,i]
        selected_items = which((complete.cases(scores) & complete.cases(selected_features)) == T)
        
        cv_fit = cv.glmnet(as.matrix(selected_features[selected_items,]),
                           scores[selected_items],
                           lambda = lambdas,nfolds = 5,keep = T)
        
        idx = cv_fit$index[1]
        y_predicted = rep(NA, nrow(selected_features))
        cur_fold_id = rep(NA, nrow(selected_features))
        y_truth = rep(NA,nrow(selected_features))
        
        y_predicted[selected_items] = cv_fit$fit.preval[,idx]
        cur_fold_id[selected_items] = cv_fit$foldid
        y_truth[selected_items] = scores[selected_items]
        
        pcc_lasso[i] = cor(y_predicted[selected_items],scores[selected_items])
        mse_lasso[i] = mean((y_predicted[selected_items]-scores[selected_items])^2)
        lambda_lasso[i] = cv_fit$lambda.min
        
        fit = glmnet(selected_features[selected_items,],scores[selected_items],
                     lambda = cv_fit$lambda.min)
        coef_lasso = cbind(coef_lasso, as.matrix(coef(fit)))
        
        pred_Y = cbind(pred_Y, y_predicted)
        fold_id = cbind(fold_id, cur_fold_id)
        truth_Y = cbind(truth_Y, y_truth)
    }
    
    lasso_performance = rbind(pcc_lasso, mse_lasso, lambda_lasso)
    rownames(lasso_performance) = c("pcc","mse","lambda")
    colnames(lasso_performance) = colnames(selected_sensory)
    #write.csv(lasso_performance, "lasso_performance.csv")
    
    colnames(coef_lasso) = colnames(selected_sensory)
    #write.csv(coef_lasso,"coef_lasso.csv")
    
    return(list(lasso_performance = lasso_performance, 
                coef_lasso = coef_lasso,
                fold_id = fold_id,
                pred_Y = pred_Y,
                truth_Y = truth_Y
    ))
}


run_svr = function(selected_sensory,selected_features){
    pcc_svr = rep(0, ncol(selected_sensory))
    mse_svr = rep(0, ncol(selected_sensory))
    costs = rep(0, ncol(selected_sensory))
    gammas = rep(0, ncol(selected_sensory))
    epsilons = rep(0, ncol(selected_sensory))
    pred_Y = NULL
    all_truth_Y = NULL
    fold_id = NULL
    
    cost = 10^seq(-1,1,0.1) 
    gamma = 10^seq(-1.8,0.2,0.1)
    epsilon = 10^seq(-2.2,-1.8,0.1)
    
    svm_tune_results = list()
    
    tc <- tune.control(cross = 5)
    for (i in 1 : ncol(selected_sensory)){
        scores = selected_sensory[,i]
        selected_items = which((complete.cases(scores) & complete.cases(selected_features)) == T)
        
        prioir_svm_radial <- tune.svm(selected_features[selected_items,], 
                                      y = scores[selected_items], 
                                      cost = cost, 
                                      gamma= gamma, epsilon = epsilon,
                                      tunecontrol = tc)
        svm_tune_results[[i]] = prioir_svm_radial
        
        
        predict_Y = rep(NA, nrow(selected_features))
        truth_Y = rep(NA, nrow(selected_features))
        cur_fold_id = rep(NA, nrow(selected_features))
        for (j in 1 : length(prioir_svm_radial$train.ind)){
            idx_train = prioir_svm_radial$train.ind[[j]]
            idx_test = setdiff(1:nrow(selected_features[selected_items,]), idx_train)
            train_x = selected_features[selected_items,][idx_train,]
            train_y = scores[selected_items][idx_train]
            test_x = selected_features[selected_items,][idx_test,]
            test_y = scores[selected_items][idx_test]
            
            svm_cur_cv = svm(x = train_x,
                             y = train_y,
                             cost = prioir_svm_radial$best.parameters$cost,
                             gamma = prioir_svm_radial$best.parameters$gamma,
                             epsilon = prioir_svm_radial$best.parameters$epsilon)
            
            predYsvm = predict(svm_cur_cv, test_x)
            
            idx_available = selected_items
            predict_Y[idx_available[idx_test]] = predYsvm
            truth_Y[idx_available[idx_test]] = test_y
            cur_fold_id[idx_available[idx_test]] = j
        }
        
        pcc_svr[i] = cor(predict_Y, truth_Y,use = "complete.obs")
        mse_svr[i] = mean((predict_Y - truth_Y)^2,na.rm=T)
        costs[i] = prioir_svm_radial$best.parameters$cost
        gammas[i] = prioir_svm_radial$best.parameters$gamma
        epsilons[i] = prioir_svm_radial$best.parameters$epsilon
        
        pred_Y = cbind(pred_Y, predict_Y)
        fold_id = cbind(fold_id, cur_fold_id)
        all_truth_Y = cbind(all_truth_Y, truth_Y)
    }
    
    svr_performance = rbind(pcc_svr, mse_svr)
    rownames(svr_performance) = c("pcc","mse")
    colnames(svr_performance) = colnames(selected_sensory)
    #write.csv(svr_performance, "svr_performance.csv")
    
    colnames(pred_Y) = colnames(selected_sensory)
    #write.csv(cbind(fold_id,pred_Y),"svr_predict_CVs.csv")
    
    return(list(svr_performance = svr_performance,
                param_svr = rbind(costs,gammas,epsilons),
                fold_id = fold_id,
                pred_Y = pred_Y,
                truth_Y = all_truth_Y))
}


run_rf = function(selected_sensory,selected_features){
    tc <- tune.control(cross = 5)
    n_features = ncol(selected_features)
    mtry = round(c(0.25,0.50,0.75,1) * n_features)
    nodesize = c(1,3,5,10)
    ntree = c(50,100,200,500,1000,2000)
    #mtry = c(5,10)
    #nodesize = c(1,3,5)
    #ntree = c(50,100)
    
    pcc_rf = rep(0, ncol(selected_sensory))
    mse_rf = rep(0, ncol(selected_sensory))
    mtrys = rep(0, ncol(selected_sensory))
    nodesizes = rep(0, ncol(selected_sensory))
    ntrees = rep(0, ncol(selected_sensory))
    
    pred_Y = NULL
    fold_id = NULL
    all_truth_Y = NULL
    
    rf_tune_results = list()
    library(randomForest)
    for (i in 1 : ncol(selected_sensory)){
        scores = selected_sensory[,i]
        selected_items = which((complete.cases(scores) & complete.cases(selected_features)) == T)
        rf_tune_result <- tune.randomForest(selected_features[selected_items,], 
                                            y = scores[selected_items], 
                                            mtry = mtry, 
                                            nodesize= nodesize, ntree = ntree,
                                            tunecontrol = tc)
        rf_tune_results[[i]] = rf_tune_result
        
        
        predict_Y = rep(NA, nrow(selected_features))
        truth_Y = rep(NA, nrow(selected_features))
        cur_fold_id = rep(NA, nrow(selected_features))
        for (j in 1 : length(rf_tune_result$train.ind)){
            idx_train = rf_tune_result$train.ind[[j]]
            idx_test = setdiff(1:nrow(selected_features[selected_items,]), idx_train)
            train_x = selected_features[selected_items,][idx_train,]
            train_y = scores[selected_items][idx_train]
            test_x = selected_features[selected_items,][idx_test,]
            test_y = scores[selected_items][idx_test]
            
            rf_cur_cv = randomForest(x = train_x,
                                     y = train_y,
                                     nodesize = rf_tune_result$best.parameters$nodesize,
                                     mtry = rf_tune_result$best.parameters$mtry,
                                     ntree = rf_tune_result$best.parameters$ntree)
            
            predYsvm = predict(rf_cur_cv, test_x)
            
            idx_available = selected_items
            predict_Y[idx_available[idx_test]] = predYsvm
            truth_Y[idx_available[idx_test]] = test_y
            cur_fold_id[idx_available[idx_test]] = j
        }
        
        pcc_rf[i] = cor(predict_Y, truth_Y,use = "complete.obs")
        mse_rf[i] = mean((predict_Y - truth_Y)^2,na.rm=T)
        mtrys[i] = rf_tune_result$best.parameters$mtry
        nodesizes[i] = rf_tune_result$best.parameters$nodesize
        ntrees[i] = rf_tune_result$best.parameters$ntree
        
        pred_Y = cbind(pred_Y, predict_Y)
        fold_id = cbind(fold_id, cur_fold_id)
        all_truth_Y = cbind(all_truth_Y, truth_Y)
    }
    
    rf_performance = rbind(pcc_rf, mse_rf)
    rownames(rf_performance) = c("pcc","mse")
    colnames(rf_performance) = colnames(selected_sensory)
    #write.csv(rf_performance, "rf_performance.csv")
    
    colnames(pred_Y) = colnames(selected_sensory)
    #write.csv(cbind(fold_id,pred_Y),"rf_predict_CVs.csv")
    
    return(list(rf_performance = rf_performance,
                param_rf = rbind(mtrys,nodesizes,ntrees),
                fold_id = fold_id,
                pred_Y = pred_Y,
                truth_Y = all_truth_Y))
}
