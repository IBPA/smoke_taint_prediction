#Smoke Taint Index Prediction Project

#Please set your working directory here.
WORKING_DIR = "C:\\Users\\Bigghost\\Documents\\GitHub\\smoke_taint_prediction"
setwd(WORKING_DIR)

#Create the output directory for storing the output figures
dir.create("output")

#Clean the working space and import necessary libraries
rm(list=ls())

#Library for data analysis
library(tsne)

#Library for prediction model
library(glmnet)             #(Lasso)
library(e1071)              #(SVM)


#Library for plotting
library(RColorBrewer)
library(circlize)           #(Color palette generator)
library(ComplexHeatmap)
library(ggplot2)


#Loading data
data = read.csv("data/dataset_for_analysis.csv",header=T,check.names = F,encoding="UTF-8",stringsAsFactors = F)
data$`Sensory_ash (smoke taint out of 100)` = as.numeric(data$`Sensory_ash (smoke taint out of 100)`)
grape_data = data[data$`Sample type` == "Grapes",]
wine_data = data[data$`Sample type` != "Grapes",]
colidx_voc = 20 : (ncol(data) - 1)
colidx_voc_with_binary_indicators = 9 : (ncol(data) - 1)


#Merge Grape Data: Median
merged_data_grape = apply(Grape_Data[,Feature_Col_Idx],2,function(x){
    sample_name = unique(Grape_Data$`Sample name`)
    results = numeric(length(sample_name))
    for (i in 1 : length(sample_name)){
        select_idx = Grape_Data$`Sample name` == sample_name[i]
        results[i] = median(x[select_idx])
    }
    names(results)= sample_name
    return(results)
})

#Merge Winen Data: Mean
merged_data_wine = apply(Wine_Data[,c(Feature_Col_Idx,N_Col)],2,function(x){
    sample_name = unique(paste0(Wine_Data$`Sample name`,"|",Wine_Data$`Sample type`))
    results = numeric(length(sample_name))
    for (i in 1 : length(sample_name)){
        select_idx = paste0(Wine_Data$`Sample name`,"|",Wine_Data$`Sample type`) == sample_name[i]
        results[i] = mean(x[select_idx])
    }
    names(results)= sample_name
    return(results)
})
write.csv(merged_data_grape,"grape_V5.csv")
write.csv(merged_data_wine,"wine_V5.csv")


#Merge grape data into wine
data_grape_merge_to_wine = lapply(1 : nrow(merged_data_wine),function(x){
    cur_grape_type = strsplit(rownames(merged_data_wine)[x],"|",fixed=T)[[1]][1]
    idx = which(cur_grape_type == rownames(merged_data_grape))
    if (length(idx) != 1){
        idx = which(startsWith(cur_grape_type, rownames(merged_data_grape)))
    }
    if (length(idx) == 0) return(NA)
    if (length(idx) != 1) stop()
    return(merged_data_grape[idx,,drop=F])
})

merged_data_both = NULL
for (i in 1 : nrow(merged_data_wine)){
    if (!is.na(data_grape_merge_to_wine[[i]])) {
        tmp = matrix((c(data_grape_merge_to_wine[[i]],merged_data_wine[i,,drop=F])),nrow=1)
        rownames(tmp) = rownames(merged_data_wine)[i]
        colnames(tmp) = c(paste0("(grape)",colnames(data_grape_merge_to_wine[[1]])),colnames(merged_data_wine))
    }else{
        next
    }
    merged_data_both = rbind(merged_data_both,tmp)
}

write.csv(merged_data_both,"both_v5.csv")



#Hist of wines
svg("output/histogram_index.svg",width=6,height=6)
dd<-data.frame(x=merged_data_both[,"Sensory_ash (smoke taint out of 100)"])
ggplot(dd, aes(x=x))+ geom_histogram(breaks=seq(0,70,10),fill="white",color="black",boundary=0) + theme_classic() +
    ggtitle(paste0("Distribution of Smoke Taint Index (N = ", nrow(merged_data_both), ")")) + 
    #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    #scale_fill_manual(values=ColorSet)+
    labs(y="Count",x="Smoke Taint Index") + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 20))  + 
    scale_x_continuous(breaks=seq(0,70,10)) +
    stat_bin(breaks=seq(0,70,10), geom="text", color='black', aes(label=..count..),
             position=position_stack(vjust = 0.5), boundary=0) 
dev.off()


#PCA
pca_result = prcomp(merged_data_both[,1:(ncol(merged_data_both)-1)],center = T,scale. = T)
pc1 = pca_result$x[,1]
pc2 = pca_result$x[,2]
pc3 = pca_result$x[,3]

pc1_label = paste0("PC 1 (", format(100*pca_result$sdev[1]/sum(pca_result$sdev), digits=2), "% of variance explained)")
pc2_label = paste0("PC 2 (", format(100*pca_result$sdev[2]/sum(pca_result$sdev), digits=2), "% of variance explained)")
abs_pc2_label = paste0("Abs(PC 2) (", format(100*pca_result$sdev[2]/sum(pca_result$sdev), digits=2), "% of variance explained)")
pc3_label = paste0("PC 3 (", format(100*pca_result$sdev[3]/sum(pca_result$sdev), digits=2), "% of variance explained)")

score = sapply(merged_data_both[,"Sensory_ash (smoke taint out of 100)"], function(x){
    if (x <= 10) return("<= 10")
    if (x <= 20) return("10 < x <= 20")
    if (x <= 30) return("20 < x <= 30")
    if (x <= 40) return("30 < x <= 40")
    if (x <= 50) return("40 < x <= 50")
    if (x <= 60) return("50 < x <= 60")
    if (x <= 70) return("60 < x <= 70")
})


point_size = 3


#SVM boundary (on PCA)


df = data.frame(pc1 = pc1, pc2 = pc2)
high_low_index = rep("Low",nrow(df))
high_low_index[merged_data_both[,"Sensory_ash (smoke taint out of 100)"] > 30] = "High"
df = cbind(df, high_low_index)
svm_res = svm(high_low_index~pc1+pc2, data = df)

grid_x = seq(-4, 18, length.out = 50)
grid_y = seq(-7, 7, length.out = 50)
grids = lapply(grid_x,function(x){
    cbind(rep(x,length(grid_y)),grid_y)
})
grids = do.call(rbind,grids)
df_grids = data.frame(pc1 = grids[,1], pc2 = grids[,2])
predict_grids = predict(svm_res, df_grids)
predict_grids = as.character(predict_grids)
predict_grids[predict_grids == "High"] = 1
predict_grids[predict_grids == "Low"] = 0
predict_grids = as.numeric(predict_grids)

df = data.frame(pc1 = pc1, pc2 = pc2, pc3 = pc3, 
                score = score, 
                cont_score = merged_data_both[,"Sensory_ash (smoke taint out of 100)"],
                abs_pc2 = abs(pc2))
pcc_val = cor(df$pc1,df$cont_score)
pcc_pval = pval_pcc(pcc_val, length(df$pc1))

pcc_val = cor(abs(df$pc2),df$cont_score)
pcc_pval = pval_pcc(pcc_val, length(df$pc1))



df_val = data.frame(x = grids[,1], y = grids[,2], zi = as.factor(predict_grids), z = predict_grids)

svg("output/pca_pc1_pc2.svg",width=6,height=4)
ggplot(df, aes(x = pc1, y = pc2, color = score)) + geom_point(size=point_size) + 
    scale_size(guide="none")  +
    theme_classic() +
    ggtitle(paste0("PCA Plot (PC1 VS PC2) (N = ", nrow(merged_data_both), ")")) + 
    #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    labs(y=pc2_label,x=pc1_label,color="Smoke Taint Index") +
    scale_color_brewer(palette="RdBu", direction=-1) + 
    stat_contour(data = df_val, aes(x=x,y=y,z=z), color="black", linetype=1,linewidth=1, breaks=c(0,1)) + 
    geom_point(data = df[svm_res$index,], aes(x = pc1, y = pc2), color = "black", shape=4)
dev.off()

svg("output/pca_pc1_score.svg",width=4,height=4)
ggplot(df, aes(x = pc1, y = cont_score)) + geom_point(size=point_size) + 
    scale_size(guide="none")  +
    theme_classic() +
    ggtitle(paste0("PC1 VS Smoke Taint Index (N = ", nrow(merged_data_both), ")", "\n(PCC = ",format(cor(df$pc1,df$cont_score),digits=2) ,")")) + 
    #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    labs(y="Smoke Taint Index",x=pc1_label) +
    scale_color_brewer(palette="RdBu", direction=-1)
dev.off()

svg("output/pca_pc2_score.svg",width=4,height=4)
ggplot(df, aes(x = abs_pc2, y = cont_score)) + geom_point(size=point_size) + 
    scale_size(guide="none")  +
    theme_classic() +
    ggtitle(paste0("Abs(PC2) VS Smoke Taint Index (N = ", nrow(merged_data_both), ")", "\n(PCC = ",format(cor(df$abs_pc2,df$cont_score),digits=2) ,")")) + 
    #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    labs(y="Smoke Taint Index",x=abs_pc2_label) +
    scale_color_brewer(palette="RdBu", direction=-1)
dev.off()

svg("output/pca_pc1_pc3.svg",width=6,height=4)
ggplot(df, aes(x = pc1, y = pc3, color = score)) + geom_point(size=point_size) + 
    scale_size(guide="none")  +
    theme_classic() +
    ggtitle(paste0("PCA Plot (PC1 VS PC3) (N = ", nrow(merged_data_both), ")")) + 
    #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    labs(y=pc3_label,x=pc1_label,color="Smoke Taint Index")+
    scale_color_brewer(palette="RdBu", direction=-1)
dev.off()

svg("output/pca_pc2_pc3.svg",width=6,height=4)
ggplot(df, aes(x = pc2, y = pc3, color = score)) + geom_point(size=point_size) + 
    scale_size(guide="none")  +
    theme_classic() +
    ggtitle(paste0("PCA Plot (PC2 VS PC3) (N = ", nrow(merged_data_both), ")")) + 
    #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    labs(y=pc3_label,x=pc2_label,color="Smoke Taint Index")+
    scale_color_brewer(palette="RdBu", direction=-1)
dev.off()

write.csv(pca_result$rotation, "output/pca_rotation.csv")



plot_pc_loading(paste0("output/pc",1,"_loading.svg"), 1, c(-0.05,0.25))
plot_pc_loading(paste0("output/pc",2,"_loading.svg"), 2, c(-0.35,0.25))



#PCA (Grape)
pca_result = prcomp(merged_data_both[,1:20],center = T,scale. = T)
pc1 = pca_result$x[,1]
pc2 = pca_result$x[,2]
pc3 = pca_result$x[,3]

pc1_label = paste0("PC 1 (", format(100*pca_result$sdev[1]/sum(pca_result$sdev), digits=2), "% of variance explained)")
pc2_label = paste0("PC 2 (", format(100*pca_result$sdev[2]/sum(pca_result$sdev), digits=2), "% of variance explained)")
abs_pc2_label = paste0("Abs(PC 2) (", format(100*pca_result$sdev[2]/sum(pca_result$sdev), digits=2), "% of variance explained)")
pc3_label = paste0("PC 3 (", format(100*pca_result$sdev[3]/sum(pca_result$sdev), digits=2), "% of variance explained)")

score = sapply(merged_data_both[,"Sensory_ash (smoke taint out of 100)"], function(x){
    if (x <= 10) return("<= 10")
    if (x <= 20) return("10 < x <= 20")
    if (x <= 30) return("20 < x <= 30")
    if (x <= 40) return("30 < x <= 40")
    if (x <= 50) return("40 < x <= 50")
    if (x <= 60) return("50 < x <= 60")
    if (x <= 70) return("60 < x <= 70")
})


point_size = 3


#SVM boundary (on PCA)
df = data.frame(pc1 = pc1, pc2 = pc2)
high_low_index = rep("Low",nrow(df))
high_low_index[merged_data_both[,"Sensory_ash (smoke taint out of 100)"] > 30] = "High"
df = cbind(df, high_low_index)
svm_res = svm(high_low_index~pc1+pc2, data = df)

grid_x = seq(-4, 18, length.out = 50)
grid_y = seq(-7, 7, length.out = 50)
grids = lapply(grid_x,function(x){
    cbind(rep(x,length(grid_y)),grid_y)
})
grids = do.call(rbind,grids)
df_grids = data.frame(pc1 = grids[,1], pc2 = grids[,2])
predict_grids = predict(svm_res, df_grids)
predict_grids = as.character(predict_grids)
predict_grids[predict_grids == "High"] = 1
predict_grids[predict_grids == "Low"] = 0
predict_grids = as.numeric(predict_grids)

df = data.frame(pc1 = pc1, pc2 = pc2, pc3 = pc3, 
                score = score, 
                cont_score = merged_data_both[,"Sensory_ash (smoke taint out of 100)"],
                abs_pc2 = abs(pc2))
pcc_val = cor(df$pc1,df$cont_score)
pcc_pval = pval_pcc(pcc_val, length(df$pc1))

pcc_val = cor(abs(df$pc2),df$cont_score)
pcc_pval = pval_pcc(pcc_val, length(df$pc1))

df_val = data.frame(x = grids[,1], y = grids[,2], zi = as.factor(predict_grids), z = predict_grids)

svg("output/pca_pc1_pc2_grape.svg",width=6,height=4)
ggplot(df, aes(x = pc1, y = pc2, color = score)) + geom_point(size=point_size) + 
    scale_size(guide="none")  +
    theme_classic() +
    ggtitle(paste0("PCA Plot (PC1 VS PC2) (N = ", nrow(merged_data_both), ")")) + 
    #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    labs(y=pc2_label,x=pc1_label,color="Smoke Taint Index") +
    scale_color_brewer(palette="RdBu", direction=-1) + 
    stat_contour(data = df_val, aes(x=x,y=y,z=z), color="black", linetype=1,linewidth=1, breaks=c(0,1)) + 
    geom_point(data = df[svm_res$index,], aes(x = pc1, y = pc2), color = "black", shape=4)
dev.off()

svg("output/pca_pc1_score_grape.svg",width=4,height=4)
ggplot(df, aes(x = pc1, y = cont_score)) + geom_point(size=point_size) + 
    scale_size(guide="none")  +
    theme_classic() +
    ggtitle(paste0("PC1 VS Smoke Taint Index (N = ", nrow(merged_data_both), ")", "\n(PCC = ",format(cor(df$pc1,df$cont_score),digits=2) ,")")) + 
    #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    labs(y="Smoke Taint Index",x=pc1_label) +
    scale_color_brewer(palette="RdBu", direction=-1)
dev.off()

svg("output/pca_pc2_score_grape.svg",width=4,height=4)
ggplot(df, aes(x = abs_pc2, y = cont_score)) + geom_point(size=point_size) + 
    scale_size(guide="none")  +
    theme_classic() +
    ggtitle(paste0("Abs(PC2) VS Smoke Taint Index (N = ", nrow(merged_data_both), ")", "\n(PCC = ",format(cor(df$abs_pc2,df$cont_score),digits=2) ,")")) + 
    #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    labs(y="Smoke Taint Index",x=abs_pc2_label) +
    scale_color_brewer(palette="RdBu", direction=-1)
dev.off()

plot_pc_loading(paste0("output/pc",1,"_loading_grape.svg"), 1, c(0,0.4))
plot_pc_loading(paste0("output/pc",2,"_loading_grape.svg"), 2, c(-0.40,0.40))



#TSNE
find_region = function(x){
    x = strsplit(x,"|",fixed=T)[[1]][1]
    idx = which(x == Data$`Sample name`)
    if (length(idx) == 0) stop()
    return(Data$`Region (AVA)`[idx[1]])
}

find_county = function(x){
    x = strsplit(x,"|",fixed=T)[[1]][1]
    idx = which(x == Data$`Sample name`)
    if (length(idx) == 0) stop()
    return(Data$County[idx[1]])
}


approaches = c("none", "min_max","z","quantile")
n_neighbors = c(2,5,10,20)
point_size = 3

for (cur_approach in approaches){
    dir.create(paste0("output/tsne_",cur_approach))
    
    if (cur_approach == "none"){
        data_tsne = merged_data_both[,1:(ncol(merged_data_both)-1)]
        normalize_approach_descriptor = "(Not Normalizaed)"
    }else if (cur_approach == "min_max"){
        data_tsne = apply(merged_data_both[,1:(ncol(merged_data_both)-1)],2,function(x){
            return((x - min(x))/(max(x) - min(x)))  
        })
        normalize_approach_descriptor = "(Min-Max Normalizaed)"
    }else if (cur_approach == "z"){
        data_tsne = apply(merged_data_both[,1:(ncol(merged_data_both)-1)],2,function(x){
            return((x - mean(x))/sd(x))  
        })
        normalize_approach_descriptor = "(Z Normalizaed)"
    }else{
        new_data = t(normalize.quantiles(t(merged_data_both[,1:(ncol(merged_data_both)-1)])))
        data_tsne = merged_data_both[,1:(ncol(merged_data_both)-1)]
        data_tsne[,] = new_data[,]
        normalize_approach_descriptor = "(Quantile Normalizaed)"
    }
    
    
    for (cur_n_neighbor in n_neighbors){
        cur_dir = paste0("output/tsne_",cur_approach,"/n_neighbor_",cur_n_neighbor)
        dir.create(cur_dir)
        
        n_neighbor_descriptor = paste0("(# of neighbor = ", cur_n_neighbor,")")
        
        tsne_title = paste0("TSNE (N = ", nrow(merged_data_both), ")", n_neighbor_descriptor, "\n",normalize_approach_descriptor)
        
        tsne_res = tsne(data_tsne,perplexity = cur_n_neighbor,max_iter = 10000)

        df_tsne = data.frame(dim1 = tsne_res[,1], dim2 = tsne_res[,2], score = score)
        point_size = 3
        
        svg(paste0(cur_dir,"/tsne_index.svg"),width=6,height=4)
        h = ggplot(df_tsne, aes(x = dim1, y = dim2, color = score)) + geom_point(size=point_size) + 
            scale_size(guide="none")  +
            theme_classic() +
            ggtitle(tsne_title) + 
            #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(size=12),
                  axis.title.x = element_text(size=12),
                  axis.text.y = element_text(size=12),
                  axis.title.y = element_text(size=12)) +
            labs(y="Dim 2",x="Dim 1",color="Smoke Taint Index")+
            scale_color_brewer(palette="RdBu", direction=-1)
        print(h)
        dev.off()
        
        
        
        
        df_tsne = data.frame(dim1 = tsne_res[,1], dim2 = tsne_res[,2], 
                             type = sapply(strsplit(rownames(merged_data_both),"|",fixed=T),function(x){x[2]}))
        svg(paste0(cur_dir,"/tsne_Scale.svg"),width=6,height=4)
        h = ggplot(df_tsne, aes(x = dim1, y = dim2, color = type)) + geom_point(size=point_size) + 
            scale_size(guide="none")  +
            theme_classic() +
            ggtitle(tsne_title) + 
            #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(size=12),
                  axis.title.x = element_text(size=12),
                  axis.text.y = element_text(size=12),
                  axis.title.y = element_text(size=12)) +
            labs(y="Dim 2",x="Dim 1",color="Scale")+
            scale_color_brewer(palette="Set3")
        print(h)
        dev.off()
        
        df_tsne = data.frame(dim1 = tsne_res[,1], dim2 = tsne_res[,2], 
                             type = substr(rownames(merged_data_both),3,4))
        svg(paste0(cur_dir,"/tsne_variety.svg"),width=6,height=4)
        h = ggplot(df_tsne, aes(x = dim1, y = dim2, color = type)) + geom_point(size=point_size) + 
            scale_size(guide="none")  +
            theme_classic() +
            ggtitle(tsne_title) + 
            #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(size=12),
                  axis.title.x = element_text(size=12),
                  axis.text.y = element_text(size=12),
                  axis.title.y = element_text(size=12)) +
            labs(y="Dim 2",x="Dim 1",color="Variety")
        print(h)
        dev.off()
        
        
        df_tsne = data.frame(dim1 = tsne_res[,1], dim2 = tsne_res[,2], 
                             type = sapply(rownames(merged_data_both),find_region))
        svg(paste0(cur_dir,"/tsne_region.svg"),width=8,height=4)
        h = ggplot(df_tsne, aes(x = dim1, y = dim2, color = type)) + geom_point(size=point_size) + 
            scale_size(guide="none")  +
            theme_classic() +
            ggtitle(tsne_title) + 
            #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(size=12),
                  axis.title.x = element_text(size=12),
                  axis.text.y = element_text(size=12),
                  axis.title.y = element_text(size=12)) +
            labs(y="Dim 2",x="Dim 1",color="Region")
        print(h)
        dev.off()
        
        
        df_tsne = data.frame(dim1 = tsne_res[,1], dim2 = tsne_res[,2], 
                             type = sapply(rownames(merged_data_both),find_county))
        svg(paste0(cur_dir,"/tsne_county.svg"),width=6,height=4)
        h = ggplot(df_tsne, aes(x = dim1, y = dim2, color = type)) + geom_point(size=point_size) + 
            scale_size(guide="none")  +
            theme_classic() +
            ggtitle(tsne_title) + 
            #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(size=12),
                  axis.title.x = element_text(size=12),
                  axis.text.y = element_text(size=12),
                  axis.title.y = element_text(size=12)) +
            labs(y="Dim 2",x="Dim 1",color="County")+
            scale_color_brewer(palette="Set3")
        print(h)
        dev.off()
    }
}



#Heatmap
prepareAnnotation = function(x, colorMapName, reverse=F, order_by_name = F){
    tab = table(x)
    tab = tab[order(tab,decreasing = T)]
    colors = brewer.pal(n = length(tab), colorMapName)
    if (reverse) colors = colors[length(colors):1]
    if (length(tab) > length(colors)){
        new_tab = tab[1:length(colors)]
        new_tab[length(new_tab)] = sum(tab) - sum(tab[1:(length(colors)-1)])
        names(new_tab)[length(new_tab)] = "others"
        
        new_x = x
        new_x[which(!(x %in% names(tab)[1:(length(colors)-1)]))] = "others"
    }else{
        new_tab = tab
        new_x = x
    }
    
    if (order_by_name){
        names(colors) = names(new_tab)[order(names(new_tab))]
    }else{
        names(colors) = names(new_tab)
    }
   
    return(list(x = new_x, colors = colors))
}


scale_annotation =  prepareAnnotation(sapply(strsplit(rownames(merged_data_both),"|",fixed=T),function(x){x[2]}),"Pastel1")
variety_annotation = prepareAnnotation(substr(rownames(merged_data_both),3,4),"Set3")
county_annotation = prepareAnnotation(sapply(rownames(merged_data_both),find_county),"Pastel2")
score_annotation = prepareAnnotation( sapply(merged_data_both[,"Sensory_ash (smoke taint out of 100)"], function(x){
    if (x <= 10) return("<= 10")
    if (x <= 20) return("10 < x <= 20")
    if (x <= 30) return("20 < x <= 30")
    if (x <= 40) return("30 < x <= 40")
    if (x <= 50) return("40 < x <= 50")
    if (x <= 60) return("50 < x <= 60")
    if (x <= 70) return("60 < x <= 70")
}), "RdBu", reverse = T, order_by_name = T)

ha = HeatmapAnnotation(
    "Scale" = scale_annotation$x,
    "Variety" = variety_annotation$x,
    "County" = county_annotation$x,
    "Smoke Taint Index" = score_annotation$x,
    col = list("Scale" = scale_annotation$colors,
               "Variety" = variety_annotation$colors,
               "County" = county_annotation$colors,
               "Smoke Taint Index" = score_annotation$colors),
    border = F
)

###(TEST REGION)
title=paste0("Heatmap (N = ", nrow(merged_data_both), ")")
testProfiles = apply(merged_data_both[,1:(ncol(merged_data_both)-1)],2,function(x){(x-min(x))/(max(x)-min(x))})
svg(filename = "output/Heatmap.svg", width = 10, height = 8)
ht_list = Heatmap(t(testProfiles), 
                  col = colorRamp2(c(0, 0.5, 1), c("green", "white", "red")),
                  name = "Normalized Compound Concentration", clustering_distance_rows = "euclidean",
                  clustering_method_rows = "complete", column_gap = unit(1, "mm"), border=F,
                  column_title = title,show_column_names = F,show_row_names = T, top_annotation = ha, 
                  show_row_dend = F,
                  heatmap_legend_param = list(direction="horizontal",title_position="leftcenter",x=unit(0,"npc"), just = c("left", "bottom"))
)
draw(ht_list, merge_legend = F, heatmap_legend_side = "bottom", legend_grouping="original",
     annotation_legend_side = "right")
dev.off()

#Analysis
avg_minmax = apply(t(testProfiles),2,mean)
pcc_val = cor(avg_minmax, merged_data_both[,ncol(merged_data_both)])
t_score = (pcc_val * sqrt(length(avg_minmax)-2)) / sqrt(1 - pcc_val^2)
p_val = pt(t_score, length(avg_minmax)-2, lower.tail = F)
#
smoke_taint_index = merged_data_both[,ncol(merged_data_both)]
county_info = sapply(rownames(merged_data_both),find_county)
t.test(smoke_taint_index[county_info=="Napa"],smoke_taint_index[county_info!="Napa"],alternative = "greater")
t.test(smoke_taint_index[county_info=="Sonoma"],smoke_taint_index[county_info!="Sonoma"])

#Uni-variate Analysis
pcc_smoke_idx = apply(merged_data_both, 2, function(x){cor(x,merged_data_both[,ncol(merged_data_both)])})
spearman_smoke_idx = apply(merged_data_both, 2, function(x){cor(x,merged_data_both[,ncol(merged_data_both)],method = "spearman")})
log_pcc_smoke_idx = apply(log(merged_data_both+1), 2, function(x){cor(x,merged_data_both[,ncol(merged_data_both)])})

correlation_matrix = cbind(pcc_smoke_idx, log_pcc_smoke_idx, spearman_smoke_idx)
colnames(correlation_matrix) = c("PCC","PCC (Log value)","Spearman")
feature_name = rownames(correlation_matrix)
feature_name = sapply(feature_name, function(x){gsub("_"," ",x)})
feature_name = sapply(feature_name, function(x){gsub(")",") ",x)})
rownames(correlation_matrix) = feature_name
correlation_matrix = correlation_matrix[-nrow(correlation_matrix),]


svg("output/heatmap_univariate.svg",width=7,height=8)
Heatmap(correlation_matrix, 
        name = "correlation", 
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        column_title = "Correlation (Features VS Smoke Taint Index)",
        column_title_gp = gpar(fontsize = 14),
        cluster_columns = F,
        cluster_rows = T,
        show_row_dend = F,
        row_names_side = "left",
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(formatC(correlation_matrix[i,j], digits=2, format="f"), x, y)
        }, 
        show_column_names = FALSE, 
        top_annotation = HeatmapAnnotation(
            text = anno_text(colnames(correlation_matrix), rot = 0, just = "center"),
            annotation_height = max_text_height(colnames(correlation_matrix))
        ), 
        show_heatmap_legend = FALSE)
dev.off()



plot_cor("output/cor_free_m-cresol.svg",
         5,5,
         x = merged_data_both[,"free_m-cresol"], 
         y = merged_data_both[,ncol(merged_data_both)],
         paste0("Free m-cresol VS Smoke Taint Index"),
         "Concentration (ug/L))"
         )
plot_cor_log("output/cor_free_m-cresol_log.svg",
         2,2,
         x = merged_data_both[,"free_m-cresol"], 
         y = merged_data_both[,ncol(merged_data_both)]
)


plot_cor("output/cor_total_syringol_grape.svg",
         5,5,
         x = merged_data_both[,"(grape)total_syringol"], 
         y = merged_data_both[,ncol(merged_data_both)],
         paste0("Total syringol (Grape) VS Smoke Taint Index"),
         "Concentration (ug/Kg))"
)
plot_cor_log("output/cor_total_syringol_grape_log.svg",
         2,2,
         x = merged_data_both[,"(grape)total_syringol"], 
         y = merged_data_both[,ncol(merged_data_both)]
)


plot_cor("output/cor_total_4-methylsyringol_grape.svg",
         5,5,
         x = merged_data_both[,"(grape)total_4-methylsyringol"], 
         y = merged_data_both[,ncol(merged_data_both)],
         paste0("Total 4-Methylsyringol (Grape) VS Smoke Taint Index"),
         "Concentration (ug/Kg))"
)
plot_cor_log("output/cor_total_4-methylsyringol_grape_log.svg",
         2,2,
         x = merged_data_both[,"(grape)total_4-methylsyringol"], 
         y = merged_data_both[,ncol(merged_data_both)]
)





raw_grape_data_with_index = lapply(1 : nrow(merged_data_both), function(x){
    cur_grape_type = strsplit(rownames(merged_data_both)[x],"|",fixed=T)[[1]][1]
    idx = which(cur_grape_type == Grape_Data[,1])
    if (length(idx) != 3){
        idx = which(startsWith(cur_grape_type, Grape_Data[,1]))
    }
    if (length(idx) != 3) stop()
    
    grape_type = rep(cur_grape_type, length(idx))
    smokeTaintIndex = rep(merged_data_both[x,ncol(merged_data_both)],length(idx))
    return((cbind(grape_type, Grape_Data[idx,Feature_Col_Idx,drop=F],smokeTaintIndex)))
})

raw_grape_data_with_index = do.call(rbind, raw_grape_data_with_index)

write.csv(cor_smoke_idx, "output/cor_features.csv")

svg("output/total_syringol_grape.svg",width=6,height=6)
plot(raw_grape_data_with_index$total_syringol, 
     raw_grape_data_with_index$smokeTaintIndex,
     main = paste0("Total Syringol (Grape) VS Smoke Taint Index\n(N = 55 (Wines) * 3 (replicates)) (PCC = ", format(cor_features["(grape)total_syringol"],digits=2),")"),
     xlab = "Total Syringol (Grape) (ug/Kg)", ylab = "Smoke Taint Index",pch=19)
dev.off()

svg("output/free_guaiacol_grape.svg",width=6,height=6)
plot(raw_grape_data_with_index$free_guaiacol, 
     raw_grape_data_with_index$smokeTaintIndex,
     main = paste0("Free Guaiacol (Grape) VS Smoke Taint Index\n(N = 55 (Wines) * 3 (replicates)) (PCC = ", format(cor_smoke_idx["(grape)free_guaiacol"],digits=2),")"),
     xlab = "Free Guaiacol (Grape) (ug/Kg)", ylab = "Smoke Taint Index",pch=19)
dev.off()

svg("output/free_4-ethylphenol_grape.svg",width=6,height=6)
plot(raw_grape_data_with_index$`free_4-ethylphenol`, 
     raw_grape_data_with_index$smokeTaintIndex,
     main = paste0("Free 4-Ethylphenol (Grape) VS Smoke Taint Index\n(N = 55 (Wines) * 3 (replicates)) (PCC = ", format(cor_smoke_idx["(grape)free_4-ethylphenol"],digits=2),")"),
     xlab = "Free 4-Ethylphenol (Grape) (ug/Kg)", ylab = "Smoke Taint Index",pch=19)
dev.off()

svg("output/total_o-cresol_grape.svg",width=6,height=6)
plot(raw_grape_data_with_index$`total_o-cresol`, 
     raw_grape_data_with_index$smokeTaintIndex,
     main = paste0("Total o-Cresol (Grape) VS Smoke Taint Index\n(N = 55 (Wines) * 3 (replicates)) (PCC = ", format(cor_smoke_idx["(grape)total_o-cresol"],digits=2),")"),
     xlab = "Total o-Cresol (Grape) (ug/Kg)", ylab = "Smoke Taint Index",pch=19)
dev.off()

svg("output/total_4-methylsyringol_grape.svg",width=6,height=6)
plot(raw_grape_data_with_index$`total_4-methylsyringol`, 
     raw_grape_data_with_index$smokeTaintIndex,
     main = paste0("Total 4-Methylsyringol (Grape) VS Smoke Taint Index\n(N = 55 (Wines) * 3 (replicates)) (PCC = ", format(cor_smoke_idx["(grape)total_4-methylsyringol"],digits=2),")"),
     xlab = "Total 4-Methylsyringol (Grape) (ug/Kg)", ylab = "Smoke Taint Index",pch=19)
dev.off()



raw_wine_data_with_index = lapply(1 : nrow(merged_data_both), function(x){
    cur_wine_type = rownames(merged_data_both)[x]
    idx = which(cur_wine_type == paste0(Wine_Data[,1],"|",Wine_Data[,2]))
    if (length(idx) != 2) stop()
    wine_type = rep(cur_wine_type, length(idx))
    smokeTaintIndex = rep(merged_data_both[x,ncol(merged_data_both)],length(idx))
    return((cbind(wine_type, Wine_Data[idx,Feature_Col_Idx,drop=F],smokeTaintIndex)))
})

raw_wine_data_with_index = do.call(rbind, raw_wine_data_with_index)


svg("output/free_m-cresol_wine.svg",width=6,height=6)
plot(raw_wine_data_with_index$`free_m-cresol`, 
     raw_wine_data_with_index$smokeTaintIndex,
     main = paste0("Free m-Cresol (Wine) VS Smoke Taint Index\n(N = 55 (Wines) * 2 (replicates)) (PCC = ", format(cor_smoke_idx["free_m-cresol"],digits=2),")"),
     xlab = "Free m-Cresol (Wine) (ug/L)", ylab = "Smoke Taint Index",pch=19)
dev.off()

svg("output/free_o-cresol_wine.svg",width=6,height=6)
plot(raw_wine_data_with_index$`free_o-cresol`, 
     raw_wine_data_with_index$smokeTaintIndex,
     main = paste0("Free o-Cresol (Wine) VS Smoke Taint Index\n(N = 55 (Wines) * 2 (replicates)) (PCC = ", format(cor_smoke_idx["free_o-cresol"],digits=2),")"),
     xlab = "Free o-Cresol (Wine) (ug/L)", ylab = "Smoke Taint Index",pch=19)
dev.off()

svg("output/total_syringol_wine.svg",width=6,height=6)
plot(raw_wine_data_with_index$`total_syringol`, 
     raw_wine_data_with_index$smokeTaintIndex,
     main = paste0("Total Syringol (Wine) VS Smoke Taint Index\n(N = 55 (Wines) * 2 (replicates)) (PCC = ", format(cor_smoke_idx["total_syringol"],digits=2),")"),
     xlab = "Total Syringol (Wine) (ug/L)", ylab = "Smoke Taint Index",pch=19)
dev.off()

svg("output/total_o-cresol_wine.svg",width=6,height=6)
plot(raw_wine_data_with_index$`total_o-cresol`, 
     raw_wine_data_with_index$smokeTaintIndex,
     main = paste0("Total o-Cresol (Wine) VS Smoke Taint Index\n(N = 55 (Wines) * 2 (replicates)) (PCC = ", format(cor_smoke_idx["total_o-cresol"],digits=2),")"),
     xlab = "Total o-Cresol (Wine) (ug/L)", ylab = "Smoke Taint Index",pch=19)
dev.off()

svg("output/free_phenol_wine.svg",width=6,height=6)
plot(raw_wine_data_with_index$`total_m-cresol`, 
     raw_wine_data_with_index$smokeTaintIndex,
     main = paste0("Free Phenol (Wine) VS Smoke Taint Index\n(N = 55 (Wines) * 2 (replicates)) (PCC = ", format(cor_smoke_idx["free_phenol"],digits=2),")"),
     xlab = "Free Phenol (Wine) (ug/L)", ylab = "Smoke Taint Index",pch=19)
dev.off()

svg("output/total_m-cresol_wine.svg",width=6,height=6)
plot(raw_wine_data_with_index$`total_m-cresol`, 
     raw_wine_data_with_index$smokeTaintIndex,
     main = paste0("Total m-Cresol (Wine) VS Smoke Taint Index\n(N = 55 (Wines) * 2 (replicates)) (PCC = ", format(cor_smoke_idx["total_m-cresol"],digits=2),")"),
     xlab = "Total m-Cresol (Wine) (ug/L)", ylab = "Smoke Taint Index",pch=19)
dev.off()





#Linear regression

features = merged_data_both[,1:(ncol(merged_data_both)-1)]
scores = merged_data_both[,ncol(merged_data_both)]
linear_model = lm(scores~features)

lm_predict = predict(linear_model)
cor(lm_predict,scores)

#Lasso regression
lambdas = 10^seq(3,-3,-0.0001)

source("FunctionSet.r")
N_Trial = 10
lasso_result_list = list()
lasso_result_grape_list = list()
svr_result_list = list()
svr_result_grape_list = list()
rf_result_list = list()
rf_result_grape_list = list()

log_lasso_result_list = list()
log_lasso_result_grape_list = list()
log_svr_result_list = list()
log_svr_result_grape_list = list()
log_rf_result_list = list()
log_rf_result_grape_list = list()

for (i in 1 : N_Trial){
    print(i)
    lasso_result_list[[i]] = run_lasso(matrix(scores,ncol=1),features) 
    lasso_result_grape_list[[i]] = run_lasso(matrix(scores,ncol=1),features[,1:20]) 
    svr_result_list[[i]] = run_svr(matrix(scores,ncol=1),features)
    svr_result_grape_list[[i]] = run_svr(matrix(scores,ncol=1),features[,1:20])
    rf_result_list[[i]] = run_rf(matrix(scores,ncol=1),features)
    rf_result_grape_list[[i]] = run_rf(matrix(scores,ncol=1),features[,1:20])
    
    log_lasso_result_list[[i]] = run_lasso(matrix(scores,ncol=1),log(features+1)) 
    log_lasso_result_grape_list[[i]] = run_lasso(matrix(scores,ncol=1),log(1+features[,1:20])) 
    log_svr_result_list[[i]] = run_svr(matrix(scores,ncol=1),log(1+features))
    log_svr_result_grape_list[[i]] = run_svr(matrix(scores,ncol=1),log(1+features[,1:20]))
    log_rf_result_list[[i]] = run_rf(matrix(scores,ncol=1),log(1+features))
    log_rf_result_grape_list[[i]] = run_rf(matrix(scores,ncol=1),log(1+features[,1:20]))
}

pcc_lasso = sapply(lasso_result_list,function(x){x$lasso_performance[1]})
pcc_lasso_grape = sapply(lasso_result_grape_list,function(x){x$lasso_performance[1]})
pcc_svr = sapply(svr_result_list,function(x){x$svr_performance[1]})
pcc_svr_grape = sapply(svr_result_grape_list,function(x){x$svr_performance[1]})
pcc_rf = sapply(rf_result_list,function(x){x$rf_performance[1]})
pcc_rf_grape = sapply(rf_result_grape_list,function(x){x$rf_performance[1]})

pcc_log_lasso = sapply(log_lasso_result_list,function(x){x$lasso_performance[1]})
pcc_log_lasso_grape = sapply(log_lasso_result_grape_list,function(x){x$lasso_performance[1]})
pcc_log_svr = sapply(log_svr_result_list,function(x){x$svr_performance[1]})
pcc_log_svr_grape = sapply(log_svr_result_grape_list,function(x){x$svr_performance[1]})
pcc_log_rf = sapply(log_rf_result_list,function(x){x$rf_performance[1]})
pcc_log_rf_grape = sapply(log_rf_result_grape_list,function(x){x$rf_performance[1]})

get_spearman = function(input){
    raw_cor_val = sapply(input,function(x){cor(x$pred_Y,x$truth_Y,method="spearman")})
    return(raw_cor_val)
}


results_pcc = list(pcc_lasso = pcc_lasso,
                   pcc_lasso_grape = pcc_lasso_grape,
                   pcc_svr = pcc_svr,
                   pcc_svr_grape = pcc_svr_grape,
                   pcc_rf = pcc_rf,
                   pcc_rf_grape = pcc_rf_grape,
                   pcc_log_lasso = pcc_log_lasso,
                   pcc_log_lasso_grape = pcc_log_lasso_grape,
                   pcc_log_svr = pcc_log_svr,
                   pcc_log_svr_grape = pcc_log_svr_grape,
                   pcc_log_rf = pcc_log_rf,
                   pcc_log_rf_grape = pcc_log_rf_grape)

heatmap_matrix = matrix(0,nrow=2,ncol=6)
heatmap_matrix[1,1] = mean(results_pcc$pcc_lasso)
heatmap_matrix[1,2] = mean(results_pcc$pcc_svr)
heatmap_matrix[1,3] = mean(results_pcc$pcc_rf)

heatmap_matrix[1,4] = mean(results_pcc$pcc_log_lasso)
heatmap_matrix[1,5] = mean(results_pcc$pcc_log_svr)
heatmap_matrix[1,6] = mean(results_pcc$pcc_log_rf)

heatmap_matrix[2,1] = mean(results_pcc$pcc_lasso_grape)
heatmap_matrix[2,2] = mean(results_pcc$pcc_svr_grape)
heatmap_matrix[2,3] = mean(results_pcc$pcc_rf_grape)

heatmap_matrix[2,4] = mean(results_pcc$pcc_log_lasso_grape)
heatmap_matrix[2,5] = mean(results_pcc$pcc_log_svr_grape)
heatmap_matrix[2,6] = mean(results_pcc$pcc_log_rf_grape)

method_name = c("Lasso","SVR","RF","Lasso","SVR","RF")
feature_selection_method = c("(Use All Features)","(Grape Features Only)")
svg("output/model_performance.svg",width=8,height=1.5)
h = Heatmap(heatmap_matrix,
        #column_title = "Model Performance (PCC) (5-Fold CV, 10 times)",
        col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
        column_split = c(rep("(Conc)",3),rep("(Log(Conc+1))",3)),
        row_split = c(1,2),
        show_row_dend = F,
        show_column_dend = F,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        column_title_gp = gpar(fontsize=12),
        column_names_gp = gpar(fontsize=10),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(formatC(heatmap_matrix[i,j], digits=2, format="f"), x, y)
        }, 
        top_annotation = HeatmapAnnotation(
            text = anno_text(method_name, rot = 0, just = "center"),
            annotation_height = max_text_height(method_name)
        ), 
        left_annotation = rowAnnotation(
            text = anno_text(feature_selection_method, rot = 0, just = "left"),
            annotation_width = max_text_width(feature_selection_method)
        ),
        row_title = NULL,
        column_title_side = "bottom",
        show_heatmap_legend = FALSE)
draw(h)
dev.off()


spearman_lasso = get_spearman(lasso_result_list)
spearman_lasso_grape = get_spearman(lasso_result_grape_list)
spearman_svr = get_spearman(svr_result_list)
spearman_svr_grape = get_spearman(svr_result_grape_list)
spearman_rf = get_spearman(rf_result_list)
spearman_rf_grape = get_spearman(rf_result_grape_list)

spearman_log_lasso = get_spearman(log_lasso_result_list)
spearman_log_lasso_grape = get_spearman(log_lasso_result_grape_list)
spearman_log_svr = get_spearman(log_svr_result_list)
spearman_log_svr_grape = get_spearman(log_svr_result_grape_list)
spearman_log_rf = get_spearman(log_rf_result_list)
spearman_log_rf_grape = get_spearman(log_rf_result_grape_list)

results_spearman = list(spearman_lasso = spearman_lasso,
                   spearman_lasso_grape = spearman_lasso_grape,
                   spearman_svr = spearman_svr,
                   spearman_svr_grape = spearman_svr_grape,
                   spearman_rf = spearman_rf,
                   spearman_rf_grape = spearman_rf_grape,
                   spearman_log_lasso = spearman_log_lasso,
                   spearman_log_lasso_grape = spearman_log_lasso_grape,
                   spearman_log_svr = spearman_log_svr,
                   spearman_log_svr_grape = spearman_log_svr_grape,
                   spearman_log_rf = spearman_log_rf,
                   spearman_log_rf_grape = spearman_log_rf_grape)

heatmap_matrix = matrix(0,nrow=2,ncol=6)
heatmap_matrix[1,1] = mean(results_spearman$spearman_lasso)
heatmap_matrix[1,2] = mean(results_spearman$spearman_svr)
heatmap_matrix[1,3] = mean(results_spearman$spearman_rf)

heatmap_matrix[1,4] = mean(results_spearman$spearman_log_lasso)
heatmap_matrix[1,5] = mean(results_spearman$spearman_log_svr)
heatmap_matrix[1,6] = mean(results_spearman$spearman_log_rf)

heatmap_matrix[2,1] = mean(results_spearman$spearman_lasso_grape)
heatmap_matrix[2,2] = mean(results_spearman$spearman_svr_grape)
heatmap_matrix[2,3] = mean(results_spearman$spearman_rf_grape)

heatmap_matrix[2,4] = mean(results_spearman$spearman_log_lasso_grape)
heatmap_matrix[2,5] = mean(results_spearman$spearman_log_svr_grape)
heatmap_matrix[2,6] = mean(results_spearman$spearman_log_rf_grape)

method_name = c("Lasso","SVR","RF","Lasso","SVR","RF")
feature_selection_method = c("(Use All Features)","(Grape Features Only)")
svg("output/model_performance_spearman.svg",width=8,height=1.5)
h = Heatmap(heatmap_matrix,
            #column_title = "Model Performance (Spearman) (5-Fold CV, 10 times)",
            col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
            column_split = c(rep("(Conc)",3),rep("(Log(Conc+1))",3)),
            row_split = c(1,2),
            show_row_dend = F,
            show_column_dend = F,
            cluster_rows = F,
            cluster_columns = F,
            show_column_names = F,
            show_row_names = F,
            column_title_gp = gpar(fontsize=12),
            column_names_gp = gpar(fontsize=10),
            cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(formatC(heatmap_matrix[i,j], digits=2, format="f"), x, y)
            }, 
            top_annotation = HeatmapAnnotation(
                text = anno_text(method_name, rot = 0, just = "center"),
                annotation_height = max_text_height(method_name)
            ), 
            left_annotation = rowAnnotation(
                text = anno_text(feature_selection_method, rot = 0, just = "left"),
                annotation_width = max_text_width(feature_selection_method)
            ),
            row_title = NULL,
            column_title_side = "bottom",
            show_heatmap_legend = FALSE)
draw(h)
dev.off()




plot_cv_result("output/lasso_regression_results.svg",width=5,height=5,lasso_result_list[[1]],"Lasso")
plot_cv_result("output/lasso_regression_results_grape_feature_only.svg",width=5,height=5,lasso_result_grape_list[[1]],"Lasso (Grape Features Only)")
plot_cv_result("output/svr_regression_results.svg",width=5,height=5,svr_result_list[[1]],"SVR (Radial)")
plot_cv_result("output/svr_regression_results_grape_feature_only.svg",width=5,height=5,svr_result_grape_list[[1]],"SVR (Radial) (Grape Features Only)")
plot_cv_result("output/rf_regression_results.svg",width=5,height=5,rf_result_list[[1]],"RF")
plot_cv_result("output/rf_regression_results_grape_feature_only.svg",width=5,height=5,rf_result_grape_list[[1]],"RF (Grape Features Only)")

plot_cv_result("output/log_lasso_regression_results.svg",width=5,height=5,log_lasso_result_list[[1]],"Lasso (Log)")
plot_cv_result("output/log_lasso_regression_results_grape_feature_only.svg",width=5,height=5,log_lasso_result_grape_list[[1]],"Lasso (Grape Features Only) (Log)")
plot_cv_result("output/log_svr_regression_results.svg",width=5,height=5,log_svr_result_list[[1]],"SVR (Radial) (Log)")
plot_cv_result("output/log_svr_regression_results_grape_feature_only.svg",width=5,height=5,log_svr_result_grape_list[[1]],"SVR (Radial) (Grape Features Only) (Log)")
plot_cv_result("output/log_rf_regression_results.svg",width=5,height=5,log_rf_result_list[[1]],"RF (Log)")
plot_cv_result("output/log_rf_regression_results_grape_feature_only.svg",width=5,height=5,log_rf_result_grape_list[[1]],"RF (Grape Features Only) (Log)")


plot_lasso_loading("output/lasso_coef_log.svg",log_lasso_result_list[[1]]$coef_lasso,0,12,2)
plot_lasso_loading("output/lasso_coef_log_grape.svg",log_lasso_result_grape_list[[1]]$coef_lasso,0,12,2)

save.image("ModelResult.rdata")






#Saturation model
lambdas = 10^seq(1,-2,-0.01)

permuted_dataset = merged_data_both[sample(nrow(merged_data_both)),]
n_cv = 5

mse_performance = rep(0,length(lambdas))
for (h in 1 : length(lambdas)){
    cur_lambda = lambdas[h]
    print(cur_lambda)
    predicted_scores = c()
    truth = c()
    for (i in 1 : n_cv){
        idx = which(1 : nrow(permuted_dataset) %% n_cv == (i-1))
        
        features = permuted_dataset[-idx,1:(ncol(permuted_dataset)-1)]
        scores = permuted_dataset[-idx,ncol(permuted_dataset)]
        
        features_test = permuted_dataset[idx,1:(ncol(permuted_dataset)-1)]
        scores_test = permuted_dataset[idx,ncol(permuted_dataset)]
        
        converted_features = NULL
        converted_features_test = NULL
        for (j in 1 : ncol(features)){
            df = data.frame(v = scores, S = features[,j])
            df_test = data.frame(v = scores_test, S = features_test[,j])
            model.drm = NULL
            try(model.drm <- drm (v ~ S, data = df, fct = MM.2()),silent = T)
            
            if (!is.null(model.drm)){
                converted_features = cbind(converted_features, predict(model.drm, newdata = df))
                converted_features_test =cbind(converted_features_test, predict(model.drm, newdata = df_test))
            }
        }
        
        fit = glmnet(converted_features,scores,lambda = cur_lambda)
        y_predicted <- predict(fit, s = cur_lambda, newx = converted_features_test)
        
        predicted_scores = c(predicted_scores, y_predicted)
        truth = c(truth, scores_test)
    }
    
    mse = mean((predicted_scores - truth)^2)
    mse_performance[h] = mse
    
    print(paste0("mse = ", mse))
}

converted_features = NULL
converted_features_test = NULL
features = merged_data_both[,1:(ncol(merged_data_both)-1)]
scores = merged_data_both[,ncol(merged_data_both)]
for (j in 1 : ncol(features)){
    df = data.frame(v = scores, S = features[,j])
    model.drm = NULL
    try(model.drm <- drm (v ~ S, data = df, fct = MM.2()),silent = T)
    
    if (!is.null(model.drm)){
        converted_features = cbind(converted_features, predict(model.drm, newdata = df))
    }
}
best_lambda = lambdas[which.min(mse_performance)]
fit = glmnet(converted_features,scores,lambda = best_lambda)
y_predicted <- predict(fit, s = best_lambda, newx = converted_features)

plot(y_predicted, scores)









#Case Study
idx_RMI_R61 = which(startsWith(as.character(raw_wine_data_with_index[,1]), "20CS_RMI_R61"))
wine_data_with_index_RMI_R61 = raw_wine_data_with_index[idx_RMI_R61,]
wine_data_with_index_RMI_R61[,c("wine_type","free_o-cresol","free_m-cresol","total_o-cresol","total_syringol","smokeTaintIndex")]

idx_RMI_R61 = which(startsWith(as.character(raw_grape_data_with_index[,1]), "20CS_RMI_R61"))
grape_data_with_index_RMI_R61 = raw_grape_data_with_index[idx_RMI_R61,]
grape_data_with_index_RMI_R61[,c("grape_type","free_guaiacol","free_4-ethylphenol","total_o-cresol","total_syringol","smokeTaintIndex")]


idx_SB_HB = which(startsWith(as.character(raw_wine_data_with_index[,1]), "20SB_HB"))
wine_data_with_index_SB_HB = raw_wine_data_with_index[idx_SB_HB,]
wine_data_with_index_SB_HB[,c("wine_type","free_o-cresol","free_m-cresol","total_o-cresol","total_syringol","smokeTaintIndex")]

idx_SB_HB = which(startsWith(as.character(raw_grape_data_with_index[,1]), "20SB_HB"))
grape_data_with_index_SB_HB = raw_grape_data_with_index[idx_SB_HB,]
grape_data_with_index_SB_HB[,c("grape_type","free_guaiacol","free_4-ethylphenol","total_o-cresol","total_syringol","smokeTaintIndex")]



idx_RMI_R24_25_58 = which(startsWith(as.character(raw_wine_data_with_index[,1]), "20CS_RMI_R24 25 58"))
wine_data_with_index_RMI_R24_25_58 = raw_wine_data_with_index[idx_RMI_R24_25_58,]
wine_data_with_index_RMI_R24_25_58[,c("wine_type","free_o-cresol","free_m-cresol","total_o-cresol","total_syringol","smokeTaintIndex")]

idx_RMI_R24_25_58 = which(startsWith(as.character(raw_grape_data_with_index[,1]), "20CS_RMI_R24 25 58"))
grape_data_with_index_RMI_R24_25_58 = raw_grape_data_with_index[idx_RMI_R24_25_58,]
grape_data_with_index_RMI_R24_25_58[,c("grape_type","free_guaiacol","free_4-ethylphenol","total_o-cresol","total_syringol","smokeTaintIndex")]


idx_TO_Tyree = which(startsWith(as.character(raw_wine_data_with_index[,1]), "20TO_Tyree"))
wine_data_with_index_TO_Tyree = raw_wine_data_with_index[idx_TO_Tyree,]
wine_data_with_index_TO_Tyree[,c("wine_type","free_o-cresol","free_m-cresol","total_o-cresol","total_syringol","smokeTaintIndex")]

idx_TO_Tyree = which(startsWith(as.character(raw_grape_data_with_index[,1]), "20TO_Tyree"))
grape_data_with_index_TO_Tyree = raw_grape_data_with_index[idx_TO_Tyree,]
grape_data_with_index_TO_Tyree[,c("grape_type","free_guaiacol","free_4-ethylphenol","total_o-cresol","total_syringol","smokeTaintIndex")]


#Cor (between grape feature and wine feature)
cor_matrix = cor(merged_data_both[,1:20],merged_data_both[,21:40])
cor_matrix = t(cor_matrix)
colnames(cor_matrix) = gsub("(grape)","",colnames(cor_matrix),fixed=T)
colnames(cor_matrix) = gsub("_"," ",colnames(cor_matrix),fixed=T)
rownames(cor_matrix) = gsub("_"," ",rownames(cor_matrix),fixed=T)
svg("output/cor_grape_wine_feature.svg",width=10,height=8)
h = Heatmap(cor_matrix,
            column_title = "(Grape VOCs)",
            row_title = "(Wine VOCs)",
            col = colorRamp2(c(-0.5, 0.5, 1.5), c("blue", "white", "red")),
            show_row_dend = F,
            show_column_dend = F,
            cluster_rows = F,
            cluster_columns = F,
            show_column_names = T,
            show_row_names = T,
            cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(formatC(cor_matrix[i,j], digits=2, format="f"), x, y)
            }, 
            column_title_side = "bottom",
            show_heatmap_legend = FALSE)
draw(h)
dev.off()

svg("output/grape_vs_wine_free_phenol.svg",width=5,height=5)
df = data.frame(x=merged_data_both[,"(grape)free_phenol"],y=merged_data_both[,"free_phenol"])
h = ggplot(df, aes(x=x,y=y)) +
    geom_point(size=3) + 
    theme_classic() + 
    ggtitle(paste0("Free Phenol Concentration\n(Grape VOCs VS Wine VOCs)(PCC = ",round(cor(df$x,df$y)*100)/100,")")) + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    labs(y="Wine VOC concetration (ug/L)",x="Grape VOC concentration (ug/Kg)")
print(h)
dev.off()

svg("output/grape_vs_wine_free_guaiacol_o-cresol.svg",width=5,height=5)
df = data.frame(x=merged_data_both[,"(grape)free_guaiacol"],y=merged_data_both[,"free_o-cresol"])
h = ggplot(df, aes(x=x,y=y)) +
    geom_point(size=3) + 
    theme_classic() + 
    ggtitle(paste0("Concentration\n(Grape VOCs VS Wine VOCs)(PCC = ",round(cor(df$x,df$y)*100)/100,")")) + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    labs(y="Free o-cresol (Wine) concetration (ug/L)",x="Free guaiacol (Grape) concentration (ug/Kg)")
print(h)
dev.off()