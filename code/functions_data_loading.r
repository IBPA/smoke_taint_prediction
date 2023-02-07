#Functions for loading data

get_median_grape_data = function(grape_data){
    median_grape_data = apply(grape_data[,colidx_voc],2,function(x){
        sample_name = unique(grape_data[,colidx_sample_name])
        results = numeric(length(sample_name))
        for (i in 1 : length(sample_name)){
            select_idx = grape_data[,colidx_sample_name] == sample_name[i]
            results[i] = median(x[select_idx])
        }
        names(results)= sample_name
        return(results)
    })
    return(median_grape_data)
}

get_mean_wine_data = function(wine_data){
    mean_wine_data = apply(wine_data[,c(colidx_voc,colidx_smoke_taint_idx)],2,function(x){
        sample_name = unique(paste0(wine_data[,colidx_sample_name],"|",wine_data$`Sample type`))
        results = numeric(length(sample_name))
        for (i in 1 : length(sample_name)){
            select_idx = paste0(wine_data[,colidx_sample_name],"|",wine_data$`Sample type`) == sample_name[i]
            results[i] = mean(x[select_idx])
        }
        names(results)= sample_name
        return(results)
    })
    return(mean_wine_data)
}

get_integrated_grape_wine_data = function(median_grape_data, mean_wine_data){
    data_grape_merge_to_wine = lapply(1 : nrow(mean_wine_data),function(x){
        cur_grape_type = strsplit(rownames(mean_wine_data)[x],"|",fixed=T)[[1]][1]
        idx = which(cur_grape_type == rownames(median_grape_data))
        if (length(idx) != 1){
            idx = which(startsWith(cur_grape_type, rownames(median_grape_data)))
        }
        if (length(idx) != 1) stop()
        return(median_grape_data[idx,,drop=F])
    })
    
    integrated_grape_wine_data = NULL
    for (i in 1 : nrow(mean_wine_data)){
        tmp = matrix((c(data_grape_merge_to_wine[[i]],mean_wine_data[i,,drop=F])),nrow=1)
        rownames(tmp) = rownames(mean_wine_data)[i]
        colnames(tmp) = c(paste0("(grape)",colnames(data_grape_merge_to_wine[[1]])),colnames(mean_wine_data))
        integrated_grape_wine_data = rbind(integrated_grape_wine_data,tmp)
    }
    return(integrated_grape_wine_data)
}