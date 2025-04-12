# -------------------------------GBLUP_cross------------------------------------
# 用于检验参考群体估计准确性，参考群内所有个体的基因型和表型都是已知的
library(BGLR)
library(plyr)
library(data.table)
install.packages("BGLR")
# 计算程序的运行时间
t1<-proc.time()

# 加载数据集查看数据维度
geno <- fread("./430geno_10w_GBLUP_BBH.csv")
ind_ID <- geno[,1]
geno <- as.matrix(geno[,-1])
phen <- fread("./430pheno_GBLUP.csv")
fix <- fread("./430fix_GBLUP.csv")
true_ebv <- fread("./430ebv_GBLUP.csv")
seed <- sample(1:10000, 1)  # 从1到10000中随机选择种子

X <- geno
Y <- phen$t1  # 有多个表型时可以直接更改列名进行选择
fix <- as.matrix(fix[,-1])

# 计算G矩阵
# G<-tcrossprod(X)
# G <- G/mean(diag(G))
# G<-tcrossprod(X)/ncol(X)

#----------------------------------参数设置-------------------------------------

# 交叉检验设置
folds<-5  # 折数
y<-Y
n <- length(y)
set.seed(seed)
sets<-rep(1:5,round(length(y)/folds,digits = 0))
sets<-sets[order(runif(nrow(X)))]

# 定义评价指标(可能用不上，但还是算一下)
COR.CV<-rep(NA,times=(folds+1))
MSE.CV<-rep(NA,times=(folds+1))
MAE.CV<-rep(NA,times=(folds+1))
RMSE.CV<-rep(NA,times=(folds+1))
yHatCV<-numeric() #定义一个变量来存储交叉验证结果

#-------------------------------BRR/BL/BayesA/BayesB/BayesC--------------------------------
result_list <- list()  # 创建一个空的列表来存储每一折的结果

for (fold in 1:folds) {
  yNa <- y
  fixNa <- fix
  whichNa <- which(sets == fold)
  yNa[whichNa] <- NA
  fixNa[whichNa] <- NA
  # 将个体划分为有表型的和无表型的
  has_phenotype <- !is.na(yNa)
  no_phenotype <- is.na(yNa)
  
  # 提取有表型个体和无表型个体的 ID
  ID_obs <- ind_ID[has_phenotype]
  ID_pred <- ind_ID[no_phenotype]
  
  # 分别提取有表型个体和无表型个体的固定效应矩阵和随机效应矩阵和表型
  y_obs <- y[has_phenotype] # 有表型个体的表型
  X_obs <- fix[has_phenotype, , drop = FALSE] # 有表型个体的固定效应
  X_pred <- fix[no_phenotype, , drop = FALSE] # 无表型个体的固定效应
  X_obs_geno <- geno[has_phenotype, , drop = FALSE] # 有表型个体的标记矩阵
  X_true_ebv <- true_ebv[no_phenotype, , drop = FALSE] # 无表型个体的实际育种值
  
  # 更新有表型个体的数量
  n_ind_obs <- nrow(X_obs)
  n_ind_pred <- nrow(X_pred)
  
  # 假设残差方差和随机效应方差
  sigma2_e <- 22.647  # 残差方差
  sigma2_u <- 13.236  # 基因效应方差
  
  # 1. 计算有表型个体的基因组关系矩阵 G_obs
  # 将012格式转为-101格式，便于计算
  M <- geno - 1
  # 构建P矩阵
  p_lower <- (apply(M,2,sum)+nrow(M))/(nrow(M)*2)
  p_upper <- 2*(p_lower-0.5)
  p_matrix <- matrix(p_upper,byrow=T,nrow=nrow(M),ncol=ncol(M))
  P <- p_matrix
  # 构建Z矩阵
  Z <-  M - P
  # 构建G矩阵
  d <- 2*sum(p_lower*(1 - p_lower))
  G_pred_obs <- Z %*% t(Z) / d
  # 获取已知个体的关系矩阵
  true_index <- which(has_phenotype == TRUE) # 列索引，有表型的个体
  false_index <- which(has_phenotype == FALSE) # 行索引，无表型的个体
  G_obs <- G_pred_obs[true_index, true_index]
  
  # 添加正则化常数，避免G_obs是奇异矩阵
  lambda <- 1e-5
  G_obs_inv <- solve(G_obs + lambda * diag(nrow(G_obs)))  # 求解G_obs的逆
  
  # 2. 构建有表型个体的化简后的混合模型方程 (MME)
  Z_obs <- diag(n_ind_obs)  # Z 是单位矩阵
  
  # MME 方程左侧矩阵 (无需包含 R^-1，直接使用 X'X 和 Z'Z + sigma2_e * G^-1)
  MME_left_obs <- rbind(
    cbind(t(X_obs) %*% X_obs, t(X_obs) %*% Z_obs),
    cbind(t(Z_obs) %*% X_obs, t(Z_obs) %*% Z_obs + sigma2_e * G_obs_inv)
  )
  
  # MME 方程右侧向量 (也无需包含 R^-1，直接使用 X'y 和 Z'y)
  MME_right_obs <- rbind(
    t(X_obs) %*% y_obs,
    t(Z_obs) %*% y_obs
  )
  
  # 3. 求解有表型个体的混合模型方程
  solution_obs <- solve(MME_left_obs) %*% MME_right_obs
  
  # 提取有表型个体的固定效应和随机效应估计值
  beta_hat_obs <- solution_obs[1:ncol(X_obs), ]  # 固定效应估计值
  u_hat_obs <- solution_obs[(ncol(X_obs) + 1):nrow(solution_obs), ]  # 随机效应估计值（有表型个体的基因效应）
  
  # 4. 基于有表型个体估计值，预测无表型个体的基因效应
  # 计算无表型个体与有表型个体之间的 G 矩阵
  # 将012格式转为-101格式，便于计算
  M <- geno - 1
  # 构建P矩阵
  p_lower <- (apply(M,2,sum)+nrow(M))/(nrow(M)*2)
  p_upper <- 2*(p_lower-0.5)
  p_matrix <- matrix(p_upper,byrow=T,nrow=nrow(M),ncol=ncol(M))
  P <- p_matrix
  # 构建Z矩阵
  Z <-  M - P
  # 构建G矩阵
  d <- 2*sum(p_lower*(1 - p_lower))
  G_pred_obs <- Z %*% t(Z) / d
  # 获取待预测个体和已知个体的关系矩阵
  true_index <- which(has_phenotype == TRUE) # 列索引，有表型的个体
  false_index <- which(has_phenotype == FALSE) # 行索引，无表型的个体
  G_pred_obs <- G_pred_obs[false_index, true_index]
  
  # 预测无表型个体的基因效应
  u_hat_pred <- G_pred_obs %*% G_obs_inv %*% u_hat_obs
  
  # 预测无表型个体的表型值
  y_hat_pred <- X_pred %*% beta_hat_obs + u_hat_pred
  
  COR.CV[fold] <- cor(X_true_ebv$s1, u_hat_pred)
  
  # 将每一折的结果存储到结果列表中
  result_list[[fold]] <- list(
    COR = COR.CV[fold],
    real_value = X_true_ebv$s1,
    pred_value = u_hat_pred
  )
}

# 将所有结果合并为一个数据框
combined_data <- data.frame()
for (fold in 1:folds) {
  fold_sample_IDs <- ind_ID[which(sets == fold)]
  real_value <- result_list[[fold]]$real_value
  pred_value <- result_list[[fold]]$pred_value
  fold_data <- data.frame(
    Fold = fold,
    SampleID = fold_sample_IDs,
    RealValue = real_value,
    PredValue = pred_value
  )
  combined_data <- rbind(combined_data, fold_data)
}

cor_mean <- mean(COR.CV[1:5])
cor1 <- as.data.frame(COR.CV[1:5])
cor_result <- cbind(cor1, cor_mean, seed)

# 保存结果到文件
write.table(cor_result, "./10w/GBLUP_cross_results.txt", row.names = FALSE, append = T)
write.table(combined_data, "./10w/GBLUP_fold.txt", row.names = FALSE, append = T)

#计算整个程序执行时间
t2 = proc.time()
running_time = t2 - t1
print(paste0('running time:',running_time[3][[1]],'s'))
write.table(running_time, "./10w/GBLUP_cross_time.txt", row.names = FALSE, append = T)

save.image('./10w/GBLUP.RData')
#清除当前空间
rm(list = ls())

