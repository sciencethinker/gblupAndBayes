# -------------------------------BayesA_cross-----------------------------------
# 用于检验参考群体估计准确性，参考群内所有个体的基因型和表型都是已知的
library(BGLR)
library(plyr)
library(data.table)
#
# 计算程序的运行时间
t1<-proc.time()

# 加载数据集查看数据维度
geno <- fread("./430geno_10w_GBLUP_BBH.csv")
sample_IDs <- geno[,1]
geno <- as.matrix(geno[,-1])
phen <- fread("./430pheno_GBLUP.csv")
fix <- fread("./430fix_GBLUP.csv")
true_gebv <- fread("./430ebv_GBLUP.csv")
seed <- sample(1:10000, 1)  # 从1到10000中随机选择种子

X <- geno
Y <- phen$t1  # 有多个表型时可以直接更改列名进行选择
V <- true_gebv$s1

# 计算G矩阵
# G<-tcrossprod(X)
# G <- G/mean(diag(G))
# G<-tcrossprod(X)/ncol(X)

#----------------------------------参数设置-------------------------------------

#参数设置
nIter<-12000 
burnIn<-2000

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

#-------------------------------BRR/BL/BayesA/BayesB/BayesC---------------------
result_list <- list()  # 创建一个空的列表来存储每一折的结果

for (fold in 1:folds) {
  yNa <- y
  fixNa <- fix
  whichNa <- which(sets == fold)
  yNa[whichNa] <- NA
  fixNa[whichNa] <- NA
  ETA <- list(list(~fixNa$c1+fixNa$c2+fix$c3,
                   X = X,
                   model = 'BayesA'))  # (更改参数‘model’来改变方法，包括BRR/BL/BayesA/BayesB/BayesC)
  fm <- BGLR(y = yNa,
             ETA = ETA,
             nIter = nIter,
             burnIn = burnIn,
             saveAt = './10w/BayesA/BayesA_')  # 这里对应的文件夹(BRR)要事先建好，不然会报错，一般跟方法名称统一，找起来方便
  yHatCV[whichNa] <- fm$yHat[fm$whichNa] - fm$mu
  COR.CV[fold] <- cor((fm$yHat[fm$whichNa]- fm$mu), V[whichNa])
  MAE.CV[fold] <- mean(abs(fm$yHat[fm$whichNa] - y[whichNa]))
  MSE.CV[fold] <- mean((fm$yHat[fm$whichNa] - y[whichNa])^2)
  RMSE.CV[fold] <- sqrt(mean((fm$yHat[fm$whichNa] - y[whichNa])^2))
  
  # 将每一折的结果存储到结果列表中
  result_list[[fold]] <- list(
    COR = COR.CV[fold],
    MAE = MAE.CV[fold],
    MSE = MSE.CV[fold],
    RMSE = RMSE.CV[fold],
    real_value = V[whichNa],
    pred_value = fm$yHat[fm$whichNa] - fm$mu
  )
}

# 将所有结果合并为一个数据框
combined_data <- data.frame()
for (fold in 1:folds) {
  fold_sample_IDs <- sample_IDs[which(sets == fold)]
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

# 保存结果到CSV文件
write.table(cor_result, "./10w/BayesA/BayesA_cross_results.txt", row.names = FALSE, append = T)
write.table(combined_data, "./10w/BayesA/BayesA_fold.txt", row.names = FALSE, append = T)

#计算整个程序执行时间
t2 = proc.time()
running_time = t2 - t1
print(paste0('running time:',running_time[3][[1]],'s'))
write.table(running_time, "./10w/BayesA/BayesA_cross_time.txt", row.names = FALSE, append = T)

save.image('./10w/BayesA/BayesA.RData')
#清除当前空间
rm(list = ls())