# -------------------------------BayesA_ind-------------------------------------
library(plyr)
library(data.table)
library(BGLR)

# 计算程序的运行时间
t1<-proc.time()

# 读入各类文件，删除行名
geno <- fread("./485geno_10w_GBLUP_BBH.csv")
ind_ID <- geno[,1]
geno <- as.matrix(geno[,-1])
phen <- fread("./485pheno_GBLUP.csv")
fix <- fread("./485fix_GBLUP.csv")
true_ebv <- fread("./485ebv_GBLUP.csv")
seed <- sample(1:10000, 1)  # 从1到10000中随机选择种子

set.seed(seed)
#设置表型和固定效应
y <- phen$t1
fix <- as.matrix(fix[,-1])

# 使用Bayes方法根据有表型个体计算育种值
bayes <- BGLR(y,ETA = list(list(~fix$c1+fix$c2+fix$c3,
                                X=geno,
                                model='BayesA')),
              nIter = 12000,burnIn = 2000,thin = 5,saveAt = "",
              S0 = NULL,R2 = 0.5, weights = NULL,verbose = TRUE,
              rmExistingFiles = TRUE, groups=NULL)

# 根据SNP效应计算无表型个体的GEBV
gebv <- bayes$yHat[bayes$whichNa] - bayes$mu

# 将个体划分为有表型的和无表型的
has_phenotype <- !is.na(y)
no_phenotype <- is.na(y)

# 分别提取有表型个体和无表型个体的固定效应矩阵和随机效应矩阵和表型
X_true_ebv <- true_ebv[no_phenotype, , drop = FALSE] # 无表型个体的实际育种值

# 输出结果
cor <- cor(X_true_ebv$s1, gebv)
cor_result <- cbind(as.data.frame(cor), seed)

#计算整个程序执行时间
t2 = proc.time()
running_time = t2 - t1
print(paste0('running time:',running_time[3][[1]],'s'))
write.table(running_time, "./10w/BayesA/BayesA_ind_time.txt", row.names = FALSE, append = T)
write.table(cor_result, "./10w/BayesA/BayesA_ind_result.txt", row.names = FALSE, append = T)
#清除当前空间
rm(list = ls())