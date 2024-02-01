library(randomForest)
library(Boruta)
library(caret)

source('module/modelRFcross.R')
source('module/indexcross.R')
source('module/varBorutacross.R')

# load genotype and phenotype markers of ASD group
df.asd = read.csv2('data_asd.csv') 
# load SNP markers of ASD and control groups (positive class: ASD)
df.snp.control = read.csv2('data_SNP_control.csv')

### 1a. Feature selection: ASD vs control groups
# feature selection stability
set.seed(123)
data = data.frame(class = c(rep(1,nrow(df.asd)), rep(0,nrow(df.snp.control))), rbind(df.asd[,2:11], df.snp.control[,c(3:12)]))
#data = data.frame(class = df.snp.control$sex, df.snp.control[,!colnames(df.snp.control) %in% c("sex", "id")])
list.index.cross = indexcross(y = data$class,folds = 3,iterations = 30, stratified = TRUE)
var.imp = varimpBoruta(data = data, folds = 3, niter = 30, list.index.cross, maxRuns = 300)

# SNP frequency in 90 feature subsets
SNP.count = sort(table(unlist(var.imp)),decreasing = T)
print(SNP.count)
# rs10099100  rs7521492 rs61847307  rs6803008  rs7122181   rs880446  rs9879311  rs1353545 
# 90          6          4          4          4          2          2          1 

# importance of SNP markers
result.boruta = Boruta(data[,2:11],data$class,doTrace=2)
print(result.boruta)
# Boruta performed 24 iterations in 1.633948 secs.
# 1 attributes confirmed important: rs10099100;
# 9 attributes confirmed unimportant: rs1353545, rs2828478, rs4904167, rs61847307, rs6803008 and 4 more;

plot(result.boruta ,sort=TRUE,las=2, cex.lab =1,main='ASD vs control',cex=1,
     xlab = "",
     cex.axis=0.9,
     ylab = "Importance")

### 1b. Feature selection: female ASD vs female control groups
df.asm.female = df.asd[df.asd$sex == '1', 2:11]
df.control.female = df.snp.control[df.snp.control$sex == '1', 3:12]  
data = data.frame(class = c(rep(1,nrow(df.asm.female)), rep(0,nrow(df.control.female))),
                  rbind(df.asm.female, df.control.female))

list.index.cross = indexcross(y = data$class,folds = 3,iterations = 30, stratified = TRUE)
var.imp = varimpBoruta(data = data, folds = 3, niter = 30, list.index.cross, maxRuns = 300)

# SNP frequency in 90 feature subsets
print(sort(table(unlist(var.imp)),decreasing = T))
# rs61847307 rs10099100  rs6803008  rs2828478  rs4904167  rs9879311  rs7122181 
# 34         15          3          2          2          2          1 

# importance of SNP markers
result.boruta = Boruta(data[,2:11],data$class,doTrace=2)
print(result.boruta)
# Boruta performed 99 iterations in 2.336752 secs.
# 1 attributes confirmed important: rs61847307;
# 8 attributes confirmed unimportant: rs1353545, rs2828478, rs4904167, rs6803008,
# rs7122181 and 3 more;
# 1 tentative attributes left: rs10099100;

### 1c. Feature selection: male ASD vs male control groups
df.asm.male = df.asd[df.asd$sex == '2', 2:11]
df.control.male = df.snp.control[df.snp.control$sex == '2', 3:12]  
data = data.frame(class = c(rep(1,nrow(df.asm.male)), rep(0,nrow(df.control.male))),
                  rbind(df.asm.male, df.control.male))

list.index.cross = indexcross(y = data$class,folds = 3,iterations = 30, stratified = TRUE)
var.imp = varimpBoruta(data = data, folds = 3, niter = 30, list.index.cross, maxRuns = 300)

# SNP frequency in 90 feature subsets
print(sort(table(unlist(var.imp)),decreasing = T))
# rs10099100  rs7521492  rs1353545   rs880446 
# 90          2          1          1 

# importance of SNP markers
result.boruta = Boruta(data[,2:11],data$class,doTrace=2)
print(result.boruta)
# Boruta performed 10 iterations in 0.4079111 secs.
# 1 attributes confirmed important: rs10099100;
# 9 attributes confirmed unimportant: rs1353545, rs2828478, rs4904167, rs61847307,
# rs6803008 and 4 more;

### 2. Feature selection and machine learning models: male ASD vs female ASD groups
set.seed(123) 
df.asd$sex <- replace(df.asd$sex, df.asd$sex == 2, 0)  
data = data.frame(class = df.asd$sex, df.asd[,-c(1,12)])
list.index.cross = indexcross(y = data$class,folds = 3,iterations = 30, stratified = TRUE)
var.imp = varimpBoruta(data = data, folds = 3, niter = 30, list.index.cross, maxRuns = 300)

# the most stable markers in 90 feature subsets
table.markers = sort(table(unlist(var.imp)),decreasing = T)
N.most.stable = table.markers[table.markers >= 9]
# abnormal_ctg       birth_head_circum            birth_length               scoliosis 
# 62                      55                      37                      33 
# rs9879311              birth_mass          smoking_mother social_difficulties_fam 
# 28                      25                      24                      23 
# brother           diagnosis_age         lip_dysmorphism              rs10099100 
# 19                      17                      13                      13 
# miscarriages               therapies 
# 12                      10 

# importance of the top most stable markers
result.boruta = Boruta(data[,names(N.most.stable)],data$class,doTrace=2)
print(result.boruta)

# Boruta performed 99 iterations in 2.762905 secs.
# 5 attributes confirmed important: abnormal_ctg, birth_head_circum, birth_length,
# rs9879311, smoking_mother;
# 6 attributes confirmed unimportant: brother, diagnosis_age, lip_dysmorphism,
# rs10099100, social_difficulties_fam and 1 more;
# 3 tentative attributes left: birth_mass, miscarriages, scoliosis;

rank.markers = data.frame(var.imp = names(result.boruta$finalDecision), importance = result.boruta$finalDecision)
Ntop.markers = rank.markers[rank.markers$importance != 'Rejected',]
result = modelRFcross(data, var.imp = Ntop.markers,
                      niter = 30,
                      folds =3,
                      list.index.cross)
print(result)
# OOB.mean 0.2540
# OOB.sd   0.0043
# mcc.mean 0.2147
# mcc.sd   0.0140
# auc.mean 0.6430
# auc.sd   0.0087
# f1.mean  0.7472
# f1.sd    0.0130