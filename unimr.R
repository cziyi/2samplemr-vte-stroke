library(TwoSampleMR)
library(MRInstruments)
library(MendelianRandomization)
library(MRPRESSO)
library(httr)
library(remotes)
library(ieugwasr)
library(MVMR)
library(readr)
library(ggplot2)
library(dplyr)


#----------cal F and Isq-----------------
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

dat <- H_AIS_to_VTE
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome
BXG = abs(BetaXG) # gene--exposure estimates are positive  
# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger
F = BXG^2/seBetaXG^2
mF = mean(F)
mF
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
Isq_unweighted
Isq_weighted
#------------STROKE AS EXPOSURE----------
allstroke <- extract_instruments("ebi-a-GCST005838")
write.csv(allstroke,"allstroke.csv")
#allstroke_exposure <-extract_instruments("ebi-a-GCST005838", p1=5e-6)
LAS <- extract_instruments("ebi-a-GCST005840",p1= 5e-6)
SVS <- extract_instruments("ebi-a-GCST005841",p1= 5e-6)
write.csv(SVS,"svsexposure.csv")
CES <- extract_instruments("ebi-a-GCST005842",p1=5e-6)
write.csv(CES,"cesexposure.csv")
AIS <- extract_instruments("ebi-a-GCST005843")
write.csv(AIS,"aisexposure.csv")

#---------snp from finngen---------------
finn_vte <- extract_instruments("finn-b-I9_VTE")
write.csv(finn_vte,"finn_vte.csv")
#-----------snp from research------------
data <- read_csv("data.csv")
formatted <- format_data(data)
VTE <- clump_data(formatted)
write.csv(VTE,"vte.csv")
#-----------merge two sets---------------
combined_vte <-read_csv("merged_vte.csv")

#-----------EXTRACT OUTCOME--------------
VTE_AS <- extract_outcome_data(snps = combined_vte$SNP, outcomes ="ebi-a-GCST005838" )
VTE_LAS <- extract_outcome_data(snps = combined_vte$SNP, outcomes = "ebi-a-GCST005840")
VTE_SVS <- extract_outcome_data(snps = combined_vte$SNP, outcomes = "ebi-a-GCST005841")
VTE_CES <- extract_outcome_data(snps = combined_vte$SNP, outcomes = "ebi-a-GCST005842")
VTE_AIS <- extract_outcome_data(snps = combined_vte$SNP, outcomes = "ebi-a-GCST005843")
#-----------FORWARD: VTE TO AS-------------
H_VTE_AS <- harmonise_data(combined_vte, VTE_AS)
result1 <- mr(H_VTE_AS)

ivw_r1<- mr(H_VTE_AS, method_list = "mr_ivw_mre")
odds1<-generate_odds_ratios(result1)

scatter1 <- mr_scatter_plot(result1, H_VTE_AS)
ggsave(filename = "scatter_VTE_AS.png", plot = scatter1[[1]], width = 10, height = 8, dpi = 300)
scatter1

single1<- mr_singlesnp(H_VTE_AS)
duplicates <- duplicated(single1$SNP)
print(single1[duplicates, ])
single1 <- single1[!duplicates, ]
forest_plot1 <- mr_forest_plot(single1)
loo1 <-mr_leaveoneout(H_VTE_AS)
mr_leaveoneout_plot(loo1)

duplicated_SNPs <- duplicated(loo1$SNP)
print(loo1[duplicated_SNPs, ])
loo1 <- loo1[!duplicated_SNPs, ]

funnel1 <- mr_funnel_plot(single1)
ggsave(filename = "funnel_VTE_AS.png", plot = funnel1[[1]], width = 10, height = 8, dpi = 300)
hetero_1 <- mr_heterogeneity(H_VTE_AS)
hetero_1
pleio_1 <- mr_pleiotropy_test(H_VTE_AS)
pleio_1
result <- pleio_1
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

run_mr_presso(H_VTE_AS)
run_mr_presso(H_VTE_AS, NbDistribution = 5000)  # 或更多，根据需要调整

#-----------FORWARD: VTE TO LAS-------------
H_VTE_LAS <- harmonise_data(combined_vte, VTE_LAS)
result2 <- mr(H_VTE_LAS)
write.csv(result2, "forward_VTE_LAS.csv")
odds2<-generate_odds_ratios(result2)

scatter2 <- mr_scatter_plot(result2, H_VTE_LAS)
scatter2
ggsave(filename = "scatter_VTE_LAS.png", plot = scatter2[[1]], width = 10, height = 8, dpi = 300)
single2<- mr_singlesnp(H_VTE_LAS)

single2 <- single2 %>%
  mutate(original_order = row_number())
# 根据 p 值和 b 值选择每个 SNP 的最佳行
single2_filtered <- single2 %>%
  group_by(SNP) %>%
  filter(p == min(p) & b == max(b)) %>%
  ungroup()
# 按原始顺序排序
single2_filtered <- single2_filtered %>%
  arrange(original_order) %>%
  select(-original_order)  # 删除原始顺序索引列

mr_forest_plot(single2_filtered)

loo2 <-mr_leaveoneout(H_VTE_LAS)
loo2 <- loo2 %>%
  mutate(original_order = row_number())

# 根据 p 值和 b 值选择每个 SNP 的最佳行
loo2_filtered <- loo2 %>%
  group_by(SNP) %>%
  filter(p == min(p) & b == max(b)) %>%
  ungroup()
# 按原始顺序排序
loo2_filtered <- loo2_filtered %>%
  arrange(original_order) %>%
  select(-original_order)  # 删除原始顺序索引列
mr_leaveoneout_plot(loo2_filtered)
funnel2 <- mr_funnel_plot(single2_filtered)
ggsave(filename = "funnel_VTE_LAS.png", plot = funnel2[[1]], width = 10, height = 8, dpi = 300)
hetero_2 <- mr_heterogeneity(H_VTE_LAS)
result <- hetero_2
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
pleio_2 <- mr_pleiotropy_test(H_VTE_LAS)
result <- pleio_2
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

run_mr_presso(H_VTE_LAS,NbDistribution = 5000)
#-----------FORWARD: VTE TO SVS-------------
H_VTE_SVS <- harmonise_data(combined_vte, VTE_SVS)
result3 <- mr(H_VTE_SVS)
odds3 <- generate_odds_ratios(result3)
scatter3 <- mr_scatter_plot(result3, H_VTE_SVS)
ggsave(filename = "scatter_VTE_SVS.png", plot = scatter3[[1]], width = 10, height = 8, dpi = 300)
scatter3
single3<- mr_singlesnp(H_VTE_SVS)
single3 <- single3 %>%
  mutate(original_order = row_number())

# 根据 p 值和 b 值选择每个 SNP 的最佳行
single3_filtered <- single3 %>%
  group_by(SNP) %>%
  filter(p == min(p) & b == max(b)) %>%
  ungroup()

# 按原始顺序排序
single3_filtered <- single3_filtered %>%
  arrange(original_order) %>%
  select(-original_order)  # 删除原始顺序索引列

mr_forest_plot(single3_filtered)

loo3 <-mr_leaveoneout(H_VTE_SVS)
loo3 <- loo3 %>%
  mutate(original_order = row_number())

# 根据 p 值和 b 值选择每个 SNP 的最佳行
loo3_filtered <- loo3 %>%
  group_by(SNP) %>%
  filter(p == min(p) & b == max(b)) %>%
  ungroup()

# 按原始顺序排序
loo3_filtered <- loo3_filtered %>%
  arrange(original_order) %>%
  select(-original_order)  # 删除原始顺序索引列
mr_leaveoneout_plot(loo3_filtered)

funnel3 <- mr_funnel_plot(single3)
ggsave(filename = "funnel_VTE_SVS.png", plot = funnel3[[1]], width = 10, height = 8, dpi = 300)

hetero_3 <- mr_heterogeneity(H_VTE_SVS)
hetero_3
result <- hetero_3
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

pleio_3 <- mr_pleiotropy_test(H_VTE_SVS)
result <- pleio_3
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

run_mr_presso(H_VTE_SVS,NbDistribution = 5000)
#-----------FORWARD: VTE TO CES-------------
H_VTE_CES <- harmonise_data(combined_vte, VTE_CES)
result4 <- mr(H_VTE_CES)
odds4 <- generate_odds_ratios(result4)

scatter4 <- mr_scatter_plot(result4, H_VTE_CES)
ggsave(filename = "scatter_VTE_CES.png", plot = scatter4[[1]], width = 10, height = 8, dpi = 300)

scatter4

single4 <- mr_singlesnp(H_VTE_CES)
single4 <- single4 %>%
  mutate(original_order = row_number())

# 根据 p 值和 b 值选择每个 SNP 的最佳行
single4_filtered <- single4 %>%
  group_by(SNP) %>%
  filter(p == min(p) & b == max(b)) %>%
  ungroup()

# 按原始顺序排序
single4_filtered <- single4_filtered %>%
  arrange(original_order) %>%
  select(-original_order)  # 删除原始顺序索引列
mr_forest_plot(single4_filtered)

loo4 <- mr_leaveoneout(H_VTE_CES)
loo4 <- loo4 %>%
  mutate(original_order = row_number())

loo4_filtered <- loo4 %>%
  group_by(SNP) %>%
  filter(p == min(p) & b == max(b)) %>%
  ungroup()

loo4_filtered <- loo4_filtered %>%
  arrange(original_order) %>%
  select(-original_order)

mr_leaveoneout_plot(loo4_filtered)

funnel4<- mr_funnel_plot(single4_filtered)
ggsave(filename = "funnel_VTE_CES.png", plot = funnel4[[1]], width = 10, height = 8, dpi = 300)

hetero_4 <- mr_heterogeneity(H_VTE_CES)
result <- hetero_4
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

pleio_4 <- mr_pleiotropy_test(H_VTE_CES)
result <- pleio_4
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

run_mr_presso(H_VTE_CES,NbDistribution = 2000)
#-----------FORWARD: VTE TO AIS-------------
H_VTE_AIS <- harmonise_data(combined_vte, VTE_AIS)
result5 <- mr(H_VTE_AIS)
write.csv(result5, "forward_VTE_AIS.csv")
odds5<-generate_odds_ratios(result5)

scatter5 <- mr_scatter_plot(result5, H_VTE_AIS)
ggsave(filename = "scatter_VTE_AIS.png", plot = scatter5[[1]], width = 10, height = 8, dpi = 300)

scatter5

single5 <- mr_singlesnp(H_VTE_AIS)
single5 <- single5 %>%
  mutate(original_order = row_number())

single5_filtered <- single5 %>%
  group_by(SNP) %>%
  filter(p == min(p) & b == max(b)) %>%
  ungroup()

single5_filtered <- single5_filtered %>%
  arrange(original_order) %>%
  select(-original_order)

mr_forest_plot(single5_filtered)

loo5 <- mr_leaveoneout(H_VTE_AIS)
loo5 <- loo5 %>%
  mutate(original_order = row_number())

loo5_filtered <- loo5 %>%
  group_by(SNP) %>%
  filter(p == min(p) & b == max(b)) %>%
  ungroup()

loo5_filtered <- loo5_filtered %>%
  arrange(original_order) %>%
  select(-original_order)
mr_leaveoneout_plot(loo5_filtered)

mr_funnel_plot(single5_filtered)

hetero_5 <- mr_heterogeneity(H_VTE_AIS)
result <- hetero_5
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

pleio_5 <- mr_pleiotropy_test(H_VTE_AIS)
result <- pleio_5
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

run_mr_presso(H_VTE_AIS,NbDistribution = 2000)

#-----------OUTPUT FORWARD RESULTS------------
write.csv(result1,"forward_VTE_AS.csv")
write.csv(odds1,"odds_VTE_AS.csv")

write.csv(result2, "forward_VTE_LAS.csv")
write.csv(odds2, "odds_VTE_LAS.csv")

write.csv(result3, "forward_VTE_SVS.csv")
write.csv(odds3, "odds_VTE_SVS.csv")

write.csv(result4, "forward_VTE_CES.csv")
write.csv(odds4, "odds_VTE_CES.csv")

write.csv(result5, "forward_VTE_AIS.csv")
write.csv(odds5, "odds_VTE_AIS.csv")


#------------reverse: AS to VTE------------

AS_to_VTE_1 <- extract_outcome_data(snps=allstroke$SNP, outcomes = "finn-b-I9_VTE")
write.csv(AS_to_VTE_1,"ASVTE1.csv")
AS_to_VTE_2 <- read_outcome_data(snps=allstroke$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                                 effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                                 samplesize_col = "samplesize")
write.csv(AS_to_VTE_2,"ASVTE2.csv")
colnames(AS_to_VTE_1)
colnames(AS_to_VTE_2)
colnames(AS_to_VTE_2)[colnames(AS_to_VTE_2) == 'chr.outcome'] <- 'chr'
colnames(AS_to_VTE_2)[colnames(AS_to_VTE_2) == 'pos.outcome'] <- 'pos'
final_colnames <- colnames(AS_to_VTE_1)
AS_to_VTE <- AS_to_VTE[, c(final_colnames, setdiff(colnames(AS_to_VTE), final_colnames))]
AS_to_VTE <- merge(AS_to_VTE_1, AS_to_VTE_2, by = "SNP", all = TRUE)
colnames(AS_to_VTE)
H_AS_to_VTE <- harmonise_data(allstroke, AS_to_VTE_1)
result6 <-mr(H_AS_to_VTE)
odds6<-generate_odds_ratios(result6)
write.csv(odds6,"odds_AS_VTE.csv")
scatter6 <- mr_scatter_plot(result6, H_AS_to_VTE)
ggsave(filename = "scatter_AS_VTE.png", plot = scatter6[[1]], width = 10, height = 8, dpi = 300)

single6 <- mr_singlesnp(H_AS_to_VTE)
mr_forest_plot(single6)
loo6 <- mr_leaveoneout(H_AS_to_VTE)
mr_leaveoneout_plot(loo6)
mr_funnel_plot(single6)

hetero_6 <- mr_heterogeneity(H_AS_to_VTE)
hetero_6
result <- hetero_6
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
pleio_6 <- mr_pleiotropy_test(H_AS_to_VTE)
result <- pleio_6
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
run_mr_presso(H_AS_to_VTE)

#------------reverse: LAS to VTE------------
# 提取结果数据
LAS_to_VTE_1 <- extract_outcome_data(snps = LAS$SNP, outcomes = "finn-b-I9_VTE")
write.csv(LAS_to_VTE_1, "LASVTE1.csv")

# 读取结果数据 （没有读取到任何相同的SNP）
LAS_to_VTE_2 <- read_outcome_data(
  snps = LAS$SNP, 
  filename = "data.csv", 
  sep = ",", 
  snp_col = "SNP", 
  beta_col = "beta", 
  se_col = "se",
  effect_allele_col = "effect_allele", 
  other_allele_col = "other_allele", 
  eaf_col = "eaf", 
  pval_col = "pval",
  samplesize_col = "samplesize"
)
LAS_to_VTE <- LAS_to_VTE_1 

# 协调数据
H_LAS_to_VTE <- harmonise_data(LAS, LAS_to_VTE)

# 进行MR分析
result7 <- mr(H_LAS_to_VTE)
write.csv(result7, "reverse_LAS_VTE.csv")

# 生成odds ratios
odds7<- generate_odds_ratios(result7)
write.csv(odds7,"odds_LAS_VTE.csv")

# 生成散点图
scatter7 <- mr_scatter_plot(result7, H_LAS_to_VTE)
ggsave(filename = "scatter_LAS_VTE.png", plot = scatter7[[1]], width = 10, height = 8, dpi = 300)
print(scatter7)

# 单SNP分析
single7 <- mr_singlesnp(H_LAS_to_VTE)
mr_forest_plot(single7)

# leave-one-out分析
loo7 <- mr_leaveoneout(H_LAS_to_VTE)
mr_leaveoneout_plot(loo7)

# 漏斗图
mr_funnel_plot(single7)

# 异质性检验
hetero_7 <- mr_heterogeneity(H_LAS_to_VTE)
result <- hetero_7
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

# 多效性检验
pleio_7 <- mr_pleiotropy_test(H_LAS_to_VTE)
result <- pleio_7
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

# 运行MR-PRESSO
run_mr_presso(H_LAS_to_VTE)


#------------reverse: SVS to VTE------------
SVS_to_VTE_1 <- extract_outcome_data(snps = SVS$SNP, outcomes = "finn-b-I9_VTE")
write.csv(SVS_to_VTE_1, "SVSVTE1.csv")

SVS_to_VTE_2 <- read_outcome_data(
  snps = SVS$SNP, 
  filename = "data.csv", 
  sep = ",", 
  snp_col = "SNP", 
  beta_col = "beta", 
  se_col = "se",
  effect_allele_col = "effect_allele", 
  other_allele_col = "other_allele", 
  eaf_col = "eaf", 
  pval_col = "pval",
  samplesize_col = "samplesize"
)
SVS_to_VTE <- SVS_to_VTE_1 

H_SVS_to_VTE <- harmonise_data(SVS, SVS_to_VTE)

result8 <- mr(H_SVS_to_VTE)
write.csv(result8, "reverse_SVS_VTE.csv")

odds8<-generate_odds_ratios(result8)
write.csv(odds8, "odds_SVS_VTE.csv")

scatter8 <- mr_scatter_plot(result8, H_SVS_to_VTE)
ggsave(filename = "scatter_SVS_VTE.png", plot = scatter8[[1]], width = 10, height = 8, dpi = 300)
print(scatter8)

single8 <- mr_singlesnp(H_SVS_to_VTE)
mr_forest_plot(single8)

loo8 <- mr_leaveoneout(H_SVS_to_VTE)
mr_leaveoneout_plot(loo8)

mr_funnel_plot(single8)

hetero_8 <- mr_heterogeneity(H_SVS_to_VTE)
result <- hetero_8
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
pleio_8 <- mr_pleiotropy_test(H_SVS_to_VTE)
result <- pleio_8
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
run_mr_presso(H_SVS_to_VTE)

#---------------------reverse:  CES to VTE ---------
CES_to_VTE_1 <- extract_outcome_data(snps = CES$SNP, outcomes = "finn-b-I9_VTE")
write.csv(CES_to_VTE_1, "CESVTE1.csv")

CES_to_VTE_2 <- read_outcome_data(
  snps = CES$SNP, 
  filename = "data.csv", 
  sep = ",", 
  snp_col = "SNP", 
  beta_col = "beta", 
  se_col = "se",
  effect_allele_col = "effect_allele", 
  other_allele_col = "other_allele", 
  eaf_col = "eaf", 
  pval_col = "pval",
  samplesize_col = "samplesize"
)
CES_to_VTE <- CES_to_VTE_1 

H_CES_to_VTE <- harmonise_data(CES, CES_to_VTE)

result9 <- mr(H_CES_to_VTE)
write.csv(result9, "reverse_CES_VTE.csv")

odds9 <- generate_odds_ratios(result9)
write.csv(odds9, "odds_CES_VTE.csv")

scatter9 <- mr_scatter_plot(result9, H_CES_to_VTE)
ggsave(filename = "scatter_CES_VTE.png", plot = scatter9[[1]], width = 10, height = 8, dpi = 300)
print(scatter9)

single9 <- mr_singlesnp(H_CES_to_VTE)
mr_forest_plot(single9)

loo9 <- mr_leaveoneout(H_CES_to_VTE)
mr_leaveoneout_plot(loo9)

mr_funnel_plot(single9)

hetero_9 <- mr_heterogeneity(H_CES_to_VTE)
result <- hetero_9
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
pleio_9 <- mr_pleiotropy_test(H_CES_to_VTE)
result <- pleio_9
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
run_mr_presso(H_CES_to_VTE)

#---------------------reverse:  AIS to VTE ---------
AIS_to_VTE_1 <- extract_outcome_data(snps = AIS$SNP, outcomes = "finn-b-I9_VTE")
write.csv(AIS_to_VTE_1, "AISVTE1.csv")

AIS_to_VTE_2 <- read_outcome_data(
  snps = AIS$SNP, 
  filename = "data.csv", 
  sep = ",", 
  snp_col = "SNP", 
  beta_col = "beta", 
  se_col = "se",
  effect_allele_col = "effect_allele", 
  other_allele_col = "other_allele", 
  eaf_col = "eaf", 
  pval_col = "pval",
  samplesize_col = "samplesize"
)
AIS_to_VTE <- AIS_to_VTE_1 

H_AIS_to_VTE <- harmonise_data(AIS, AIS_to_VTE)

result10 <- mr(H_AIS_to_VTE)
write.csv(result10, "reverse_AIS_VTE.csv")

odds10<-generate_odds_ratios(result10)
write.csv(odds10, "odds_AIS_VTE.csv")

scatter10 <- mr_scatter_plot(result10, H_AIS_to_VTE)
ggsave(filename = "scatter_AIS_VTE.png", plot = scatter10[[1]], width = 10, height = 8, dpi = 300)
print(scatter10)

single10 <- mr_singlesnp(H_AIS_to_VTE)
mr_forest_plot(single10)

loo10 <- mr_leaveoneout(H_AIS_to_VTE)
mr_leaveoneout_plot(loo10)

mr_funnel_plot(single10)

hetero_10 <- mr_heterogeneity(H_AIS_to_VTE)
result <- hetero_10
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
pleio_10 <- mr_pleiotropy_test(H_AIS_to_VTE)
result <- pleio_10
  write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
run_mr_presso(H_AIS_to_VTE)

#--------------reverse:output----------
write.csv(result6,"reverse_AS_VTE.csv")
write.csv(odds6,"odds_AS_VTE.csv")
odds7 <- generate_odds_ratios(result7)
write.csv(odds7, "odds_LAS_VTE.csv")

odds8 <- generate_odds_ratios(result8)
write.csv(odds8, "odds_SVS_VTE.csv")

odds9 <- generate_odds_ratios(result9)
write.csv(odds9, "odds_CES_VTE.csv")

odds10 <- generate_odds_ratios(result10)
write.csv(odds10, "odds_AIS_VTE.csv")
