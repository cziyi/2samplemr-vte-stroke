library(TwoSampleMR)
library(MRInstruments)
library(MendelianRandomization)
library(MRPRESSO)
library(httr)
library(remotes)
library(ieugwasr)
library(MVMR)
library(readr)
ieugwasr::api_status()
user()
id_allstroke <-"ebi-a-GCST005838"
id_las <- "ebi-a-GCST005840"
id_svs <- "ebi-a-GCST005841"
id_ces <- "ebi-a-GCST005842"
id_ais <- "ebi-a-GCST005843"
id_vte <- "finn-b-I9_VTE"
#-------多变量 smoking AS-------
smoking_AS <- mv_extract_exposures(c("ieu-b-4877","ieu-b-142",id_allstroke))
#out1<- read_outcome_data(snps=smoking_AS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
#                         effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
#                         samplesize_col = "samplesize")
out1<- extract_outcome_data(smoking_AS$SNP, id_vte)
h1 <- mv_harmonise_data(smoking_AS, out1)
result <- mv_multiple(h1)
mv_multiple(h1)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

smoking_LAS <- mv_extract_exposures(c("ieu-b-4877","ieu-b-142",id_las))
#out2<- read_outcome_data(snps=smoking_LAS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
#                         effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
#                         samplesize_col = "samplesize")
out2<- extract_outcome_data(smoking_LAS$SNP, id_vte)
h2 <- mv_harmonise_data(smoking_LAS, out2)
result <- mv_multiple(h2)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

smoking_SVS <- mv_extract_exposures(c("ieu-b-4877","ieu-b-142",id_svs))
#out3<- read_outcome_data(snps=smoking_SVS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
#                         effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
 #                        samplesize_col = "samplesize")
out3 <- extract_outcome_data(smoking_SVS$SNP, id_vte)
h3 <- mv_harmonise_data(smoking_SVS, out3)
result <- mv_multiple(h3)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

smoking_CES <- mv_extract_exposures(c("ieu-b-4877","ieu-b-142",id_ces))
#out4<- read_outcome_data(snps=smoking_CES$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
#                         effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
#                         samplesize_col = "samplesize")
out4 <- extract_outcome_data(smoking_CES$SNP, id_vte)
h4 <- mv_harmonise_data(smoking_CES, out4)
result <- mv_multiple(h4)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

smoking_AIS <- mv_extract_exposures(c("ieu-b-4877","ieu-b-142",id_ais))
out5<- read_outcome_data(snps=smoking_AIS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                         effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                         samplesize_col = "samplesize")
out5 <- extract_outcome_data(smoking_AIS$SNP, id_vte)
h5 <- mv_harmonise_data(smoking_AIS, out5)
result <- mv_multiple(h5)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

#-------多变量 obesity and stroke---------
obesity_AS <- mv_extract_exposures(c("finn-b-E4_OBESITY",id_allstroke))
out6 <-  read_outcome_data(snps=obesity_AS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                                effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                                samplesize_col = "samplesize")
out6 <- extract_outcome_data(obesity_AS$SNP, id_vte)
h6 <-mv_harmonise_data(obesity_AS, out6)
result <- mv_multiple(h6)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)


#---------------------------------
obesity_LAS <- mv_extract_exposures(c("finn-b-E4_OBESITY",id_las))
out7 <-  read_outcome_data(snps=obesity_LAS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                           effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                           samplesize_col = "samplesize")
out7 <- extract_outcome_data(obesity_LAS$SNP, id_vte)
h7 <- mv_harmonise_data(obesity_LAS, out7)
result <- mv_multiple(h7)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
#----------------------------------

obesity_SVS <- mv_extract_exposures(c("finn-b-E4_OBESITY",id_svs))
out8 <-  read_outcome_data(snps=obesity_SVS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                           effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                           samplesize_col = "samplesize")
out8 <- extract_outcome_data(obesity_SVS$SNP, id_vte)
h8 <- mv_harmonise_data(obesity_SVS, out8)
result <- mv_multiple(h8)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
#----------------------------------

obesity_CES <- mv_extract_exposures(c("finn-b-E4_OBESITY",id_ces))
out9 <-  read_outcome_data(snps=obesity_CES$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                           effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                           samplesize_col = "samplesize")
out9 <- extract_outcome_data(obesity_CES$SNP, id_vte)
h9 <- mv_harmonise_data(obesity_CES, out9)
result <- mv_multiple(h9)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
#----------------------------------
obesity_AIS <- mv_extract_exposures(c("finn-b-E4_OBESITY",id_ais))
out10 <-  read_outcome_data(snps=obesity_AIS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out10 <- extract_outcome_data(obesity_AIS$SNP, id_vte)
h10 <- mv_harmonise_data(obesity_AIS, out10)
result <- mv_multiple(h10)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)


#-------多变量 hyperlipidemia and stroke---------
hyperlip_AS <- mv_extract_exposures(c("finn-b-E4_HYPERLIPNAS",id_allstroke))
out11 <-  read_outcome_data(snps=hyperlip_AS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out11 <- extract_outcome_data(hyperlip_AS$SNP, id_vte)
h11 <- mv_harmonise_data(hyperlip_AS, out11)
result <- mv_multiple(h11)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

#----------------------------------
hyperlip_LAS <- mv_extract_exposures(c("finn-b-E4_HYPERLIPNAS", id_las))
out12 <-  read_outcome_data(snps=hyperlip_LAS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out12 <- extract_outcome_data(hyperlip_LAS$SNP, id_vte)
h12 <- mv_harmonise_data(hyperlip_LAS, out12)
result <- mv_multiple(h12)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

#----------------------------------
hyperlip_SVS <- mv_extract_exposures(c("finn-b-E4_HYPERLIPNAS", id_svs))
out13 <-  read_outcome_data(snps=hyperlip_SVS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out13 <- extract_outcome_data(hyperlip_SVS$SNP, id_vte)
h13 <- mv_harmonise_data(hyperlip_SVS, out13)
result <- mv_multiple(h13)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

#----------------------------------
hyperlip_CES <- mv_extract_exposures(c("finn-b-E4_HYPERLIPNAS", id_ces))
out14 <-  read_outcome_data(snps=hyperlip_CES$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out14 <- extract_outcome_data(hyperlip_CES$SNP, id_vte)
h14 <- mv_harmonise_data(hyperlip_CES, out14)
result <- mv_multiple(h14)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
#----------------------------------
hyperlip_AIS <- mv_extract_exposures(c("finn-b-E4_HYPERLIPNAS", id_ais))
out15 <-  read_outcome_data(snps=hyperlip_AIS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out15 <- extract_outcome_data(hyperlip_AIS$SNP, id_vte)
h15 <- mv_harmonise_data(hyperlip_AIS, out15)
result <- mv_multiple(h15)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

#----多变量 hyperTENS and stroke------------------------------
hypertens_AS <- mv_extract_exposures(c("finn-b-I9_HYPTENS", id_allstroke))
out16 <-  read_outcome_data(snps=hypertens_AS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out16
out16 <- extract_outcome_data(hypertens_AS$SNP, id_vte)
h16 <- mv_harmonise_data(hypertens_AS, out16)
result <- mv_multiple(h16)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

#----------------------------------

hypertens_LAS <- mv_extract_exposures(c("finn-b-I9_HYPTENS", id_las))
out16 <-  read_outcome_data(snps=hypertens_LAS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out16 <- extract_outcome_data(hypertens_LAS$SNP, id_vte)
h16 <- mv_harmonise_data(hypertens_LAS, out16)
result <- mv_multiple(h16)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

#----------------------------------
hypertens_SVS <- mv_extract_exposures(c("finn-b-I9_HYPTENS", id_svs))
out17 <-  read_outcome_data(snps=hypertens_SVS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out17 <- extract_outcome_data(hypertens_SVS$SNP, id_vte)
h17 <- mv_harmonise_data(hypertens_SVS, out17)
result <- mv_multiple(h17)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

#----------------------------------
hypertens_CES <- mv_extract_exposures(c("finn-b-I9_HYPTENS", id_ces))
out18 <-  read_outcome_data(snps=hypertens_CES$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out18 <- extract_outcome_data(hypertens_CES$SNP, id_vte)
h18 <- mv_harmonise_data(hypertens_CES, out18)
result <- mv_multiple(h18)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
#----------------------------------
hypertens_AIS <- mv_extract_exposures(c("finn-b-I9_HYPTENS", id_ais))
out19 <-  read_outcome_data(snps=hypertens_AIS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out19 <- extract_outcome_data(hypertens_AIS$SNP, id_vte)
h19 <- mv_harmonise_data(hypertens_AIS, out19)
result <- mv_multiple(h19)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)



#---------------------多变量 hemiple stroke-------------
hemiple_AS <- mv_extract_exposures(c("finn-b-G6_HEMIPLE", id_allstroke))
out21 <-  read_outcome_data(snps=hemiple_AS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out21
out21 <- extract_outcome_data(hemiple_AS$SNP, id_vte)
h21 <- mv_harmonise_data(hemiple_AS, out21)
result <- mv_multiple(h21)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
#----------------------------------
hemiple_LAS <- mv_extract_exposures(c("finn-b-G6_HEMIPLE", id_las))
out22 <-  read_outcome_data(snps=hemiple_LAS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out22 <- extract_outcome_data(hemiple_LAS$SNP, id_vte)
h22 <- mv_harmonise_data(hemiple_LAS, out22)
result <- mv_multiple(h22)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
#----------------------------------
hemiple_SVS <- mv_extract_exposures(c("finn-b-G6_HEMIPLE", id_svs))
out23 <-  read_outcome_data(snps=hemiple_SVS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out23 <- extract_outcome_data(hemiple_SVS$SNP, id_vte)
h23 <- mv_harmonise_data(hemiple_SVS, out23)
result <- mv_multiple(h23)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
#----------------------------------
hemiple_CES <- mv_extract_exposures(c("finn-b-G6_HEMIPLE", id_ces))
out24 <-  read_outcome_data(snps=hemiple_CES$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out24 <- extract_outcome_data(hemiple_CES$SNP, id_vte)
h24 <- mv_harmonise_data(hemiple_CES, out24)
result <- mv_multiple(h24)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
#----------------------------------
hemiple_AIS <- mv_extract_exposures(c("finn-b-G6_HEMIPLE", id_ais))
out25 <-  read_outcome_data(snps=hemiple_AIS$SNP, filename ="data.csv",sep = ",", snp_col = "SNP", beta_col = "beta", se_col ="se",
                            effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col ="eaf", pval_col = "pval",
                            samplesize_col = "samplesize")
out25
out25 <- extract_outcome_data(hemiple_AIS$SNP, id_vte)
h25 <- mv_harmonise_data(hemiple_AIS, out25)
result <- mv_multiple(h25)
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)


#---general output code---
#
result <- 
write.table(result, "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)
#