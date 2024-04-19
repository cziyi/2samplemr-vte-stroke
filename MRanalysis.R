ieugwasr::api_status()
library(TwoSampleMR)
library(MRInstruments)
library(MendelianRandomization)
library(MRPRESSO)
data("gwas_catalog")



#-----------------------#
#import VTE database:   #
#-----------------------#

formatted <- format_data(data)
VTE_new<- clump_data(formatted)

#-----------------------------------#
#extract outcome, VTE as exposure:  #
#-----------------------------------#

VTEnew_to_AS <- extract_outcome_data(snps = VTE_new$SNP, outcomes ="ebi-a-GCST005838" )
VTEnew_to_LAS <- extract_outcome_data(snps = VTE_new$SNP, outcomes = "ebi-a-GCST005840")
VTEnew_to_SVS <- extract_outcome_data(snps = VTE_new$SNP, outcomes = "ebi-a-GCST005841")
VTEnew_to_CES <- extract_outcome_data(snps = VTE_new$SNP, outcomes = "ebi-a-GCST005842")
VTEnew_to_AIS <- extract_outcome_data(snps = VTE_new$SNP, outcomes = "ebi-a-GCST005843")

#-------------------
# analysis: 
# VTE to All Stroke
#
#-------------------

H_VTEnew_to_AS <- harmonise_data(VTE_new, VTEnew_to_AS)
mr_results_1 <- mr(H_VTEnew_to_AS)
mr_randomivw1<- mr(H_VTEnew_to_AS, method_list = "mr_ivw_mre")
scatter1 <- mr_scatter_plot(mr_results_1, H_VTEnew_to_AS)
scatter1
scatter1[1]$`XTcQ9L.ebi-a-GCST005838`$labels$x <- "SNP effect on VTE (Thibord et al.)"
scatter1[1]$`XTcQ9L.ebi-a-GCST005838`$labels$y <- "SNP effect on AS"
single1<- mr_singlesnp(H_VTEnew_to_AS)
forest1<-mr_forest_plot(single1)
forest1$`XTcQ9L.ebi-a-GCST005838`$labels$x<-"MR effect size for VTE on AS"
forest1

loo1 <-mr_leaveoneout(H_VTEnew_to_AS)
looplot1<- mr_leaveoneout_plot(loo1)
looplot1$`XTcQ9L.ebi-a-GCST005838`$labels$x <-"MR leave-one-out sensitivity analysis for VTE on AS" 
looplot1

funnel1 <- mr_funnel_plot(single1)
funnel1

hetero_1 <- mr_heterogeneity(H_VTEnew_to_AS)
pleio_1 <- mr_pleiotropy_test(H_VTEnew_to_AS)
run_mr_presso(H_VTEnew_to_AS, NbDistribution = 2000)
generate_odds_ratios(mr_results_1)
#---------
# analysis: 
# VTE to LAS
#
#----------
H_VTEnew_to_LAS <-harmonise_data(VTE_new, VTEnew_to_LAS)
mr_results_2 <- mr(H_VTEnew_to_LAS)
mr_randomivw2<- mr(H_VTEnew_to_LAS, method_list = "mr_ivw_mre")

scatter2 <- mr_scatter_plot(mr_results_2, H_VTEnew_to_LAS)
scatter2$`XTcQ9L.ebi-a-GCST005840`$labels$x <- "SNP effect on VTE (Thibord et al.)"
scatter2$`XTcQ9L.ebi-a-GCST005840`$labels$y <-"SNP effect on LAS"
scatter2

single2<- mr_singlesnp(H_VTEnew_to_LAS)
forest2<- mr_forest_plot(single2)
forest2$`XTcQ9L.ebi-a-GCST005840`$labels$x <-"MR effect size for VTE on LAS"
forest2

loo2 <-mr_leaveoneout(H_VTEnew_to_LAS)
looplot2 <-mr_leaveoneout_plot(loo2)
looplot2$`XTcQ9L.ebi-a-GCST005840`$labels$x <-"MR leave-one-out sensitivity analysis for VTE on LAS"
looplot2

mr_leaveoneout_plot(loo2)
mr_funnel_plot(single2)
hetero_2 <- mr_heterogeneity(H_VTEnew_to_LAS)
pleio_2 <- mr_pleiotropy_test(H_VTEnew_to_LAS)
run_mr_presso(H_VTEnew_to_LAS, NbDistribution = 2000)
generate_odds_ratios(mr_results_2)
#---------
# analysis: 
# VTE to SVS
#
#----------

H_VTEnew_to_SVS <-harmonise_data(VTE_new, VTEnew_to_SVS)
mr_results_3 <- mr(H_VTEnew_to_SVS)
mr_randomivw3<- mr(H_VTEnew_to_SVS, method_list = "mr_ivw_mre")

scatter3 <- mr_scatter_plot(mr_results_3, H_VTEnew_to_SVS)
scatter3$`x8p9tV.ebi-a-GCST005841`$labels$x <-"SNP effect on VTE (Thibord et al.)"
scatter3$`x8p9tV.ebi-a-GCST005841`$labels$y <-"SNP effect on SVS"
scatter3

single3<- mr_singlesnp(H_VTEnew_to_SVS)
forest3 <-mr_forest_plot(single3)
forest3$`x8p9tV.ebi-a-GCST005841`$labels$x<-"MR effect size for VTE on SVS"
forest3

loo3 <-mr_leaveoneout(H_VTEnew_to_SVS)
looplot3<-mr_leaveoneout_plot(loo3)
looplot3$`x8p9tV.ebi-a-GCST005841`$labels$x <-"MR leave-one-out sensitivity analysis for VTE on SVS"
looplot3

mr_funnel_plot(single3)

hetero_3<-mr_heterogeneity(H_VTEnew_to_SVS)
pleio_3<-mr_pleiotropy_test(H_VTEnew_to_SVS)
run_mr_presso(H_VTEnew_to_SVS, NbDistribution = 2000)
generate_odds_ratios(mr_results_3)
#---------
# analysis: 
# VTE to CES
#
#----------

H_VTEnew_to_CES <-harmonise_data(VTE_new, VTEnew_to_CES)
mr_results_4<-mr(H_VTEnew_to_CES)
mr_randomivw4<- mr(H_VTEnew_to_CES, method_list = "mr_ivw_mre")

hetero4<-mr_heterogeneity(H_VTEnew_to_CES)
pleio4<-mr_pleiotropy_test(H_VTEnew_to_CES)

scatter4 <- mr_scatter_plot(mr_results_4, H_VTEnew_to_CES)
scatter4$`x8p9tV.ebi-a-GCST005842`$labels$x<-"SNP effect on VTE (Thibord et al.)"
scatter4$`x8p9tV.ebi-a-GCST005842`$labels$y <-"SNP effect on CES"
scatter4

single4<-mr_singlesnp(H_VTEnew_to_CES)
forest4<-mr_forest_plot(single4)
forest4$`x8p9tV.ebi-a-GCST005842`$labels$x<-"MR effect size for VTE on CES"
forest4

loo4 <-mr_leaveoneout(H_VTEnew_to_CES)
looplot4<-mr_leaveoneout_plot(loo4)
looplot4$`x8p9tV.ebi-a-GCST005842`$labels$x <-"MR leave-one-out sensitivity analysis for VTE on CES"
looplot4

mr_funnel_plot(single4)

run_mr_presso(H_VTEnew_to_CES, NbDistribution = 2000)
generate_odds_ratios(mr_results_4)
#---------
#---------
# analysis: 
# VTE to AIS
#
#----------

H_VTEnew_to_AIS <-harmonise_data(VTE_new, VTEnew_to_AIS)
mrresults5<-mr(H_VTEnew_to_AIS)
mr_randomivw5<- mr(H_VTEnew_to_AIS, method_list = "mr_ivw_mre")
scatter5 <- mr_scatter_plot(mrresults5, H_VTEnew_to_AIS)
scatter5$`x8p9tV.ebi-a-GCST005843`$labels$x<-"SNP effect on VTE (Thibord et al.)"
scatter5$`x8p9tV.ebi-a-GCST005843`$labels$y <-"SNP effect on AIS"
scatter5

single5<-mr_singlesnp(H_VTEnew_to_AIS)
forest5<-mr_forest_plot(single5)
forest5$`x8p9tV.ebi-a-GCST005843`$labels$x<-"MR effect size for VTE on AIS"
forest5

loo5 <-mr_leaveoneout(H_VTEnew_to_AIS)
looplot5<-mr_leaveoneout_plot(loo5)
looplot5$`x8p9tV.ebi-a-GCST005843`$labels$x <-"MR leave-one-out sensitivity analysis for VTE on AIS"
looplot5

mr_funnel_plot(single5)

hetero5<-mr_heterogeneity(H_VTEnew_to_AIS)
pleio5<-mr_pleiotropy_test(H_VTEnew_to_AIS)
directionality_test(H_VTEnew_to_AIS)
run_mr_presso(H_VTEnew_to_AIS, NbDistribution = 2000)
generate_odds_ratios(mrresults5)

#-------------------------------------
# extract instrument SNPs for strokes:  
#-------------------------------------

allstroke <- extract_instruments("ebi-a-GCST005838")
#allstroke_exposure <-extract_instruments("ebi-a-GCST005838", p1=5e-6)
LAS <- extract_instruments("ebi-a-GCST005840",p1= 5e-6)
SVS <- extract_instruments("ebi-a-GCST005841",p1= 5e-6)
CES <- extract_instruments("ebi-a-GCST005842",p1=5e-6)
AIS <- extract_instruments("ebi-a-GCST005843")

##-----------------------
## AS to VTE:
##-----------------------

AS_to_VTE <- extract_outcome_data(snps=allstroke$SNP, outcomes = "finn-b-I9_VTE")
H_AS_to_VTE <-harmonise_data(allstroke, AS_to_VTE)
mr_results6 <- mr(H_AS_to_VTE)
mr_randomivw6<- mr(H_AS_to_VTE, method_list = "mr_ivw_mre")
mr(H_AS_to_VTE, method_list = "mr_ivw_mre")

scatter6 <- mr_scatter_plot(mr_results6, H_AS_to_VTE)
scatter6$`ebi-a-GCST005838.finn-b-I9_VTE`$labels$x <-"SNP effect on AS"
scatter6$`ebi-a-GCST005838.finn-b-I9_VTE`$labels$y <-"SNP effect on VTE (FinnGen)"
scatter6

single6<-mr_singlesnp(H_AS_to_VTE)
forest6<-mr_forest_plot(single6)
forest6$`ebi-a-GCST005838.finn-b-I9_VTE`$labels$x<-"MR effect size for AS on VTE"
forest6

loo6 <-mr_leaveoneout(H_AS_to_VTE)
looplot6<-mr_leaveoneout_plot(loo6)
looplot6$`ebi-a-GCST005838.finn-b-I9_VTE`$labels$x <-"MR leave-one-out sensitivity analysis for AS on VTE"
looplot6

mr_funnel_plot(single6)

hetero6 <-mr_heterogeneity(H_AS_to_VTE)
pleio6<-mr_pleiotropy_test(H_AS_to_VTE)
directionality_test(H_AS_to_VTE)
run_mr_presso(H_AS_to_VTE, NbDistribution = 2000)
odds6<-generate_odds_ratios(mr_results6)

##-----------------------
## LAS to VTE:
##-----------------------
LAS_to_VTE <- extract_outcome_data(snps = LAS$SNP, outcomes = "finn-b-I9_VTE")
H_LAS_to_VTE <-harmonise_data(LAS, LAS_to_VTE)
mr_results7 <-mr(H_LAS_to_VTE)
mr_randomivw7<- mr(H_LAS_to_VTE, method_list = "mr_ivw_mre")

scatter7 <- mr_scatter_plot(mr_results7, H_LAS_to_VTE)
scatter7$`ebi-a-GCST005840.finn-b-I9_VTE`$labels$x <-"SNP effect on LAS"
scatter7$`ebi-a-GCST005840.finn-b-I9_VTE`$labels$y <-"SNP effect on VTE (FinnGen)"
scatter7

single7<-mr_singlesnp(H_LAS_to_VTE)
forest7<-mr_forest_plot(single7)
forest7$`ebi-a-GCST005840.finn-b-I9_VTE`$labels$x<-"MR effect size for LAS on VTE"
forest7

loo7 <-mr_leaveoneout(H_LAS_to_VTE)
looplot7<-mr_leaveoneout_plot(loo7)
looplot7$`ebi-a-GCST005840.finn-b-I9_VTE`$labels$x <-"MR leave-one-out sensitivity analysis for LAS on VTE"
looplot7

mr_funnel_plot(single7)

hetero7 <- mr_heterogeneity(H_LAS_to_VTE)
pleio7<-mr_pleiotropy_test(H_LAS_to_VTE)
directionality_test(H_LAS_to_VTE)
run_mr_presso(H_LAS_to_VTE, NbDistribution = 2000)
generate_odds_ratios(mr_results7)
##-----------------------
## SVS to VTE:
##-----------------------
SVS_to_VTE <- extract_outcome_data(snps = SVS$SNP, outcomes = "finn-b-I9_VTE")
H_SVS_to_VTE <-harmonise_data(SVS, SVS_to_VTE)
mr_results8 <-mr(H_SVS_to_VTE)
mr_randomivw8<- mr(H_SVS_to_VTE, method_list = "mr_ivw_mre")


scatter8 <- mr_scatter_plot(mr_results8, H_SVS_to_VTE)
scatter8$`ebi-a-GCST005841.finn-b-I9_VTE`$labels$x <-"SNP effect on SVS"
scatter8$`ebi-a-GCST005841.finn-b-I9_VTE`$labels$y <-"SNP effect on VTE (FinnGen)"
scatter8

single8<-mr_singlesnp(H_SVS_to_VTE)
forest8<-mr_forest_plot(single8)
forest8$`ebi-a-GCST005841.finn-b-I9_VTE`$labels$x<-"MR effect size for SVS on VTE"
forest8

loo8 <-mr_leaveoneout(H_SVS_to_VTE)
looplot8<-mr_leaveoneout_plot(loo8)
looplot8$`ebi-a-GCST005841.finn-b-I9_VTE`$labels$x <-"MR leave-one-out sensitivity analysis for SVS on VTE"
looplot8

mr_funnel_plot(single8)

oo8<-generate_odds_ratios(mr_randomivw8)
hetero8 <- mr_heterogeneity(H_SVS_to_VTE)
pleio8<-mr_pleiotropy_test(H_SVS_to_VTE)
directionality_test(H_SVS_to_VTE)
run_mr_presso(H_SVS_to_VTE, NbDistribution = 2000)
generate_odds_ratios(mr_results8)

##-----------------------
## CES to VTE:
##-----------------------
CES_to_VTE <- extract_outcome_data(snps = CES$SNP, outcomes = "finn-b-I9_VTE")
H_CES_to_VTE <-harmonise_data(CES, CES_to_VTE)
mr_results9 <-mr(H_CES_to_VTE)
mr_randomivw9<- mr(H_CES_to_VTE, method_list = "mr_ivw_mre")

scatter9 <- mr_scatter_plot(mr_results9, H_CES_to_VTE)
scatter9$`ebi-a-GCST005842.finn-b-I9_VTE`$labels$x <-"SNP effect on CES"
scatter9$`ebi-a-GCST005842.finn-b-I9_VTE`$labels$y <-"SNP effect on VTE (FinnGen)"
scatter9

single9<-mr_singlesnp(H_CES_to_VTE)
forest9<-mr_forest_plot(single9)
forest9$`ebi-a-GCST005842.finn-b-I9_VTE`$labels$x<-"MR effect size for CES on VTE"
forest9

loo9 <-mr_leaveoneout(H_CES_to_VTE)
looplot9<-mr_leaveoneout_plot(loo9)
looplot9$`ebi-a-GCST005842.finn-b-I9_VTE`$labels$x <-"MR leave-one-out sensitivity analysis for CES on VTE"
looplot9

mr_funnel_plot(single9)

hetero9 <- mr_heterogeneity(H_CES_to_VTE)
pleio9 <- mr_pleiotropy_test(H_CES_to_VTE)
directionality_test(H_CES_to_VTE)

run_mr_presso(H_CES_to_VTE, NbDistribution = 2000)
generate_odds_ratios(mr_results9)
##-----------------------
## AIS to VTE:
##-----------------------
AIS_to_VTE <- extract_outcome_data(snps = AIS$SNP, outcomes = "finn-b-I9_VTE")
H_AIS_to_VTE <-harmonise_data(AIS, AIS_to_VTE)
mr_results10 <-mr(H_AIS_to_VTE)
mr_randomivw10<- mr(H_AIS_to_VTE, method_list = "mr_ivw_mre")

scatter10 <- mr_scatter_plot(mr_results10, H_AIS_to_VTE)
scatter10$`ebi-a-GCST005843.finn-b-I9_VTE`$labels$x <-"SNP effect on AIS"
scatter10$`ebi-a-GCST005843.finn-b-I9_VTE`$labels$y <-"SNP effect on VTE (FinnGen)"
scatter10

single10<-mr_singlesnp(H_AIS_to_VTE)
forest10<-mr_forest_plot(single10)
forest10$`ebi-a-GCST005843.finn-b-I9_VTE`$labels$x<-"MR effect size for AIS on VTE"
forest10

loo10 <-mr_leaveoneout(H_AIS_to_VTE)
looplot10<-mr_leaveoneout_plot(loo10)
looplot10$`ebi-a-GCST005843.finn-b-I9_VTE`$labels$x <-"MR leave-one-out sensitivity analysis for AIS on VTE"
looplot10

mr_funnel_plot(single10)

hetero10 <- mr_heterogeneity(H_AIS_to_VTE)
pleio10 <- mr_pleiotropy_test(H_AIS_to_VTE)
directionality_test(H_AIS_to_VTE)
run_mr_presso(H_AIS_to_VTE, NbDistribution = 2000)
generate_odds_ratios(mr_results10)

