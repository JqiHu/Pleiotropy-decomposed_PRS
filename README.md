# README file of descriptions of scripts
Robust pleiotropy-decomposed polygenic scores(PD-PRSs) identify distinct contributions to elevated coronary artery disease (CAD) polygenic risk
##
# Step 1. CAD-related regions identification
Use GNOVA to calculate genetic correlations between CAD and 2,128 GWAS summary statistics. \
Folder: ./1_1_gnova_2128 \
all_cad_gnova.R # generate all_cad_gnova.job \
all_cad_gnova.job # jobs to run gnova one by one for 2,128 sum stats \
all_cad_gnova.sh # create all_cad_gnova.pbs and submit it \
gnova_res.R # organize results from GNOVA and output FDR control results 
##
Then organize traits with similar patterns together and define biologically meaningful pathways manually. \
In total, 43 summary stats have been selected and 8 clusters have been identified. \
Use GNOVA to check the genetic correlations between 43 selected traits/diseases. \
Folder: ./1_2_gnova_43 \
gnova_43.R # generate gnova_43.job \
gnova_43.job # jobs to run pair-wise GNOVA \
gnova_43.sh # create gnova_43.pbs and submit it \
gnova_43_res.R # organize results from GNOVA and output correlation coefficient/p-value matrix 
##
# Step 2. CAD PRS construction
Use different PRS methods to define the best PRS for CAD. \
The best we use is from AnnoPred. \
The weights for CAD PRS can be obtained upon request to Prof.Hongyu Zhao (hongyu.zhao@yale.edu). 
##
# Step 3. Decomposition of CAD PRS into PD-PRSs
Partition chromosomes into 2,353 LD blocks based on the partition file provided by SUPERGNOVA. \
Conduct SUPERGNOVA to calculate the local genetic correlations between 43 traits/diseases from step 1 and CAD. \
Folder: ./3_1_supergnova \
job_gene.R # generate su43.job \
su43.job # jobs to run SUPERGNOVA between CAD and 43 traits/diseases \
su43.sh # create su43.pbs and submit it \
su43.pbs # job to be submitted 
## 
Define pleiotropy-decomposed regions as regions where the traits/disease are significantly correlated with CAD. \
Select pleiotropy-decomposed SNPs out of these regions. \
Folder: ./3_2_snp_select \
select_snps.R # Extract SNPs within CAD-correlated regions for input traits \
job_gene.R # generate select_snps.job \
select_snps.job # jobs to run select_snps.R for 43 traits \
select_snps.sh # create select_snps.pbs and submit it \
select_snps.pbs # job to be submitted \
select_final.R # combine traits to form clusters and overlapping regions were moved to the clusters with larger absolute values of local genetic covariance \
select_ns_snps.R # form non-specific SNP subset/SNPs not in any of the pleiotropy-decomposed regions. 
##
Partition the CAD PRS into PD-PRSs \
Folder: ./3_3_prs_calculation \
ps_prs.job # parallel jobs to calculate PD-PRSs \
prs_prepare.R # merge all PD-PRS files into one table 
##
# Step 4. Statistical analysis
# 4.1 Phenome-wide association analysis of PD-PRSs
Correlations between PD-PRSs and corresponding CAD-related traits/diseases \
Folder: ./4_1_prs_analysis \
prs_corr.R # correlations between 9 PD-PRSs \
prs_pheno_corr.R # correlations between 9 PD-PRSs and 59 phenotypes \
hr.R # prediction of 9 PD-PRSs for CAD 

# 4.2 Identify pleiotropy subgroups among people with high CAD genetic risk
Genetic subgroups for CAD identification based on PD-PRSs.\
Folder: ./4_2_subgroup \
classification.R # classify high-risk (top 5% CAD PRS) into 9 subgroups (allow overlapping) \
pattern_size.R # identify subgroup patterns for high-risk subjects \
rela_change.R # relative changes on specific traits between certain subgroup and remaining high-risk subjects \
compare.R # statistical tests on the differences in specific traits between certain subgroups and remaining high-risk subjects 

# 4.3 Interaction analyses
Interactions between PD-PRSs/subgroups and certain traits. \
Folder: ./4_3_interaction \
step0_data_prep_INT.R # Transform phenotype values by adjusting confounding factors and PD-PRS and apply inverse normal transformation on residuals \ 
step1_inter_test.R # Test the interactions between transformed traits/diseases and PRSs \
step2_inter_BH.R # B-H FDR control for multiple comparisons \
step3_inter_subgroup_hr.R # Re-define a categorical variable with 4 levels using subgroup and dichotomized traits and calculate the HRs \
step4_inter_subgroup_arr.R # Comparing the absolute risk (proportion of CAD) in four levels of the categorical variable 

##
# Other folder information
## 1. ./trait_prs/
Construct PRS for representative traits in each pleiotropy-decomposed region to compare with PD-PRSs. \
beta_prep.R # format gwas summary stats for PRS-CS \
beta_job.R # generate beta.job \
beta.job #  calculate PRScs betas for 8 traits \
beta.sh # generate beta.pbs and submit for job running \
beta.pbs # job to be submitted \
score.job # parallele prs calculations 
  
hr_trait_prs.R # prediction performance of 9 trait prs for CAD 

## 2. ./vis/
Contains scripts for figure generations. \
The name corresponds to the figure tag in the article. 
#### Main figures
Figure2_corr_traits_prs.R # Figure 2. Correlations between phenotypes and PD-PRSs. \
Figure3_subgroup_stack.R # Figure 3. PD-PRS profiles in different pleiotropy subgroups. \
Figure4_subgroup_relative_change.R # Figure 4. Relative changes of quantitative traits in pleiotropy subgroups. \
Figure5_SuppFigure6_interaction_subgroup.R # Figure 5. Interactions between PD-PRS subgroup and smoking status. And also Supplementary Figure 6. \ 
#### Supplementary figures
SuppFigure1_corr_traits.R # Supplementary Figure 1. Genetic correlations between traits. \
SuppFigure3_HR.R # Supplementary Figure 3. Hazards ratios for CAD of overall CAD PRS, 8 PD-PRSs, and 1 NS-PRS. \
SuppFigure4_pattern.R # Supplementary Figure 4. Number of subjects in the top subgroup patterns. \
SuppFigure8_corr_pdprs.R # Supplementary Figure 8. Correlations between 8 PD-PRS and 1 NS-PRS. 


## 3. ./sensitivity_analysis
Scripts for sensitivity analyses. \

#### Hierarchical clustering for 43 selected traits
step1_hc.R # Perform hierarchical clustering on genetic correlation matrix; generate heatmaps for number of clusters = 5, 7, 9, 15. \
step2_snp_select.R # Assign SNPs to newly defined clusters. \
step3_pd_prs.sh # Calculate PD-PRSs based on newly defined clusters. \
step4_pd_prs_organize.R # Organize PD-PRSs into one table. \
step5_subgroup_classification.R # Partition high-risk subjects into genetic subgroups. \
step6_subgroup_rela_change.R # Calculate relative changes of traits of interest in subgroups. \
step7_res_organize.R # Plot Supplemental Figure 2. \

#### Subgroup thresholds at 1% or 10%
step1_subgroup_classification.R # Partition high-risk subjects into genetic subgroups based on certain threshold. \
step2_rela_change.R # Calculate relative changes of traits of interest in subgroups. \
step3_vis.R # Plot Supplemental Figure 5. \

#### Sensitivity analyses for interactions
##### Analysis 1. Adjust for baseline covariates only
analysis1_step0_data_prep_INT.R # Transform phenotype values by adjusting baseline covariates and PD-PRS and apply inverse normal transformation on residuals \
analysis1_step1_inter_test.R # Test the interactions between transformed traits/diseases and PRSs \

##### Analysis 2. Directly use residuals without INT
analysis2_step0_data_prep_Residual.R # Take Residuals after adjusting confounding and PD-PRS \
analysis2_step1_inter_test.R # Test the interactions between transformed traits/diseases and PRSs \

step2_inter_BH.R # B-H FDR control for multiple comparisons \

## 4. ./simulations
Scripts for simulations of 1 target trait and 3 quantitative traits. \

##### PD-PRS framework
step1_beta_simulation.R # Simulate beta region by region. \
step2_SNP_subset.R # Assign SNPs to clusters. \
step3_PD_PRS.sh # Calculate PD-PRSs based on classified SNP subsets. \
step3_PD_PRS_organize.R # Organize PD-PRSs. \
step4_1_Y_genetic.sh # Simulate the genetic phenotype. \
step4_2_Y_simulate.R # Simulate phenotypes. \
step5_subgroup_classification.R # Partition high-risk subjects into subgroups. \
step5_subgroup_stack.R # Generate stack plots for compositions of PD-PRSs across subgroups. \
step5_subgroup_relachange.R # Calculate relative changes of traits of interest in subgroups. \
step6_inter_prep_INT.R # Transform simulated traits by inverse normal transfomration. \
step6_inter_test.R # Test the interactions between transformed traits and PRSs. \
step6_inter_test_subgroup.R # Test the interactions between transformed traits and subgroup membership. \

##### Method by Chasman et al., 2019
./method_chasman/step1_LD_prune.sh # LD pruning on simulated GWAS. \
./method_chasman/step2_SNP_select.R # Select causal SNPs for trait 1. \
./method_chasman/step3_scale_beta.R # Scale Beta. \
./method_chasman/step4_beta_decomposition.R # Decompose Phenotype-SNP matrix. \
./method_chasman/step5_prs.sh # Calculate component PRSs and overall PRS. \
./method_chasman/step5_prs_organize.R # Organize PRSs. \
./method_chasman/step6_assoc.R # Associations between component PRSs and simulated traits. \
./method_chasman/step7_subgroup_classification.R # Partition high-risk subjects into subgroups. \
./method_chasman/step7_subgroup_stack.R # Generate stack plots for compositions of PD-PRSs across subgroups. \
./method_chasman/step7_subgroup_relachange.R # Calculate relative changes of traits of interest in subgroups. \

##### Visualize results
step7_res_vis.R # Generate Supplemental Figure 7. \
