# README file of decriptions of scripts
Robust pleiotropy-decomposed polygenic scores(PD-PRSs) identify distinct contributions to elevated coronary artery disease(CAD) polygenic risk

# Step 1. CAD-related regions identification
Use GNOVA to calculate genetic correlations between CAD and 2,128 GWAS summary statistics.
Folder: ./1_1_gnova_2128
## Scripts
all_cad_gnova.R # generate all_cad_gnova.job
all_cad_gnova.job # jobs to run gnova one by one for 2,128 sum stats
all_cad_gnova.sh # create all_cad_gnova.pbs and submit it
gnova_res.R # organize results from GNOVA and output FDR control results

Then organize traits with similar patterns together and define biologically meaninful pathways manually.
In total, 43 summary stats have been selected and 8 clusters have been identified.
Use GNOVA to check the genetic correlations between 43 selected traits/diseases.
Folder: ./1_2_gnova_43
## Scripts
gnova_43.R # generate gnova_43.job
gnova_43.job # jobs to run pair-wise GNOVA
gnova_43.sh # create gnova_43.pbs and submit it
gnova_43_res.R # organize results from GNOVA and output correlation coefficient/p-value matrix

# Step 2. CAD PRS construction
Use different PRS methods to define the best PRS for CAD.
The best we use is from AnnoPred.
The weights for CAD PRS can be obtained upon request to Prof.Hongyu Zhao (hongyu.zhao@yale.edu).

# Step 3. Decomposition CAD PRS into PD-PRSs
Partition chromosomes into 2,353 LD blocks based on the partition file provided by SUPERGNOVA.
Conduct SUPERGNOVA to calculate the local genetic correlations between 43 traits/diseases from step 1 and CAD.
Folder: ./3_1_supergnova
## Scripts
job_gene.R # generate su43.job
su43.job # jobs to run SUPERGNOVA between CAD and 43 traits/diseases
su43.sh # create su43.pbs and submit it
su43.pbs # job to be submitted

Define pleiotropy-decomposed regions as regions where the traits/disease are significantly correlated with CAD.
Select pleiotropy-decomposed SNPs out of these regions.
Folder: ./3_2_snp_select
## Scripts
select_snps.R # Extract SNPs within CAD-correlated regions for input traits
job_gene.R # generate select_snps.job
select_snps.job # jobs to run select_snps.R for 43 traits
select_snps.sh # create select_snps.pbs and submit it
select_snps.pbs # job to be submitted

select_final.R # combine traits to form clusters and overlapping regions were moved to the clusters with larger absolute value of local genetic covariance
select_ns_snps.R # form non-specific SNP subset/SNPs not in any of the pleiotropy-decomposed regions.

Partition the CAD PRS into PD-PRSs
Folder: ./3_3_prs_calculation
## Scripts
ps_prs.job # parallel jobs to calculate PD-PRSs
prs_prepare.R # merge all PD-PRS files into one table

# Step 4. Statistical analysis
# 4.1 Phenome-wide association analysis of PD-PRSs
Correlations between PD-PRSs and corresponding CAD-related traits/diseases
Folder: ./4_1_prs_analysis
## Scripts
prs_corr.R # correlations between 9 PD-PRSs
prs_pheno_corr.R # correlations between 9 PD-PRSs and 59 phenotypes
hr.R # prediction of 9 PD-PRSs for CAD

# 4.2 Identify subgroups among people with high CAD genetic risk
Genetic subgroups for CAD identification based on PD-PRSs.
Folder: ./4_2_subgroup
## Scripts
clasiffication.R # classify high-risk (top 5% CAD PRS) into 9 subgroups (allow overlapping)
rela_change.R # relative changes on specific traits between certain subgroup and remaining high-risk subjects
compare.R # statistical tests on the differences in specific traits between certain subgroup and remaining high-risk subjects

# 4.3 Interaction analyses
Interactions between PD-PRSs/subgroups and certain traits.
Folder: ./4_3_interaction
## Scripts
### psPRS * traits
step1_inter_hr.R # test for interactions between 29 traits and 9 PD-PRSs + 1 overall CAD PRS
step2_inter_hr_sensitivity.R # test for interactions between traits [adjusting confounding factors and PD-PRS] and 9 PD-PRSs + 1 overall CAD PRS
step3_inter_subgroup_hr.R # Re-define a categorical variables with 4 levels using subgroup and dichotomized traits and calculated the HRs
step4_inter_subgroup_arr.R # Comparing the absolute risk (proportion of CAD) in four levels of the categorical variable

# Step 5. Additional analysis: integrated PRS
Combine the PD-PRSs together to see if any improvement in prediction performance.
Folder: ./integrate_prs
## Scripts
integrate.R # combine all PD-PRSs together and compare the HR per sd and HR of top 10% v.s. bottom 10%


# Other folder information
1. ./trait_prs/
Construct PRS for representative traits in each pleiotropy-decomposed region to compare with PD-PRSs.
## Scripts
beta_prep.R # format gwas summary stats for PRS-CS
beta_job.R # generate beta.job 
beta.job #  calculate PRScs betas for 8 traits
beta.sh # generate beta.pbs and submit for job running
beta.pbs # job to be submitted
score.job # parallele prs calculations

hr_trait_prs.R # prediction performance of 9 trait prs for CAD

2. vis
Contains scripts for figure generations.
The name is correspondant to figure tag in the article.
## Scripts
Figure2_corr_traits_prs.R # Figure 2. Correlations between phenotypes and PD-PRSs. 
Figure3_subgroup_stack.R # Figure 3. PD-PRS profiles in different pleiotropy subgroups.
Figure4_subgroup_relative_change.R # Figure 4. Relative changes of quantitative traits in pleiotropy subgroups.
Figure5_SuppFigure3_interaction_subgroup.R # Figure 5. Interactions between PD-PRS subgroup and blood pressures, smoking status.
# and Supplementary Figure 3.  Interactions between PD-PRS subgroup and triglycerides, sleep duration, Vitamin D, FVC, and FEV1
#### Supplementary figures
SuppFigure1_corr_traits.R # Supplementary Figure 1. Genetic correlations between traits. 
SuppFigure2_HR.R # Supplementary Figure 2. Hazards ratios for CAD of overall CAD PRS, 8 PD-PRSs, and 1 NS-PRS.
SuppFigure4_corr_pdprs.R # Supplementary Figure 4. Correlations between 8 PD-PRS and 1 NS-PRS.









