Scripts for 2018-2020 paired eDNA/Beach Seine data
Note: all scripts must be run independantly (i.e., restart to clear environment and packages loaded)

01_data_cleaning:
		- formatting colnames
		- pool PCR replicates
		- remove samples not used in this analyses
		- Check that within site dissimilarity is <0.49 (it checks out - 
		none >0.49) (Kelly 2018)
		- identify minimum sequence depth  (23988) - had to drop one 
		BS paired site: GOG1, sept 2019
			-  discarded sites/samples that aren't used in analyses
		- rarify reads (don't use McKnight et al 2018, reasons why) before 1726 ASVs, 
		24629880 reads; after 954 ASVs 3310344
		- create site-level data from sample data - because it is all duplicated
	IN: 	asv_matrix - sequence_table.12S.merged.w_ASV_names.length_var.merged_dataset.txt
			sample data - Calvert_12S_metadata.csv
			taxa assignment - 12S_ASV_sequences.length_var.blast.out.txt
	OUT:	site data - sitesurvey_data.rds
			sample data - sample_data.rds
			pooled ASVs - pooled_asvMatrix.rds

02a_occupancy_model_formatting:
	Format ASV matrix for occupancy model
		- make data binary; remove 0s
		- make a list of dataframes for each ASV that will each be run through occupancy model	
		- if we only removed singletons 215 ASVs, 3221992
	IN:		pooled ASVs - pooled_asvMatrix.rds
	OUT:	ASVs by sample binary (matrix) - ASV_by_sample.rds
			detection histories - ASVlist.rds
			ordered list of ASVs - ASV.rds
			
02b_occ_model:
	Run Royle-Link occupancy model
	IN:		detection histories - ASVlist.rds
			ordered list of ASVs - ASV.rds
	OUT:	occupancy probabilities by ASV - occProb_royallink.rds
	
02c_occ_formatting:
	Format read count-by-site with occupancy thresholds
		-remove observations with <80% probability of occupancy: 110 ASVs remain, 3078676 reads
	IN:		*requires data from occ_wrangling02a that is not saved as a temporary file - spec (dataframe)
			occupancy probabilities - occProb_royallink.rds
			ASVs by sample (matrix) - ASV_by_sample.rds
	OUT:	occupancy probability by sample - occprob_by_sample.rds
			read count by sample, low occupancy removed - data12se_asvmatrix_lor_12s.rds

02d_occ_plotting:
	Plotting occupancy probability
			- look at distribution of occupancy probabilities	
	IN: 	*requires DF from occ_wrangling02a
			occupancy probabilities by ASV - occProb_royallink.rds
	OUT:	scatter plot of occupancy probabilities positive detections 
			distogram of distribution of occupancy probabilities 

02e_occ_diagnosis:
	This code generates traceplots and density plots calculates as well as summary statistics 
	for a set of ASVs ranging in their number of positive detections (2 - 243) in entire detection histories
	IN:		detection histories - ASVlist.rds
	OUT:	plots of probability of occupancy for differing numbers of detections

03_BS_taxized:
	resolve taxonomy for beach seining species names
	IN:		Hakai fish taxa and codes - fishcodes.csv
	OUT:	Taxonomy resolved - taxa_BS_resolved.rds

04_eDNA_taxa_assignment:
	NOTE: this will now replace assignment_corrections04.R and eDNA_taxized05b_new.R
	- use the top 10 NCBI hits to determine taxa clusters that are indestinguishable at 12S mifish
	Steps:1. Identifies if taxa are native to the Northeastern Pacific
			a. this is done manually at the moment because we are also detecting freshwater 
			species from rivers and estuaries. If we compile a list of probably species, this
			could be replaced. Probably have to be updated every time because of novel results.
		2. Identifies all ASVs and taxa with multiple equal top hits
           		a. If this can be resolved by removing non-natives (from step 1), it is.
           		b. If not, taxa are grouped, within family or genus - note: groups are assigned within a taxonomic level, not across that level, allowing for other distinct species or sub groups within that level to maintain their unique identity
		3. Identifies ASVs and taxa assigned to a single species
           		a. If a species is non-native, the next closest native species is assigned (only within genus) - this only happened on one occasion - I think that if the problem was more prevalent another approach might be needed
           		b. All  taxa that were assigned to groups in Step 2 are assigned to those groups.

	IN: 	ASV by site from occupancy model - data12se_asvmatrix_lor_12s.rds
		site survey data - sitesurvey_data.rds
		top 10 blast hits - 12S_ASV_sequences.length_var.blast.out20221101.txt
	OUT:  	Table of species/species groups (for supplementary figure) - taxonomy_groups_12s_eDNA20221103.csv
	   	ASV taxonomy revised - ASV_taxonomy_eDNA20221103.rds
		LCT by site - LCT_by_site.rds	NOTE: should probably rewrite this when reconciled with BS
		taxonomy_eDNA.rds taxonomy table of of species and their group assignments - taxonomy_eDNA.rds
			-use this to reconcile BS taxonomy, must be done after biomass estimates are made

05_bs_data_formatting:
	Correct for mis-ids in beach seine surveys (clustering species into LITs (in code "LCT") for indistinguishable species.
	Resolve taxonomy for eDNA and beach seining - so both share same LITs
	IN:	beach seine taxa abundance-by-site - FABSMasterData.xlsx(sheet = "abundance")
		Taxonomy resolved - taxa_BS_resolved.rds
		environmental and habitat data - FABSMasterData.xlsx(sheet = "habitat")
	OUT:	taxonomy table for BS - BS_taxonomy_reconciled.rds
		taxonomy table for eDNA - ASV_taxonomy_reconciled.rds
		LCT-by-site matrix (with replicates for Occ. modelling) - bs_specmatrix_retainreplicates.rds
		LCT-by-site matrix  - bs_specmatrix.rds
		resolved taxonomy table for both methods - taxonomy_combined.rds
		
06_length_weight:
	Turn count data into biomass.
	-Using bayesian estimates of mass:length ratios from higher taxonomy, Froese, R., J. Thorson and R.B. Reyes Jr., 2014. A Bayesian approach for estimating length-weight relationships in fishes. J. Appl. Ichthyol. 30(1):78-85.
	-In the field a subsample of fish were measured when >20 were identified.  In these cases, the weight 
	distribution was sampled with replacement equal to the number of observations of that species. This was repeated 1000
	times and the mean and SD taken. Occaisionally, fish could not be measured, but were only counted. In these cases,
	samples were drawn from the global survey-wide distribution for that species.
	Taxa are measured at sites by abundance, weight (grams), as well as abundance/m^2 (seafloor area), abundance/m^3 (seawater volume), weight/m^2 (seafloor area), and weight/m^3 (seawater volume)
	IN:	specandgen_traits_mod.csv
		bs_lengths.csv
		bs_abundance.csv
		names_for_length_1.csv - a manually modified CSV where taxa (Query) is assigned a species to extract traits data (NOTE: do not use as taxonomy table)
		eDNA_long.rds - for filtering to shared surveys
	OUT: 	bs_long.rds - contains all measures of observations (count, weight, and their area/volume densities)

07_environmental_data_prep:
	Format environmental data, calculate temperature gradients and date difference, PCA of seawater mixing
	
	OUT: environmental.rds - formated environmental data
	
10_species_detections_compared:
	Format data to compare read abundance and biomass and biomass densities, plotting
	Plots of richness differences across environmental gradients
	
11_Trees:
	Build a dendrogram from higher taxonomy and colour code by method. 
	Plot frequency of positive detections by method.

12_richness_analyses:
	Plot differences in richness. run models to predict differences in richness

13_spatial_turnover:
	Calculations of turnover, exploring trends, plotting, generalized dissimilarity models

14_test_statistics:
	test statistices presented in paper
	
15_patchiness 
	summarize spatial data of vegetation and exposure. Plot density-histogram showing difference in scales of patchiness.

16_SADs



