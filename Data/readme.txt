Readme
Primary preprocessed and processed datasets.

Data:

	Preprocessed data:
	
	eDNA
	
		Sequences and bioinformatics:
	
			Sequences amplified from water samples (fastq files) are available at NCBI under BioProject ID PRJNA1051585: BioSamples: 39222712 - 39222995 
			Data supplied by the author: Bioinformatics scripts are available at: https://zenodo.org/records/10183046.

		post-bioinformatics:
	
			asv_matrix - Data\2022_10_31\sequence_table.12S.merged.w_ASV_names.length_var.merged_dataset.txt
			sample data - Data\2022_10_31\Calvert_12S_metadata.csv
			taxa assignment - Data\2022_10_31\12S_ASV_sequences.length_var.blast.out.txt
		
	
		beach seine
		
			beach seine capture data (long format) - Data\2022_10_31\bs_data\bs_abundance.csv
			habitat data - Data\2022_10_31\bs_data\bs_habitat.csv
		
	
	Processed data (used in analysis):
	
		beach seine taxa observations cleaned and processed (long format) - Data/2022_10_31/derived_data/bs_long.rds
		eDNA taxa observations cleaned and processed (long format) - Data/2022_10_31/derived_data/eDNA_long.rds
		eDNA and beach seine taxa observations (long format) - Data/2022_10_31/derived_data/shared_long.rds
		
		resolved taxonomy table for both methods - Data/2022_10_31/derived_data/taxonomy_combined.rds
		
		survey dates,and environmental data - Data/2022_10_31/derived_data/environmental.rds
		
		habitat preferences of detected fish taxa - Figures_tables/LCT_table_habitat_preferences.csv
		
		seagrass patchs for BC coast - Data/2022_10_31/spatial/seagrass_dims.csv
		giant kelp patchs for BC coast - Data/2022_10_31/spatial/macro_dims.csv
		exposure patchs for BC coast - Data/2022_10_31/spatial/trial.csv