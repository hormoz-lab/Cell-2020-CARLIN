# CARLIN Paper (Cell 2020)

This repository provides instructions and code to reproduce all results, numerics and figures from the [CARLIN paper](https://doi.org/10.1016/j.cell.2020.04.048):

> S. Bowling, D. Sritharan, F. G. Osorio, M. Nguyen, P. Cheung, 
A. Rodiguez-Fraticelli, S. Patel, W-C. Yuan, Y. Fujiwara, B. E. Li, S. H. Orkin, 
S. Hormoz, F. D. Camargo. "An engineered CRISPR/Cas9 mouse line for 
simultaneous readout of lineage histories and gene expression profiles 
in single cells." Cell (2020), https://doi.org/10.1016/j.cell.2020.04.048 
>

The instructions will walk you through how to (i) download relevant software and data, (ii) constitute processed files from raw data, (iii) invoke the CARLIN pipeline and other software to produce results, and finally (iv) regenerate figures and reproduce numerics from the paper.

For convenience, there are six main paths referred to below that appear in multiple steps:

	CODE_PATH	- directory where the CARLIN pipeline code is cloned
	PAPER_PATH	- directory where this repository is cloned
	RAW_PATH  	- location of unpaired Illumina sequencing data
	PROCESSED_PATH	- location of data processed by PEAR and CellRanger
	ANALYSIS_PATH	- output from running the CARLIN pipeline on amplicon data, and Seurat on transcriptome data
	RESULTS_PATH	- output from custom scripts for manuscript preparation

These directories mirror my own workflow when working on this paper, and make it easy to separate data from code, and checkpoint results.

## Download Software

1. If you haven't done so, first [download and install the CARLIN pipeline](https://gitlab.com/hormozlab/carlin) according to the instructions, making sure to update the MATLAB path as instructed. The results in the manuscript match commit `f20a032f41d57a7a96950c86d71eea0250540894` so be sure to use that version (via `git checkout`), if you want to match to reproduce numerics exactly.

2. Download the code from this repository into the directory given by PAPER_PATH:

	```bash
	$ git clone https://gitlab.com/hormozlab/cell_2020_carlin.git ${PAPER_PATH}
	```

3. [Download SRA-Tools](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit) (to extract sequencing data from NCBI SRA) into SRA_PATH. Follow the [installation](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) and [configuration](https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration) steps.

4. [Download PEAR V0.9.11](https://www.h-its.org/downloads/pear-academic) (the paired-end read merger we used in the paper) into PEAR_PATH. 

5. [Download CatchAll 4.0](http://www.northeastern.edu/catchall/downloads.html) for fitting allele abundance and making the bank, and copy the full path to where you downloaded CatchAll (including binary name) into:

		CODE_PATH/@Bank/CatchAllPath.txt

6. [Download Seurat V3.1.2](https://satijalab.org/seurat/install.html) for integrating single-cell transcriptomic datasets and performing differential gene expression.

## Download Data

The data for this study corresponds to NCBI [GEO GSE146792](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146972)/[Bioproject PRJNA612835](https://www.ncbi.nlm.nih.gov//bioproject/PRJNA612835)/[SRA SRP252955](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP252955). 

7. Download all single-cell count matrix files into PROCESSED_PATH. These files were generated using CellRanger (`cellranger mkfastq` and `cellranger count`).

	```bash
	$ mkdir -p ${PROCESSED_PATH}
	$ cd ${PROCESSED_PATH}
	$ wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146972/suppl/GSE146972_RAW.tar
	$ tar -xvf GSE146972_RAW.tar
	$ rm GSE146972_RAW.tar
	```

8. Retrieve the amplicon sequencing data corresponding to the single-cell runs, using the SRA-tools fastq-dump utility. These files were generated using CellRanger (`cellranger mkfastq`),

	```bash
	$ cd ${PROCESSED_PATH}
	$ grep sc10x ${PAPER_PATH}/sample_sheet.csv | cut -d ',' -f11 | grep -v ';' > SampleList.txt	
	$ ./${SRA_PATH}/bin/fastq-dump $(<SampleList.txt) -O . --split-files --gzip	
	```

9. Archiving the data in NCBI flattens the descriptively named sub-directory structure and renames some files. The ugly commands below are just to recreate the expected hierarchy and recover the original file names, which the downstream code expects.

	```bash
	$ cd ${PROCESSED_PATH}
	$ find -name "*.h5" | sed "s/^.*\(Processed.*Transcriptome\).*$/\1/" | sed "s:_:/:g" | sed "s@Processed/@@g" | xargs -n1 mkdir -p
	$ find -name "*.h5" | sed "s/^.*\(Processed.*Transcriptome\).*$/\1/" | sed "s:_:/:g" | sed "s@Processed/@@g" | sed "s@Transcriptome@Amplicon@g" | xargs -n1 mkdir -p
	$ find -name "*.h5" | xargs -I % sh -c 'echo %; echo % | sed "s/^.*\(Processed.*\).*$/\1/" | sed -e :1 -e "s@\(.*\)_\(.*filtered\)@\1/\2@;t1" | sed 's@Processed/@@g' ' | xargs -n2 mv
	$ grep -f SampleList.txt ${PAPER_PATH}/sample_sheet.csv | grep -v ';' | cut -d ',' -f6,8,11 | sed 's%/Replicate.*,%,%' | xargs -I % sh -c ' echo % | sed "s@.*,@@" | sed "s/$/_1.fastq.gz/"; echo % | sed "s@\([^,]*\),\([^,]*\),\([^,]*\)@\2/\1_R1_001.fastq.gz@g" ' | xargs -n2 mv
	$ grep -f SampleList.txt ${PAPER_PATH}/sample_sheet.csv | grep -v ';' | cut -d ',' -f6,8,11 | sed 's%/Replicate.*,%,%' | xargs -I % sh -c ' echo % | sed "s@.*,@@" | sed "s/$/_2.fastq.gz/"; echo % | sed "s@\([^,]*\),\([^,]*\),\([^,]*\)@\2/\1_R2_001.fastq.gz@g" ' | xargs -n2 mv
	$ rm SampleList.txt
	```

10. Retrieve the unpaired single-end bulk amplicon sequencing data into RAW_PATH using the SRA-tools fastq-dump utility:

    ```bash
	$ mkdir -p ${RAW_PATH}
	$ cd ${RAW_PATH}
	$ grep Bulk ${PAPER_PATH}/sample_sheet.csv | cut -d ',' -f11 | grep -v ';' > SampleList.txt
	$ ./${SRA_PATH}/bin/fastq-dump $(<SampleList.txt) -O . --split-files --gzip
	$ rm SampleList.txt
	```

## Preprocess Data	

11. Merge the single-end reads from the bulk amplicon sequencing runs into paired-end reads using PEAR:

    ```bash
	$ ${PAPER_PATH}/process_paired_end_reads.sh ${PEAR_PATH}/bin/pear ${PAPER_PATH}/sample_sheet.csv ${RAW_PATH} ${PROCESSED_PATH} -SRA
    ```

## Run CARLIN Pipeline	

12. Temporarily add PAPER_PATH to MATLAB's search path.

	```MATLAB
	>> addpath(genpath(PAPER_PATH));
	```

13. Process all the bulk amplicon data from the various experiments:

    ```MATLAB
	>> run_paper_experiments('ChronicInduction',    PROCESSED_PATH, ANALYSIS_PATH);			% Related to Figure 1
	>> run_paper_experiments('PulsedInduction',     PROCESSED_PATH, ANALYSIS_PATH);			% Related to Figure 2
	>> run_paper_experiments('Sanger',              PROCESSED_PATH, ANALYSIS_PATH);			% Related to Figure 2
	>> run_paper_experiments('InVitroPhylogeny',    PROCESSED_PATH, ANALYSIS_PATH);			% Related to Figure 2
	>> run_paper_experiments('TissuePanel',         PROCESSED_PATH, ANALYSIS_PATH);			% Related to Figure 3
	>> run_paper_experiments('Diversity', 	        PROCESSED_PATH, ANALYSIS_PATH);			% Related to Figure 3
	>> run_paper_experiments('Inducible',           PROCESSED_PATH, ANALYSIS_PATH);			% Related to Figure 3
	>> run_paper_experiments('InVivoPhylogeny',     PROCESSED_PATH, ANALYSIS_PATH);			% Related to Figure 4
    ```

14. To process single-cell data, first filter transcriptome data to generate reference barcode lists:

    ```MATLAB
	>> run_paper_experiments('Transcriptome', PROCESSED_PATH, ANALYSIS_PATH);
    ```

15. Then process single-cell amplicon data:

    ```MATLAB
	>> run_paper_experiments('EB',	        PROCESSED_PATH, ANALYSIS_PATH);					% Related to Figure 5
	>> run_paper_experiments('5FU',	        PROCESSED_PATH, ANALYSIS_PATH);					% Related to Figure 6
	>> run_paper_experiments('SCReplicate', PROCESSED_PATH, ANALYSIS_PATH);					% Related to Figure 6
    ```

## Custom Analysis for Paper

16. There are some single-cell samples that need to be pooled together (e.g. mouse bones, replicates)

    ```MATLAB
	>> pool_samples(ANALYSIS_PATH);
    ```

17. Create the CARLIN allele bank (and banks for negative controls):

    ```MATLAB
	>> create_allele_banks(ANALYSIS_PATH, RESULTS_PATH);
    ```

18. Run some simulations pertaining to CARLIN allele diversity:

    ```MATLAB
	>> simulate_multimapped_alleles(RESULTS_PATH);
	>> simulate_diversity_replicates(ANALYSIS_PATH, RESULTS_PATH);
    ```

19. Run simulations for in vitro and in vivo lineage reconstruction:

    ```MATLAB
	>> simulate_invitro_phylogeny(ANALYSIS_PATH, RESULTS_PATH);
	>> simulate_invivo_phylogeny(ANALYSIS_PATH, RESULTS_PATH);
    ```

20. Align the datasets from single-cell experiments, generate clusters and look at marker genes using Seurat:
	
    ```R
	> setwd(PAPER_PATH)
    > source('Analysis_EB.R')
	> Analysis_EB(PROCESSED_PATH, ANALYSIS_PATH, RESULTS_PATH)
	> source('Analysis_5FU.R')
	> Analysis_5FU(PROCESSED_PATH, ANALYSIS_PATH, RESULTS_PATH)
    ```

21. Overlay CARLIN alleles on to transcriptome for single-cell experiments:
	
    ```MATLAB
	>> analyze_EB_experiment(ANALYSIS_PATH, RESULTS_PATH);
	>> analyze_5FU_experiment(ANALYSIS_PATH, RESULTS_PATH);
	```

22. Do differential gene expression for 5-FU experiment using Seurat:

	```R
	> setwd(PAPER_PATH)
    > source('DGE_5FU.R')
	> DGE_5FU(ANALYSIS_PATH, RESULTS_PATH)
	```

23. Compute statistics for single-cell experiments:

	```MATLAB
	>> generate_SC_stats(RESULTS_PATH);
    ```

## Generate Outputs

24. Prepare main and supplementary figures:

    ```MATLAB
	>> make_schematic_subplots(RESULTS_PATH);							% Figure 1B
	>> make_chronic_induction_subplots(ANALYSIS_PATH, RESULTS_PATH);	% Figures 1C-H, S1C-F
	>> make_pulsed_induction_subplots(ANALYSIS_PATH, RESULTS_PATH);		% Figures 2A, S1G
	>> make_invitro_phylogeny_subplots(RESULTS_PATH);					% Figures 2CD, S4C
	>> make_tissue_panel_subplots(ANALYSIS_PATH, RESULTS_PATH); 		% Figures 3B, S3AB
	>> make_diversity_subplots(RESULTS_PATH);							% Figures 3C-I, S3C-E
	>> make_inducible_subplots(ANALYSIS_PATH, RESULTS_PATH);			% Figure S3F
	>> make_invivo_phylogeny_subplots(RESULTS_PATH);					% Figures 4, S4D
	>> make_EB_subplots(RESULTS_PATH);									% Figures 5, S5
	>> make_5FU_subplots(RESULTS_PATH);									% Figures 6, S6
	>> make_algorithm_subplots(RESULTS_PATH);							% Figure S2A-C
	>> make_alignment_subplots(ANALYSIS_PATH, RESULTS_PATH);			% Figure S2D
	>> make_toy_model_subplots(RESULTS_PATH);							% Figure S4B
    ```

25. Make supplementary tables:

    ```MATLAB
	>> generate_supplemental_tables(ANALYSIS_PATH, RESULTS_PATH);
    ```

26. Generate all the statistics and numerical quantities found in the manuscript text:

    ```MATLAB
	>> generate_manuscript_stats(PROCESSED_PATH, ANALYSIS_PATH, RESULTS_PATH);
    ```

## Miscellany

- Data that doesn't have the right format to be included in NCBI databases but is needed to reproduce some results (Sanger sequencing data used in Figure 2CD and fragment analyzer data used in Figure S2D) can be found [here](https://gitlab.com/hormozlab/cell_2020_carlin/-/tree/master/data).
- The figures' schematic panels and graphical abstract can be found [here](https://gitlab.com/hormozlab/cell_2020_carlin/-/tree/master/figures/schematics). 
- The source files for the fragment analyzer traces shown in Figure S1A can be found [here](https://gitlab.com/hormozlab/cell_2020_carlin/-/tree/master/figures/fa_trace).
- The source files for the FACS plots shown in Figure S5AB can be found [here](https://gitlab.com/hormozlab/cell_2020_carlin/-/tree/master/figures/facs).

#### Prepared By: Duluxan Sritharan
