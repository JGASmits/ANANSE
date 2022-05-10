import ananse
import glob
import os
import subprocess as sp
from gimmemotifs.scanner import scan_regionfile_to_table
from loguru import logger
import pandas as pd
from ananse.peakpredictor import predict_peaks
from ananse.utils import check_path
import ananse.network


work_dir = '/ceph/rimlsfnwi/data/moldevbio/zhou/jsmits/Publications/scANANSE/'
output_dir = f'{work_dir}analysis'
#lets loop over the single cell directory
data_dir = f'{work_dir}ANANSE/tests/data/GRCh38.p13_chr11/single_cell/'
Bamfiles = glob.glob(f"{data_dir}*.bam")
TPMfiles = glob.glob(f"{data_dir}*.tsv")
genome = f'{work_dir}ANANSE/tests/data/GRCh38.p13_chr11/GRCh38.p13/GRCh38.p13.fa'
ncpus = 12
pfm_file = f'{work_dir}ANANSE/tests/data/GRCh38.p13_chr11/gimme.vertebrate.v5.0_chrom11.pfm'
factors_file = f'{work_dir}ANANSE/tests/data/GRCh38.p13_chr11/gimme.vertebrate.v5.0.factors_chrom11.txt'
region_file = f'{work_dir}ANANSE/tests/data/GRCh38.p13_chr11/regions.txt'
region_bed_file = f'{work_dir}ANANSE/tests/data/GRCh38.p13_chr11/regions.bed'


logger.info("Scanning regions for motifs.")
gimmescanfile = f'{output_dir}/all_peaks_scanned.pfmscores'


if not os.path.exists(f"{gimmescanfile}"): 
    motif_df = scan_regionfile_to_table(
        region_file, genome, "score", ncpus=ncpus,
    )
    motif_df.to_csv(gimmescanfile,index=True, header=True,sep = '\t')
else:
    motif_df = pd.read_table(gimmescanfile,
     sep="\t", header = 0, comment="#", index_col = 0)

#motif_df_working = pd.read_table("/ceph/rimlsfnwi/data/moldevbio/zhou/jsmits/Ananse_test_data/analysis/peakpred_V10/LSC_peaks_scanned.pfmscores",
#     sep="\t", header = 0, comment="#", index_col = 0)

for bamfile in Bamfiles:
    sample_name = os.path.basename(bamfile)[:-4]
    logger.info(sample_name)
    ATAC_list = []
    ATAC_list.append(bamfile)
    if not os.path.exists( f'{output_dir}/{sample_name}_binding.h5'): 
        logger.info(f'Running binding for {sample_name}')
        peak_output = predict_peaks(output_dir, 
        atac_bams=ATAC_list, 
        regions=region_file,
        genome=genome,
        #factors = factors_file,
        #pfmfile=pfm_file,
        pfmscorefile=gimmescanfile, 
        jaccard_cutoff=0.0, 
        ncore=4)
        os.rename(f'{output_dir}/binding.h5', f'{output_dir}/{sample_name}_binding.h5')
    else:
        logger.info(f'h5 file found for {sample_name}')
    
    logger.info(f'running ananse network {sample_name}')
    ananse.network.Network(
        genome=genome,  # checked in CLI
        gene_bed=check_path(args.annotation),

        include_promoter=True,
        include_enhancer=True,
        full_output=True,
    ).run_network(
        binding=f'{output_dir}/{sample_name}_binding.h5',
        fin_expression=f"{data_dir}{sample_name}.tsv",
        outfile=f'{output_dir}/{sample_name}_network.tsv',
    )
#all = 7 min for file C all TFs
#all = 7 min for file C all TFs
