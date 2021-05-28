# conda activate clonalh

# qmap submit chvariants_filter_vep.qmap

import pandas as pd
import os
from io import StringIO
import gzip
import click


@click.command()

@click.option('--input_maf',
              '-i',
              required = True,
              help="path to gvcf file from uk biobank")
@click.option('--output_folder',
              '-o',
              required = True,
              help="path to output folder")

def filter_chvariants(input_maf, output_folder):
    """
    Filter MAF files (with VEP annotation) from UK Biobank:
    filter by ALT, select canonical transcript, rearrange
    """

    ### 1. GET PANDAS DATA FRAME
    # Transform file to pandas dataframe eliminating comment lines
    lines = ''.join([line for line in gzip.open(input_maf, 'rt') if not line.startswith("##")])
    maf = pd.read_csv(StringIO(lines), sep= '\t')


    ### 2. FILTER

    ### Filter by ALT
    # Eliminate variants with >2 ALT
    maf = maf[~(maf["ALT_num"] > 2)]
    # Eliminate variants with 2 ALT that are indels
    maf = maf[~((maf["ALT_num"] == 2) & (maf["var_type"] == "Indel"))]
    # Eliminate variants with AD ALT = 0
    maf = maf[maf["AD_alt"] > 0]

    ### Eliminate variants w/o VEP
    maf = maf[maf['VEP'] != "."]

    ### Reset indexes
    maf.reset_index(drop=True, inplace=True)


    ### 3. SELECT AND FILTER CANONICAL TRANSCRIPT (VEP)

    ## 3A. Get CH canonical transcript codes
    # Read file
    ch_ctrans_coord = pd.read_csv("/workspace/projects/clonalhemato_ukb/analysis_50k_202103/CH_genes_canonical_transcripts_202103/results/CH_genes_canonical_transcripts_ENST.tsv",
                                  sep="\t", header = None)
    # Transform to list
    ENST = list(ch_ctrans_coord[0])

    ## 3B. Select VEP info corresponing to the canonical transcripts
    # Define function
    def get_canonical(vep):
        """
        From the VEP annotation, get the info corresponing to the canonical transcript
        (if there is any VEP annotation)
        """
        VEP_annotations = vep.split(',')
        # Read all VEP annotations (different transcripts)
        for annotation in VEP_annotations:
            # Chech that there is VEP annotation
            if len(annotation) > 1:
                # Select by canonical transcript code
                if annotation.split('|')[6] in ENST:
                    return annotation
    # Create new column with canonical transcript VEP info
    maf['VEP_cano'] = maf['VEP'].apply(lambda x: get_canonical(x))

    ## 3C. Eliminate variants w/o canonical transcript
    maf = maf.dropna(axis=0, subset=['VEP_cano'])


    ### 4. EXTRACT VEP INFO AND REARRANGE DATAFRAME

    ## Get VEP info
    maf['Consequence'] = maf['VEP_cano'].apply(lambda x: x.split('|')[1])
    maf['IMPACT'] = maf['VEP_cano'].apply(lambda x: x.split('|')[2])
    maf['SYMBOL'] = maf['VEP_cano'].apply(lambda x: x.split('|')[3])
    maf['ENST'] = maf['VEP_cano'].apply(lambda x: x.split('|')[6])
    maf['CDS_pos'] = maf['VEP_cano'].apply(lambda x: x.split('|')[13])
    maf['Prot_pos'] = maf['VEP_cano'].apply(lambda x: x.split('|')[14])
    maf['AA_change'] = maf['VEP_cano'].apply(lambda x: x.split('|')[15])
    maf['rs_ID'] = maf['VEP_cano'].apply(lambda x: x.split('|')[17])
    maf['ENSP'] = maf['VEP_cano'].apply(lambda x: x.split('|')[24])
    maf['gnomAD_AF'] = maf['VEP_cano'].apply(lambda x: x.split('|')[25])
    maf['gnomAD_AF_NFE'] = maf['VEP_cano'].apply(lambda x: x.split('|')[31])

    ## Rearrange dataframe
    maf = maf[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'AD_alt',
               'DP', 'VAF_alt', 'VAF_ref', 'ALT_num', 'var_type', 'Consequence', 'IMPACT',
               'SYMBOL', 'ENST', 'CDS_pos', 'Prot_pos', 'AA_change', 'rs_ID', 'ENSP', 
               'gnomAD_AF', 'gnomAD_AF_NFE', 'GT', 'AD', 'GQ',
               'PGT', 'PID', 'PL', 'SB', 'INFO', 'VEP_cano', 'FILTER', 'VEP']]
    #maf = maf.rename(columns={"#CHROM" : "CHROM"})
    
    ### 5. Save file
    # Get file name without extension form initial path and create temp output folder+file
    file = input_maf.split("/")[-1].split(".")[0]
    file_folder = output_folder + file + ".maf.vep2.gz"
    # Save
    maf.to_csv(file_folder, sep="\t", index=False, compression='gzip')

if __name__ == "__main__":
    filter_chvariants()
