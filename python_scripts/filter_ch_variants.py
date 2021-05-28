# conda activate clonalh

# qmap submit chvariants_filter_vep.qmap

import pandas as pd
import os
import copy
import click
import subprocess

@click.command()

@click.option('--input_maf',
              '-i',
              required = True,
              help="path to gvcf file from uk biobank")
@click.option('--output_folder',
              '-o',
              required = True,
              help="path to output folder")

def filter_chgenes(input_maf, output_folder):
    """
    Starting from MAF files from UKB 50K (gVCF filterd by CH gene coordinates)
    filters and arranges the data: select real variants, split variants with 2+ ALTs,
    annotate the data, etc AND run VEP
    """

    ### 1. Open MAF file
    maf = pd.read_csv(input_maf, sep= '\t')


    ### 2. Select only rows with real variants from MAFs
    # Select only variants NOT having "<NON_REF>" in "ALT"
    maf_var = maf[~maf["ALT"].str.startswith("<NON_REF>")]


    ### 3. Rearrange info from last column
    # Split column and add column names
    split_ukbcol = maf_var[maf_var.columns[-1]].str.split(":",expand=True)
    split_ukbcol.columns = ["GT","AD","DP","GQ","PGT","PID","PL","SB"]
    # Reorder columns 4-7
    for i in range(split_ukbcol.shape[0]):
        if split_ukbcol.iloc[i, 6] is None:
            split_ukbcol.iloc[i, 6] = split_ukbcol.iloc[i, 4]
            split_ukbcol.iloc[i, 4] = None
        if split_ukbcol.iloc[i, 7] is None:
            split_ukbcol.iloc[i, 7] = split_ukbcol.iloc[i, 5]
            split_ukbcol.iloc[i, 5] = None
    # Select columns from original maf and merge with splitted columns
    maf_var_spl = pd.concat([maf_var.iloc[:,[0,1,3,4,5,6,7]],
                             split_ukbcol, 
                             maf_var.iloc[:,8]], axis=1)


    ### 4. Eliminate "<NON_REF>" from ALT column
    # Define function
    def eliminate_nonref(alt):
        """"
        Eliminates "<NON_REF>" from ALT
        Initially checks that "<NON_REF>" is present
        """
        if alt.split(",")[-1] == "<NON_REF>":
            wo_nonref = ",".join(alt.split(",")[:-1])
            return wo_nonref
        else:
            return alt
    # Eliminate "<NON_REF>"
    maf_var_spl["ALT"] = maf_var_spl["ALT"].apply(lambda x: eliminate_nonref(x))
    
    
    ### 5. Eliminate variants with AD=None, DP=None, or DP=0
    maf_var_spl = maf_var_spl[~((maf_var_spl['AD'].isnull()) |
                            (maf_var_spl['DP'].isnull()) |
                            (maf_var_spl['AD'].str.split(',', expand=True)[0].isnull()) |
                            (maf_var_spl['AD'].str.split(',', expand=True)[1].isnull()) |
                            (maf_var_spl['DP']=="0"))]


    ### 6. Split variants with >1 ALT & Calculate VAF

    ### 6A. Variants 1 ALT
    # Select variants with 1 ALT
    df1 = maf_var_spl[~maf_var_spl['ALT'].str.contains(',')]
    # Calculate VAF (& add number of ALT=1)
    df1['AD_alt'] = df1['AD'].str.split(',', expand=True)[1].astype(int)
    df1['VAF_alt'] = df1['AD'].str.split(',', expand=True)[1].astype(int) / df1['DP'].astype(int)
    df1['VAF_ref'] = df1['AD'].str.split(',', expand=True)[0].astype(int) / df1['DP'].astype(int)
    df1['ALT_num'] = 1

    ### 6B. Variants 2+ ALT
    # Select variants with 2+ ALT & transform df to list
    df2 = maf_var_spl[maf_var_spl['ALT'].str.contains(',')]
    df2_list = df2.values.tolist()
    # Divide variants in 1 line per ALT & calculate VAF
    df2_newlist = []
    # Loop through all variants
    for row in df2_list:
        # Extract info of ALTs and VAF
        ALTs = row[4].split(',')
        ADs = row[8].split(',')
        DP = row[9]
        # Loop to create 1 line per ALT
        for i in range(0,len(ALTs)):
            newrow = copy.deepcopy(row)
            # Take ALT and substitute column        
            newrow[4] = ALTs[i]
            # Take AD corresponding ALT in new column
            newrow.append([])
            newrow[16] = int(ADs[i+1])
            # Calculate VAF from corresponding ALT in new column
            newrow.append([])
            newrow[17] = int(ADs[i+1]) / int(DP)
            # Calculate VAF from REF in new column
            newrow.append([])
            newrow[18] = int(ADs[0]) / int(DP)
            # Annotate number of total ALT in new column
            newrow.append([])
            newrow[19] = len(ALTs)
            # Append variant to new list
            df2_newlist.append(newrow)       
    # Transform to dataframe and change column names
    df2_newdf = pd.DataFrame(df2_newlist)
    df2_newdf.columns = df1.columns

    ### 6C. Concatenate variants 1ALT and 2+ALT
    maf_var_spl_1alt = pd.concat([df1, df2_newdf], ignore_index=True)


    ### 7. Annotate variant type
    # Define function
    def variant_type(var):
        if (len(var[3]) == len(var[4])) & (len(var[3]) == 1):
            return 'SNV'
        elif len(var[3]) != len(var[4]):
            return 'Indel'
        elif (len(var[3]) == len(var[4])) & (len(var[3]) > 1):
            return 'MNV'
        else:
            return 'Unknown'
    # Add type of variant
    maf_var_spl_1alt['var_type'] = maf_var_spl_1alt.apply(lambda x: variant_type(x), axis=1)


    ### 8. Reorder columns and create new column for VEP output
    maf_var_spl_1alt["VEP"] = ""
    maf_var_spl_1alt = maf_var_spl_1alt[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'VEP',
                                         'AD_alt', 'DP', 'VAF_alt', 'VAF_ref', 'ALT_num', 'var_type',
                                         'GT', 'AD', 'GQ', 'PGT', 'PID', 'PL', 'SB', 'INFO']]
    
    
    ### 9. Save temporary file
    # Get file name without extension form initial path and create temp output folder+file
    file = input_maf.split("/")[-1].split(".")[0]
    temp_file = output_folder + "temp/" + file + ".maf.filt"
    # Save
    maf_var_spl_1alt.to_csv(temp_file, sep="\t", index=False, compression='gzip')
    
    
    ### 10. Run VEP in the terminal
    # Get output file name (path+file)
    final_file = output_folder + file + ".maf.vep"
    # Run VEP with singularity
    subprocess.run(
        'singularity exec /home/sdemajo/vep.simg vep -i ' + temp_file + ' -o ' + final_file + ' --assembly GRCh38 --no_stats --cache --offline --symbol --protein --vcf --canonical --af_gnomad --dir /workspace/datasets/vep --compress_output gzip', shell=True, executable='/bin/bash')
    
if __name__ == "__main__":
    filter_chgenes()

