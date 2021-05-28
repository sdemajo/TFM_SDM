import pandas as pd
import os
from io import StringIO
import gzip
from tqdm import tqdm
import click
import warnings
warnings.filterwarnings('ignore')


@click.command()

@click.option('--first_file',
              '-first',
              required = True,
              help="first file to include in the loop")
@click.option('--last_file',
              '-last',
              required = True,
              help="last file to include in the loop")

def identify_chvariants(first_file, last_file):

    # List MAF fiels
    maf_files = os.listdir('/workspace/projects/clonalhemato_ukb/ch_variants_50k_mafvep')
    maf_files = [file for file in maf_files if '.maf.vep2.gz' in file]

    # Create empty dataframes with the same columns as the MAFs
    maf = pd.read_csv("/workspace/projects/clonalhemato_ukb/ch_variants_50k_mafvep/4033478_23161_0_0_ch.maf.vep2.gz", sep= '\t')
    CH_mutations = maf[maf['SYMBOL']== 'xxxxxx']
    CH_mutations['patient'] = [] # add column to incorporate patient code

    # Loop the selected files
    for file in tqdm(maf_files[int(first_file):int(last_file)]):
        filename = '/workspace/projects/clonalhemato_ukb/ch_variants_50k_mafvep/' + file
        maf = pd.read_csv(filename, sep= '\t')

        # Filters 1,2,3: Select CH mutations by AD, VAF, & coding region
        maf_CH = maf[(maf['AD_alt'] >= 2) & (maf['VAF_alt'] <= 0.3333) & (maf['Prot_pos'].notnull())]

        # Filter 4: Select by consequence severity
        consequence_select = "transcript_ablation|splice_acceptor_variant|splice_donor_variant|stop_gained|frameshift_variant|stop_lost|start_lost|transcript_amplification|inframe_insertion|inframe_deletion|missense_variant|protein_altering_variant|splice_region_variant|incomplete_terminal_codon_variant|start_retained_variant|stop_retained_variant"
        maf_CH = maf_CH[maf_CH['Consequence'].str.contains(consequence_select, na=False)]

        # Concatenate CH mutations with patient code
        if len(maf_CH) > 0:
            patient_code = file.split('_')[0]
            maf_CH['patient'] = patient_code
            CH_mutations = pd.concat([CH_mutations, maf_CH])


    # Save
    output_file = "CHmutations_" + first_file + "_" + last_file + ".txt.gz"
    CH_mutations.to_csv('/workspace/projects/clonalhemato_ukb/analysis_50k_202103/identify_ch_202104/results/' + output_file,
                        sep="\t", index=False, compression='gzip')

    
if __name__ == "__main__":
    identify_chvariants()
      