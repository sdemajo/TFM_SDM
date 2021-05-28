import gzip
from io import StringIO
import pandas as pd
import pybedtools
import click


@click.command()

@click.option('--input_gvcf',
              '-i',
              required = True,
              help="path to gvcf file from uk biobank")
@click.option('--input_coordinates',
              '-i_coor',
              required = True,
              help="path to bed file with gene coordinates to extract")
@click.option('--output_folder',
              '-o',
              required = True,
              help="path to output folder")

def extract_chgenes(input_gvcf, input_coordinates, output_folder):
    """
    Processes gVCF files from UK Biobank and extracts variants corresponding to 
    the coordinates of CH genes
    """
       
    ## Open gVCF and transform to pandas dataframe eliminating comment lines
    # Extract lines that are not comments lines (##)
    lines = ''.join([line for line in gzip.open(input_gvcf, 'rt') if not line.startswith("##")])
    # Transform to pandas dataframe
    gvcf = pd.read_csv(StringIO(lines), sep= '\t')
    
    ## Transform to BedTool object using pybedtools
    # Duplicate genomic position column
    gvcf.insert(loc = 2, column = 'POS2', value = gvcf["POS"])
    # Transform to bedtool object
    gvcf_bed = pybedtools.BedTool.from_dataframe(gvcf)

    ## Open CH gene coordinate data (bed file)
    ch_ctrans_coord_bed = pybedtools.BedTool(input_coordinates)

    ## Intersect using "u=True" to obtain all rows in gvcf_bed that overlap with ch_ctrans_coord_bed
    vcf_chgenes = gvcf_bed.intersect(ch_ctrans_coord_bed, u=True)

    ## Transform to pandas data frame and recover column names
    vcf_chgenes_df = vcf_chgenes.to_dataframe()
    vcf_chgenes_df.columns = gvcf.columns

    ## Save file
    # Get file name without extension form initial path
    file = input_gvcf.split("/")[-1].split(".")[0]
    # Save
    vcf_chgenes_df.to_csv(output_folder + "/" + file + "_ch.maf.gz", sep="\t", index=False, compression='gzip')
    
if __name__ == "__main__":
    extract_chgenes()
