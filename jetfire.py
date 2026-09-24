### Project Jetfire ###
# Creator - L. Kraesing

# Import needed libraries
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import numpy as np

# Set seaborn configurations
sns.set_theme()
sns.set_style("white")
sns.set_palette("tab20")

# Functions for running

def change_dir():
    """ Change directory to current working directory """
    os.chdir(os.getcwd())

#os.chdir(r"Y:\BSSH_project\DRAGEN TruSight Oncology 500 v2_6_2_4_ds.ee24f1e982424c32861d4815a86b17ec\Logs_Intermediates\Gis\FFPE_colon_tissue1_rep2_dna")

#os.chdir(r"Y:\BSSH_project\DRAGEN TruSight Oncology 500 v2_6_2_4_ds.ee24f1e982424c32861d4815a86b17ec\Logs_Intermediates")

def create_path_files():
    """ Get full path to b-allele files for data ingestion
    :return List of b-allele file paths"""

    b_allele_files = []

    for root, dirs, files in os.walk(os.getcwd()):
        for file in files:
            if file.endswith("_bAllele.tsv"):

                path_file_name = root + os.sep + file
                b_allele_files.append(path_file_name)

    return b_allele_files

def create_df_for_plotting(filename, id):
    """ Creates dataframe with desired formatting for plotting
    :param filename: full path to b-allele file
    :param id: id of sample
    :return: dataframe with desired formatting
    """
    df = pd.read_csv(filename, sep="\t", names=["Chromosome", "Position", "RS_id", "BAF"])
    df["Global_position"] = df["Chromosome"] + ":" + df["Position"].astype("str")
    df["sample"] = id
    return df

def labeling_details(df):
    """ Labeling details for chromosome positions in the plot
    :param df: df created from create_df_for_plotting
    :return: dict with labeling details

    """
    list_with_chromosomes = list(df["Chromosome"].unique())
    list_of_positions = []
    label_positions = []

    for chr in list_with_chromosomes:
        len_ = df.loc[df["Chromosome"] == chr].__len__()
        list_of_positions.append(len_)

    for idx in range(len(list_of_positions)):
            label_pos = list_of_positions[idx] // 2 + sum(list_of_positions[0:idx])
            label_positions.append(label_pos)

    list_of_lengths = np.cumsum(list_of_positions)

    return list_with_chromosomes, label_positions, list_of_lengths


def b_allele_plotting(df, list_with_chromosomes, label_positions, sample_name, list_of_lengths=None):
    """ b-allele plotting
    :param df: df created from create_df_for_plotting
    :param list_with_chromosomes: list of chromosomes
    :param label_positions: labeling positions
    :param list_of_lengths: list of lengths
    :param sample_name: sample name
    :return: SVG figure plot saved to directory containing bAllelel.tsv files"""

    fig, axs = plt.subplots(figsize=(15, 8), nrows=1)

    # Scatter plot creation
    for chr in list_with_chromosomes:
        sns.scatterplot(data=df.loc[df["Chromosome"] == chr], x="Global_position", y="BAF")

    # Title fitting
    axs.set_title("Sample: " + df["sample"].unique()[0], fontsize=16)

    # xticks fitting
    axs.set_xlabel('Chromosome', fontsize=15)
    axs.set_xticks(label_positions)
    axs.set_xticklabels([i.title() for i in list_with_chromosomes], rotation=80)

    # yticks fitting
    axs.set_ylabel('B-Allele Frequency', fontsize=15)
    axs.set_yticks(np.arange(0, 1.1, 0.1))

    # overall fitting
    axs.tick_params(axis='both', which='major', labelsize=12)
    axs.set_xlim(-75, df.__len__() + 75)

    # Introduce a horizontal line for every 0.1 step on the y-axis.
    for i in np.arange(0.1,1,0.1):
        axs.axhline(y=i, linestyle="solid", color="grey", alpha=0.1)

    plt.tight_layout()

    fig.savefig(f'{sample_name}_b-allele.svg', dpi=300, bbox_inches='tight', format="svg")


## Rolling
def rolling():

    change_dir()
    file_path = create_path_files()

    # Create dataframes for the b-allele files.
    dfs_b_allele = []
    for b_allele_file, i in zip(file_path, range(len(file_path))):
       sample_name = b_allele_file.split("\\")[-1].split("_bAllele")[0]
       #print(f"Processing {sample_name}")
       sample_df = create_df_for_plotting(b_allele_file, sample_name)
       list_with_chromosomes, label_positions, list_of_lengths = labeling_details(sample_df)
       b_allele_plotting(sample_df, list_with_chromosomes=list_with_chromosomes, label_positions=label_positions, list_of_lengths=list_of_lengths, sample_name=sample_name)

def main():
    rolling()
if __name__ == "__main__":
    main()










