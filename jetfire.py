### Project Jetfire ###
# Creator - L. Kraesing

# Import needed libraries
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import numpy as np
import itertools

# Set seaborn configurations
sns.set_theme()
sns.set_style("white")
sns.set_palette("tab20")

# Functions for running

def change_dir():
    """ Change directory to current working directory """
    os.chdir(os.getcwd())

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

    fig, axs = plt.subplots(nrows=1, ncols=len(list_with_chromosomes), figsize=(10, 2))
    fig.subplots_adjust(wspace=0)

    palette = itertools.cycle(sns.color_palette())

    for idx, chr in enumerate(list_with_chromosomes):

        df = df_test.loc[df_test["Chromosome"] == chr]

        sns.scatterplot(data=df, x="Position", y="BAF",
                        ax=axs[idx],
                        color=palette.__next__(),
                        s=20)

        axs[idx].set_yticks([])
        axs[idx].spines['right'].set_visible(False)
        axs[idx].spines['left'].set_visible(False)
        #axs[idx].set_xlim(0, df["Position"].max()+10000)

        # Y-axis and y-label fitting
        fig.supylabel('BAF', fontsize=15)
        #axs[idx].get_yaxis().set_visible(False)
        axs[idx].set_ylim(-0.02, 1.02)
        axs[idx].set(ylabel=None)
        for i in np.arange(0, 1.1, 0.1):
            axs[idx].axhline(y=i, linestyle="solid", color="grey", alpha=0.1)

        # X-axis and x-label fitting
        fig.supxlabel('Chromosome', fontsize=15)
        axs[idx].set(xlabel=None)
        axs[idx].set_xticks([np.median(df["Position"])])
        axs[idx].set_xticklabels([chr.title()], rotation=80)

    axs[0].set_yticks(np.arange(0, 1.1, 0.1))
    axs[0].spines['left'].set_visible(True)
    axs[-1].spines['right'].set_visible(True)
    fig.suptitle("Sample: " + df["sample"].unique()[0], fontsize=16, fontweight="bold", y=0.90)

    fig.savefig(f'{sample_name}_b-allele.svg', dpi=300, bbox_inches='tight', format="svg")


## Rolling
def rolling():

    change_dir()
    file_path = create_path_files()

    # Create dataframes for the b-allele files.
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










