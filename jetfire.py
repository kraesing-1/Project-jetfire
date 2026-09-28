### Project Jetfire ###
# Creator - L. Kraesing

# Import needed libraries
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import numpy as np
import itertools
from scipy.stats import iqr

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


def b_allele_plotting(df, list_with_chromosomes, sample_name, label_positions=None, list_of_lengths=None):
    """ b-allele plotting
    :param df: df created from create_df_for_plotting
    :param list_with_chromosomes: list of chromosomes
    :param label_positions: labeling positions
    :param list_of_lengths: list of lengths
    :param sample_name: sample name
    :return: SVG figure plot saved to directory containing bAllelel.tsv files"""

    fig, axs = plt.subplots(nrows=1, ncols=len(list_with_chromosomes), figsize=(10, 2.5))
    fig.subplots_adjust(wspace=0)

    palette = itertools.cycle(sns.color_palette())

    for idx, chr in enumerate(list_with_chromosomes):

        df = df.loc[df["Chromosome"] == chr]

        sns.scatterplot(data=df, x="Position", y="BAF",
                        ax=axs[idx],
                        color=palette.__next__(),
                        s=10)

        axs[idx].set_yticks([])
        axs[idx].spines['right'].set_visible(False)
        axs[idx].spines['left'].set_visible(False)

        # Y-axis and y-label fitting
        fig.supylabel('BAF', fontsize=10, x=0.075)
        #axs[idx].get_yaxis().set_visible(False)
        axs[idx].set_ylim(-0.02, 1.02)
        axs[idx].set(ylabel=None)
        for i in np.arange(0, 1.1, 0.1):
            axs[idx].axhline(y=i, linestyle="solid", color="grey", alpha=0.1)

        # X-axis and x-label fitting
        fig.supxlabel('Chromosome', fontsize=10, y=-0.25)
        axs[idx].set(xlabel=None)
        axs[idx].set_xticks([np.median(df["Position"])])
        axs[idx].set_xticklabels([chr.title()], rotation=80, fontsize=10)

    axs[0].set_yticks(np.arange(0, 1.1, 0.1))
    axs[0].spines['left'].set_visible(True)
    axs[-1].spines['right'].set_visible(True)
    fig.suptitle("Sample: " + df["sample"].unique()[0], fontsize=10, fontweight="bold", y=1)

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


class OutlierFilter:

    def __init__(self, df):
        self.df = df

    def filter(self, chromosome):
        """ Calculate upper and lower thresholds for outliers.
            :return DataFrame with outliers. """

        self.chromosome = chromosome

        chr_df = self.df[self.df["Chromosome"] == self.chromosome]
        quartile_1 = np.quantile(chr_df["LogR"], 0.25)
        quartile_3 = np.quantile(chr_df["LogR"], 0.75)
        interquartile_range = q3 - q1

        lower_threshold = quartile_1 - (1.5 * interquartile_range)
        upper_threshold = quartile_3 + (1.5 * interquartile_range)

        return chr_df.loc[(chr_df["LogR"] >= lower_threshold) & (chr_df["LogR"] <= upper_threshold)]





os.getcwd()
tt = pd.read_csv(r"FFPE_colon_tissue1_rep2_dna_logRatio.tsv", sep="\t", names=["Chromosome", "Start", "Stop", "RS_id", "LogR"])

fig, axs = plt.subplots(nrows=1, ncols=len(list_with_chromosomes), figsize=(10, 2.5))
fig.subplots_adjust(wspace=0)

palette = itertools.cycle(sns.color_palette())

for idx, chr in enumerate(list_with_chromosomes):

    df = OutlierFilter(tt).filter(chr)


    sns.scatterplot(data=df, x="Start", y="LogR",
                    ax=axs[idx],
                    color=palette.__next__(),
                    s=10)


    axs[idx].set_yticks([])
    axs[idx].spines['right'].set_visible(False)
    axs[idx].spines['left'].set_visible(False)

    # Y-axis and y-label fitting
    fig.supylabel('LogR', fontsize=10, x=0.075)
    # axs[idx].get_yaxis().set_visible(False)
    axs[idx].set_ylim(-2, 2)
    axs[idx].set(ylabel=None)
    #for i in np.arange(0, 1.1, 0.1):
    #    axs[idx].axhline(y=i, linestyle="solid", color="grey", alpha=0.1)

    # X-axis and x-label fitting
    fig.supxlabel('Chromosome', fontsize=10, y=-0.25)
    axs[idx].set(xlabel=None)
    axs[idx].set_xticks([np.median(df["Start"])])
    axs[idx].set_xticklabels([chr.title()], rotation=80, fontsize=10)

axs[0].set_yticks(np.arange(-2, 3, 1))
axs[0].spines['left'].set_visible(True)
axs[-1].spines['right'].set_visible(True)
fig.suptitle("Sample: " + df["sample"].unique()[0], fontsize=10, fontweight="bold", y=1)