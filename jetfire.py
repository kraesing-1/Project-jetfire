### Project Jetfire ###
# Creator - L. Kraesing

# Import needed libraries
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import numpy as np
import itertools

# Variables for plot naming
B_ALLELE = "_bAllele.svg"
LOGR = "_logr.svg"


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

class ListOfChromosomes:

    def __init__(self, df):
        self.df = df

    def column(self, column):
        self.column = column

        list_of_chromosomes = list(self.df[column].unique())

        return list_of_chromosomes


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


def b_allele_plotting(df, sample_name, label_positions=None, list_of_lengths=None):
    """ b-allele plotting
    :param df: df created from create_df_for_plotting
    :param list_with_chromosomes: list of chromosomes
    :param label_positions: labeling positions
    :param list_of_lengths: list of lengths
    :param sample_name: sample name
    :return: SVG figure plot saved to directory containing bAllelel.tsv files"""

    list_with_chromosomes = ListOfChromosomes(df).column("Chromosome")

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
        interquartile_range = quartile_3 - quartile_1

        lower_threshold = quartile_1 - (1.5 * interquartile_range)
        upper_threshold = quartile_3 + (1.5 * interquartile_range)

        return chr_df.loc[(chr_df["LogR"] >= lower_threshold) & (chr_df["LogR"] <= upper_threshold)]


class Visualize:

    def __init__(self, LogRatio_df, bAllele_df, coverage_df):
        self.LogRatio_df = LogRatio_df
        self.bAllele_df = bAllele_df
        self.coverage_df = coverage_df
        self.chromosomes = (ListOfChromosomes(self.LogRatio_df).column("Chromosome"))

        self.fig, self.axs = plt.subplots(nrows=3, ncols=len(self.chromosomes), figsize=(20, 7))
        self.fig.subplots_adjust(wspace=0)

    def plot_coverage(self):

        self.palette = itertools.cycle(sns.color_palette())

        for idx, chr in enumerate(self.chromosomes):

            coverage_df = self.coverage_df.loc[self.coverage_df["Chromosome"] == chr]

            ax_0 = self.axs[0, idx]
            sns.scatterplot(data=coverage_df,
                            x="Position",
                            y="Coverage",
                            ax=ax_0,
                            color=self.palette.__next__(),
                            s=10)

            ax_0.spines["right"].set_visible(False)
            ax_0.spines["left"].set_visible(False)
            ax_0.set_ylim(-100, 1500)
            ax_0.set_yticks([])
            ax_0.set_xticks([])
            ax_0.set(ylabel=None)
            ax_0.set(xlabel=None)
            ax_0.axhline(y=np.mean(coverage_df["Coverage"]), linestyle="solid", color="black", alpha=0.5)
            #ax_0.axhline(y=np.mean(coverage_df["Coverage"].groupby(np.arange(len(coverage_df)) // 20)), linestyle="solid", color="black", alpha=0.5)

            for i in np.arange(0, 1500, 250):
                ax_0.axhline(y=i, linestyle="solid", color="grey", alpha=0.1)

        return self

    def plot_logr(self):

        self.palette = itertools.cycle(sns.color_palette())

        for idx, chr in enumerate(self.chromosomes):

            df = OutlierFilter(self.LogRatio_df).filter(chr)
            ax_1 = self.axs[1,idx]
            sns.scatterplot(data=df,
                            x="Start",
                            y="LogR",
                            ax=ax_1,
                            color=self.palette.__next__(),
                            s=10)

            ax_1.spines["right"].set_visible(False)
            ax_1.spines["left"].set_visible(False)
            ax_1.set_ylim(-2, 2)
            ax_1.set_yticks([])
            ax_1.set(ylabel=None)
            ax_1.set(xlabel=None)
            ax_1.set_xticks([])

            for i in np.arange(-2, 3, 1):
                ax_1.axhline(y=i, linestyle="solid", color="grey", alpha=0.1)

        return self

    def plot_bAllele(self):

        self.palette = itertools.cycle(sns.color_palette())

        for idx, chr in enumerate(self.chromosomes):

            bAllele_df = self.bAllele_df.loc[self.bAllele_df["Chromosome"] == chr]

            ax_2 = self.axs[2,idx]
            sns.scatterplot(data=bAllele_df,
                            x="Position",
                            y="BAF",
                            ax=ax_2,
                            color=self.palette.__next__(),
                            s=10)

            ax_2.spines["right"].set_visible(False)
            ax_2.spines["left"].set_visible(False)
            ax_2.set_ylim(-0.02, 1.02)
            ax_2.set_yticks([])
            ax_2.set(ylabel=None)
            ax_2.set(xlabel=None)
            ax_2.set_xticks([np.median(bAllele_df["Position"])])
            ax_2.set_xticklabels([chr.title()], rotation=80, fontsize=10)

            for i in np.arange(0, 1.1, 0.1):
                ax_2.axhline(y=i, linestyle="solid", color="grey", alpha=0.1)

        return self


    def format_figure(self):
        # y-axis coverage
        self.axs[0, 0].set_ylabel("Coverage")
        self.axs[0, 0].set_yticks(np.arange(0, 1500, 250))
        self.axs[0, 0].spines['left'].set_visible(True)
        self.axs[0, -1].spines['right'].set_visible(True)

        # y-axis LogR
        self.axs[1,0].set_ylabel("LogR")
        self.axs[1,0].set_yticks(np.arange(-2, 3, 1))
        self.axs[1,0].spines['left'].set_visible(True)
        self.axs[1,-1].spines['right'].set_visible(True)

        # y-axis B-Allele Frequency
        self.axs[2,0].set_ylabel("B-Allele Frequency")
        self.axs[2,0].set_yticks(np.arange(0, 1.1, 0.1))
        self.axs[2,0].spines['left'].set_visible(True)
        self.axs[2,-1].spines['right'].set_visible(True)
        #sample = self.df["sample"].unique()[0]

        # x-axis for chromosome
        self.fig.supxlabel('Chromosome', y=-0.1)

        #self.fig.suptitle(f"Sample: {sample}", fontsize=10, fontweight="bold", y=1)
        return self

    def plot(self):
        self.plot_coverage()
        self.plot_logr()
        self.plot_bAllele()

        self.format_figure()

        return self.fig




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





os.chdir(r"Y:\BSSH_project\DRAGEN TruSight Oncology 500 v2_6_2_4_ds.ee24f1e982424c32861d4815a86b17ec\Logs_Intermediates\Gis\FFPE_colon_tissue1_rep2_dna")
os.getcwd()
logR = pd.read_csv(r"FFPE_colon_tissue1_rep2_dna_logRatio.tsv", sep="\t", names=["Chromosome", "Start", "Stop", "RS_id", "LogR"])
b_allele = pd.read_csv(r"FFPE_colon_tissue1_rep2_dna_bAllele.tsv", sep="\t", names=["Chromosome", "Position", "RS_id", "BAF"])
coverage = pd.read_csv(r"FFPE_colon_tissue1_rep2_dna_raw.tsv", sep="\t", usecols=[0,1,3], names=["Chromosome", "Position", "Coverage"])

coverage = coverage.loc[coverage["Chromosome"] != "chrY"]
coverage["norm_cov"] = coverage["Coverage"] / np.median(coverage["Coverage"])
coverage["log2_ratio"] = coverage["norm_cov"].apply(lambda x: np.log2(x))

b_allele2 = b_allele.loc[b_allele["Chromosome"] == "chr2"]
b_allele_above = b_allele2.loc[(b_allele["BAF"] >= 0.5) & (b_allele["BAF"] <= 0.95)]