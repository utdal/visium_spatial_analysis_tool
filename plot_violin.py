# Imports
import constants as constants, warnings, random, string, click
import matplotlib.pyplot as plt, pandas as pd, scanpy as sc, matplotlib as mpl, os
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)


@click.command()
@click.option("--file_path", prompt="Enter file-path where the Processed_files directory is present", type=str,
              help="file-path where the Processed_files directory is located")
def plot_violin(file_path: str):
    """
    This method generates all the violin plots for the "Single", "Multiple", "Surrounding" and "Other" dataframes
    in a violin plot.

    :param file_path: file-path where the Processed_files directory is located, this is the very folder where the
    directories "Single" and "Multiple" are present from where this function reads the respective dataframes to
    generate the histogram plots.
    :return: This method returns a dict{} - response_template with STATUS: SUCCESS/FAILED and RESPONSE: File-path
    where the violin plots are stored for the data["Multiple", "Single", "Surrounding", "Other"] for each dataset.
    """
    # sc.settings.verbosity = 3
    # sc.logging.print_header()
    # sc.settings.set_figure_params(dpi=80, facecolor="white")
    violin_plot_data = file_path+os.sep+constants.PROCESSED_FILES+os.sep+constants.SEURAT_INPUT
    plots_dir = file_path+os.sep+constants.PROCESSED_FILES+os.sep+constants.PLOTS
    for file in sorted(os.listdir(violin_plot_data)):
        print("Running for {}".format(file))
        temp_data_frame = pd.read_csv(violin_plot_data+os.sep+file)
        temp_data_frame = temp_data_frame.set_index([constants.GENENAME])
        adata = sc.AnnData(temp_data_frame)
        # Checking and generating random variable of length 10
        check_list = list()
        for plot in os.listdir(plots_dir):
            file_name, extension = os.path.splitext(plot)
            check_list.append(file_name.split('_')[-1])
        # Random char-num generator of length 10
        while True:
            random_var = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
            if random_var not in check_list:
                break
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
        sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts'], jitter=0.4, rotation=90, show=False)
        temp_plot_path = plots_dir + os.sep + 'violin_plot_' + str(file.split('_')[0]) + constants.PNG_FILE
        plt.savefig(temp_plot_path)


if __name__ == '__main__':
    plot_violin()














