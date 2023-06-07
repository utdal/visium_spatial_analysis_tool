# import libraries
import random, string, click
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import os, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt, seaborn as sns,constants as constants


@click.command()
@click.option("--file_path", prompt="Enter file-path where the Processed_files directory is present", type=str,
              help="file-path where the Processed_files directory is located")
@click.option("--plot_type", prompt="Enter plot_type of the tSNE plot to be generated(pdf/png)", type=str,
              help="plot_type describes the type of plot to be saved (.png/.pdf)", default='png')
def quality_control_for_visium_data(file_path, plot_type='png'):
    """
    This method generates all the tSNE plots for the "Multiple", "Surrounding" and "Other" dataframes in a same
    plot. To achieve this, the function also processes the data through PCA() to reduce the dimensions and transposes
    the dataframe.

    :param file_path: file-path where the Processed_files directory is located, this is the very folder
    where the directories "Surrounding", "Other" and "Multiple" are present from where this function reads the
    respective dataframes to generate the tSNE-plots.

    :param plot_type: plot_type describes the type of plot to be saved (.png/.pdf) and the values plot_type takes
    are PDF/pdf (or) PNG/png

    :return: This method returns a dict{} - response_template with STATUS: SUCCESS/FAILED and RESPONSE: File-path
    where the tSNE-plots are stored for the data["Multiple", "Surrounding" and "Other"] for each dataset.
    """
    try:
        plots_list = list()  # empty list to return all the plot names
        response_template = constants.response_template.copy()  # copying response template
        # Directory variables
        processed_file_dir = file_path + os.sep + constants.PROCESSED_FILES
        multiple_dir = file_path + os.sep + constants.PROCESSED_FILES + os.sep + constants.MULTIPLE_NEUR
        surrounding_dir = file_path + os.sep + constants.PROCESSED_FILES + os.sep + constants.SURROUNDING_NEUR
        other_neur_dir = file_path + os.sep + constants.PROCESSED_FILES + os.sep + constants.OTHER_NEUR
        plots_dir = processed_file_dir + os.sep + constants.PLOTS

        # Custom parms
        label_size = 30
        mpl.rcParams['xtick.labelsize'] = label_size
        mpl.rcParams['ytick.labelsize'] = label_size
        mpl.rcParams["font.sans-serif"] = ["Arial"]
        mpl.rcParams["font.family"] = "Arial"
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42

        # Exception Handling
        if not os.path.exists(file_path+os.sep+constants.PROCESSED_FILES):
            raise Exception("The file path: {} provided does not have Processed_files directory.".format(file_path))

        # Checking and creating the directories(Plots) if that do not exist
        if not os.path.exists(plots_dir):
            os.mkdir(plots_dir)
        # Checking if the no. of files in Final_matrix and Neuronal_barcodes folders are same
        if not len(os.listdir(multiple_dir)) == len(os.listdir(surrounding_dir)) == len(os.listdir(other_neur_dir)):
            raise Exception(constants.FILES_NOT_MATCHING_2)

        # zip, os.listdir & sorted functions help us get the respective Final matrix associated with Neural barcodes
        for multiple_neur, surrounding_neur, other_neur in zip(sorted(os.listdir(multiple_dir)),
                                                               sorted(os.listdir(surrounding_dir)),
                                                               sorted(os.listdir(other_neur_dir))):
            print("Running for {}, {} and {}".format(multiple_neur, surrounding_neur, other_neur))
            # Checking and generating random variable of length 10
            check_list = list()
            for file in os.listdir(plots_dir):
                file_name, extension = os.path.splitext(file)
                check_list.append(file_name.split('_')[-1])
            # random char-num generator
            while True:
                random_var = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
                if random_var not in check_list:
                    break

            # Loading the dataframes from Multple, Surrounding and Other datasets
            df_neurons = pd.read_csv(multiple_dir + os.sep + multiple_neur)
            df_surrounding = pd.read_csv(surrounding_dir + os.sep + surrounding_neur)
            df_other = pd.read_csv(other_neur_dir + os.sep + other_neur)
            # Merging all the 3 datasets on GENE-ID and GENE-NAME
            df_neurons_surr = pd.merge(df_neurons, df_surrounding,
                                       on=[constants.GENEID, constants.GENENAME],
                                       how='outer')
            df_all = pd.merge(df_neurons_surr, df_other,
                              on=[constants.GENEID,constants.GENENAME],
                              how='outer')
            # Dropping GENE-ID and GENE-NAME columns
            df_all_values = df_all.drop([constants.GENEID, constants.GENENAME], axis=1)
            trans_df_all_values = np.transpose(df_all_values)  # creating a transpose matrix
            # Applying the PCA() - dimensionality reduction and selecting only 50 components
            pca_50_dim = PCA(n_components=50).fit_transform(trans_df_all_values)
            tsne_2d = TSNE(n_components=2).fit_transform(pca_50_dim)  # applying the TSNE function
            tsne_df = pd.DataFrame(tsne_2d, columns=[constants.TSNE1, constants.TSNE2])
            full_df = pd.concat([tsne_df, trans_df_all_values.reset_index()], axis=1)  # concatenating the above vals
            # Setting-up a column to identify if the row is from Multiple, Surrounding or Other
            full_df[constants.IDENTIFIER] = [constants.MULTIPLE_NEUR_QC] * (len(df_neurons.columns)-2) + \
                                    [constants.SURROUNDING_NEUR] * (len(df_surrounding.columns)-2) + \
                                    [constants.OTHER_NEUR] * (len(df_other.columns)-2)
            fig, ax = plt.subplots(figsize=(18, 16))  # defining the size of the axes
            sns.scatterplot(data=full_df, x=constants.TSNE1, y=constants.TSNE2,
                            hue=constants.IDENTIFIER, ax=ax, s=50)
            fig.suptitle('TSNE Plot for Single Neuron data - {}'.format(multiple_neur.split('.')[0].split('_')[-1]), fontsize=28)
            ax.legend(fontsize=22)
            plt.xlabel('t-SNE_1', fontsize=20)  # setting the x-label
            plt.ylabel('t-SNE_2', fontsize=20)  # setting the y-label
            if plot_type.upper() == "PDF":
                plot_name = plots_dir + os.sep + 'tSNE_sample_' + str(multiple_neur.split('.')[0].split('_')[-1]) + '_' + str(random_var) + constants.PDF_FILE
            else:
                plot_name = plots_dir + os.sep + 'tSNE_sample_' + str(multiple_neur.split('.')[0].split('_')[-1]) + '_' + str(random_var) + constants.PNG_FILE
            plt.savefig(plot_name)  # saving the plot as a .pdf file to plots_name
            plots_list.append(plot_name)  # appending the plot_name to a list, to output in the response_template

        response_template[constants.STATUS] = constants.SUCCESS
        response_template[constants.RESPONSE] = "The plot is saved here: {}".format(plots_list)

    except Exception as error_msg:
        response_template[constants.STATUS] = constants.FAILED
        response_template[constants.RESPONSE] = error_msg

    finally:
        print(response_template)
        return response_template


if __name__ == '__main__':
    quality_control_for_visium_data()

# Testing
# ----------------------------------------------------------
# file_path = r"C:\Users\NXI220005\Desktop\Visium_near_single_neuron_workflow\Visium_near_single_neuron_workflow"
# file_path = "/Users/nikhil/Downloads/Visium_near_single_neuron_workflow/"  # replace this path
# make sure that the folders "Final_matrix" and "Neuronal_barcodes" are in the file_path above
# quality_control_for_visium_data(file_path, 'jhf')
# quality_control_for_visium_data(file_path, 'Pdf')

# TSNE using scanpy
# import scanpy as sc
# import matplotlib.pyplot as plt
# sc.pp.filter_genes(adata, min_cells=3)
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# sc.pp.scale(adata)
# sc.tl.tsne(adata)
# ----------------------------------------------------------
