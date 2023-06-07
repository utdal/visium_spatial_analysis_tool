# Imports
import constants as constants
import os, pandas as pd, numpy as np, click


@click.command()
@click.option("--final_matrix_fp", prompt="Enter file-path for Final_matrix directory", type=str,
              help='Directory path to final matrix - Final_matrix')
@click.option("--neuronal_barcodes_fp", prompt="Enter file-path for Neuronal_barcode directory", type=str,
              help='Directory path to neuronal barcodes - Neuronal Barcodes')
@click.option("--processed_files", prompt="Enter file-path where Processed_files directory needs to be generated",
              type=str, help='Directory path where processed files need to be saved')
# @click.option("--run_param", required=True, prompt="run_param(Ex.:single/multiple", type=str,
#               help='run_param can be single or multiple')
def fetch_gene_expr_selected_neuronal_barcodes(final_matrix_fp, neuronal_barcodes_fp, processed_files):
    """
    This function processes two sets of files based on Single or Multiple neurons in the run_param parameter, the
    only discretion here is if the neurons in the dataset is "single" or NOT.
    After processing, the files are saved to the directories under /Processed_files; "Single" and "Multiple"

    :param final_matrix_fp: file-path to the Final_matrix datasets
     [Ex.: /Users/nikhil/Downloads/Visium_near_single_neuron_workflow/Final_matrix]

    :param neuronal_barcodes_fp: file-path to the Neuronal_barcode datasets
     [Ex.: /Users/nikhil/Downloads/Visium_near_single_neuron_workflow/Neuronal_barcodes]

    :param processed_files: file-path where the Processed_files directory needs to be saved, this is the very folder
    where the directories "Single" and "Multiple" are created accordingly.

    :param run_param: Run parameter says which run to be triggered, is it Single or Multiple. The user might not
    need to run for both Single and Multiple neurons always.

    :return: This method returns a dict{} - response_template with STATUS: SUCCESS/FAILED and RESPONSE: File-path
    where the processed files are stored for "Single" and "Multiple" neurons.
     """
    try:
        counter = constants.COUNTER
        col_names = constants.COLNAMES_VISIUM_SINGLE_NEURONS  # column headers to be used in the vsl_matrix
        response_template = constants.response_template.copy()  # copying response template
        # single_neur_counter, all_neur_counter = constants.COUNTER, constants.COUNTER
        processed_files = processed_files + os.sep + constants.PROCESSED_FILES  # Processed_files folder variable

        # Exceptions for the Final_matrix and Neuronal_barcodes directories
        if not os.path.exists(final_matrix_fp):
            raise Exception("Final_matrix path: {} does not exist".format(final_matrix_fp))
        if not os.path.exists(neuronal_barcodes_fp):
            raise Exception("Neuronal_barcodes path: {} does not exist".format(neuronal_barcodes_fp))
        # Checking if the no. of files in Final_matrix and Neuronal_barcodes folders are same
        if not len(os.listdir(final_matrix_fp)) == len(os.listdir(neuronal_barcodes_fp)):
            raise Exception(constants.FILES_NOT_MATCHING_1)

        # Checking and creating the directories(Single/Multiple) if that do not exist
        if not os.path.exists(processed_files):
            os.mkdir(processed_files)
        if not os.path.exists(processed_files+os.sep+constants.SINGLE_NEUR):
            os.mkdir(processed_files+os.sep+constants.SINGLE_NEUR)
        if not os.path.exists(processed_files+os.sep+constants.MULTIPLE_NEUR):
            os.mkdir(processed_files+os.sep+constants.MULTIPLE_NEUR)

        # zip, os.listdir & sorted functions help us get the respective Final matrix associated with Neural barcodes
        for final_matrix, neural_barcode in zip(sorted(os.listdir(final_matrix_fp)),
                                                sorted(os.listdir(neuronal_barcodes_fp))):
            print("Running: {} -- {}".format(final_matrix, neural_barcode))
            # Loading final matrix into a pandas dataframe, i.e. output of SpaceRanger
            vsl_matrix = pd.read_csv(final_matrix_fp + os.sep + final_matrix, names=col_names)
            # converting the vsl_matrix df to unstacked format, i.e. a feature_counts matrix
            all_barcodes = vsl_matrix.pivot_table(index=[constants.GENEID, constants.GENENAME],
                                                  columns=constants.BARCODE, values=constants.GENEEXPR).fillna(0)
            # Loading neuronal barcodes of interest that were manually selected in the Loupe browser
            neurons_df = pd.read_csv(neuronal_barcodes_fp + os.sep + neural_barcode)

            # Checking the run_param from the SINGLE_LIST and MULTIPLE_LIST
            # if run_param in constants.SINGLE_LIST:
            #     var, var_neur = constants.SINGLE_NEUR, constants.SINGLE_NEURON_ID
            #     # Filtering single neurons
            #     single_neurons = neurons_df[neurons_df[constants.NEURONS].isin(constants.SINGLE_LIST)]
            #     # Fetching the barcodes for the respective neurons
            #     barcodes = single_neurons[constants.BARCODE].to_list()
            #     temp_file_name = processed_files + os.sep + constants.SINGLE_NEUR + os.sep + \
            #                      constants.SINGLE_NEURON_ID + str(single_neur_counter) + constants.CSV_FILE
            #     single_neur_counter = single_neur_counter + 1
            # else:
            #     var, var_neur = constants.MULTIPLE_NEUR, constants.MULTIPLE_NEURON_ID
            #     # Fetching the barcodes for the respective neurons
            #     barcodes = neurons_df[constants.BARCODE].to_list()
            #     temp_file_name = processed_files + os.sep + constants.MULTIPLE_NEUR + os.sep + \
            #                      constants.MULTIPLE_NEURON_ID + str(all_neur_counter) + constants.CSV_FILE
            #     all_neur_counter = all_neur_counter + 1

            # Single neurons
            single_neurons = neurons_df[neurons_df[constants.NEURONS].isin(constants.SINGLE_LIST)]
            # Fetching the barcodes for the respective neurons
            barcodes = single_neurons[constants.BARCODE].to_list()
            temp_file_name = processed_files + os.sep + constants.SINGLE_NEUR + os.sep + \
                             constants.SINGLE_NEURON_ID + str(final_matrix.split(".")[0].split('_')[-1]) + constants.CSV_FILE
            return_df = all_barcodes[np.intersect1d(all_barcodes.columns, barcodes)]
            # return_df = return_df.add_prefix("id{}-".format(final_matrix.split('.')[0][-1]))
            return_df = return_df.add_prefix("id{}-".format(final_matrix.split('.')[0].split('_')[-1]))
            return_df.to_csv(temp_file_name)  # saving the dataframe to the temp_file_name path
            # Fetching the barcodes for the respective neurons - All/Multiple neurons
            barcodes = neurons_df[constants.BARCODE].to_list()
            temp_file_name = processed_files + os.sep + constants.MULTIPLE_NEUR + os.sep + \
                             constants.MULTIPLE_NEURON_ID + str(final_matrix.split(".")[0].split('_')[-1]) + constants.CSV_FILE
            return_df = all_barcodes[np.intersect1d(all_barcodes.columns, barcodes)]
            # Change ID to sampleID, so we can always associate these barcodes with the respective sample
            return_df = return_df.add_prefix("id{}-".format(final_matrix.split('.')[0].split('_')[-1]))
            return_df.to_csv(temp_file_name)  # saving the dataframe to the temp_file_name path
            counter = counter + 1

        response_template[constants.STATUS] = constants.SUCCESS
        response_template[constants.RESPONSE] = "Files generated are saved in the directories: {} & {}".format(
            processed_files + os.sep + constants.SINGLE_NEUR,
            processed_files + os.sep + constants.MULTIPLE_NEUR)

    except Exception as error_msg:
        response_template[constants.STATUS] = constants.FAILED
        response_template[constants.RESPONSE] = error_msg

    finally:
        print(response_template)  # returns a response template
        return response_template


if __name__ == '__main__':
    fetch_gene_expr_selected_neuronal_barcodes()

# Testing
# ----------------------------------------------------------
# file_path = "/Users/nikhil/Downloads/Visium_near_single_neuron_workflow/"  # replace this path
# make sure that the folders "Final_matrix" and "Neuronal_barcodes" are in the file_path above
# file_path = r"C:\Users\NXI220005\Desktop\Visium_near_single_neuron_workflow\Visium_near_single_neuron_workflow"
# Single
# response_template = fetch_gene_expr_selected_neuronal_barcodes(file_path + os.sep + "Final_matrix",
#                                                                file_path + os.sep + "Neuronal_barcodes",
#                                                                file_path, 'single')
# Multiple
# response_template = fetch_gene_expr_selected_neuronal_barcodes(file_path + os.sep + "Final_matrix",
#                                                                file_path + os.sep + "Neuronal_barcodes",
#                                                                file_path, 'multiple')
# ----------------------------------------------------------
