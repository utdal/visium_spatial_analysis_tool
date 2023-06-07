# imports
from functools import reduce
import constants as constants
import os, pandas as pd, numpy as np, click, dask.dataframe as dd


def df_merger(df_file_loc):
    """
    This method merges the dataframes to a single dataframe
    :param df_file_loc: path where the .csv files exist
    :return: merged dataframe
    """
    try:
        df_list = list()
        for file in os.listdir(df_file_loc):
            temp_df = pd.read_csv(df_file_loc+os.sep+file)
            try:
                temp_df.drop(constants.GENEID, axis=1, inplace=True)
                df_list.append(temp_df)
            except:
                df_list.append(temp_df)
                continue
        merged_df = reduce(lambda left, right: pd.merge(left, right, on=constants.GENENAME, how='outer'), df_list)
        merged_df.fillna(value=0, inplace=True)
        suffixes = merged_df.groupby(constants.GENENAME).cumcount().astype(str)
        merged_df[constants.GENENAME] += '_' + suffixes
        merged_df[constants.GENENAME] = merged_df[constants.GENENAME].str.replace('_0', '')
        return merged_df

    except Exception as error_msg:
        return error_msg


def df_merger_dask(df_file_loc):
    """
    This method merges the dataframes to a single dataframe
    :param df_file_loc: path where the .csv files exist
    :return: merged dataframe
    """
    try:
        df_list = list()
        for file in os.listdir(df_file_loc):
            temp_df = pd.read_csv(df_file_loc+os.sep+file)
            try:
                temp_df.drop(constants.GENEID, axis=1, inplace=True)
                df_list.append(temp_df)
            except:
                df_list.append(temp_df)
                continue
        merged_df = reduce(lambda left, right: dd.merge(left, right, on=constants.GENENAME, how='outer'), df_list)
        merged_df = merged_df.fillna(value=0)
        suffixes = merged_df.groupby(constants.GENENAME).cumcount().astype(str)
        merged_df[constants.GENENAME] += '_' + suffixes
        merged_df[constants.GENENAME] = merged_df[constants.GENENAME].str.replace('_0', '')
        return merged_df

    except Exception as error_msg:
        return error_msg


def tally_cols_function(act_file_path, merge_file_path, comb_length):
    const = 0
    for act_file in os.listdir(act_file_path):
        const = const + pd.read_csv(act_file_path + os.sep + act_file).shape[1]
    if (const - 2 * comb_length + 1) == (pd.read_csv(merge_file_path).shape[1]):
        return True
    else:
        return False


@click.command()
@click.option("--final_matrix_fp", prompt="Enter file-path for Final_matrix directory", type=str,
              help='Directory path to final matrix - Final_matrix')
@click.option("--neuronal_barcodes_fp", prompt="Enter file-path for Neuronal_barcode directory", type=str,
              help='Directory path to neuronal barcodes - Neuronal Barcodes')
@click.option("--processed_files", prompt="Enter file-path where Processed_files directory needs to be generated",
              type=str, help='Directory path where processed files need to be saved')
def fetch_gene_expr_surrounding_and_other_barcodes(final_matrix_fp, neuronal_barcodes_fp, processed_files):
    """
    This function processes two sets of files based on Surrounding and Other neurons, here the surrounding neurons
    are saved to /Processed_files/Surrounding directory and other neurons are saved to /Processed_files/Others
    directory.

    :param final_matrix_fp: file-path to the Final_matrix datasets
     [Ex.: /Users/nikhil/Downloads/Visium_near_single_neuron_workflow/Final_matrix]

    :param neuronal_barcodes_fp: file-path to the Neuronal_barcode datasets
     [Ex.: /Users/nikhil/Downloads/Visium_near_single_neuron_workflow/Neuronal_barcodes]

    :param processed_files: file-path where the Processed_files directory needs to be saved, this is the very folder
    where the directories "Single" and "Multiple" are created accordingly.

    :return: This method returns a dict{} -  response_template with STATUS: SUCCESS/FAILED and RESPONSE: File-path
    where the processed files are stored for "Surrounding" and "Other" neurons.
    """
    try:
        counter = constants.COUNTER  # defining counter
        col_names = constants.COLNAMES_VISIUM_SINGLE_NEURONS  # column-name headers
        response_template = constants.response_template.copy()
        processed_files = processed_files + os.sep + constants.PROCESSED_FILES  # Processed_files folder variable

        # Exceptions for the Final_matrix and Neuronal_barcodes directories
        if not os.path.exists(final_matrix_fp):
            raise Exception("Final_matrix path: {} does not exist".format(final_matrix_fp))
        if not os.path.exists(neuronal_barcodes_fp):
            raise Exception("Neuronal_barcodes path: {} does not exist".format(neuronal_barcodes_fp))
        # Checking if the no. of files in Final_matrix and Neuronal_barcodes folders are same
        if not len(os.listdir(final_matrix_fp)) == len(os.listdir(neuronal_barcodes_fp)):
            raise Exception(constants.FILES_NOT_MATCHING_1)

        # Checking and creating the directories(Surrounding & Other) if that do not exist
        if not os.path.exists(processed_files):
            os.mkdir(processed_files)
        if not os.path.exists(processed_files + os.sep + constants.SURROUNDING_NEUR):
            os.mkdir(processed_files + os.sep + constants.SURROUNDING_NEUR)
        if not os.path.exists(processed_files + os.sep + constants.OTHER_NEUR):
            os.mkdir(processed_files + os.sep + constants.OTHER_NEUR)

        # zip, os.listdir & sorted functions help us get the respective Final matrix associated with Neural barcodes
        for final_matrix, neural_barcode in zip(sorted(os.listdir(final_matrix_fp)),
                                                sorted(os.listdir(neuronal_barcodes_fp))):
            print("Running: {} -- {}".format(final_matrix, neural_barcode))
            # Loading final matrix into a pandas dataframe, i.e. output of SpaceRanger
            vsl_matrix = pd.read_csv(final_matrix_fp + os.sep + final_matrix,
                                     names=col_names)
            # Converting the vsl_matrix df to unstacked format, i.e. a feature_counts matrix
            all_barcodes = vsl_matrix.pivot_table(index=[constants.GENEID, constants.GENENAME],
                                                  columns=constants.BARCODE, values=constants.GENEEXPR)
            all_barcodes = all_barcodes.fillna(0)  # or inplace=True can also be used
            neurons_df = pd.read_csv(neuronal_barcodes_fp + os.sep + neural_barcode)  # neuronal barcodes dataframe
            barcodes_list = neurons_df[constants.BARCODE].to_list()  # barcodes
            # Fetching surrounding barcodes
            single_nuer = vsl_matrix[vsl_matrix[constants.BARCODE].isin(barcodes_list)]
            neur_cord = single_nuer.loc[:, [constants.BARCODE, 'X', 'Y']].drop_duplicates()
            # Fetching 6 coordinates around each neuronal barcode
            surround_df = pd.concat([pd.DataFrame((neur_cord['X'] + 1, neur_cord['Y'] + 1)).T,
                                     pd.DataFrame((neur_cord['X'] + 0, neur_cord['Y'] + 2)).T,
                                     pd.DataFrame((neur_cord['X'] - 1, neur_cord['Y'] + 1)).T,
                                     pd.DataFrame((neur_cord['X'] - 1, neur_cord['Y'] - 1)).T,
                                     pd.DataFrame((neur_cord['X'] + 0, neur_cord['Y'] - 2)).T,
                                     pd.DataFrame((neur_cord['X'] + 1, neur_cord['Y'] - 1)).T])
            surround_df.columns = ['X', 'Y']
            barcodes_ids = vsl_matrix.loc[:, [constants.BARCODE, 'X', 'Y']]
            # Fetching barcodes based on the coordinates
            surround_coord = pd.merge(surround_df, barcodes_ids, on=['X', 'Y'], how='inner')
            # Remove duplicates from neuronal barcodes #if a surrounding barcode overlaps a neuron, exclude it
            temp_surround_df = surround_coord[
                surround_coord[constants.BARCODE].isin(barcodes_list) == False].drop_duplicates()
            barcodes_interest2 = list(temp_surround_df[constants.BARCODE])
            # Change id to corresponding id
            surround = all_barcodes[np.intersect1d(all_barcodes.columns, barcodes_interest2)]
            surround = surround.add_prefix("id{}-".format(final_matrix.split('.')[0].split('_')[-1]))
            temp_file_name = processed_files + os.sep + constants.SURROUNDING_NEUR + os.sep \
                             + constants.SURROUNDING_NEURON_ID + str(final_matrix.split(".")[0].split('_')[-1]) + constants.CSV_FILE
            surround.to_csv(temp_file_name)
            # Get other barcodes
            other_barcodes = list(set(all_barcodes.columns) - set(barcodes_list + barcodes_interest2))
            other_barcodes_df = all_barcodes[np.intersect1d(all_barcodes.columns, other_barcodes)]
            # Change id to respective sample
            other_barcodes_df = other_barcodes_df.add_prefix("id{}-".format(final_matrix.split('.')[0].split('_')[-1]))
            # Save file, change id accordingly
            temp_file_name = processed_files + os.sep + constants.OTHER_NEUR + os.sep \
                             + constants.OTHER_NEURON_ID + str(final_matrix.split(".")[0].split('_')[-1]) + constants.CSV_FILE
            other_barcodes_df.to_csv(temp_file_name)
            counter = counter + 1

        all_neur_fp = processed_files + os.sep + constants.MULTIPLE_NEUR
        other_neur_fp = processed_files + os.sep + constants.OTHER_NEUR
        surr_neur_fp = processed_files + os.sep + constants.SURROUNDING_NEUR
        if not os.path.exists(processed_files + os.sep + constants.SEURAT_INPUT):
            os.mkdir(processed_files + os.sep + constants.SEURAT_INPUT)
        print("Running: merge on Multiple neuron files, results are saved in {} directory.".format(
            processed_files + os.sep + constants.SEURAT_INPUT))
        multiple_sample_merge = processed_files + os.sep + constants.SEURAT_INPUT + os.sep + 'multiple_sample_merge.csv'
        df_merger_dask(all_neur_fp).to_csv(multiple_sample_merge, index=False, chunksize=10000)
        print("Running: merge on Other files, results are saved in {} directory.".format(
            processed_files + os.sep + constants.SEURAT_INPUT))
        other_sample_merge = processed_files + os.sep + constants.SEURAT_INPUT + os.sep + 'other_sample_merge.csv'
        df_merger_dask(other_neur_fp).to_csv(other_sample_merge, index=False, chunksize=10000)
        print("Running: merge on Surrounding files, results are saved in {} directory.".format(
            processed_files + os.sep + constants.SEURAT_INPUT))
        surrounding_sample_merge = processed_files + os.sep + constants.SEURAT_INPUT + os.sep + 'surrounding_sample_merge.csv'
        df_merger_dask(surr_neur_fp).to_csv(surrounding_sample_merge, index=False, chunksize=10000)

        # # Additional QC Checks
        # multi_constant, surr_constant, other_constant = 0, 0, 0
        # multiple_dir = processed_files + os.sep + constants.MULTIPLE_NEUR
        surrounding_dir = processed_files + os.sep + constants.SURROUNDING_NEUR
        other_neur_dir = processed_files + os.sep + constants.OTHER_NEUR
        # combined_length = len(os.listdir(multiple_dir) + os.listdir(surrounding_dir) + os.listdir(other_neur_dir))
        #
        # for direc, mer_file in zip([multiple_dir, surrounding_dir, other_constant],
        #                            [multiple_sample_merge, surrounding_sample_merge, other_sample_merge]):
        #     print(direc, mer_file, combined_length)
        #     if tally_cols_function(direc, mer_file, combined_length) is True: continue
        #     else: raise Exception("There is a merge conflict in the files, please check with the administrator.")

        response_template[constants.STATUS] = constants.SUCCESS
        response_template[constants.RESPONSE] = "Files created in {} & {}".format(other_neur_dir, surrounding_dir)

    except Exception as error_msg:
        response_template[constants.STATUS] = constants.FAILED
        response_template[constants.RESPONSE] = error_msg

    finally:
        print(response_template)
        return response_template


if __name__ == '__main__':
    fetch_gene_expr_surrounding_and_other_barcodes()

# Testing
# ----------------------------------------------------------
# file_path = r"C:\Users\NXI220005\Desktop\Visium_near_single_neuron_workflow\Visium_near_single_neuron_workflow"
# file_path = "/Users/nikhil/Downloads/Visium_near_single_neuron_workflow/"  # replace this path
# make sure that the folders "Final_matrix" and "Neuronal_barcodes" are in the file_path above
# fetch_gene_expr_surrounding_and_other_barcodes(file_path + os.sep + "Final_matrix",
#                                                file_path + os.sep + "Neuronal_barcodes",
#                                                file_path)
# ----------------------------------------------------------
