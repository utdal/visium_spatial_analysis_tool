# imports
from functools import reduce
import pandas as pd, os, warnings
import constants as constants, click
# from visium_surrounding_others_neurons import df_merger
warnings.simplefilter(action='ignore', category=FutureWarning)


def cleanup_unique_and_gene_names(dataframe):
    """
    This function removes all the duplicate(s) from "Gene Name" column in the given dataframe, and replaces all the
    gene mappings from the gene map dictionary.
    :param dataframe: pandas dataframe
    :return: cleaned pandas dataframe free from duplicates and gene_map keys in the dataframe
    """
    suffixes = dataframe.groupby(constants.GENENAME).cumcount().astype(str)
    dataframe[constants.GENENAME] += '_' + suffixes
    dataframe[constants.GENENAME] = dataframe[constants.GENENAME].str.replace('_0', '')
    return dataframe.replace({constants.GENENAME: constants.GENE_MAP_DICT})


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
        return merged_df

    except Exception as error_msg:
        return error_msg


@click.command()
@click.option("--file_path", prompt="Enter file-path where Processed_files directory is present", type=str,
              help="file-path where the Processed_files directory is located")
@click.option("--threshold", prompt="threshold", type=int, default=0)
@click.option("--filters", prompt="filters(Ex. SNAP25,RLE1)",
              type=click.STRING, default='None')
@click.option("--filter_threshold", prompt="filter_threshold", type=int, default=0)
def remove_barcodes_with_low_expression(file_path, threshold=None, filters=None, filter_threshold=0):
    """
    This function removes the barcodes that are expressed low, and also apply two different filters based on
    filters and threshold

    :param file_path: file-path where the /Processed_files/Single is located

    :param threshold: by default threshold is set to None, but a threshold value will be used to filter the
    single neuron dataframe

    :param filter_threshold: by default filter_threshold is 0, but a filter_threshold is only applied when filters
    param is not None

    :param filters: a string of Genes separated by a comma ","
        [Ex. "SNAP25,SNAP52"]

    :return: The function creates a folder "Cluster" and saves all the files generated post filtration
    """
    try:
        response_template = constants.response_template.copy()   # copying response template

        # Directory variables
        single_neur_fp = file_path + os.sep + constants.PROCESSED_FILES + os.sep + constants.SINGLE_NEUR
        single_clust_fp = file_path + os.sep + constants.PROCESSED_FILES + os.sep + constants.CLUSTER

        # Exception Handling
        if not os.path.exists(file_path + os.sep + constants.PROCESSED_FILES):
            raise Exception("The file path: {} provided does not have Processed_files directory.".format(file_path))

        # zip, os.listdir & sorted functions help us get the respective Final matrix associated with Neural barcodes
        for single_neur_file in sorted(os.listdir(single_neur_fp)):
            print("Running: {}".format(single_neur_file))
            if not os.path.exists(single_clust_fp):
                os.mkdir(single_clust_fp)
            # Loading the single neurons dataframe
            single_neur_df = pd.read_csv(single_neur_fp + os.sep + single_neur_file)
            # Here is the cleanup step just before saving the df's for further analysis
            single_neur_df = cleanup_unique_and_gene_names(single_neur_df)  # cleaning the dataframe
            single_neur_df.set_index(constants.GENENAME, inplace=True)  # set gene name as index
            single_neur_df = single_neur_df.drop([constants.GENEID], axis=1)  # remove gene ID column
            single_neur_df = single_neur_df.T  # transposing the dataframe
            # Filter conditions
            if filters != 'None':
                # subset = single_neur_df[single_neur_df[filters] > filter_threshold].T # subset filtering condition
                for filter_field in filters.split(','):
                    try:
                        subset = single_neur_df.loc[single_neur_df[filter_field] > filter_threshold]
                    except Exception as error_msg:
                        raise Exception("Filters: {} does not exist in the data".format(filter_field))
                subset = subset.T
                if threshold is not None:
                    cleaned_df = subset.drop([col for col, val in subset.sum().iteritems() if val < threshold], axis=1)
                else: cleaned_df = subset
                cleaned_df.to_csv(single_clust_fp + os.sep + constants.SINGLE_NEUR_CLEANUP_THRESHOLD +
                                  str(single_neur_file.split('.')[0].split('_')[-1]) + constants.CSV_FILE)
            else:
                subset = single_neur_df.T
                if threshold is not None:
                    cleaned_df = subset.drop([col for col, val in subset.sum().iteritems() if val < threshold], axis=1)
                else: cleaned_df = subset
                cleaned_df.to_csv(single_clust_fp + os.sep + constants.SINGLE_NEUR_CLEANUP_THRESHOLD +
                                  str(single_neur_file.split('.')[0].split('_')[-1]) + constants.CSV_FILE)
        processed_files = file_path + os.sep + constants.PROCESSED_FILES
        print("Running: merge on Single neuron files, results are saved in {} directory.".format(
            processed_files + os.sep + constants.SEURAT_INPUT))
        df_merger(single_clust_fp).to_csv(
            processed_files + os.sep + constants.SEURAT_INPUT + os.sep + 'single_sample_merge.csv', index=False)
        response_template[constants.STATUS] = constants.SUCCESS
        response_template[constants.RESPONSE] = constants.BARCODES_RESPONSE

    except Exception as error_msg:
        response_template[constants.STATUS] = constants.FAILED
        response_template[constants.RESPONSE] = error_msg

    finally:
        print(response_template)
        return response_template


if __name__ == '__main__':
    remove_barcodes_with_low_expression()

# Testing
# ----------------------------------------------------------
# file_path = r"C:\Users\NXI220005\Desktop\Visium_near_single_neuron_workflow\Visium_near_single_neuron_workflow"
# file_path = "/Users/nikhil/Downloads/Visium_near_single_neuron_workflow/"  # replace this path
# make sure that the folders "Final_matrix" and "Neuronal_barcodes" are in the file_path above
# remove_barcodes_with_low_expression(file_path=file_path, threshold=1000,
#                                     filters='SNAP25,SCYL3', filter_threshold=1)
# ----------------------------------------------------------
0