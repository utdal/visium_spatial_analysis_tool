# import libraries
import typer
import constants as constants, click
import os, random, string, pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns


@click.command()
@click.option("--file_path", prompt="Enter file-path where the Processed_files directory is present", type=str,
              help="file-path where the Processed_files directory is located")
@click.option("--upper_threshold", prompt="upper_threshold", type=int,
              help="upper_threshold is to filter the ceiling counts")
@click.option("--lower_threshold", prompt="lower_threshold", type=int, default=0,
              help="lower_threshold is to filter the floor counts")
def generate_hist_plot(file_path, upper_threshold, lower_threshold=0):
    """
    This method generates all the histogram plots for the "Single", "Multiple", "Surrounding" and "Other" dataframes
    in a histogram plot.

    :param file_path: file-path where the Processed_files directory is located, this is the very folder where the
    directories "Single" and "Multiple" are present from where this function reads the respective dataframes to
    generate the histogram plots.

    :param upper_threshold: upper_threshold is to filter the counts data towards the ceiling

    :param lower_threshold: lower_threshold is to filter the counts data towards the floor

    :return: This method returns a dict{} - response_template with STATUS: SUCCESS/FAILED and RESPONSE: File-path
    where the histogram plots are stored for the data["Multiple", "Single", "Surrounding", "Other"] for each dataset.
    """
    try:
        plots_list = list()  # empty list to return all the plot names
        counter = constants.COUNTER
        response_template = constants.response_template.copy()  # copying response template
        # Directory variables
        processed_file_dir = file_path + os.sep + constants.PROCESSED_FILES
        single_neur_fp = processed_file_dir + os.sep + constants.SINGLE_NEUR
        all_neur_fp = processed_file_dir + os.sep + constants.MULTIPLE_NEUR
        other_neur_fp = processed_file_dir + os.sep + constants.OTHER_NEUR
        surr_neur_fp = processed_file_dir + os.sep + constants.SURROUNDING_NEUR
        plots_dir = processed_file_dir + os.sep + constants.PLOTS
        # Exception Handling
        if lower_threshold >= upper_threshold:
            raise Exception("Lower threshold cannot be greater that the upper threshold, please correct it.")
        # Checking and creating the directories(Plots) if that do not exist
        if not os.path.exists(processed_file_dir):
            os.mkdir(processed_file_dir)
        if not os.path.exists(plots_dir):
            os.mkdir(plots_dir)

        # Exception Handling
        if not os.path.exists(file_path + os.sep + constants.PROCESSED_FILES):
            raise Exception("The file path: {} provided does not have Processed_files directory.".format(file_path))

        # Checking if the no. of files in Final_matrix and Neuronal_barcodes folders are same
        if not len(os.listdir(single_neur_fp)) == len(os.listdir(all_neur_fp)):
            raise Exception(constants.FILES_NOT_MATCHING_3)

        # zip, os.listdir & sorted functions help us get the respective Final matrix associated with Neural barcodes
        for single_neuron, multiple_neuron, surround_neuron, other_neuron  in zip(sorted(os.listdir(single_neur_fp)),
                                                                                  sorted(os.listdir(all_neur_fp)),
                                                                                  sorted(os.listdir(surr_neur_fp)),
                                                                                  sorted(os.listdir(other_neur_fp))):
            print("Running for {}, {}, {} and, {}".format(single_neuron, multiple_neuron, surround_neuron, other_neuron))
            # Checking and generating random variable of length 10
            check_list = list()
            sample_ext = str(single_neuron.split('.')[0].split('_')[-1])
            for file in os.listdir(plots_dir):
                file_name, extension = os.path.splitext(file)
                check_list.append(file_name.split('_')[-1])
            # Random char-num generator of length 10
            while True:
                random_var = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
                if random_var not in check_list:
                    break
            # Loading the data into sample_df
            sample_df1 = pd.read_csv(single_neur_fp + os.sep + single_neuron)
            # Sum of reads per barcode(i.e. in this file it is  a column)
            barcode_count_per_seq1 = sample_df1.sum(axis=0)
            sample_df2 = pd.read_csv(all_neur_fp + os.sep + multiple_neuron)
            barcode_count_per_seq2 = sample_df2.sum(axis=0)
            sample_df3 = pd.read_csv(surr_neur_fp + os.sep + surround_neuron)
            barcode_count_per_seq3 = sample_df3.sum(axis=0)
            sample_df4 = pd.read_csv(other_neur_fp + os.sep + other_neuron)
            barcode_count_per_seq4 = sample_df4.sum(axis=0)
            barcode_count_per_seq1, barcode_count_per_seq2 = barcode_count_per_seq1[2:], barcode_count_per_seq2[2:]
            barcode_count_per_seq3, barcode_count_per_seq4 = barcode_count_per_seq3[2:], barcode_count_per_seq4[2:]
            fig, axes = plt.subplots(2, 4, figsize=(18, 10))  # defining the size of the axes
            # Plot 1 - single neurons
            sns.histplot(x=barcode_count_per_seq1, ax=axes[0, 1], kde=True, bins=30)
            axes[0, 1].set_xlabel("Single neurons")
            axes[0, 1].set_ylabel(constants.COUNTS)
            axes[0, 1].set_title("Distribution of Single Neurons - {}".format(sample_ext))
            # Plot 2 - multiple neurons
            sns.histplot(x=barcode_count_per_seq2, ax=axes[0, 0], kde=True, bins=30)
            axes[0, 0].set_xlabel("All neurons")
            axes[0, 0].set_ylabel(constants.COUNTS)
            axes[0, 0].set_title("Distribution of all Neurons - {}".format(sample_ext))
            # Plot 3 - single neurons - upon applying threshold

            sns.histplot(x=barcode_count_per_seq1[
                (barcode_count_per_seq1 >= lower_threshold) & (upper_threshold >= barcode_count_per_seq1)],
                         ax=axes[1, 1], kde=True, bins=30)
            axes[1, 1].set_xlabel("Single neurons")
            axes[1, 1].set_ylabel(constants.COUNTS)
            axes[1, 1].set_title("Distribution of Single Neurons - {}".format(sample_ext))
            # Plot 4 - multiple neurons - upon applying threshold
            sns.histplot(x=barcode_count_per_seq2[
                (barcode_count_per_seq2 >= lower_threshold) & (upper_threshold >= barcode_count_per_seq2)],
                         ax=axes[1, 0], kde=True, bins=30)
            axes[1, 0].set_xlabel("All neurons")
            axes[1, 0].set_ylabel(constants.COUNTS)
            axes[1, 0].set_title("Distribution of all Neurons - {}".format(sample_ext))
            # Plot 5 - surrounding neurons
            sns.histplot(x=barcode_count_per_seq3, ax=axes[0, 2], kde=True, bins=30)
            axes[0, 2].set_xlabel("Surrounding")
            axes[0, 2].set_ylabel(constants.COUNTS)
            axes[0, 2].set_title("Distribution of Surrounding - {}".format(sample_ext))
            # Plot 6 - surrounding neurons - upon applying threshold
            sns.histplot(x=barcode_count_per_seq3[
                (barcode_count_per_seq3 >= lower_threshold) & (upper_threshold >= barcode_count_per_seq3)],
                         ax=axes[1, 2], kde=True, bins=30)
            axes[1, 2].set_xlabel("Surrounding")
            axes[1, 2].set_ylabel(constants.COUNTS)
            axes[1, 2].set_title("Distribution of Surrounding - {}".format(sample_ext))
            # Plot 7 - other neurons
            sns.histplot(x=barcode_count_per_seq4, ax=axes[0, 3], kde=True, bins=30)
            axes[0, 3].set_xlabel("Other")
            axes[0, 3].set_ylabel(constants.COUNTS)
            axes[0, 3].set_title("Distribution of Other - {}".format(sample_ext))
            # Plot 8 - other neurons - upon applying threshold
            sns.histplot(x=barcode_count_per_seq4[
                (barcode_count_per_seq4 >= lower_threshold) & (upper_threshold >= barcode_count_per_seq4)],
                         ax=axes[1, 3], kde=True, bins=30)
            axes[1, 3].set_xlabel("Other")
            axes[1, 3].set_ylabel(constants.COUNTS)
            axes[1, 3].set_title("Distribution of Other - {}".format(sample_ext))
            plt.tight_layout()
            plt.savefig(plots_dir + os.sep + 'hist_plot_' + sample_ext + '_' + str(random_var) + constants.PNG_FILE)
            plots_list.append(plots_dir + os.sep + 'hist_plot_' + sample_ext + '_' + str(random_var) + constants.PNG_FILE)
            counter = counter + 1
            response_template[constants.STATUS] = constants.SUCCESS
            response_template[constants.RESPONSE] = "The plot is saved here: {}".format(plots_list)

    except Exception as error_msg:
        response_template[constants.STATUS] = constants.FAILED
        response_template[constants.RESPONSE] = error_msg

    finally:
        print(response_template)
        return response_template


if __name__ == '__main__':
    generate_hist_plot()

# Testing
# ----------------------------------------------------------
# file_path = r"C:\Users\NXI220005\Desktop\Visium_near_single_neuron_workflow\Visium_near_single_neuron_workflow"
# file_path = "/Users/nikhil/Downloads/Visium_near_single_neuron_workflow/"  # replace this path
# make sure that the folders "Final_matrix" and "Neuronal_barcodes" are in the file_path above
# generate_hist_plot(file_path, 5000, 2000)  # check the upper and lower threshold(s) if need be
# ----------------------------------------------------------

