import click
from plot_histograms import generate_hist_plot
from quality_control import quality_control_for_visium_data
from remove_barcodes_low_expr import remove_barcodes_with_low_expression
from visium_single_multiple_neurons import fetch_gene_expr_selected_neuronal_barcodes
from visium_surrounding_others_neurons import fetch_gene_expr_surrounding_and_other_barcodes


@click.command()
@click.option("--function_call", prompt="Function call(select a number)", type=int)
def visium_func_call(function_call):
    try:
        if function_call == 1:
            print("\n")
            print("1. Fetching single and multiple neurons;")
            return fetch_gene_expr_selected_neuronal_barcodes()
        elif function_call == 2:
            print("\n")
            print("2. Fetching surrounding and other data;")
            return fetch_gene_expr_surrounding_and_other_barcodes()
        elif function_call == 3:
            print("\n")
            print("3. QC and plotting tSNE plots;")
            return quality_control_for_visium_data()
        elif function_call == 4:
            print("\n")
            print("4. Plotting histograms based on ceiling and floor thresholds;")
            return generate_hist_plot()
        elif function_call == 5:
            print("\n")
            print("5. Removing the low-expressed barcodes;")
            return remove_barcodes_with_low_expression()
        else:
            print("Not Implemented Error")

    except Exception as error_msg:
        return error_msg


if __name__ == '__main__':
    print(""" \n
Welcome to visium spatial transcriptomics application ...
-----------------------------------------------\n
1. Fetch single and multiple neuronal data
2. Fetch surrounding and other neuronal data
3. QC and plot tSNE plots
4. Plot histograms based on thresholds
5. Remove barcodes with low expression
-----------------------------------------------\n""")
    visium_func_call()
