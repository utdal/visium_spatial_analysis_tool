# constants file
COLNAMES_VISIUM_SINGLE_NEURONS = ['Barcode', 'X', 'Y',
                                  'Gene ID', 'Gene Name',
                                  'Text', 'Gene Expression']
SINGLE_LIST = ['single', 'Single', 'SINGLE', 'neuron']
MULTIPLE_LIST = ['multiple', 'Multiple', 'MULTIPLE', 'multi']

# constants
GENEID = 'Gene ID'
GENENAME = 'Gene Name'
GENEEXPR='Gene Expression'
IDENTIFIER = 'Identifier'
COUNTER = 1
NEURONS = 'neurons'
BARCODE = 'Barcode'
ID_PREFIX = 'id1-'
CSV_FILE = ".csv"
PNG_FILE = ".png"
PDF = "PDF"
PNG = "PNG"
PDF_FILE = ".pdf"
TSNE1 = "tSNE1"
TSNE2 = "tSNE2"

SINGLE_NEURON_ID = "single_neurons_sample_"
MULTIPLE_NEURON_ID = "multiple_neurons_sample_"
ALL_NEURON_ID = "single_multiple_neurons_sample_"
OTHER_NEURON_ID = "other_sample_"
SURROUNDING_NEURON_ID = 'surrounding_sample_'

# response_template
response_template = {'Status': None,
                     'Response': None}
SUCCESS = 'Success'
FAILED = 'Failed'
STATUS = 'Status'
RESPONSE = 'Response'

SINGLE_NEUR = 'Single'
CLUSTER = 'Seurat_input_to_cluster'
MULTIPLE_NEUR = 'Multiple'
ALL_NEUR = 'Single_and_Multiple'
SURROUNDING_NEUR = "Surrounding"
OTHER_NEUR = "Other"
SEURAT_INPUT = 'Seurat_input_to_violin_plots'
MULTIPLE_NEUR_QC = "Neurons"
SINGLE_NEUR_CLEANUP_FILTER = "Single_neur_cleanup_filter_"
SINGLE_NEUR_CLEANUP = "Single_neur_cleanup_"
SINGLE_NEUR_CLEANUP_THRESHOLD = "Single_neur_cleanup_"
SINGLE_NEUR_CLEANUP_THRESHOLD = "Single_neur_cleanup_threshold_"
PROCESSED_FILES = 'Processed_files'
PLOTS = 'Plots'
COUNTS = "Counts"
TYPER_MAIN_CALL_ERROR = "Error occurred while running the function call"
FILES_NOT_MATCHING_1 = "The files in Final_matrix & Neuronal_barcodes are not a match, please check and try again."
FILES_NOT_MATCHING_2 = "The files in Multiple, Surrounding and Other are not a match, please check and try again."
FILES_NOT_MATCHING_3 = "The files in Single and Multiple are not a match, please check and try again."
BARCODES_RESPONSE = "Barcodes with low expression are removed successfully."

GENE_MAP_DICT = {
    'MARCH1': 'MARCHF1', 'MARC1': 'MTARC1', 'MARCH2': 'MARCHF2', 'MARC2': 'MTARC2', 'MARCH3': 'MARCHF3',
    'MARCH4': 'MARCHF4', 'MARCH5': 'MARCHF5', 'MARCH6': 'MARCHF6', 'MARCH7': 'MARCHF7', 'MARCH8': 'MARCHF8',
    'MARCH9': 'MARCHF9', 'MARCH10': 'MARCHF10', 'MARCH11': 'MARCHF11', 'SEPT1': 'SEPTIN1', 'SEPT2': 'SEPTIN2',
    'SEPT3': 'SEPTIN3', 'SEPT4': 'SEPTIN4', 'SEPT5': 'SEPTIN5', 'SEPT6': 'SEPTIN6', 'SEPT7': 'SEPTIN7',
    'SEPT8': 'SEPTIN8', 'SEPT9': 'SEPTIN9', 'SEPT10': 'SEPTIN10', 'SEPT11': 'SEPTIN11', 'SEPT12': 'SEPTIN12',
    'SEPT14': 'SEPTIN14', 'DEC1': 'DELEC1'
}
