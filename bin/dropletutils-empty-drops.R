#!/usr/bin/env Rscript 

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options

option_list = list(
  make_option(
    c("-i", "--input-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which a serialized R SingleCellExperiment object can be found."
  ),
  make_option(
    c("-l", "--lower"),
    action = "store",
    default = 100,
    type = 'numeric',
    help = "A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets."
  ),
  make_option(
    c("-n", "--niters"),
    action = "store",
    default = 10000,
    type = 'numeric',
    help = "An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations."
  ),
  make_option(
    c("-m", "--test-ambient"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "A logical scalar indicating whether results should be returned for barcodes with totals less than or equal to lower."
  ),
  make_option(
    c("-g", "--ignore"),
    action = "store",
    default = NULL,
    type = 'numeric',
    help = "A numeric scalar specifying the lower bound on the total UMI count, at or below which barcodes will be ignored."
  ),
  make_option(
    c("-r", "--retain"),
    action = "store",
    default = NULL,
    type = 'numeric',
    help = "A numeric scalar specifying the threshold for the total UMI count above which all barcodes are assumed to contain cells."
  ),
  make_option(
    c("-f", "--filter-empty"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Should barcodes estimated to have no cells be removed from the output object?"
  ),
  make_option(
    c("-d", "--filter-fdr"),
    action = "store",
    default = 0.01,
    type = 'numeric',
    help = "FDR filter for removal of barcodes with no cells"
  ),
  make_option(
    c("-t", "--output-text-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name of text file in which to store output data frame."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized SingleCellExperiment object."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_text_file', 'output_object_file'))

# Now we're hapy with the arguments, load DropletUtils and do the work

suppressPackageStartupMessages(require(DropletUtils))

# Input from serialized R object

single_cell_experiment <- readRDS(opt$input_object_file)

# Run the function. Functions used internally by emptyDrops (see e.g.
# edgeR::goodTuringProportions() expect integers, and emptyDrops will no work
# correctly with non-integers from e.g. Alevin et al, so we test with integers
# only. 

testmat <- counts(single_cell_experiment)
testmat@x <- round(testmat@x)

empty <- emptyDrops(testmat, lower=opt$lower, niters=opt$niters, test.ambient = opt$test_ambient, ignore = opt$ignore, retain = opt$retain)

# Report to STDOUT on likely empty cells

# Some parameters aren't interesting for reporting purposes (e.g. file
# locations), so hide from the summary

nonreport_params <- c('input_object_file', 'output_object_file', 'help', 'output_text_file')
opt_table <- data.frame(value=unlist(opt), stringsAsFactors = FALSE)
opt_table <- opt_table[! rownames(opt_table) %in% nonreport_params, , drop = FALSE]

is.cell <- empty$FDR <= opt$filter_fdr & ! is.na(empty$FDR)
cat(c(
  paste0(
    'At an FDR of ', opt$filter_fdr,', estimate that ',
    sum(is.cell, na.rm = TRUE),
    ' barcodes have cells.'
  ),
  ifelse( opt$filter_empty, paste("Will filter to",  sum(is.cell, na.rm = TRUE), 'barcodes.'), ''),
  '\nParameter values:',
  capture.output(print(opt_table))
), sep = '\n')     

# Output a tsv

write.table(cbind(Barcode=rownames(empty), empty), file = opt$output_text_file, sep="\t", quote = FALSE, row.names = FALSE)

# Stash the results in the object metadata

colnames(empty) <- paste0('empty', colnames(empty))
colData(single_cell_experiment)[,colnames(empty)] <- empty

# Filter empty cells if specified

if (opt$filter_empty){
    nonempty_cells <- rownames(empty)[is.cell]
    single_cell_experiment <- single_cell_experiment[, nonempty_cells]
}

# Re-save the object

saveRDS(single_cell_experiment, file = opt$output_object_file)
