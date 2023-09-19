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
    help = "File name in which a serialized R SingleCellExperiment object can be found"
  ),
  make_option(
    c("-p", "--prop"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A numeric scalar or, if bycol=TRUE, a file with ncol(x) lines specifying a value for each cell. All values should lie in [0, 1] specifying the downsampling proportion for the matrix or for each cell."
  ),
  make_option(
    c("-c", "--bycol"),
    action = "store",
    default = TRUE,
    type = 'logical',
    help = "A logical scalar indicating whether downsampling should be performed on a column-by-column basis."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized SingleCellExperiment object."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'prop', 'output_object_file'))

# Now we're hapy with the arguments, load DropletUtils and do the work

suppressPackageStartupMessages(require(DropletUtils))

# Input from serialized R object

single_cell_experiment <- readRDS(opt$input_object_file)

# Test if we have a single or multiple lines for 'prop', and derive the value

if ( !is.na(as.numeric(opt$prop)) ){
  prop <- as.numeric(opt$prop)
} else if ( file.exists(opt$prop) ){
  prop <- readLines(opt$prop)
  if ( all(!is.na(as.numeric(opt$prop))) ){
    if ( length(prop) == length( single_cell_experiment ) ){
      prop <- as.numeric(prop)
    }else{
      stop(paste0(length(prop), ' proportions supplied, but object has ', ncol(single_cell_experiment), 'cells.'))
    }
  }else{
    stop(paste('Proportions supplied in', opt$prop, 'are not numeric.'))
  }
} else{
  stop(paste(prop, 'is neither a fixed proportion nor a valid file'))
}

# Keep in in single_cell_experiment format

counts(single_cell_experiment) <- downsampleMatrix(counts(single_cell_experiment), prop = prop, bycol = opt$bycol)

# Output to a serialized R object

saveRDS(single_cell_experiment, file = opt$output_object_file)