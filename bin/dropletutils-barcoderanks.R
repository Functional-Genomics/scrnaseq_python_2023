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
    c("-l", "--lower"),
    action = "store",
    default = 100,
    type = 'numeric',
    help = "A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets."
  ),
  make_option(
    c("-f", "--fit-bounds"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "A string, '<lower>,<upper>', specifying the lower and upper bouunds on the total UMI count for spline fitting."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized SingleCellExperiment object."
  ),
  make_option(
    c("-p", "--output-png-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized SingleCellExperiment object."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_png_file', 'output_object_file'))

# Now we're hapy with the arguments, load DropletUtils and do the work

suppressPackageStartupMessages(require(DropletUtils))

# Input from serialized R object

single_cell_experiment <- readRDS(opt$input_object_file)

# Calculate the barcode ranks

fitbounds <- opt$fit_bounds
if ( ! is.null(fitbounds) ){
  fitbounds <- wsc_split_string(fitbounds)
}

br.out <- barcodeRanks(counts(single_cell_experiment), lower = opt$lower, fit.bounds = fitbounds)

# The output object. Note: suspect a labelling issue in the output vectors, so
# assuming they're ordered correctly: https://support.bioconductor.org/p/116153/

colData(single_cell_experiment)$barcodeRank <-  br.out$rank
colData(single_cell_experiment)$barcodeTotal <-  br.out$total
colData(single_cell_experiment)$barcodeFitted <-  br.out$fitted

# Re-save the object

saveRDS(single_cell_experiment, file = opt$output_object_file)

# Also make an image

png(filename = opt$output_png_file, width = 1000, height = 1000)

# Making a plot.
plot(
  br.out$rank,
  br.out$total,
  log = "xy",
  xlab = "Rank",
  ylab = "Total"
)
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col = "red")
abline(h = metadata(br.out)$knee, col = "dodgerblue", lty = 2)
abline(h = metadata(br.out)$inflection,
       col = "forestgreen",
       lty = 2)
legend(
  "bottomleft",
  lty = 2,
  col = c("dodgerblue", "forestgreen"),
  legend = c("knee", "inflection")
)

dev.off()
