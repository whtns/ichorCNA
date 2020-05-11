# file:   ichorCNA.R
# authors: Gavin Ha, Ph.D.
#          Fred Hutch
# contact: <gha@fredhutch.org>
#
#         Justin Rhoades
#          Broad Institute
# contact: <rhoades@broadinstitute.org>

# ichorCNA: https://github.com/broadinstitute/ichorCNA
# date:   July 24, 2019
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This script is the main script to run the HMM.

# library(optparse)
library(HMMcopy)
library(GenomicRanges)
library(GenomeInfoDb)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')



## load ichorCNA library or source R scripts
if (!is.null(libdir) && libdir != "None"){
	source(paste0(libdir,"/R/utils.R"))
	source(paste0(libdir,"/R/segmentation.R"))
	source(paste0(libdir,"/R/EM.R"))
	source(paste0(libdir,"/R/output.R"))
	source(paste0(libdir,"/R/plotting.R"))
} else {
    library(ichorCNA)
}

runIchorCNA(WIG = "inst/extdata/MBC_315.ctDNA.reads.wig", gcWig = "inst/extdata/gc_hg19_1000kb.wig", mapWig = "inst/extdata/map_hg19_1000kb.wig")
