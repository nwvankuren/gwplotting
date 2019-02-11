#' gwplotting: A package to simplify plotting genome-wide statistics
#'
#' This package contains a set of small functions that perform three main tasks:
#' -Loading files in common data formats
#' -Reordering constituent scaffolds according to chromosome assignments or
#' lengths
#' -Plotting common statistics
#'
#' @section Required Packages:
#'
#' This package utilizes a variety of functions from Hadley Wickham's tidyverse
#' package. Update your R and RStudio to the latest versions, then install the
#' latest version of tidyverse.
#'
#' @section Common Format:
#'
#' I have tried to make each loading function return a standard format tibble,
#' containing four columns with the following headers:
#'   scaf - scaffold name
#'   ps   - the position on the scaffold (this may be a single site, or the start
#'          or midpoint of a window)
#'   stat - the statistic value at that position.
#'   chr  - the chromosome (this is initialized with a unique number for each
#'          scaffold for most of the load functions, but is replaced when doing
#'          the actual reordering)
#'
#' @section Loading Functions:
#'
#' As of Feb 2019, I have written in the following loading functions. Except
#' where noted, they output the Common Format described above. See their help
#' messages for more detail. The function names are fairly explanatory.
#'
#' load_gemma_gwas
#' load_vcftools_stats
#' load_plink_gwas
#' load_abbababa
#'
#' @section Reordering Functions:
#'
#' There are two ways I usually reorder scaffolds: by length or by assignments
#' to another species' chromosomes. See these function help messages for more
#' details.
#'
#' reorder_scaffolds
#' reorder_by_scaf_len
#' get_cumulative_positions
#'
#' @section Plotting Functions:
#'
#' There is one main function, but I will add as I go. This function wraps in
#' get_cumulative_positions before plotting.
#'
#' plot_genomewide_data
#' plot_region_data
#'
#' @section Linkage Disequilibrium Functions:
#'
#' There are various things that you can do with LD measurements.
#'
#' load_plink_ld
#' calculate_ld_decay
#' calculate_windowed_ld
#'
#' @section Example:
#'
#' x <- load_gemma_gwas( 'file.assoc.txt.gz', pval = 'p_wald' )
#' y <- reorder_by_scaf_len( x, 'scaffolds.chromSizes' )
#' z <- plot_genomewide_data( y, type = 'gwas', 'scaffolds.chromSizes',
#'      plotting_column = 'stat' )
#' tiff( 'myplot.tiff', width = 4, height = 2, units = 'in', res = 600 )
#' z
#' dev.off()
#'
#'
#' @docType package
#' @name gwplotting
NULL
