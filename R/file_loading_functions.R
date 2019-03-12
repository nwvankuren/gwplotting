# Loading functions

#' Load a GEMMA association results file.
#'
#' This function will read a (possibly gzipped) GEMMA results file (.assoc.txt)
#' and return a tibble containing just the scaffold, position, and statistic.
#'
#' @param file GEMMA file (.assoc.txt). Can be gzipped.
#' @param pval Statistic column to read. This must match the column name
#'     exactly and can be one of the following: p_wald, p_lrt, or p_score.
#'
#' @return A three-column tibble containing scaffold, position, and stat
#' @export
#'
#' @examples
#' a <- system.file("extdata", "test.gemma_gwas.txt.gz",
#'                  package = "gwplotting")
#'
#' b <- load_gemma_gwas( a, pval = 'p_wald' )
#' b
#'
load_gemma_gwas <- function( file, pval = 'p_wald' ){

  pp <- c( 'p_wald', 'p_lrt', 'p_score')

  if( pval %in% pp ){

    tf <- readr::read_table2( file, col_names = T )
    tf <- dplyr::select( tf, chr, ps, pval )
    colnames( tf ) <- c( 'scaf', 'ps', 'stat' )
    tf <- dplyr::mutate( tf, scaf = stringr::str_replace( scaf, "[|]size[:digit:]*$", "") )

    # Add a false chr column with a number for each chromosome
    scaf_order <- tibble::tibble( scaf = rle( tf$scaf )$values ,
                          num = as.numeric(seq( 1:length(unique( tf$scaf )))))

    tf$chr <- scaf_order$num[ match( tf$scaf, scaf_order$scaf ) ]

    return( tf )

  } else {

    stop("pval can only take values p_wald, p_lrt, or p_score")

  }
}

#' Load a VCFtools results file.
#'
#' This function will read a (possibly gzipped) VCFtools results file
#' from windowed Fst analysis (.windowed.weir.fst), nucleotide diversity
#' analysis (.windowed.pi), or Tajima's D (.Tajima.D) and return a tibble
#' containing just the scaffold, position, and statistic.
#'
#' @param file VCFtools results file. Can be gzipped.
#' @param stat Name of the column to use as statistic. This will be mean_fst,
#'     weighted_fst, pi, tajimas_d, ...
#' @param position Coordinate to use as position: start, end, or midpoint
#' @param min_sites Minimum number of sites required to keep a window
#'
#' @return A three-column tibble containing scaffold, bin position, and
#'     statistic
#' @export
#'
#' @examples
#' a <- system.file("extdata", "test.vcftools_pi.txt.gz",
#'                  package = "gwplotting")
#' b <- load_vcftools_stats( a, stat = 'pi' )
#' b
#'
load_vcftools_stats <- function( file, stat = 'mean_fst', min_sites = 0,
                                 position = 'midpoint' ){
  # Correct column name
  if( stat == 'mean_fst' ){
    stat <- 'MEAN_FST'
  } else if( stat == 'weighted_fst' ){
    stat <- 'WEIGHTED_FST'
  } else if( stat == 'pi' ) {
    stat <- 'PI'
  } else if( stat == 'tajimas_d'){
    stat <- 'TajimaD'
  } else {
    stop("stat can currently be only: mean_fst, weighted_fst, pi, or tajimas_d")
  }

  # Read in table
  tf <- readr::read_table2( file, col_names = T )
  tf <- dplyr::mutate( tf, chr = stringr::str_replace( CHROM, "[|]size[:digit:]*$", "") )

  if( 'BIN_START' %in% colnames( tf ) ){

    if( position == "midpoint" ){
      tf <- dplyr::mutate( tf, mp = ((BIN_END + BIN_START) / 2) )
    } else if( position == "start" ){
      tf <- dplyr::mutate( tf, mp = BIN_START )
    } else {
      tf <- dplyr::mutate( tf, mp = BIN_END )
    }

  } else {

    if( position == "start" ){
      warning("No BIN_START in input file, setting position to BIN_END")
    }

    tf <- dplyr::select( tf, CHROM, BIN_END, stat )

  }

  # Fst correction
  if( grepl( "FST", stat ) ){
    tf <- mutate( tf, stat = ifelse( stat < 0, yes = 0, no = stat ))
  }

  # Min sites
  if( 'N_VARIANTS' %in% colnames( tf ) ){
    tf <- dplyr::filter( tf, N_VARIANTS >= min_sites )
  }

  tf <- dplyr::select(tf, CHROM, mp, stat )
  colnames(tf) <- c( 'scaf', 'ps', 'stat' )

  # Add a false chr column with a number for each chromosome, in case you don't
  # reorder you scaffolds in any way.
  scaf_order <- tibble::tibble( scaf = rle( tf$scaf )$values ,
                                num = as.numeric(seq( 1:length(unique( tf$scaf )))))

  tf$chr <- scaf_order$num[ match( tf$scaf, scaf_order$scaf ) ]

  return( tf )

}

#' Load a PLINK association results file.
#'
#' This function will read a (possibly gzipped) PLINK results file (.assoc.txt)
#' and return a tibble containing just the scaffold, position, and statistic.
#'
#' @param file PLINK association file
#'
#' @return A three-column tibble containing scaffold, bin midpoint, and
#'     statistic value.
#' @export
#'
#' @examples
#' a <- system.file("extdata", "test.plink_gwas.qassoc.gz",
#'                  package = "gwplotting")
#' b <- load_plink_gwas( a )
#' b
load_plink_gwas <- function( file ){

  tf <- readr::read_table2( file, col_names = T )
  tf <- dplyr::mutate( tf, CHR = stringr::str_replace( CHR, "[|]size[:digit:]*$", "") )

  if( colnames( tf )[3] == "UNADJ" ){ # It has been adjusted
    # tf <- select( tf, CHR, BP, FDR_BH )
    stop("Cannot handle adjusted p-values yet. You can paste your unadjusted
         and ajusted files together and cut the proper columns if you want -
         you will need CHR, BP, and P (FDR_BH if you want the adjusted values)")
  } else { # It hasn't been adjusted
    tf <- dplyr::select( tf, CHR, BP, P )
  }

  colnames( tf ) <- c( 'scaf', 'ps', 'stat' )

  # Add a false chr column with a number for each chromosome
  scaf_order <- tibble::tibble( scaf = rle( tf$scaf )$values ,
                        num = as.numeric(seq( 1:length(unique( tf$scaf )))))

  tf$chr <- scaf_order$num[ match( tf$scaf, scaf_order$scaf ) ]

  return( tf )

}

#' Load a results file from Simon Martin's ABBABABAwindows.py script.
#'
#' This function was developed strictly following his GitHub description of the
#' output file. The columns in this output file are: scaffold, start, end,
#' sites, sitesUsed, ABBA, BABA, D, fd, fdM. You can choose to plot any of the
#' statistics, and filter windows by the number of sites used.
#'
#' @param file Input file name
#' @param stat Statistic(s) to keep. This can be a vector of column names.
#' @param position Use the start, end, or midpoint as the position?
#' @param min_sites Minimum number of sitesUsed required for window inclusion.
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' a <- system.file("extdata", "test.abbababa.csv.gz",
#'                  package = "gwplotting")
#' b <- load_abbababa( a, stat = 'D' )
#' b
load_abbababa <- function( file, stat = 'D', position = 'midpoint',
                           min_sites = 0 ){

  tf <- readr::read_csv( file, col_names = T )
  tf <- dplyr::mutate( tf, scaffold = stringr::str_replace( scaffold, "[|]size[:digit:]*$", "") )

  # Filter out windows with small numbers of sites
  tf <- dplyr::filter( tf, sitesUsed >= min_sites )

  # Position
  if( position == 'midpoint' ){
    if( 'mid' %in% colnames(tf) ){
      tf <-dplyr::select( tf, scaffold, mid, tidyselect::one_of( stat ))
    } else{
      tf <- dplyr::mutate( tf, ps = ((end + start) / 2) )
      tf <- dplyr::select( tf, scaffold, ps, tidyselect::one_of(stat))
    }
  } else if( position == 'start' ){
    tf <- dplyr::select( tf, scaffold, start, tidyselect::one_of(stat) )
  } else if( position == 'end' | ! position %in% c('midpoint', 'start', 'end') ){
    tf <- dplyr::select( tf, scaffold, end, tidyselect::one_of(stat) )
  }

  # Standardize column names
  colnames( tf )[1:2] <- c('scaf','ps')
  if( length( stat ) == 1 ){ colnames(tf)[3] <- 'stat' }

  # Add a false chr column with a number for each chromosome
  scaf_order <- tibble::tibble( scaf = rle( tf$scaf )$values ,
                        num = as.numeric(seq( 1:length(unique( tf$scaf )))))

  tf$chr <- scaf_order$num[ match( tf$scaf, scaf_order$scaf ) ]

  return( tf )

}


#' Load a plink LD file
#'
#' This function loads a file containing pairwise LD estimates generated by
#' plink. A typical command to generate this file would be:
#' \emph{plink --bfile input --allow-extra-chr --r2 inter-chr gz \
#' --ld-window-r2 0.0 --out test.plink --chr scaffold_1000}
#' This will only keep SNPs on the same chromosome, because it assumes that you
#' do not know the distance between different scaffolds, even those located to
#' the same chromosome according to reordering information.
#'
#' @param file Name of a file containing LD values calculated using
#' @param max_dist Maximum allowed distance between two SNPs, in bp
#' @param keep_ids Keep SNP IDs in returned tibble? T / F
#' @param stat Which statistic to keep? r2, d, and/or dp (D' and D'-signed)
#'
#' @return A tibble with pairwise LD measures
#' @export
#'
#' @examples
#' a <- system.file("extdata", "test.plink.ld.gz", package = "gwplotting" )
#' b <- load_plink_ld( a )
#' b
load_plink_ld <- function( file, max_dist = 50000, keep_ids = FALSE,
                           stat = c('r2', 'dp') ){

  # cols to keep at the end
  keep <- c('scaf','ps','dist',stat)

  if( keep_ids ){ keep <- c( keep, 'snpA' ) }

  tf <- readr::read_table2( file, col_names = T ) %>% select( -X9 )

  cnames <- c('scaf','ps','snpA','scaf2','ps2','snpB','r2','dp')
  colnames(tf) <- cnames[ 1:ncol( tf ) ]

  tf <- dplyr::filter( tf, scaf == scaf2 ) %>%
    dplyr::mutate( dist = ps2 - ps ) %>%
    dplyr::select( tidyselect::one_of( keep ) )

  # Add a false chr column with a number for each chromosome
  scaf_order <- tibble::tibble( scaf = rle( tf$scaf )$values ,
                                num = as.numeric(seq( 1:length(unique( tf$scaf )))))

  tf$chr <- scaf_order$num[ match( tf$scaf, scaf_order$scaf ) ]

  return( tf )

}

#' Load a results file from Simon Martin's popgenWindows.py script.
#'
#' This function was developed strictly following his GitHub description of the
#' output file. The columns in this output file will vary depending on how many
#' populations you specify. It is up to you to determine what are the necessary
#' columns to keep. Minimally, the output from that script is: scaffold, start,
#' end, mid, sites, pi_popA, pi_popB, dxy_popA_popB, Fst_popA_popB. But you may
#' have all the pairwise comparisons for each specified populations. This
#' function will return all of the statistic columns.
#'
#' @param file Input file name
#' @param position Use the start, end, or midpoint as the position?
#' @param min_sites Minimum number of sites required for window inclusion.
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' a <- system.file("extdata", "test.popgen.csv.gz",
#'                  package = "gwplotting")
#' b <- load_pgw( a )
#' b
load_pgw <- function( file, position = 'midpoint',
                           min_sites = 0 ){

  tf <- readr::read_csv( file, col_names = T )
  tf <- dplyr::mutate( tf, scaffold = stringr::str_replace( scaffold, "[|]size[:digit:]*$", "") )

  # Filter out windows with small numbers of sites
  tf <- dplyr::filter( tf, sites >= min_sites ) %>% dplyr::select( -sites )

  # Position
  if( position == 'midpoint' ){
    tf <-dplyr::select( tf, -start, -end )
  } else if( position == 'start' ){
    tf <- dplyr::select( tf, scaffold, -end, -mid )
  } else if( position == 'end' | ! position %in% c('midpoint', 'start', 'end') ){
    tf <- dplyr::select( tf, -start, -mid )
  }

  # Standardize column names
  colnames( tf )[1:2] <- c('scaf','ps')

  # Add a false chr column with a number for each chromosome
  scaf_order <- tibble::tibble( scaf = rle( tf$scaf )$values ,
                                num = as.numeric(seq( 1:length(unique( tf$scaf )))))

  tf$chr <- scaf_order$num[ match( tf$scaf, scaf_order$scaf ) ]

  return( tf )

}
