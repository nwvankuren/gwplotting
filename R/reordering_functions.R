#' Reorder scaffolds based on assignment to chromosomes.
#'
#' This function will reorder the results of a genome-wide statistical
#' analysis based on the contents of a file that assigns each scaffold to a
#' reference genome assembly. The assignments can be generated using, e.g.
#' OrderScaffoldsByBlatingProteins.pl, but must have the following columns: \cr\cr
#' scaf, scafLen, chr, strand, median_pos \cr\cr
#' where scaf is scaffold name, scafLen is the scaffold length, chr is the
#' chromosome assignment, strand is the direction of the scaffold relative to
#' the chromosome, and median pos is the position of the scaffold along the
#' chromosome. This last bit is used to relatively order the scaffolds and must
#' be increasing along the chromosome length. \cr\cr
#' The main difficulty is handling the variety of different species' chromosome
#' level assemblies. Some, like \emph{Melitaea cinxia} and
#' \emph{Papilio xuthus} follow the convention of "chr1", etc. Others, however,
#' like the Heliconius genomes, haven't specified complete chromosomes and have
#' a unique naming system, typically e.g. Hmel201003o. These need to be handled
#' differently. I can only guarantee that reordering information from the
#' following genomes can currently be handled:\cr\cr
#' \emph{Heliconius melpomene} v2.5 (hmel)\cr
#' \emph{H. erato demophoon} v1.0 (herd)\cr
#' \emph{Bombyx mori} chromosome (bmor)\cr
#' \emph{Papilio xuthus} chromosome (pxut)\cr
#' \emph{Melitaea cinxia} (mcin)\cr\cr
#' The appropriate abbreviation should be supplied to the "species" argument
#' to the reordering function.
#'
#' @param input A tibble containing scaffold, ps,
#'     stat, and chr as the first three columns.
#' @param assignments Name of the file containing the reordering information.
#' @param species The abbreviation for the species assembly you ordered to.
#'
#' @return Returns the input tibble with chr colun detailing which
#'     reference genome chromosomes the draft genome scaffolds map to. Only
#'     input rows that reside on mapped scaffolds are included.
#' @export
#'
#' @examples
#' a1 <- system.file("extdata", "test.gemma_gwas.txt.gz",
#'                   package = "gwplotting")
#'
#' a2 <- system.file("extdata", "test.reordering.txt.gz",
#'                   package = "gwplotting")
#'
#' b <- load_gemma_gwas( a1, pval = 'p_wald' )
#' b <- reorder_scaffolds( input = b, assignments = a2, species = 'pxut' )
#' b
reorder_scaffolds <- function( input , assignments, species ){

  # Currently handles:
  # pxut, bmor = chr1
  # mcin = 1
  # hmel = Hmel201003o
  # herd = Herato2001

  chr_species <- c('pxut','bmor','mcin')
  ok_species <- c('hmel','herd',chr_species)

  if( ! species %in% ok_species ){
    cat(paste0("I can't yet handle your specified species ",species))
    stop(paste0("The only acceptable species are: ",paste(ok_species,sep=",")))
  }

  reordered <- readr::read_table2( assignments , col_names = T  )

  # OLD VERSION of OrderScaffoldsByBlatProteins.pl and CURRENT (4-4-19) version
  # of ConvertRagooToMapping.pl
  #`#scaf`     scafLen scafMin scafMax chr    chrMin  chrMax strand num_proteins median_pos
  #<chr>         <int>   <int>   <int> <chr>   <int>   <int> <chr>         <int>
  #1 scaffold49  1564542   53542  901564 chr7  6339484 6621196 +                19
  #2 scaffold460  191485   28352  184843 chr2  3270890 3397095 +                 8
  #3 scaffold580  104787   14511   68652 chr3  8062310 8083970 +                 3

  # 4-3-2019 version of OrderScaffoldsByBlatProteins.pl
  # `#scaf`         scafLen	 scafMin scafMax	scafMid	   chr	chrMin	chrMax	 chrMedian strand	num_proteins
  # NW_020662793.1	10603640 66045	 10488078	5111737.75  4	 7641530	16690293 12731507	 NA	    312
  # NW_020662794.1	 3247479 49402    3244253	1903466.75 24	 5034737   9231567  5843544	 NA	     68
  # NW_020662795.1	 1923411  1723	  1923242	 827769.5	  2 11870059	15057231	3105893	  -	    111

  colnames( reordered )[1] <- 'scaf'

  # Handle new version of output
  if( 'chrMedian' %in% colnames(reordered) ){
    colnames(reordered)[ colnames( reordered ) == "chrMedian" ] <- 'median_pos'
  }

  # Strip platanus scaffold sizes if they're there
  reordered <- dplyr::mutate( reordered,
                  scaf = stringr::str_replace( scaf, "[|]size[:digit:]*$", "") ) %>%
    dplyr::mutate( orichr = paste0(chr,":",chrMin,"-",chrMax,"(",strand,")"))

  # cat("scaffold names are like so:",as.character(reordered[1,1]),"\n")
  # Handle the different species -----------------------------------------------

  if( species %in% chr_species ){

    # Strip 'chr' from chromosome names, if they're there. This simplifies -
    # will plot in numerically increasing chromosome order.
    if( grepl( "chr", reordered[1,5] ) ){
      reordered$chr <- as.numeric( unlist( purrr::map( reordered$chr,
                                          stringr::str_replace, "chr", "")))
    }
  }

  # Order it. This will order Hmel and Herd scaffolds properly
  reordered <- reordered[ order( reordered$chr, reordered$median_pos ), ]

  # Keep only those SNPs on reordered scaffolds
  input <- dplyr::filter( input, scaf %in% reordered$scaf )

  # Some scaffolds hit the reference on the reverse strand
  for_rev <- as.character( reordered$scaf[ reordered$strand == "-" &
                                             ! is.na( reordered$strand )] )

  for( rscaf in 1:length(for_rev) ){

    # Get the length of the scaffold that needs to be reversed, get the
    # positions
    scaf_len <- reordered$scafLen[ reordered$scaf == for_rev[rscaf] ]
    pos1 <- input$ps[ input$scaf == for_rev[rscaf] ]

    # The new position is length - pos1. Check that there are sites on that
    # scaffold first.
    if( length( pos1 ) > 0 ){
      pos2 <- scaf_len - pos1
      input$ps[ input$scaf == for_rev[rscaf] ] <- pos2
    }

  }


  # Actually do the reordering - by chrom then by position. Add a new column to
  # input that tells you which reference scaffold / chromosome the scaffolds
  # match
  input$chr <- reordered$chr[ match( input$scaf, reordered$scaf ) ]
  input$mpos <- reordered$median_pos[ match( input$scaf, reordered$scaf ) ]
  input$orichr <- reordered$orichr[ match( input$scaf, reordered$scaf ) ]

  # Sanity check
  input <- filter(input, ! is.na( chr ) )

  # Order by chromosome, then relative position on the chrom, then position on
  # the scaffold.
  input <- input[ order( input$chr,
                         input$mpos,
                         input$ps ), ]

  input <- dplyr::select( input, -mpos )

  # Parse down to chromosomes for special species. For some reason str_sub won't
  # return numerics, so you need to replace first
  if( species == 'hmel' ){
    input <- input %>% filter( ! grepl('Hmel200', chr)) %>%
      dplyr::mutate( chr = stringr::str_replace(chr, "Hmel2","") ) %>%
      dplyr::mutate( chr = as.numeric( stringr::str_sub( chr, 1,2 )))

  } else if( species == 'herd' ){
    # make a temporary chr_scaf column for sorting
    input <- input %>% dplyr::mutate( chr = stringr::str_replace(chr, "Herato","") ) %>%
      dplyr::mutate( chr = as.numeric( stringr::str_sub( chr, 1,2 )))
  }

  return( input )

}

#' Reorder scaffolds based on their size.
#'
#'This function will order a genome stats tibble according to increasing
#'scaffold size.
#'
#' @param input A three-column tibble containing (at minimum) scaffold, ps, and
#'     stat as the first three columns.
#' @param scaffold_lengths Name of the two-column tab-delimited file containing
#'     scaffold name and length.
#' @param min_length Minimum length of scaffolds to keep
#'
#' @return A four-column tibble with scaffolds ordered in increasing size.
#' @export
#'
#' @examples
#' a1 <- system.file("extdata", "test.gemma_gwas.txt.gz",
#'                  package = "gwplotting")
#'
#' a2 <- system.file("extdata", "test.chromSizes.txt.gz",
#'                   package = "gwplotting" )
#'
#' b <- load_gemma_gwas( a1, pval = 'p_wald' )
#' c <- reorder_by_scaf_len( b, a2 )
#' c
reorder_by_scaf_len <- function( input, scaffold_lengths, min_length = 0 ){

  # Get the scaffolds and lengths
  lens <- readr::read_table2( scaffold_lengths, col_names = F )
  colnames( lens ) <- c( 'scaf', 'length' )
  lens <- dplyr::mutate( lens,
                  scaf = stringr::str_replace( scaf, "[|]size[:digit:]*$", ""))

  lens <- lens[ order( lens$length, decreasing = T ) ,]

  #numeric indices
  lids <- 1:nrow(lens)

  # What index is shortest desired length?
  cutoff <- min( lids[ lens$length < min_length ] )

  # Assign "chromosome" numbers according to size order
  input$chr <- lids[ match( input$scaf, lens$scaf ) ]

  # Remove too short
  if( length( lens$length[lens$length < min_length] ) > 0 ){
    cutoff <- min( lids[ lens$length < min_length ], na.rm = T )
    input <- filter(input, chr < cutoff )
  }

  # Order
  input <- input[ order( as.numeric( input$chr ),
                         as.numeric( input$ps ) ), ]

  return( input )

}

#' Assign cumulative positions to postions split by scaffolds.
#'
#'This function will add a new column to your tibble, bp_cum with the cumulative
#'sum of positions in the ps column. This makes it so that the stat in each row in the
#'tibble can be plotted sequentially.
#'
#' @param input A tibble with the columns scaf (scaffold), ps
#'     (position), and stat (statistic) at the minimum. This function operates
#'     on the ps column.
#' @param buffer Extra amount to add between consecutive scaffolds or chromosomes.
#' @param after Where to add buffer, either after scaffolds or chromosomes.
#' @param scaffold_lengths Path to the two-column tab-delimited file containing
#'     scaffold name and length.
#'
#' @return The original tibble, containing a new column, bp_cum, with cumulative positions.
#' @export
#'
#' @examples
#' a1 <- system.file("extdata", "test.gemma_gwas.txt.gz",
#'                   package = "gwplotting" )
#'
#' a2 <- system.file("extdata", "test.chromSizes.txt.gz",
#'                   package = "gwplotting" )
#'
#' b <- load_gemma_gwas( a1, pval = 'p_wald' )
#' b <- get_cumulative_positions( input = b, scaffold_lengths = a2,
#' buffer=1000, after = 'scaffolds' )
#' b
get_cumulative_positions <- function( input, scaffold_lengths, buffer = 0,
                                      after = 'scaffolds' ){

  if( after == "scaffolds" ){
    after <- 'scaf_num'
  } else if( after == "chromosomes" ){
    after <- 'chr'
  } else {
    stop("\"after\" must be set to either scaffolds or chromosomes")
  }

  # Each subsequent scaffold gets the cumulative length of previous scaffolds
  # added to each position. Have to make totLens separately because you can't
  # put an order(match()) command in a pipe.

  # Get lengths of each scaffold, then get cumulative positions
  tot_lens <- readr::read_table2( scaffold_lengths, col_names = F )
  colnames( tot_lens ) <- c( 'scaf', 'length' )
  tot_lens <- dplyr::mutate( tot_lens,
              scaf = stringr::str_replace( scaf, "[|]size[:digit:]*$", "") )

  # Get order of scaffolds in the input file
  scaf_order <- tibble::tibble( scaf = rle( input$scaf )$values ,
                        num = as.numeric( seq(1:length( unique( input$scaf )))))

  # For keeping order in input
  input$scaf_num <- scaf_order$num[ match( input$scaf, scaf_order$scaf ) ]

  tot_lens <- dplyr::filter( tot_lens, tot_lens$scaf %in% scaf_order$scaf )
  tot_lens <- tot_lens[ match( scaf_order$scaf, tot_lens$scaf ), ]

  # Actually shift the positions
  input <- tot_lens %>%

    # Get cumulative chromosome sizes
    dplyr::mutate( cum_sum = cumsum( length ) - length ) %>%

    # Add the matching running total to a new column in gwas
    dplyr::left_join( input, . , by = c("scaf" = "scaf")) %>%

    # Add a new column of ps + cumSum + a buffer
    dplyr::mutate( bp_cum = ps + cum_sum + buffer * (chr - 1 ))


  input <- dplyr::select( input, -cum_sum, -length )

  return( input )

}



#' Generate reordering information for a VCF file
#'
#' This function will produce two files that you can use in conjunction with
#' PLINK to reorder your VCF according to chromosome assignments.
#'
#' @param vcf The name of the VCF file that you want to reorder.
#' @param scaffold_lengths Draft genome scaffold lengths.
#' @param assignments Name of the file containing the reordering information.
#' @param species The abbreviation for the species assembly you ordered to.
#' @param out Prefix of the files to write to.
#'
#' @return Nothing. Produces two files: out.exclude.txt and out.reordering.txt for
#'     use with PLINK and your original VCF file.
#' @export
#'
#' @examples
generate_reordering_for_vcf <- function( vcf, scaffold_lengths, assignments,
                                         species, out ){

  # Prepare output filenames
  for_exclude <- paste0( out, '.exclude.txt' )
  for_reorder <- paste0( out, '.reordering.txt' )

  # Read in only the first three columns.
  x <- read_table2( file = vcf, comment = "#", col_names = F,
                    col_types = cols_only( X1 = 'c', X2 = 'i', X3 = 'c'))

  # Make consistent with helper functions
  colnames(x) <- c('scaf','ps','stat')
  x$chr <- 1

  # need bp_cum for each chromosome, not as a whole. Assign to chromosomes.
  y <- reorder_scaffolds( x, assignments = assignments, species = species ) %>%
    select( scaf, ps, stat, chr )

  # Get SNPs excluded after reordering
  setdiff( x$stat, y$stat ) -> excluded
  write.table( x= excluded, file = for_exclude, quote = F, col.names = F,
               row.names = F )
  rm(excluded)
  rm(x)

  # Get cumulative positions per chr

  chrs <- unique( y$chr )

  for( i in chrs ){

    z <- y %>% filter( chr == i ) %>%
      get_cumulative_positions( ., scaffold_lengths = scaffold_lengths ) %>%
      select( stat, chr, bp_cum )

    write_delim( x = z, path = for_reorder, delim = " ", append = T,
                 col_names = F )

    message("Finished ",i," \n")
  }

  message("Finished. Apply reordering information to your VCF using the",
          "following\nPLINK command:\n\nplink --vcf ", vcf, "--exclude ",
          for_exclude, " --allow-extra-chr --allow-no-sex --make-bed ",
          "--update-chr ", for_reorder, " 2 1 --update-map ", for_reorder,
          " 3 1 \n\nYou can then recode as VCF using PLINK, if needed.")
}
