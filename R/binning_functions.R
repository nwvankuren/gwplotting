#' Bin LD values for decay plotting
#'
#' @param input A tibble generated by load_plink_ld
#' @param bin_size Size of bins, in bp, in which to calculate LD
#' @param stat Calculate decay for r2 or dp (dprime)?
#' @param max_dist Maximum distance between SNPs to allow.
#'
#' @return A tibble
#' @export
#'
#' @examples
#' a <- system.file("extdata", "test.plink.ld.gz", package = "gwplotting" )
#' b <- load_plink_ld( a, max_dist = 10000 )
#' c <- calculate_ld_decay( b, bin_size = 200 )
#' c
calculate_ld_decay <- function( input, bin_size = 200, max_dist = 50000,
                                stat = 'r2' ){

  # Don't need all that other stuff for genome-wide decay
  input <- dplyr::select( input, dist, stat ) %>%
    dplyr::filter( dist <= max_dist )

  # What happens when we're inconsistent?
  if( max_dist > max( input$dist ) ){
    warning(paste0("Your maximum allowed distance is set to a value greater than
            the maximum distance found in your input. Setting max_dist to
            ",max(input$dist)))
    max_dist <- max( input$dist, na.rm = T )
  }

  # Bin it
  input$bin <- ggplot2::cut_width( x = input$dist, width = bin_size,
                                   boundary = 0, closed = 'left' )

  # use standard evaluation summarize_
  todo <- paste0('mean(',stat,')')

  # Group and mean
  for_return <- input %>% dplyr::group_by( bin ) %>%
    dplyr::summarize_( mean_stat = todo )

  # Add a column for the max distance between SNPs in the bin
  for_return <- dplyr::mutate( for_return, position = bin ) %>%
    dplyr::mutate( position = stringr::str_replace( position, "\\[.*,","" )) %>%
    dplyr::mutate( position = as.numeric( stringr::str_replace( position, "\\)|\\]","" ))) %>%
    dplyr::select( position, mean_stat )

}


#' Function for applying to single-scaffold LD tibbles. Used by the
#' calculate_windowed_ld function.
#'
#' @param input A tibble
#' @param window Window size, in bp
#' @param step Step size, in bp
#'
#' @return A summary tibble
#' @export
#'
#' @examples
ld_roller <- function( input, window, step ){

  start <- min( input$ps, na.rm = T )
  end <- max( input$ps, na.rm = T )

  nwin <- round(( end - start ) / step)

  if( nwin == 0 ){
    nwin <- 1
  }

  message( paste0("Calculating means in a total of ", nwin, " windows."))

  # desired output: ps, stat, nvar, chr
  results <- matrix( nrow = nwin,
                     ncol = 4,
                     dimnames = list( 1:nwin, c('ps','stat','nvar','chr') ) )

  # Window by window
  win_num <- 1

  while( start + window <= end ){

    results[ win_num, 'ps' ] <- start + 0.5 * window

    # subset
    vars <- filter( input, ps > start & ps <= start + window )

    results[ win_num, 'nvar' ] <- length( unique( vars$ps ))
    results[ win_num, 'stat' ] <- mean( vars$stat, na.rm = T )

    win_num <- win_num + 1
    start <- start + step


    if( ( win_num %% 100 ) == 0 ){
      message(sprintf( "Finished %.2f%% of windows", (win_num/nwin*100) ))
    }
  }

  results[,'chr'] <- input$chr[1]

  results <- tibble::as_tibble( results )

  return( results )

}

#' Calculate average LD values in sliding windows
#'
#' This function will take as input a tibble containing a PLINK LD file
#' loaded with load_plink_ld(), calculate mean LD values (based on the stat
#' column values), and return a new tibble.
#'
#' @param input A tibble constructed with load_plink_ld
#' @param window Size of window, in bp, to calculate mean LD in
#' @param step Step size between windows
#'
#' @return A four-column tibble containing summary statistics
#' @export
#'
#' @examples
calculate_windowed_ld <- function( input, window, step ){

  # input:
  # A tibble: 39,636,531 x 5
  #   scaf           ps  dist    dp    chr
  #   <chr>       <dbl> <dbl>  <dbl> <dbl>
  # 1 scaffold160    82   497 0.771      1
  # 2 scaffold160    82   521 0.609      1
  # 3 scaffold160    82   565 0.190      1

  # Desired output:
  # scaf ps stat nsites chr

  # Scaffold by scaffold

  output <- input %>% group_by( scaf ) %>% nest() %>%
    mutate( means = map( data, ld_roller, window = window, step = step ) ) %>%
    select( scaf, means ) %>%
    unnest() %>%
    ungroup()

  return( output )
}
