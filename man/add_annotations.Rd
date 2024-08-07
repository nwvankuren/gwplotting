% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{add_annotations}
\alias{add_annotations}
\title{Create a tibble containing information for plotting gene models using ggplot
and geom_rect().}
\usage{
add_annotations(
  input,
  ymin,
  ymax,
  scaffold_lengths,
  gff = NA,
  feature = "gene",
  mapping = "chr"
)
}
\arguments{
\item{input}{Tibble containing scaf, ps, stat, chr, and bp_cum ONLY for the
region you want to plot. That is, your data has already been reordered.}

\item{ymin}{Minimum y positions for geom_rects.}

\item{ymax}{Maximum y positions for geom_rects.}

\item{scaffold_lengths}{Path to a two-column tab-delimited file containing
scaffold names and their lengths.}

\item{gff}{Path to a GFF or BED format file, if plotting gene models. If left
NA, then elements in the output correspond to constituent scaffolds in the
plotted region.}

\item{feature}{GFF feature type to plot. Usually genes. Not extensively
tested.}

\item{mapping}{Path to AGP file, or "chr" if the assembly is already at
chromosome level and annotated.}
}
\value{
A four-column tibble containing an ID, xmins, xmaxs, ymins, and ymaxs
    for use in a geom_rect() call.
}
\description{
This function will calculate the coordinates of scaffold spans or gene models
for addition to a ggplot. It will use scaffold spans encoded in the input
tibble. GFF features of the specified type will be plotted if a GFF file is
specified.
}
\details{
The input must only contain information from the region that you want to plot.
The reason for this limitation is that you may want to plot gene models from
multiple scaffolds and these will need to be reordered and reoriented based
on the scaffold-to-chromosome mapping. You can simply filter in the function
call instead of creating another object, though.
}
\examples{
a1 <- system.file("extdata", "test.gemma_gwas.txt.gz",
                 package = "gwplotting")
a2 <- system.file("extdata", "test.chromSizes.txt.gz",
                  package="gwplotting")

# Load statistics, add cumulative positions (bp_cum)
b <- load_gemma_gwas( a1, pval = 'p_wald' )
c <- get_cumulative_positions( b, scaffold_lengths = a2 )

# Get scaffold rects
scafs <- add_annotations( c, ymin = 0, ymax = 0.5 )

scafs

# Add rects to a plot of the region data
d <- plot_region_data( c, chromosome = 1, type = 'gwas',
                       plotting_column = 'stat', scaffold_lengths = a2) +
     geom_rect( inherit.aes = F,
                data = scafs,
                mapping = aes( xmin = xmins,
                               xmax = xmaxs,
                               ymin = ymins,
                               ymax = ymaxs),
                color = 'black',
                alpha = 0.5,
                size = 0.25 )

d
}
