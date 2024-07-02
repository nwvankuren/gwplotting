library(tidyverse)

# Get some variables
gwas <- '~/Dropbox/gwas_2024/results/plow.mim.filtered_gwas.assoc.txt.gz'
scaf_lens <- '~/Dropbox/gwas_2024/info/low_mim_ragtag.genome'
gff3 <- '~/Dropbox/gwas_2024/info/low_mim_ragtag.paplowMliftoff.gff3'

# Load data
x <- load_gemma_gwas(gwas,
                     pval = 'p_lrt')

# Reorder by chromosome. This is a new function that assumes data are already
# to chromosome level.
y <- reorder_by_chromosome(x)

# Get cumulative plotting positions.
z <- get_cumulative_positions( input = y,
                               scaffold_lengths = scaf_lens,
                               buffer = 500000,
                               after = 'chromosomes' )

# Make a genome-wide plot.
plot_genomewide_data( z, type = 'gwas', scaf_lens, plotting_column = 'stat', col1 = 'black',
                      'gray')

# Chromosome must be in z. "scaffold" matches the scaf column, "chromosome"
# matches the chr column. Flexibility.
target <- "chr17"
chr <- 17
from_pos <- 5400000
to_pos <- 5800000

# Plot data within a particular region. This is a ggplot object, so you can
# always add more layers to it, e.g. theme(...) or ylim(...) to change aspects.
plot_region_data( z, chromosome = chr, from = from_pos, to = to_pos )

# Get gene model coordinates. Have to filter input because otherwise it's really
# complicated. Alternatively, make a new object and use that as input:
# zz <- z %>%
#         filter( scaf == target, ps >= from_pos,
#         ps <= to_pos)
gene_models <- add_annotations(z %>% filter( scaf == target, ps >= from_pos,
                                             ps <= to_pos), ymin = -6,
                               ymax = 0, scaffold_lengths = scaf_lens,
                               gff = gff3, feature = 'gene', mapping = "chr")

# Put it all together. Be sure to set your ylims to include your new gene
# model coordinates, otherwise you won't see them because they're outside the
# plotting area. Probably a better way to do this...
plot_region_data(z, chromosome = chr, from = from_pos, to = to_pos ) +
  ylim( -10, 100 ) +
  geom_rect( inherit.aes = F,
             data = gene_models,
             mapping = aes( xmin = xmins,
                            xmax = xmaxs,
                            ymin = ymins,
                            ymax = ymaxs),
             color = 'black',
             fill = 'blue',
             alpha = 1,
             size = 0.25 )

# Can modify just as before. Maybe want to color a particular gene model(s)
gene_colors <- rep('gray',nrow(gene_models))
gene_colors[ gene_models$ids %in% c('dsx','sir2','pros') ] <- 'red'

plot_region_data(z, chromosome = chr, from = from_pos, to = to_pos ) +
  ylim( -10, 100 ) +
  geom_rect( inherit.aes = F,
             data = gene_models,
             mapping = aes( xmin = xmins,
                            xmax = xmaxs,
                            ymin = ymins,
                            ymax = ymaxs),
             color = 'black',
             fill = gene_colors,
             alpha = 1,
             size = 0.25 )
