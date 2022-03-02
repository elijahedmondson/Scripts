#devtools::install_github("BodenmillerGroup/neighbouRhood")
#https://github.com/BodenmillerGroup/neighbouRhood/blob/master/vignettes/example_permutation_analysis.md#2-calculate-the-baseline-statistics
library(data.table)
library(dplyr)
library(magrittr)
library(dtplyr)
library(ggplot2)
library(parallel)
library(neighbouRhood)
library(gplots)
library(RColorBrewer)


# path to a (potentially modified) cellprofiller like object measurements file
fn_cells = system.file("extdata", "cell.csv", package = "neighbouRhood", mustWork = TRUE)
# path to the Object relationships
fn_relationship = system.file("extdata", "Object relationships.csv", package = "neighbouRhood", mustWork = TRUE)

n_perm = 100

# Number of cores used for multicore:
ncores=4
# 1. Load and prepare data
dat_cells = fread(fn_cells)
dat_relation = fread(fn_relationship)

# 2. Prepare Data
d = prepare_tables(dat_cells, dat_relation)


# 3. Calculate baseline stats
dat_baseline = apply_labels(d[[1]], d[[2]]) %>%
  aggregate_histo()

# 4. Calculate permutation stats
set.seed(12312)
dat_perm = rbindlist(mclapply(1:n_perm, function(x){
  dat_labels = shuffle_labels(d[[1]])
  apply_labels(dat_labels, d[[2]]) %>%
    aggregate_histo()
},mc.cores = ncores
), idcol = 'run') 

# 5. Visualize
ggplot(dat_perm %>% filter(group==1), aes(x=ct)) +
  facet_grid(FirstLabel~SecondLabel)+
  geom_histogram() +
  geom_vline(data=dat_baseline%>% filter(group==1),aes(xintercept=ct), color='red')
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

# 6. Calculate p-values
dat_p <- calc_p_vals(dat_baseline, dat_perm, n_perm = 1000, p_tresh = 0.01) 

# 7. Generate a heatmap of the number of significant interactions for the labels
pmat = dcast(dat_p, 'FirstLabel ~ SecondLabel', value.var = 'sigval', fun.aggregate = sum,
             fill=0, drop=F)

rname = pmat$FirstLabel

pmat = pmat %>%
  select(-c('FirstLabel')) %>%
  as.matrix()

row.names(pmat) <- rname

cols = rev(brewer.pal(11,'Spectral'))
cmap = colorRampPalette(cols)

# 8. Plot the heatmap

hr <- hclust(dist(pmat), method="ward.D")
heatmap.2(pmat,
          Colv = as.dendrogram(hr),
          Rowv = as.dendrogram(hr),
          trace = "none",
          col=cmap(75),
          density.info ='none'
          #comments = data.frame(names = row.names(tab_Prot))
)

