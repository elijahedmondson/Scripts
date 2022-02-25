################################################################################
# Estimate the number of crossovers in HS mice to get an estimate of initial
# founder proportions.
# Daniel Gatti
# dan.gatti@jax.org
# 2021-08-28
################################################################################
options(stringsAsFactors = FALSE)

### LIBRARIES ###

library(tidyverse)
library(qtl2convert)
library(qtl2)

### VARIABLES ###

base_dir    = '/media/dmgatti/data1/ColoState/HS'
data_dir    = file.path(base_dir, 'data')
data_file   = file.path(data_dir, 'HZEproject_qtl2.Rdata')
results_dir = file.path(base_dir, 'results')

# NOTE: Megamuga genotypes.
muga_dir    = '/media/dmgatti/data0/MUGA'
marker_file = file.path(muga_dir, 'mm_uwisc_v1.csv')

founder_names = c('A/J', 'AKR/J', 'BALB/cJ', 'C3H/HeJ', 'C57BL/6J', 'CBA/J', 'DBA/2J', 'LP/J')

HScolors = CCcolors
HScolors[1] = '#FFC800'
names(HScolors) = founder_names

### CODE ###

# Read in the big Rdata file.
load(data_file)

# Read in markers, convert to map, and synch with probs.
markers = read.csv(marker_file)[,1:3]
markers = subset(markers, chr %in% c(1:19, 'X'))
markers$bp_mm10 = markers$bp_mm10 * 1e-6
colnames(markers)[3] = 'pos'
map = map_df_to_list(markers, pos_column = 'pos')

for(i in seq_along(probs)) {
  
  map[[i]] = map[[i]][dimnames(probs[[i]])[[3]]]
  
} # for(i)

# Remove objects that we don't need.
rm(addcovar, K)

# Count crossovers.
# We can use maxmarg because it needs the 36 state probs to do this correctly.
# We'll just use the autosomes.
gt = vector('list', length(probs) - 1)
names(gt) = names(probs)[-20]

for(chr in names(probs)) {

  print(paste('CHR', chr))

  gt[[chr]] = matrix('', nrow = nrow(probs[[chr]]), ncol = dim(probs[[chr]])[3], 
                     dimnames = list(rownames(probs[[chr]]), dimnames(probs[[chr]])[[3]]))
  pr = round(2 * probs[[chr]])
  
  for(j in 1:dim(pr)[3]) {
  
    x = apply(pr[,,j] == 1, 1, which)
    het = which(sapply(x, length) == 2)
    gt[[chr]][het, j] = apply(sapply(x[het], names), 2, paste0, collapse = '')
    x = apply(pr[,,j] == 2, 1, which)
    hom = which(sapply(x, length) == 1)
    tmp = sapply(x[hom], names)
    gt[[chr]][hom, j] = paste0(tmp, tmp)
    other = which(gt[[chr]][,j] == '')
    if(length(other) > 0) {
      tmp = apply(apply(probs[[chr]][other,,j, drop = F], 1, rank, 
                        ties.method = 'first') >= 7, 2, which)
      gt[[chr]][other, j] = paste0(colnames(probs[[chr]])[tmp[1,]], 
                                   colnames(probs[[chr]])[tmp[2,]])
    } # if(length(other) > 0)
    stopifnot(all(gt[[chr]][,j] != ''))
    
  } # for(j)

} # for(chr)

num_crossovers = vector('list', length(gt))
names(num_crossovers) = names(gt)

# Count the number of crossovers in each sample.
for(chr in names(gt)) {
  
  curr_gt = data.frame(t(gt[[chr]]))
  curr_gt = lapply(curr_gt, factor)
  curr_gt = lapply(curr_gt, as.numeric)
  cx = lapply(curr_gt, diff)
  cx = lapply(cx, '!=', 0)
  cx = lapply(cx, which)
  num_crossovers[[chr]] = sapply(cx, length)

} # for(chr)

total_cx = matrix(unlist(num_crossovers), nrow = length(num_crossovers[[1]]), 
                  ncol = length(num_crossovers), 
                  dimnames = list(names(num_crossovers[[1]]), names(num_crossovers)))
saveRDS(total_cx, file = file.path(results_dir, 'num_crossovers.rds'))
total_cx = rowSums(total_cx[,-20])

hist(total_cx)

# Get the founder proportions.
fpr = apply(probs[[1]], c(2,3), mean) %>%
        t() %>% 
        as.data.frame() 

for(chr in names(probs)[-1]) {
  print(paste('CHR', chr))
  fpr = bind_rows(fpr, apply(probs[[chr]], c(2,3), mean) %>%
                    t() %>% 
                    as.data.frame())
} # for(chr)

fpr = right_join(markers, rownames_to_column(fpr, var = 'marker')) %>% 
        filter(!is.na(chr))
colnames(fpr) = c('marker', 'chr', 'pos', founder_names)

saveRDS(fpr, file.path(results_dir, 'founder_prop_genome.rds'))

# Founder proportions along chromosomes.
fpr %>%
  mutate(chr = factor(chr, levels = c(1:19, 'X'))) %>%
  pivot_longer(cols = `A/J`:`LP/J`, names_to = 'founder', values_to = 'freq') %>%
  ggplot() +
  geom_hline(aes(yintercept = 0.125), linetype = 2, color = 'grey50') +
  geom_line(aes(pos, freq, color = founder)) +
  scale_color_manual(values = HScolors) +
  facet_wrap(~chr, ncol = 2, dir = 'v') +
  labs(title = "Founder Allele Frequency by Chromosome")


# Founder proportions by chromosome.
fpr %>%
  mutate(chr = factor(chr, levels = c(1:19, 'X'))) %>%
  pivot_longer(cols = `A/J`:`LP/J`, names_to = 'founder', values_to = 'freq') %>%
  ggplot(aes(founder, freq, fill = founder)) +
    geom_boxplot() +
    scale_fill_manual(values = HScolors) +
    facet_wrap(~chr)

# Genome founder proportions.
fpr %>%
  mutate(chr = factor(chr, levels = c(1:19, 'X'))) %>%
  pivot_longer(cols = `A/J`:`LP/J`, names_to = 'founder', values_to = 'freq') %>%
  ggplot(aes(founder, freq, fill = founder)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = 0.125), linetype = 'dashed') +
  scale_fill_manual(values = HScolors)

# Get mean and median contributions from each founder across the genome.
fpr %>%
  mutate(chr = factor(chr, levels = c(1:19, 'X'))) %>%
  pivot_longer(cols = `A/J`:`LP/J`, names_to = 'founder', values_to = 'freq') %>% 
  group_by(founder) %>% 
  summarize(mean_prop   = mean(freq, na.rm = TRUE),
            median_prop = median(freq, na.rm = TRUE)) %>% 
  write_csv(file = file.path(results_dir, 'genome_founder_prop.csv'))






