################################################################################
# Select the 129S proportion and number of outcrossing generations for genail
# model.
# Run the haplotype reconstruction pipeline on several different 129S 
# proportions and generations and produce a table of:
# - number of crossovers
# - LOD at Chr 2 locus
# - median 129S proportion across genome
# - heritability of survival
# for 129S proportions of 10, 15, 20, 25, 30, 40 & 50
# for generations of 2, 10, 20, 30, 40 & 50.
# Daniel Gatti
# Dan.gatti@jax.org
# 2021-08-31
# 
# ASSUMPTIONS: The genotypes have already been formatted into qtl2 format,
#              all of the required files are in one directory and we only 
#              need to change teh cross_info file to update founder
#              proportions and generation.
# Files in directory:
#     atmko_hs.json
#     covar.csv
#     cross_info.csv
#     founder_geno.csv
#     gmap.csv
#     markers.csv
#     pheno.csv
#     pmap.csv
#     sample_geno.csv
################################################################################
options(stringsAsFactors = FALSE)

### LIBRARIES ###

library(qtl2convert)
library(qtl2)

### VARIABLES ###

base_dir = '/media/dmgatti/data1/ColoState/ATM'
hap_dir  = file.path(base_dir, 'haplo_reconstr')
qtl2_dir = file.path(hap_dir,  'qtl2')
fig_dir  = file.path(base_dir, 'figures')
results_dir = file.path(base_dir, 'results')
data_dir = file.path(base_dir, 'data')
pheno_file      = file.path(data_dir, 'phenotypes', 'atmko_hs_phenotypes_cleaned.csv')
covar_file      = file.path(qtl2_dir, 'covar.csv')
geno_file       = file.path(qtl2_dir, 'sample_geno.csv')
cross_info_file = file.path(qtl2_dir, 'cross_info.csv')
map_file        = file.path(qtl2_dir, 'pmap.csv')
control_file    = file.path(qtl2_dir, 'atmko_hs.json')

muga_dir     = '/media/dmgatti/data0/MUGA/'
hs_snp_file  = file.path(muga_dir, 'atmkohsnpt_variants_129S1_SvImJ.sqlite')

founder_names = c("AJ", "AKR", "BALB", "C3H", "C57", "CBA", "DBA", "LP", "129SvE")

hs_founder_prop_file = '/media/dmgatti/data1/ColoState/HS/results/genome_founder_prop.csv'

HScolors = CCcolors
HScolors = c(CCcolors, '129SvE' = '#BB5500')
HScolors[1] = '#FFC800'
names(HScolors) = founder_names

snp_func = qtl2::create_variant_query_func(hs_snp_file)

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

### FUNCTIONS ###

# Simplify the sample IDs by removing the genotyping batch and removing
# duplicates.
# pheno_file: string containing the full path to the phenotype file.
# covar_file: string containing the full path to the covariates file.
# cross_info_file: string containing the full path to the cross_info file.
# geno_file: string containing the full path to the genotype file.
# Overwrites these three files.
clean_sample_ids = function(pheno_file, covar_file, cross_info_file, geno_file) {
  
  pheno = read.csv(pheno_file)
  covar = read.csv(covar_file)
  cross_info = read.csv(cross_info_file)
  geno  = read.csv(geno_file)
  
  # Remove the batch IDs from the samples.
  ids = sapply(strsplit(pheno$id, split = '\\.'), '[', 2)
  pheno$id = ids
  covar$id = ids
  cross_info$id = ids
  colnames(geno)[-1] = ids
  
  # Remove duplicaes.
  dupl = which(duplicated(ids))
  pheno = pheno[-dupl,]
  covar = covar[-dupl,]
  cross_info = cross_info[-dupl,]
  geno  = geno[,-(dupl + 1)]
  
  # Write out files.
  write.csv(pheno, file = pheno_file, row.names = FALSE, quote = FALSE)
  write.csv(covar, file = covar_file, row.names = FALSE, quote = FALSE)
  write.csv(cross_info, file = cross_info_file, row.names = FALSE, quote = FALSE)
  write.csv(geno,  file = geno_file,  row.names = FALSE, quote = FALSE)
  
} # clean_sample_ids()



# Update the cross_info file to contain the current founder proportions and
# outcrossing generation.
# info_file: string containing full path to cross_info file.
# prop129: floating point number between 0 and 1, indicating the proportion
#          of 129S strain in the cross.
# gen: integer indicating the outcrossing generation.
# hs_prop: numeric vector containing the proportions of the other8 HS founders.
update_cross_info = function(info_file, prop129, gen, hs_prop) {
  
  # We use 10000 as the total contribution. (100%)
  # All values are scaled to add to 10000.
  total = 10000
  
  # Read cross_info.
  cross_info = read.csv(info_file)
  
  # Set generation.
  cross_info$ngen = gen
  
  # Set founder proportions.
  p129 = prop129 / 100
  hs_prop = hs_prop * (1.0 - p129) * total
  
  cross_info$A = hs_prop[1]
  cross_info$B = hs_prop[2]
  cross_info$C = hs_prop[3]
  cross_info$D = hs_prop[4]
  cross_info$E = hs_prop[5]
  cross_info$F = hs_prop[6]
  cross_info$G = hs_prop[7]
  cross_info$H = hs_prop[8]
  cross_info$I = p129 * total
  
  # Verify that all columns add to the total.
  stopifnot(all(abs(rowSums(cross_info[,LETTERS[1:9]]) - total) < 1e-6))
  
  write.csv(cross_info, file = info_file, row.names = FALSE, quote = FALSE)
  
} # update_cross_info()


# pr: list containing full, 36-state genoprobs.
count_crossovers = function(pr) {
  
  retval = setNames(rep(0, nrow(pr[[1]])), rownames(pr[[1]]))
  
  for(chr in seq_along(pr)) {
    
    maxgt = apply(pr[[chr]], c(1, 3), which.max)
    retval = retval + rowSums(maxgt[,-1] != maxgt[,-ncol(maxgt)])

  } # for(chr)
  
  return(retval)

} # count_crossovers()

# Get the het and homozygous proportions for 129S.
# pr: list containing full, 36-state genoprobs.
# map: list containing marker positions.
contrib_129 = function(pr, map) {
  
  n_samples = nrow(pr[[1]])
  n_markers = sum(sapply(map, length))
  
  # Keep only 129S genotypes ("I").
  keep = grep('I', colnames(pr[[1]]))
  pr   = lapply(pr, function(z) { z[,keep,] })
  Ibychr = sapply(lapply(pr, colSums), rowSums)
  Ibygen = rowSums(Ibychr) / n_markers / n_samples
  
  return(list(het = sum(Ibygen[1:8]),
              hom = Ibygen[9]))
  
} # contrib_129()


# Get the 129S homozygous block lengths.
# pr: list containing full, 36-state genoprobs.
# map: list containing marker positions.
block_len_129 = function(pr, map) {
  
  # Keep only 129S homozygote column ("II").
  pr = lapply(pr, function(z) { z[,'II',] })
  # probs elements are 2D matrices now.
  
  
} # block_len_129()


# pr: list containing 9-state allele probs.
count_129_prop = function(pr) {
  
  n_markers = sum(sapply(pr, dim)[3,])
  p129 = lapply(pr, function(z) { z[,'I',] })
  p129 = sapply(p129, rowSums)
  return(rowSums(p129) / n_markers)
  
} # count_129_prop()


# control_file: string containing full path to qtl2 control file.
reconstruct_haplotypes = function(control_file) {
  
  cross = qtl2::read_cross2(file = control_file, quiet = FALSE)
  
  # Run haplotype reconstruction.
  print('Reconstructing haplotypes...')
  return(qtl2:::calc_genoprob2(cross = cross, quiet = FALSE, cores = 8))

} # reconstruct_haplotypes()


# probs: list containing 9-state allele probs.
# pheno: data.frame containing phenotype.
est_heritability = function(probs, pheno) {
  
  print('Estimating Heritability...')
  allK = qtl2::calc_kinship(probs, type = 'overall', cores = 8)
  return(est_herit(pheno = pheno[,'rz_days',drop = F], kinship = allK))

} # est_heritability()



# Map the Chr 2 locus.
map_chr2 = function(probs, pheno, map, pct, gen) {
  
  print('Mapping on Chr 2...')
  K = qtl2::calc_kinship(probs, type = 'loco', cores = 8)
  
  gwas = scan1snps(genoprobs = probs, pheno = pheno[,'rz_days',drop = F], 
                   map = map, kinship = K, keep_all_snps = FALSE, 
                   query_func = snp_func, chr = 2, cores = 8)
  
  png(file.path(fig_dir, paste0('rz_survival_', pct,'pct_', gen,'_gen.png')), 
      width = 1000, height = 800, res = 128)
  plot_snpasso(gwas$lod, gwas$snpinfo, main = paste('Pct 129', pct, 'Gen', gen, 'RankZ Survival'))
  dev.off()
  
  return(max(gwas$lod))
  
} # map_chr2()



### MAIN ###

# Start by removing duplicates and the batch IDs.
# Only need to do once.
#clean_sample_ids(pheno_file = pheno_file, 
#                 covar_file = covar_file,
#                 cross_info_file = cross_info_file,
#                 geno_file  = geno_file)

# Use the HS HZE genotypes to estimate the HS founder proportions.
hs_founder_prop = read.csv(hs_founder_prop_file)

pheno = read.csv(pheno_file)
colnames(pheno)[1] = 'id'
rownames(pheno) = pheno$id
pheno$rz_days   = rankZ(pheno$days)

map = read.csv(map_file)
map = map_df_to_list(map, pos_column = 'pos')

smry_table = expand.grid(prop129 = seq(10, 50, 10),
                        gen = c(seq(10, 100, 10)))
smry_table$crossovers = 0
smry_table$obsprop129 = 0
smry_table$heritability = 0
smry_table$het129     = 0
smry_table$hom129     = 0
smry_table$chr2_lod   = 0

for(i in 1:nrow(smry_table)) {
  
  prop129 = smry_table$prop129[i]
  gen     = smry_table$gen[i]
  
  print(paste('p129', prop129, 'gen', gen))

  update_cross_info(info_file = cross_info_file, 
                    prop129   = prop129,
                    gen       = gen, 
                    hs_prop   = hs_founder_prop$mean_prop)

  probs = reconstruct_haplotypes(control_file = control_file)

  # Count crossovers.
  smry_table$crossovers[i] = median(count_crossovers(probs))
  
  # Synch up map & probs.
  for(j in seq_along(probs)) {
    map[[j]] = map[[j]][dimnames(probs[[j]])[[3]]]
  } # for(j)
  
  stopifnot(sapply(map, length) == sapply(probs, dim)[3,])
  
  # Get mean 129S het and homozygous contributions.
  het_homo_129 = contrib_129(probs, map)
  smry_table$het129[i] = het_homo_129$het
  smry_table$hom129[i] = het_homo_129$hom
  
  # Convert to allele probs.
  print('Converting alleleprobs...')
  aprobs = genoprob_to_alleleprob(probs = probs, quiet = FALSE, cores = 8)
  rm(probs)
  gc()

  smry_table$obsprop129[i] = median(count_129_prop(aprobs))
  
  # Synch samples between pheno & probs.
  samples = intersect(pheno$id, rownames(aprobs[[1]]))
  stopifnot(length(samples) == 676)
  
  pheno  = pheno[samples,]
  aprobs = lapply(aprobs, function(z) { z[samples,,] })
  for(j in seq_along(aprobs)) {
    map[[j]] = map[[j]][dimnames(aprobs[[j]])[[3]]]
  } # for(j)
  
  stopifnot(sapply(map, length) == sapply(aprobs, dim)[3,])
  
  # Set attributes on probs that didn't get set.
  attr(aprobs, 'crosstype') = 'genail9'
  attr(aprobs, 'is_x_chr') = setNames(c(rep(F, 19), T), c(1:19, 'X'))
  attr(aprobs, 'alleles') = LETTERS[1:9]
  attr(aprobs, 'alleleprobs') = TRUE
  attr(aprobs, 'class') = c('calc_genoprob', 'list')
  
  smry_table$heritability[i] = est_heritability(aprobs, pheno)
  smry_table$chr2_lod[i]     = map_chr2(aprobs, pheno, map, prop129, gen)
  
  write.csv(smry_table, file = file.path(results_dir, 'p129_gen_summary.csv'),
            row.names = FALSE)
  
  rm(aprobs)
  gc()

} # for(i)

