library(broom)
library(viridis)
library(tidyverse)
library(Rfast) ## has fast, multi-feature ttests function

############################################################
## function to simulate a simple high dimensional dataset
gen_xy <- function( frac_feat = 0.1, frac_samples = 0.5, 
                    mean_diff = 1, n_feat = 100, n_samples = 90){
  # percent of true features and number of positive samples
  n_true_feat = round(n_feat * frac_feat)
  n_positives = round(n_samples * frac_samples)
  
  # init the matrix w/ random normal numbers
  x = sapply(1:n_feat,function(x) rnorm(n_samples, mean = 0, sd = 1) %>% t())
  if(frac_feat > 0){
    # make some of samples and features different
    x[1:n_positives, 1:n_true_feat] = 
      sapply(1:n_true_feat,function(x) rnorm(n_positives, mean = mean_diff, sd = 1) %>% t())
  }

  # labels
  y = c(rep(1, n_positives), rep(0, n_samples-n_positives))
  return(list(x %>% as.matrix(), y))
}

###################################
## test with subset of permutations with same direction
emp_p_with_same_sign <- function( x, y, p, c, n_perm = 100, ncores = 1){
  nfeat = ncol(x); out = list()
  for(i in seq(n_perm)){
    #permute the y's
    y2 = sample(y) 
    # test over all the features, ttests is really fast, test here is 2-sided
    out = c(out, list(ttests(x[y2==0,], y = x[y2==1,])))
  }
  c_perm = sapply(out, "[", 1:nfeat, 1) # coefficient from permuted t-test
  p_perm = sapply(out, "[", 1:nfeat, 2) # pvalue from permuted t-test
  # count the number of times permutations more significant w/ same sign
  p2 = sapply(1:nfeat, function(i) 
    sum( p_perm[i,] <= p[i] & sign(c_perm[i,]) == sign(c[i])))
  totalPermSameSign = sapply(1:nfeat, function(i) sum(sign(c_perm[i,]) == sign(c[i])))
  return( p2 / (totalPermSameSign + 1))
}

set.seed(1234)

df = data.frame()
for(frac in c(0, .125, .25, .5)){
  print(paste('working on', frac))
  ## generate a random dataset w/ increasing % features different
  out1 = gen_xy(frac_feat = frac, frac_samples = 0.5, n_feat = 1000)
  x = out1[[1]]; y = out1[[2]]
  
  ## do the parametric t-test
  df_param = ttests(x = x[y==0,], y = x[y==1,])
  p = df_param[,'pvalue']
  c = df_param[,'stat']
  
  ## do the modified 1-sided test with only the permutations in the same direction
  p_modi = emp_p_with_same_sign( x, y, p, c, n_perm = 1000)
  
  ## append all of these p-values to the end of the data frame
  df = rbind(df, data.frame(sample = seq_along(p),
                            pi0 = frac, 
                            parametic.ttest = p, 
                            permuted.SignChecking = p_modi))
}

## make into a long data frame for plotting
df2 = df %>% pivot_longer(-c(sample, pi0))
pdf(paste0('perm_p_plot_0to0.5_20221225.pdf'), width = 6, height = 8)
(p1 = ggplot(df2, aes(x = value, fill = name)) + 
  geom_histogram(bins = 20) + theme_bw() + 
  facet_grid(pi0~name, scales = 'free_y') +
  scale_fill_viridis('', option = 'viridis', discrete=TRUE) +
  theme(legend.position = 'bottom'))
dev.off()

p1





