# bmediatR

A Bayesian model selection approach to mediation analysis. 

## Example

The following code services as a simple vignette for using the package.

``` r
library(bmediatR)

## Simulate Collaborative Cross data
# Balanced design matrix with 20 observations of each founder allele
# thus 160 CC strains
balanced_matrix <- sim_balanced_locus(founder_allele_reps = 20)

## Simulate based on bi-allelic SNP evenly split between the founder strains
# Matrix that maps eight founder alleles to two SNP alleles
M_single <- model_matrix_from_ID("0,0,0,0,1,1,1,1")
# SNP design matrix
SNP_X <-  balanced_matrix %*% M_single

## Set seed so simulation is replicable
set.seed(10)

## Simulate a mediator for which the SNP explains 70% of its variation
simple_m <- sim_cc_single_locus(locus_matrix = balanced_matrix, 
                                num_replicates = 1, 
                                num_sim = 1,
                                M_ID = "0,0,0,0,1,1,1,1", 
                                impute = TRUE,
                                strain_effect_size = 0,
                                qtl_effect_size = 0.7)

## Simulate an target for the mediator simulates 60% of its variation
simple_y <- sim_target(simple_m$data, mediator_effect_size = 0.6)

## Simulate a null mediator (SNP has no effect)
simple_m_null <- sim_cc_single_locus(locus_matrix = balanced_matrix, 
                                     num_replicates = 1, 
                                     num_sim = 1,
                                     M_ID = "0,0,0,0,1,1,1,1", 
                                     impute = TRUE,
                                     strain_effect_size = 0,
                                     qtl_effect_size = 0,
                                     sim_label = "sim_m_null")

## Run mediation analysis
true_med <- bmediatR(y = simple_y[,1], 
                     M = cbind(simple_m$data[,1, drop = FALSE], simple_m_null$data[,1, drop = FALSE]), 
                     X = M_single_X,
                     ln_prior_c = "reactive")
null_med <- bmediatR(y = simple_y[,1], 
                     M = cbind(simple_m$data[,1, drop = FALSE], simple_m_null$data[,1, drop = FALSE]), 
                     X = M_single_X,
                     ln_prior_c = "reactive")

## Plot posterior probabilities
plot_posterior_bar(true_med, mediator_id = "sim_m_1", relabel_x = "simulated m") + ggtitle("Mediation through M")
plot_posterior_bar(null_med, mediator_id = "sim_m_null_1", relabel_x = "simulated null m") + ggtitle("Mediation through null M")
```
