# scanonevar_2 <- function(mean_formula,
#                          var_formula,
#                          phenotype,
#                          covariate_df,
#                          mean_locus_list,
#                          var_locus_list) {
#
#   stopifnot(length(mean_locus_list) == length(var_locus_list))
#   stopifnot(all(names(mean_locus_list) %in% names(covariate_df)))
#   stopifnot(all(names(var_locus_list) %in% names(covariate_df)))
#
#   LLs <- rep(NA, length(mean_locus_list))
#   for (locux_idx in 1:length(mean_locus_list)) {
#
#     this_locus_df <- covariate_df
#     this_locus_df[,names(mean_locus_list[[locus_idx]])] <- mean_locus_list[[locus_idx]]
#     this_locus_df[,names(var_locus_list[[locus_idx]])] <- var_locus_list[[locus_idx]]
#
#     fit <- dglm::dglm(formula = mean_formula,
#                       dformula = var_formula,
#                       data = this_locus_df,
#                       method = 'ml')
#
#   }
# }