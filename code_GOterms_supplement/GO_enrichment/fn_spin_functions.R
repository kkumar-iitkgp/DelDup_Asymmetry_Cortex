#-----------------------------------------------------------------------------
## Spatial permutation functions. 
#-----------------------------------------------------------------------------

## These rotate only the reference/hierarchy map when requested.

normalize_perm <- function(perm.id, nroi) {
  perm.id <- as.matrix(perm.id)
  if (min(perm.id,na.rm=TRUE)==0) perm.id <- perm.id+1L
  stopifnot(nrow(perm.id)==nroi)
  perm.id
}

spin_cor_one_map <- function(x_static, y_map, perm.id, corr.type='spearman', two_tailed=TRUE) {
  ok <- is.finite(x_static) & is.finite(y_map)
  if (!all(ok)) stop('spin_cor_one_map currently expects complete maps in the same ROI space.')
  perm.id <- normalize_perm(perm.id,length(y_map))
  obs <- cor(x_static,y_map,method=corr.type)
  null <- apply(perm.id,2,function(idx) cor(x_static,y_map[idx],method=corr.type))
  if (two_tailed) p <- (1+sum(abs(null)>=abs(obs)))/(length(null)+1) else p <- (1+sum(null>=obs))/(length(null)+1)
  list(rho=obs,p_spin=p,null=null)
}

spin_delta_r_one_map <- function(del, dup, hierarchy, perm.id, corr.type='spearman', two_tailed=TRUE) {
  # Requested contrast: Δr = cor(DUP,H) - cor(DEL,H).
  # ONLY H is spun. DEL and DUP maps, and therefore their relative ranking, remain static.
  stopifnot(length(del)==length(dup),length(del)==length(hierarchy))
  perm.id <- normalize_perm(perm.id,length(hierarchy))
  obs <- cor(dup,hierarchy,method=corr.type) - cor(del,hierarchy,method=corr.type)
  null <- apply(perm.id,2,function(idx) {
    hp <- hierarchy[idx]
    cor(dup,hp,method=corr.type) - cor(del,hp,method=corr.type)
  })
  if (two_tailed) p <- (1+sum(abs(null)>=abs(obs)))/(length(null)+1) else p <- (1+sum(null>=obs))/(length(null)+1)
  list(delta_r=obs,p_spin=p,null=null)
}

spin_topography_bundle <- function(del, dup, hierarchy, perm.id, corr.type='spearman') {
  delta_z <- as.numeric(scale(dup)) - as.numeric(scale(del))
  a <- spin_cor_one_map(del,hierarchy,perm.id,corr.type)
  b <- spin_cor_one_map(dup,hierarchy,perm.id,corr.type)
  c <- spin_cor_one_map(delta_z,hierarchy,perm.id,corr.type)
  d <- spin_delta_r_one_map(del,dup,hierarchy,perm.id,corr.type)
  tibble::tibble(p_spin_Del_H=a$p_spin,p_spin_Dup_H=b$p_spin,p_spin_DeltaZ_H=c$p_spin,p_spin_Delta_r=d$p_spin)
}

# Optional wrapper for the user's existing both-profile function for full Del-Dup spatial correlations.
spin_deldup_with_existing_function <- function(del,dup,perm.id,corr.type='spearman') {
  if (!exists('perm.sphere.p_two_tail_both_avg',mode='function')) return(NA_real_)
  out <- perm.sphere.p_two_tail_both_avg(del,dup,perm.id,corr.type=corr.type)
  # Adapt this line if your local function returns a named list rather than a scalar p-value.
  if (length(out)==1 && is.numeric(out)) return(as.numeric(out))
  if (is.list(out) && !is.null(out$p)) return(as.numeric(out$p))
  if (is.list(out) && !is.null(out$pval)) return(as.numeric(out$pval))
  warning('Could not automatically extract p-value from perm.sphere.p_two_tail_both_avg output; returning NA.')
  NA_real_
}
