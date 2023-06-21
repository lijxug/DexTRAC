#' @import stats
NULL

#' True-binding cutoffs estimator
#' The function to estimate cutoffs for identifying true binding counts
#' from negative controls per antibody.
#'
#' @param atbd_mt Matrix of antibody binding counts from all cells
#' @param experiment_cells Vector. Which cells are in the experiment group. Cells not specified here are considered as in the ctrol group.
#' @param ctrl_p_cutoff The controled false postive rate
#'
#' @return a named list
#' @export
#'
estimate_truecounts = function(atbd_mt, experiment_cells, ctrl_p_cutoff = 0.0005){

  # Initiate ----
  # atbd_mt = gpexpr_mt
  # ctrl_vec = major_celltype_vec
  # ctrl_type = 'CD4'
  # --------------

  stopifnot(all(names(experiment_cells) %in% rownames(atbd_mt)))

  sink('extra-message-from-fitdist')
  atbd_lst = list()
  for(i_atbd in colnames(atbd_mt)){
    plt_df = data.frame(
      cellID = rownames(atbd_mt),
      x = atbd_mt[, i_atbd],
      group = ifelse(rownames(atbd_mt) %in% experiment_cells, 'target', 'ctrl')
    )

    y = plt_df[plt_df$group == 'ctrl', 'x']
    fit_res = fitdistrplus::fitdist(
      data = y,
      distr = 'nbinom',
      method = 'mme',
      start = list(size = 0.05, mu = 0.05)
    )
    # plot(fit_res)
    if (any(is.na(fit_res$estimate))) {
      t_cutoff = 0
    } else {
      t_cutoff = qnbinom(
        ctrl_p_cutoff,
        size = fit_res$estimate['size'],
        mu = fit_res$estimate['mu'],
        lower.tail = F
      )
    }
    attached_lgl = atbd_mt[, i_atbd] > t_cutoff
    atbd_lst[[i_atbd]] = list(
      plt_df = plt_df,
      fit_res = fit_res,
      cutoff = t_cutoff,
      attached_lgl = attached_lgl,
      count_vec = atbd_mt[, i_atbd]
    )
    atbd_lst[[i_atbd]]$count_trim_vec = atbd_lst[[i_atbd]]$count_vec
    atbd_lst[[i_atbd]]$count_trim_vec[!attached_lgl] = 0
  }

  names(atbd_lst) = colnames(atbd_mt)

  sink()
  return(atbd_lst)
}


#' Calculate enrichment for clonotypes and antigens
#'
#' @param atbd_lst, calculated by `estimate_truecounts`
#' @param experiment_cells Vector. Which cells are in the experiment group. Cells not specified here are considered as in the ctrol group.
#' @param cloneid_vec Named vector of the clone id for each cell
#' @param adjust.method Specify the method for `p.adjust`
#' @param n_cores By default all the cores available
#' @return A data frame
#' @export
#'
estimate_enrichment = function(atbd_lst,
                               experiment_cells,
                               cloneid_vec,
                               adjust.method = 'BH',
                               n_cores = NULL) {
  n_cores = ifelse(is.null(n_cores), parallel::detectCores(), n_cores)

  bind_mt = do.all(lapply(atbd_lst, function(x) {
    x$attached_lgl
  }), what = cbind)

  bind_mt = bind_mt[experiment_cells,]
  flt_cloneid_vec = cloneid_vec[rownames(bind_mt)]
  clone_counts = table(flt_cloneid_vec)
  clone_counts = clone_counts[clone_counts > 1]
  clonal_cloneids =  names(clone_counts)
  n_all_cells = nrow(bind_mt)
  clone_atbd_lst = pbmcapply::pbmclapply(clonal_cloneids, mc.cores = n_cores, function(t_tcr_clone) {
    clone_atbd_vec = colSums(bind_mt[names(flt_cloneid_vec)[flt_cloneid_vec == t_tcr_clone], ])
    t_clone.size = sum(flt_cloneid_vec == t_tcr_clone)
    t_res_lst = list()
    for (t_antibody in colnames(bind_mt)) {
      # t_antibody = "HBVgp2_ct_env"
      t_antibody_size = sum(bind_mt[, t_antibody])
      test_mt = matrix(
        c(
          clone_atbd_vec[t_antibody],
          t_clone.size - clone_atbd_vec[t_antibody],
          t_antibody_size - clone_atbd_vec[t_antibody],
          n_all_cells - t_clone.size - t_antibody_size + clone_atbd_vec[t_antibody]
        ),
        ncol = 2,
        dimnames = list(c('Bind', 'NotBind'), c('Clone', 'NotClone'))
      )
      fisher_res = fisher.test(test_mt, alternative = 'greater')
      p_hyper = phyper(
        clone_atbd_vec[t_antibody] - 1,
        m = t_antibody_size,
        n = n_all_cells - t_antibody_size,
        k = t_clone.size,
        lower.tail = F
      )
      t_res_lst[[t_antibody]] = list(
        test_mt = test_mt,
        fisher_res = fisher_res,
        p_hyper = p_hyper,
        data_df = data.frame(
          CloneID = t_tcr_clone,
          Antibody = t_antibody,
          Fisher.Pval = fisher_res$p.value,
          Fisher.OR = fisher_res$estimate,
          Hyper.Pval = p_hyper
        )
      )
    }
    return(t_res_lst)
  })
  names(clone_atbd_lst) = clonal_cloneids
  pval_df = bind_rows(pbmcapply::pbmclapply(clone_atbd_lst, mc.cores = 10, function(x) {
    t_tbl = bind_rows(lapply(x, function(y) {
      rownames(y$data_df) = NULL
      return(y$data_df)
    }))
    return(t_tbl)
  }))
  pval_df$Fisher.Padj = p.adjust(pval_df$Fisher.Pval, method = adjust.method)
  pval_df$Hyper.Padj = p.adjust(pval_df$Hyper.Pval, method = adjust.method)

  return(pval_df)
}



#' Distinguish high and low avidities from enriched antigens for each clone
#'
#' @param pval_df pval_df calculated by `estimate_enrichment`
#' @param target_clones User-defined clones of interest. Setting to include all the clones with at least one significantly enriched antigen is recommended.
#' @param p_norm_cutoff The controled false postive rate
#' @export
#' @return a list
#'
determine_avidities = function(pval_df, target_clones, p_norm_cutoff = 0.005) {
  avidity_df = pval_df
  avidity_df$Fisher.OR[avidity_df$Fisher.OR == Inf] = max(avidity_df$Fisher.OR[avidity_df$Fisher.OR != Inf])
  x = log10(avidity_df[avidity_df$CloneID %in% target_clones & avidity_df$Fisher.OR > 1, 'Fisher.OR'])

  # Fitting
  GMMfit_obj = mixtools::normalmixEM(x = x, k = 2)
  which_low = which(GMMfit_obj$mu == min(GMMfit_obj$mu))
  logOR_cutoff = qnorm(
    p_norm_cutoff,
    mean = GMMfit_obj$mu[which_low],
    sd = GMMfit_obj$sigma[which_low],
    lower.tail = F
  )

  avidity_df$avidity = ifelse(avidity_df$Fisher.OR > 10 ^ logOR_cutoff, 'High', 'Low')
  avidity_df$avidity[avidity_df$Fisher.OR == 0] = 'None'
  avidity_df$avidity[avidity_df$Fisher.Padj >= padj_cutoff] = 'NotSignif'


  highAff_lst = lapply(target_clones, function(t_tcr_clone) {
    tmp_df = avidity_df[avidity_df$CloneID == t_tcr_clone, c('Antibody', 'avidity')]
    y = setNames(nm = tmp_df$Antibody, tmp_df$avidity)
    y[y != 'NotSignif']
  })
  names(highAff_lst) = target_clones

  return(list(
    x = x,
    avidity_df = avidity_df,
    GMMfit_obj = GMMfit_obj,
    logOR_cutoff = logOR_cutoff,
    highAff_lst = highAff_lst
  ))
}
