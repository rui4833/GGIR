applyDFA = function(ts) {
  ## --------------------------------------------------------------------------##
  ## Detrended Fluctuation Analysis (DFA)
  ## author: Ian Meneghel Danilevicz 
  ## orcid: 0000-0003-4541-0524
  ## email: ian.meneghel-danilevicz @ inserm.fr
  ## --------------------------------------------------------------------------##
  
  ## --------------------------------------------------------------------------##
  ## Functions: SSP, DFA, DFA_aux
  ## co-author: Victor Barreto Mesquita 
  ## <victormesquita40@hotmail.com>
  ## --------------------------------------------------------------------------##
  
  ## --------------------------------------------------------------------------##
  ## Functions: ABI, SSP, DFA, DFA_aux
  ## Activity balance index (ABI)
  ## Self-similarity parameter (SSP) or alpha
  ## Detrended Fluctuation Analysis (DFA)
  ## auxiliar function DFA_aux
  ## --------------------------------------------------------------------------##
  
  # Steps of the project:
  # 1. the function SSP requires a vector, which is the individual time series as for example the ENMO time series
  # file = IMP$metashort$ENMO
  # there is an example at the end
  
  
  # Declare functions
  
  ABI = function(x){
    #' @title Activity balance index (ABI)
    #' 
    #' @description This function estimates the Activity balance index (ABI), which is a transformation of the self-similarity parameter (SSP), also known as scaling exponent or alpha.  
    #' @param x the estimated self-similarity parameter (SSP)
    #'
    #' @details ABI = exp(-abs(SSP-1)/exp(-2))
    #'
    #' @return The estimated Activity balance index (ABI) is a real number between zero and one.
    #'
    #' @authors Ian Meneghel Danilevicz  
    #' 
    #' @references C.-K. Peng, S.V. Buldyrev, S. Havlin, M. Simons, H.E. Stanley, A.L. Goldberger Phys. Rev. E, 49 (1994), p. 1685
    #' Mesquita, Victor & Filho, Florencio & Rodrigues, Paulo. (2020). Detection of crossover points in detrended fluctuation analysis: An application to EEG signals of patients with epilepsy. Bioinformatics. 10.1093/bioinformatics/btaa955. 
    #'
    #' @examples 
    #' # Estimate Activity balance index of a very known time series available on R base: the sunspot.year.
    #' 
    #' ssp = SSP(sunspot.year, scale = 1.2)
    #' abi = ABI(ssp)
    #'
    if (is.na(x)) {
      y = NA
    }else{
      y = exp(-abs(x - 1) * exp(2))  
    }
    return(y)  
  }
  
  SSP = function(file, scale = 2^(1/8), box_size = 4, m = 1){
    #' @title Estimated self-similarity parameter
    #' 
    #' @description This function estimates the self-similarity parameter (SSP), also known as scaling exponent or alpha.  
    #' @usage alpha_hat(file,scale = 2^(1/8),box_size = 4,m=1) 
    #' @param file Univariate time series (must be a vector or data frame)
    #' @param scale Specifies the ratio between successive box sizes (by default scale = 2^(1/8))
    #' @param box_size Vector of box sizes (must be used in conjunction with scale = "F")
    #' @param m An integer of the polynomial order for the detrending (by default m=1)
    #'
    #' @details The DFA fluctuation can be computed in a geometric scale or for different choices of boxes sizes.
    #'
    #' @return Estimated alpha is a real number between zero and two.
    #' 
    #' @note It is not possible estimating alpha for multiple time series at once. 
    #'
    #' @authors Ian Meneghel Danilevicz  and Victor Mesquita 
    #' 
    #' @references C.-K. Peng, S.V. Buldyrev, S. Havlin, M. Simons, H.E. Stanley, A.L. Goldberger Phys. Rev. E, 49 (1994), p. 1685
    #' Mesquita, Victor & Filho, Florencio & Rodrigues, Paulo. (2020). Detection of crossover points in detrended fluctuation analysis: An application to EEG signals of patients with epilepsy. Bioinformatics. 10.1093/bioinformatics/btaa955. 
    #'
    #' @examples 
    #' # Estimate self-similarity of a very known time series available on R base: the sunspot.year.
    #' # Then the spend time with each method is compared.
    #' 
    #' SSP(sunspot.year, scale = 2)
    #' SSP(sunspot.year, scale = 1.2)
    #'
    #' time1 = system.time(SSP(sunspot.year, scale = 1.2))
    #' time2 = system.time(SSP(sunspot.year, scale = 2))
    #'
    #' time1
    #' time2  
    #'
    if (any(is.na(file))) {
      alpha_hat = NA
    } else {
      dfa_hat = DFA(as.vector(file), scale = scale, box_size = box_size, m = m)    
    }
    est_ols = lm(log(dfa_hat[,2]) ~ log(dfa_hat[,1]))
    alpha_hat = est_ols$coefficients[[2]]    
    return(alpha_hat)
  } 
  
  DFA_aux = function(j, box_size, ninbox2, file, y_k, m, N){
    aux_j = numeric(box_size[j] * ninbox2[j])
    fit = y_k
    for (i in seq_len(box_size[j] * ninbox2[j])) {
      if (i == 1) {
        aux_j[1] = box_size[j]
        mod_i = lm(y_k[i:aux_j[i]] ~ poly(c(i:aux_j[i]), m, raw = TRUE))
        fit[i:(aux_j[i])] = mod_i$fitted.values
      }
      if (i >= 2) {
        aux_j[i] = aux_j[i - 1] + box_size[j]
        mod_i = lm(y_k[(aux_j[i - 1] + 1):(aux_j[i])] ~ poly(c((aux_j[i - 1] +1):(aux_j[i])), m, raw = TRUE))
        fit[(aux_j[i - 1] + 1):(aux_j[i])] = mod_i$fitted.values
      }
      if (i >= ninbox2[j]) {
        aux_j[i] <- 0
      }
    }
    DFA = sqrt((1 / N) * sum((y_k[1:(box_size[j] * ninbox2[j])] -
                                fit[1:(box_size[j] * ninbox2[j])]) ^ 2))
    Results = c(round(box_size[j], digits = 0), DFA)
    return(Results)
  }
  
  DFA = function(file, scale = 2^(1/8), box_size = 4, m = 1){
    #' @title Detrended Fluctuation Analysis
    #' 
    #' @description ...
    #' @usage DFA(file, scale = 2^(1/8), box_size = 4, m = 1)
    #' @param file Univariate time series (must be a vector or data frame)
    #' @param scale Specifies the ratio between successive box sizes (by default scale = 2^(1/8))
    #' @param box_size Vector of box sizes (must be used in conjunction with scale = "F")
    #' @param m An integer of the polynomial order for the detrending (by default m=1)
    #'
    #' @details The DFA fluctuation can be computed in a geometric scale or for different choices of boxes sizes.
    #'
    #' @return Estimated alpha is a real number between zero and two.
    #' 
    #' @note It is not possible estimating alpha for multiple time series at once. 
    #'
    #' @authors Ian Meneghel Danilevicz  and Victor Mesquita 
    #' 
    #' @references C.-K. Peng, S.V. Buldyrev, S. Havlin, M. Simons, H.E. Stanley, A.L. Goldberger Phys. Rev. E, 49 (1994), p. 1685
    #' Mesquita, Victor & Filho, Florencio & Rodrigues, Paulo. (2020). Detection of crossover points in detrended fluctuation analysis: An application to EEG signals of patients with epilepsy. Bioinformatics. 10.1093/bioinformatics/btaa955. 
    #'
    #' @examples 
    #' # Estimate self-similarity of a very known time series available on R base: the sunspot.year.
    #' # Then the spend time with each method is compared.
    #' 
    if (class(file) == "data.frame") {
      file = file[, 1]
    }
    N = length(file)
    if (scale != "F") {
      box_size <- NULL
      n = 1
      n_aux <- 0
      box_size[1] <- 4
      for (n in 1:N) {
        while (n_aux < round(N/4)) {
          n_aux <- box_size[1]
          n = n + 1
          box_size[n] <- ceiling(scale * box_size[n - 1])
          n_aux <- box_size[n] + 4
        }
      }
    }
    ninbox2 <- NULL
    for (j in 1:length(box_size)) {
      ninbox2[j] <- N %/% box_size[j]
    }
    aux_seq = seq_along(ninbox2)  
    aux_length = aux_seq[length(aux_seq)]  
    y_k = cumsum(file) - mean(file)
    aux_mat = matrix(nrow = aux_length, ncol = 2)
    for (j in seq_along(ninbox2)) {
      aux_mat[j,] = DFA_aux(j, box_size, ninbox2, file, y_k, m, N)
    }
    colnames(aux_mat) <- c("box", "DFA")
    aux_list = aux_mat
    return(aux_list)
  }
  
 
  #------------------------------
  # apply functions
  #------------------------------
  ssp = SSP(ts, scale = scale)
  abi = ABI(ts)
  
  invisible(list(ssp = ssp, abi = abi))
}