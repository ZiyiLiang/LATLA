#' Thresholding function of LATLA
#'
#' This function first ranks all hypotheses according to the weighted
#' p-values and then determines a threshold along the ranking to control
#' the FDR. See algorithm 1 in reference paper for more details.
#'
#' @param pvs Vector of p-values under the null hypotheses.
#' @param pis The sparsity levels.
#' @param ws Weights for the p-values.
#' @param q FDR level.
#' @return A list of three: 1.number of rejections; 2.rejection threshold;
#'         3.indication of whether each index is rejected.
#' @importFrom stats density dnorm pnorm
#' @export
latla_thres<-function(pvs, pis, ws, q)
{
  ## implementing neat
  ## Arguments
  # pvs: p-values
  # pis: conditional probabilities
  # ws: weights proposed by neat
  # q: FDR level
  ## Values
  # de: the decision
  # th: the threshold for weighted p-values

  m<-length(pvs)
  pws<-pvs/ws
  st.pws<-sort(pws)
  fdps<-sum(ws*(1-pis))*st.pws/(1:m)
  de<-rep(0, m)
  if(sum(fdps<=q)==0)
  {
    k<-0
    pwk<-1
  }
  else
  {
    k<-max(which(fdps<=q))
    pwk<-st.pws[k]
    de[which(pws<=pwk)]<-1
  }
  y<-list(nr=k, th=pwk, de=de)
  return (y)
}

#' Data-driven weights of LATLA
#'
#' The function calculates the oracle-assisted data-driven weights described
#' in the algorithm 2 in LATLA paper.
#' @param x Vector of the primary statistics.
#' @param d Distance matrix.
#' @param pis Estimated sparsity level.
#' @param mu0 Vector of means under the null hypotheses.
#' @param sd0 Vector of standard deviations under the null hypotheses.
#' @param q FDR level, default value is set to be 0.05.
#' @param h Bandwidth for kernel estimation, default is set to 'auto', where
#'          the bandwidth is automatically selected by build-in method 'ste'.
#' @param eps Determines the size of the neighborhood, only m^(1-eps) data
#'            points will be used for kernel estimation. eps should range
#'            from 0 to 1.
#' @return The oracle-assisted weights.
#' @export
############## Data-driven Weights #####################################
latla_weights <- function(x, d, pis, mu0, sd0, q=0.05, h='auto', eps=0.1){
  m <- length(x)
  # eps should be in [0,1], check validity
  if ((eps < 0) | (eps > 1)){
    stop('Invalid eps, eps should be in [0,1].')
  }
  nb <- floor(m^(1-eps))
  # check default argument of bandwidth
  if (h=='auto'){
    zds=density(x,bw="SJ-ste")
    h=zds$bw
  }

  tor.est <- rep(0, m)
  for (i in 1:m) {
    vhi <- dnorm(d[i,], 0, h)
    vhi[i] = 0
    kht <- dnorm(x-x[i], 0, h)
    # estimated conditional density evaluated at ti
    fti <- sum(vhi*kht)/sum(vhi)
    tor.est[i] <- (1-pis[i])*dnorm(x[i],mu0[i],sd0[i])/fti
  }
  st.tor <- sort(tor.est)

  # calculating the moving average
  mva <- rep(0,m)
  for (i in 1:m) {
    mva[i] <- sum(st.tor[1:i])/i
  }
  # calculating the oracle threshold
  th <- max(which(mva<=q))
  th <- st.tor[th]

  # calculating the weights
  weights <- rep(0,m)
  pt <- seq(0, max(x), 0.1)
  nt <- seq(min(x), 0, 0.1)
  t_star <- rep(0,m)
  n_p <- length(pt)
  n_n <- length(nt)
  for(i in 1:m){
    vhi <- dnorm(d[i,], 0, h)
    vhi[i]=0
    # only consider the test points in the neighborhood
    if (nb < m){
      st.nb<- sort(vhi)[m-nb]
      vhi[which(vhi<st.nb)]<-0
    }
    # adapt the asymmetricity
    if(x[i]>=0){
      pt_mat <- matrix(rep(pt,m),m,n_p,byrow = TRUE)
      x_mat <- matrix(rep(x,n_p),m,n_p)
      kht <- dnorm(x_mat - pt_mat, 0, h)
      ft <- c(vhi %*% kht)/sum(vhi)
      tor<-(1-pis[i])*dnorm(pt,mu0[i],sd0[i])/ft
      tor[which(tor>1)]<-1
      if(length(which(tor<=th))==0){
        t_star[i] <- max(x)
      }
      else{
        t_star[i] = pt[min(which(tor<=th))]
      }
      weights[i] <- 1-pnorm(t_star[i], mu0[i], sd0[i])
    }
    else{
      nt_mat <- matrix(rep(nt,m),m,n_n,byrow = TRUE)
      x_mat <- matrix(rep(x,n_n),m,n_n)
      kht <- dnorm(x_mat - nt_mat, 0, h)
      ft <- c(vhi %*% kht)/sum(vhi)
      tor<-(1-pis[i])*dnorm(nt,mu0[i],sd0[i])/ft
      tor[which(tor>1)]<-1
      if(length(which(tor<=th))==0){
        t_star[i] <- min(x)
      }
      else{
        t_star[i] = nt[max(which(tor<=th))]
      }
      weights[i] <- pnorm(t_star[i], mu0[i], sd0[i])
    }
  }
  nu<-10e-5
  weights[which(weights<nu)]<-nu # stabilization
  return(weights)
}

#' Estimated sparsity level
#'
#' The is a kernel-based function that estimates the local sparsity levels
#' which are subsequently used to construct oracle-assisted weights.
#' @param x Vector of the primary statistics.
#' @param d Distance matrix.
#' @param pval p-values of the primary statistics
#' @param tau A threshold to approximate the 'null set', roughly speaking,
#'            we deem p-values greater than tau as null points. Default is
#'            set to 0.9, see reference paper for more discussions on
#'            parameter selection.
#' @param h Bandwidth for kernel estimation, default is set to 'auto',
#'          where the bandwidth is automatically selected by build-in
#'          method 'ste'(solve the equation).
#' @param eps Determines the size of the neighborhood, only m^(1-eps) data
#'            points will be used for kernel estimation. eps should range
#'            from 0 to 1.
#' @return The estimated sparsity levels
#' @export
latla_pis <- function(x, d, pval, tau=0.9, h='auto', eps=0.1)
{
  m <- length(x)
  # eps should be in [0,1], check validity
  if ((eps < 0) | (eps > 1)){
    stop('Invalid eps, eps should be in [0,1].')
  }
  nb <- floor(m^(1-eps))
  if (h=='auto'){
    zds=density(x,bw="SJ-ste")
    h=zds$bw
  }

  p.est <- rep(0, m)
  for (i in 1:m){
    kht <- dnorm(d[i,], 0, h)
    kht[i] = 0 # considers only j!=i
    if (nb < m){
      st.nb<- sort(kht)[m-nb]
      kht[which(kht<st.nb)]<-0
    }
    p.est[i] <- sum(kht[which(pval>=tau)])/((1-tau)*sum(kht))
  }
  p.est[which(p.est>1)] <- 1
  return(1-p.est)
}
