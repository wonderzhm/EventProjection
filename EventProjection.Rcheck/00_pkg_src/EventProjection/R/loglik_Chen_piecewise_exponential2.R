#' @title Log-likelihood function for piecewise-exponential distribution with cured population
#' @description Provide log-likelihood function for piecewise-exponential distribution with cured population
#'
#' @param df The subject-level event data, including \code{time}
#'   and \code{event}.
#' @param piecewiseSurvivalTime A vector with length m that specifies the time
#'   intervals for the piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2
#'   event intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param par a vector with m elements denoting the log of the hazard rate in intervals
#' @param p cure rate
#'
#' @return
#' The negative value of the log-likelihood function given parameter
#' \code{par} and the dataset \code{df}
#'
#' @references 
#' \itemize{
#' \item Chen, Tai-Tsang. "Predicting analysis times in randomized clinical trials with cancer immunotherapy." 
#' BMC medical research methodology 16.1 (2016): 1-10.
#' }
#' 
#'
#' @export
#'
#'
#'
loglik_Chen_piecewise_exponential2<-function(par=NULL,df,piecewiseSurvivalTime, p){
  delta=df$event
  t=df$time
  u = piecewiseSurvivalTime[piecewiseSurvivalTime < max(t)]
  ucut = c(u, max(t))
  J = length(u)
  n=dim(df)[1]
  if(is.null(par)){
  par=rep(0,J)
  }
  par0=rep(NA,J+2)

  par0[1]=p
  par0[2:(J+1)]=exp(par[1:J])
  par0[J+2]=exp(par[J])
  time_interval_length=diff(ucut)

  sum_interval=c(0,cumsum((par0[2:(J+1)])*time_interval_length))
  interval_index=rep(NA,n)

  for(i in 1:n){
    interval_index[i]=findInterval(t[i],ucut)
  }

  S0=exp(-sum_interval[interval_index]-(t-ucut[interval_index])*par0[interval_index+1])
  f0=par0[interval_index+1]*S0

  part1=delta*(log(1-p)+log(f0))
  part2=(1-delta)*log(p+(1-p)*S0)
  neg_loglik=-sum(part1+part2)
  return(neg_loglik)
}
