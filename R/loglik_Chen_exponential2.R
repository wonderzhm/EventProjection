#' @title Log-likelihood function for exponential distribution with cured population
#' @description Provide log-likelihood function for exponential distribution with cured population. 
#'
#' @param df The subject-level event data, including \code{time}
#'   and \code{event}.
#' @param par the log of the hazard rate.
#' @param p cure rate
#'
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
loglik_Chen_exponential2<-function(par,df,p){
  sigma=exp(-par)
  delta=df$event
  t=df$time

  f0=1/sigma*exp(-(t/sigma))
  S0=exp(-(t/sigma))


  part1=delta*(log(1-p)+log(f0))
  part2=(1-delta)*log(p+(1-p)*S0)
  neg_loglik=-sum(part1+part2)
  return(neg_loglik)
}
