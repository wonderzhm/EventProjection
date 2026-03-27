#' @title Log-likelihood function for log-logistic distribution with cured population
#' @description Provide log-likelihood function for log-logistic distribution with cured population
#'
#' @param df The subject-level event data, including \code{time}
#'   and \code{event}.
#' @param par a vector with two elements denoting 
#'  the log of shape and the log of the scale parameter of the log-logistic distribution.
#' @param p cure rate
#'
#' @return
#' The negative value of the log-likelihood function given parameter
#' \code{par} and the dataset \code{df}
#'
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
loglik_Chen_log_logistic2<-function(par,df,p){
  a=exp(par[1])
  b=exp(par[2])
  delta=df$event
  t=df$time

  f0=a/b*(t/b)^(a-1)/((1+(t/b)^a)^2)
  S0=1/(1+(t/b)^a)

  ###
  part1=delta*(log(1-p)+log(f0))
  part2=(1-delta)*log(p+(1-p)*S0)
  neg_loglik=-sum(part1+part2)
  return(neg_loglik)
}

