#' @title Log-likelihood function for log-normal distribution with cured population
#' @description Provide log-likelihood function for log-normal distribution with cured population
#'
#' @param df The subject-level event data, including \code{time}
#'   and \code{event}.
#' @param par a vector with two elements denoting
#'  the mean and the log of standard deviation parameter of the log-normal distribution.
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
loglik_Chen_log_normal2<-function(par,df,p){
  mu=par[1]
  sigma=exp(par[2])
  delta=df$event
  t=df$time

  f0=1/(sqrt(2*pi)*sigma*t)*exp(-(log(t)-mu)^2/(2*sigma^2))
  S0=1-stats::pnorm(log(t),mean=mu,sd=sigma)


  part1=delta*(log(1-p)+log(f0))
  part2=(1-delta)*log(p+(1-p)*S0)
  neg_loglik=-sum(part1+part2)
  return(neg_loglik)
}
