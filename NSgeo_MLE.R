
NSgeo.MLE<- function(formula,coords,rand.effect.domain,
                                          data,ID.coords,
                                          kappa,fixed.rel.nugget=NULL,start.cov.pars,
                                          method="BFGS",messages=TRUE,
                                          return_se=FALSE) {
  
  start.cov.pars <- as.numeric(start.cov.pars)
  if(any(start.cov.pars<0)) stop("start.cov.pars must be positive.")
  kappa <- 0.5
  methods = "nlminb"
  if(any(method==c("BFGS","nlminb"))==FALSE) stop("method must be either BFGS or nlminb.")
  mf <- model.frame(formula,data=data)
  y <- as.numeric(model.response(mf))
  
  m <- length(y)
  
  D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))
  coords.aux <- as.matrix(model.frame(coords,data))
  rand.effect.domain.aux <- as.numeric(
    as.matrix(
      model.frame(rand.effect.domain,data))[,1])
  p <- ncol(D)
  
  coords <- unique(cbind(coords.aux,rand.effect.domain.aux))[,1:2]
  rand.eff.domain <- unique(cbind(coords.aux,rand.effect.domain.aux))[,3]
  U.s <- dist(coords)
  U.re <- dist(rand.eff.domain)
  
  n.coords <- as.numeric(tapply(ID.coords,ID.coords,length))
  DtD <- t(D)%*%D
  Dty <- as.numeric(t(D)%*%y)
  D.tilde <- sapply(1:p,function(i) as.numeric(tapply(D[,i],ID.coords,sum)))
  y.tilde <- tapply(y,ID.coords,sum)
  
  beta.start <- coef(lm(formula,data=data))
  
  mu.start <- as.numeric(D%*%beta.start)
  omega2.start <- mean((y-mu.start)^2)
  sigma2.start <- 0.5*omega2.start
  
  phi.start <- start.cov.pars[1]
  phi.re.start <- start.cov.pars[2]
  
  start.par <- log(c(phi.start,phi.re.start,
                     omega2.start/sigma2.start))
  
  profile.log.lik <- function(par) {
    phi <- exp(par[1])
    phi.re <- exp(par[2])
    omega2.star <- exp(par[3])
    
    R.s <- varcov.spatial(dists.lowertri=U.s,cov.pars=c(1,phi),
                          kappa=0.5)$varcov
    R.re <- varcov.spatial(dists.lowertri=U.re,
                           cov.pars=c(1,phi.re),
                           kappa=0.5)$varcov
    R <- R.s*R.re
    R.inv <- solve(R)
    Omega <- R.inv
    diag(Omega) <- (diag(Omega)+n.coords/omega2.star)
    
    Omega.inv <- solve(Omega)
    M.beta <- DtD/omega2.star+
      -t(D.tilde)%*%Omega.inv%*%D.tilde/(omega2.star^2)
    v.beta <- Dty/omega2.star-t(D.tilde)%*%Omega.inv%*%y.tilde/(omega2.star^2)
    beta.hat <- as.numeric(solve(M.beta)%*%v.beta)
    
    M.det <- t(R*n.coords/omega2.star)
    diag(M.det) <- diag(M.det)+1
    
    mu.hat <- as.numeric(D%*%beta.hat)
    diff.y <- y-mu.hat
    diff.y.tilde <- tapply(diff.y,ID.coords,sum)
    sigma2.hat <- as.numeric((sum(diff.y^2)/omega2.star-
                                t(diff.y.tilde)%*%Omega.inv%*%diff.y.tilde/(omega2.star^2))/
                               m)
    
    out <- -0.5*(m*log(sigma2.hat)+m*log(omega2.star)+
                   as.numeric(determinant(M.det)$modulus))
    return(out)
  }
  
  est.profile <- nlminb(start.par,
                        function(x) -profile.log.lik(x),
                        control=list(trace=1*messages))
  
  phi <- exp(est.profile$par[1])
  phi.re <- exp(est.profile$par[2])
  omega2.star <- exp(est.profile$par[3])
  
  R.s <- varcov.spatial(dists.lowertri=U.s,cov.pars=c(1,phi),
                        kappa=0.5)$varcov
  R.re <- varcov.spatial(dists.lowertri=U.re,
                         cov.pars=c(1,phi.re),
                         kappa=0.5)$varcov
  R <- R.s*R.re
  R.inv <- solve(R)
  Omega <- R.inv
  diag(Omega) <- (diag(Omega)+n.coords/omega2.star)
  
  Omega.inv <- solve(Omega)
  M.beta <- DtD/omega2.star+
    -t(D.tilde)%*%Omega.inv%*%D.tilde/(omega2.star^2)
  v.beta <- Dty/omega2.star-t(D.tilde)%*%Omega.inv%*%y.tilde/(omega2.star^2)
  beta.hat <- as.numeric(solve(M.beta)%*%v.beta)
  
  M.det <- t(R*n.coords/omega2.star)
  diag(M.det) <- diag(M.det)+1
  
  mu.hat <- as.numeric(D%*%beta.hat)
  diff.y <- y-mu.hat
  diff.y.tilde <- tapply(diff.y,ID.coords,sum)
  sigma2.hat <- as.numeric((sum(diff.y^2)/omega2.star-
                              t(diff.y.tilde)%*%Omega.inv%*%diff.y.tilde/(omega2.star^2))/
                             m)
  
  est.profile$par <- c(beta.hat,sigma2=sigma2.hat,phi=phi,phi=phi.re,
                       omega2=omega2.star*sigma2.hat)
  if(return_se) {
    log.lik <- function(par) {
      beta <- par[1:p]
      sigma2 <- exp(par[p+1])
      phi <- exp(par[p+2])
      phi.re <- exp(par[p+3])
      omega2.star <- exp(par[p+4])
      
      R.s <- varcov.spatial(dists.lowertri=U.s,cov.pars=c(1,phi),
                            kappa=0.5)$varcov
      R.re <- varcov.spatial(dists.lowertri=U.re,
                             cov.pars=c(1,phi.re),
                             kappa=0.5)$varcov
      R <- R.s*R.re
      R.inv <- solve(R)
      Omega <- R.inv
      diag(Omega) <- (diag(Omega)+n.coords/omega2.star)
      
      Omega.inv <- solve(Omega)
      M.beta <- DtD/omega2.star+
        -t(D.tilde)%*%Omega.inv%*%D.tilde/(omega2.star^2)
      v.beta <- Dty/omega2.star-t(D.tilde)%*%Omega.inv%*%y.tilde/(omega2.star^2)
      beta.hat <- as.numeric(solve(M.beta)%*%v.beta)
      
      M.det <- t(R*n.coords/omega2.star)
      diag(M.det) <- diag(M.det)+1
      
      mu <- as.numeric(D%*%beta)
      diff.y <- y-mu
      diff.y.tilde <- tapply(diff.y,ID.coords,sum)
      sigma2.hat <- as.numeric((sum(diff.y^2)/omega2.star-
                                  t(diff.y.tilde)%*%Omega.inv%*%diff.y.tilde/(omega2.star^2)))
      
      out <- -0.5*(m*log(sigma2)+m*log(omega2.star)+
                     as.numeric(determinant(M.det)$modulus)+
                     sigma2.hat/sigma2)
      return(as.numeric(out))
    }
    library(numDeriv)
    cat("Computing standard errors (this may be time demanding) \n")
    H <- numDeriv::hessian(function(x) log.lik(x),c(beta.hat,log(sigma2.hat),
                                                    log(phi),log(phi.re),
                                                    log(omega2.star)))
    est.profile$se <- sqrt(diag(solve(-H)))
  }
  return(est.profile)
} 