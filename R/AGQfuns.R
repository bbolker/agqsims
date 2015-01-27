##' simulation code, based approximately on Culcita data set
##' not longitudinal; approximates toenail if n.ttt=2
##' simulates covariate x, but not currently used

##' @param n.blocks number of levels of RE grouping variable
##' @param n.ttt number of levels of (categorical) fixed effect
##' @param n.rep number of replicates per ttt*block combination
##' @param N.binom number of trials per binomial sample
##' @param x.range range of continuous covariate (unused)
##' @param beta fixed-effects parameter vector
##' @param theta RE parameter vector (Cholesky factor)
##' @param seed random-number seed
##' @export
simfun_culc <- function(n.blocks=10,
                        n.ttt=4,
                        n.rep=3,N.binom=1,
                        x.range=c(-1,1),
                        ## intercept, ttt effects, slope
                        ## default gives homogeneous log-odds of 0 -> prob=0.5
                        beta=c(5,-3.75,-4.4,-5.5),
                        theta=3.5,
                        seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    dd <- expand.grid(block=factor(seq(n.blocks)),
                      ttt=factor(letters[seq(n.ttt)]),
                      rep=seq(n.rep))
    ## note: not used
    dd$x <- runif(nrow(dd),min=x.range[1],max=x.range[2])
    dd$N <- N.binom
    ## the messages about "theta/beta parameter vector not named"
    ##  are annoying ...
    sim.ok <- FALSE
    while (!sim.ok) {
        dd$y <- suppressMessages(
            simulate(~ttt+(1|block),family="binomial",
                     size=N.binom,
                     newdata=dd,
                     newparams=list(beta=beta,theta=theta))[[1]])
        dtab <- with(dd,tapply(y,list(block,ttt),FUN=sum))
        sim.ok <- !any(apply(dtab,2,function(x) (all(x==n.rep) || all(x==0))))
    }
    return(dd)
}

##' Wald variance-covariance matrix of lme4 RE parameters
##'
##' @param fit an lme4 fit
##' @param \dots required for method compatibility
##' @details by brute force, finds the Wald var-cov matrix
##' of lme4 RE parameters (on the 'user', i.e. sd-corr-sigma scale)
##' @export
##' @importFrom numDeriv hessian
##' @importFrom lme4 refitML
##' @importFrom lme4 isREML
##' @importFrom lme4 isGLMM
vcov.VarCorr.merMod <- function(fit,...) {
    ## unfortunately, since we need access to the original fit,
    ## we can't *actually* make this a method that applies to
    ## VarCorr.merMod objects ...
    object <- VarCorr(fit)
    if (isREML(fit)) {
        warning("refitting model with ML")
        fit <- refitML(fit)
    }
    useSc <- attr(object,"useSc")
    dd <- lme4:::devfun2(fit,useSc=useSc,signames=FALSE)
    vdd <- as.data.frame(object,order="lower.tri")
    pars <- vdd[,"sdcor"]
    npar0 <- length(pars)
    if (isGLMM(fit)) {
        pars <- c(pars,fixef(fit))
    }
    hh1 <- hessian(dd,pars)
    vv2 <- 2*solve(hh1)
    if (isGLMM(fit)) {
        vv2 <- vv2[1:npar0,1:npar0,drop=FALSE]
    }
    nms <- apply(vdd[,1:3],1,
                 function(x) paste(na.omit(x),collapse="."))
    dimnames(vv2) <- list(nms,nms)
    return(vv2)
}

## OBSOLETE
sumfun <- function(fitted) {
  if (class(fitted)=="glmerMod") {
    c(fixef(fitted),REvar=sqrt(unlist(VarCorr(fitted))),
      logLik=logLik(fitted))
  } else if (class(fitted)=="glmmML") {
    c(fitted$coefficients,
      REvar=fitted$sigma,logLik=-fitted$deviance/2)
  }  else  {
    c(fixef(fitted),REvar=as.numeric(VarCorr(fitted)[1,2]),
      logLik <- NA)
  }
}


##' computes summary statistics for
##' 
##' @param fitted a fitted GLMM model (glmer, glmmML, or glmmPQL)
##' @param truevals optional vector of true parameter values (beta, RE sd)
##' @export
sumfun2 <- function(fitted,truevals=NULL) {
  if (is(fitted,"glmerMod")) {
    pars <- c(fixef(fitted),REvar=sqrt(unlist(VarCorr(fitted))))
    ss <- summary(fitted)
    ## if (is.function(family)) { ## HACK
    ## family <- fitted@resp$family$family
    ## }
    stders <- c(coef(ss)[,"Std. Error"],
                sqrt(diag(vcov.VarCorr.merMod(fitted))))
    logLik <- logLik(fitted)
  } else if (is(fitted,"glmmML")) {
    pars <- c(fitted$coefficients,
              REvar=fitted$sigma)
    stders <- c(rep(fitted$coef.sd,
                    length.out=length(fitted$coefficients)),
                ## rep() is minor hack for cases where fitted$coef.sd
                ##   is NA ...
                fitted$sigma.sd)
    logLik <- -fitted$deviance/2
  } else if (is(fitted,"glmmPQL"))   {
      pars <- c(fixef(fitted),REvar=as.numeric(VarCorr(fitted)[1,2]))  
      stders <- c(summary(fitted)$tTable[,"Std.Error"],
                  ## transform Wald SD of log-sd to Wald SD of sd ...
                  sqrt(fitted$apVar[1,1])*pars["REvar"])
      logLik <- NA  ## not available for PQL fit
  }
  r <- cbind(c(pars,logLik=logLik),c(stders,NA))
  colnames(r) <- c("est","se")
  if (!is.null(truevals)) {
    r <- cbind(r,c(truevals,NA))
    colnames(r) <- c("est","se","true")
  }
  return(r)
}


################# fit_lme4 ########

sseq <- function(n) {
    if (n==0) 0 else 1:n
}

##' fit GLMM for varying values
##'
##' @param data data
##' @param formula model formula (fixed-effect or full, depending on pkg)
##' @param family model family
##' @param pkg package to use for fitting
##' @param cluster cluster variable specification (for glmmPQL or glmmML)
##' @param maxAGQ maximum AGQ value (0=PQL)
##' @param AGQvec vector of specified AGQ orders
##' @param verbose print verbose output
##' @param truevals vector of true parameter values
##'
##' @export
##' @importFrom MASS glmmPQL
fit_gen <- function(data,formula,family,
                    cluster=NULL,
                    pkg=c("lme4","glmmML","glmmPQL"),
                    maxAGQ=10,
                    AGQvec=sseq(maxAGQ),
                    verbose=TRUE,
                    truevals = NULL) {
  pkg <- match.arg(pkg)

  ## Preliminary fit to get baseline shape of results
  ## could do better by figuring out what sumfun2() SHOULD return
  if (pkg=="glmmML") {
    ## silly hacks required to deal with the way that
    ## glmmML evaluates its arguments
    assign("data",data,envir=environment(formula))
    assign("formula",formula,envir=environment(formula))
    assign("cluster",cluster,envir=environment(formula))
    fit0 <- glmmML(formula,
                   family= family,
                   data = data,
                   cluster=eval(parse(text=paste0("data$",
                                                  as.name(cluster)))),
                   method="ghq" ,n.points = 1)
  } else if (pkg=="lme4") {
    ### ?? why necessary ???
    ## assign("data",data,envir=environment(formula))
    ## assign("formula",formula,envir=environment(formula))
    fit0 <- glmer(formula,
                  family= family, data = data,
                  nAGQ = 1)
    ## hack for vcov.VarCorr.merMod
    fit0@call$family <- family
    fit0@call$data <- data
  } else {
    fit0<- glmmPQL(formula,
                   family=family,
                   data=data,
                   random=formula(paste0("~1|",as.name(cluster))))
  }
  
  fit0_sum <- sumfun2(fit0,truevals=truevals)
  nvar <- nrow(fit0_sum)
  res1 <- array(NA,dim=c(length(AGQvec),nvar+2,ncol(fit0_sum)),
                dimnames=list(AGQ=AGQvec,
                              var=c(rownames(fit0_sum),
                                    c("t.user","t.elapsed")),
                              type=colnames(fit0_sum)))
  for(j in seq_along(AGQvec)) {
      if (verbose) cat(j,AGQvec[j],"\n")
      if (pkg=="lme4") {
          st <- system.time(fit1 <-
              try(glmer(formula,
                        family= family, data = data,
                        nAGQ = AGQvec[j])))
          fit1@call$family <- family  ## hack (see above)
          fit1@call$nAGQ <- AGQvec[j]  ## hack (see above)
          fit1@call$data <- data
      } else if (pkg=="glmmML"){
          st <- system.time(fit1 <-
              try(glmmML(formula,
                         family= family, data = data,
                         cluster=eval(parse(text=
                                                paste0("data$",as.name(cluster)))), 
                         method="ghq" ,n.points = AGQvec[j])))
      } else  {
          if (AGQvec[j]>0) stop("shouldn't be trying glmmPQL with AGQ>0")
          st <- system.time(fit1 <-
              try(glmmPQL(formula,
                          family=family,
                          data=data,
                          random=as.formula(paste("~1|",as.name(cluster))))))
      }
    
      ## slightly hacky way to avoid NA Hessian issues
      if (is(fit1,"try-error") || 
          (pkg=="lme4" && any(is.na(fit1@optinfo$derivs$Hessian)))) {
          res1[j,,] <- NA                                         
      } else {
          ss2 <- sumfun2(fit1,truevals=truevals)
          st2 <- matrix(c(st[c(1,3)],rep(NA,2*(ncol(ss2)-1))),
                        nrow=2)
          res1[j,,] <- rbind(ss2,st2)
      }
  }
  return(res1)
}




########### Graph 1 ############

##' @importFrom reshape2 melt
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 facet_wrap
##' @importFrom ggplot2 labs
graph1 <- function( d) {
    t_d <- t(d[,-which(colnames(d) %in% c("logLik","t.user")),
               type="est"])
    melt_d <- melt(t_d)
    p1 <- ggplot(melt_d,aes(x=AGQ,y=value))+
        geom_line()+
            geom_point()+facet_wrap(~var,scale="free",nrow=1)+
                
                labs(x="Number of Adaptive Gauss-Hermite quadrature points")
    return(p1)
}

graph_SAS <- function(d) {
    melt_d <- melt(d,id.vars="AGQ")   
    p1 <- ggplot(melt_d,aes(x=AGQ,y=value))+
        geom_line()+
            geom_point()+facet_wrap(~variable,scale="free",nrow=1)+
                
                labs(x="Number of Adaptive Gauss-Hermite quadrature points")
    return(p1)
    
}

##'  center and scale by 'correct' answer (max AGQ)
##' @param d a 3-dimensional numeric array: { AGQ * var * type }
##' @param v name of random-effects variance
##' @param skipcols variables to skip
##' @param std how to standardize
##' @param ret what summary to return
##' @export
standardize_all <- function(d,v,
                            skipcols=c("(Intercept)",
                            "logLik","t.user","t.elapsed"),
                            std=c("max","se","true"),
                            ret=c("rmse","bias","var")) {
    std <- match.arg(std)
    ret <- match.arg(ret)
    d2 <- d[ , -which(dimnames(d)$var %in% skipcols), ]
    mel_d <- dcast(melt(d2),AGQ+var~type)
    mel2_d <- transform(mel_d ,
                        varcat=ifelse(as.character(var) %in% v,
                        "RE","FE"))
    sumf <- switch(std,
                   max =function(x) transform(x,stdvalue=(est-tail(est,1))/tail(est,1)),
                   se  =function(x) transform(x,stdvalue=(est-tail(est,1))/tail(se,1)),
                   true=function(x) transform(x,stdvalue=(est-tail(true,1))/tail(true,1)))
    mel3_d <- ddply(mel2_d,"var",sumf)
    ## browser()
    sumf2 <- switch(ret,
                    rmse=function(x) summarise(x,sumvalue=sqrt(mean(stdvalue^2))),
                    bias=function(x) summarise(x,sumvalue=mean(stdvalue)))
    mel3_sum_d <- ddply(mel3_d,c("varcat","AGQ"),sumf2)
    return(mel3_sum_d)
}
## FIXME: should split standardization -- what to subtract vs
##  what to divide by
## want to be able to use this for real data *or* sim data,


##' compute FE and RE summaries
##' @param x a 3-dimensional numeric array: { AGQ * vars * type }
##' @param v name of RE variance
##' @param std how to standardize
##' @param ret what to return
sim_to_rmsetab <- function(x,v,std=NULL,ret="rmse") {
    if (is.null(std)) {
        tmp0 <- standardize_maxval(x,v=v)
    } else {
        tmp0 <- standardize_all(x,v=v,std=std,ret=ret)
    }
    ## rearrange to make FE and RE separate columns
    tmp1 <- dcast(tmp0,AGQ~varcat,value.var=names(tmp0)[3])
    ## convert to a matrix (drop first column)
    tmp2 <- as.matrix(tmp1[,-1])
    ## ... and restore dimnames
    dimnames(tmp2) <- list(AGQ=tmp1$AGQ,varcat=colnames(tmp1)[-1])
    return(tmp2)
}

##' compute RMSE and return a reduced table
##' @importFrom abind abind
##' @param S a 4-dimensional numeric array: { sim * AGQ * vars * type }
##' @param v ""
##' @param std ""
##' @param ret ""
##' @param dims ""
##' @param calc.range ""
rmsetab_to_combdf <- function(S,v,std=NULL,ret="rmse",dims=2:3,
                              calc.range=TRUE) {
    rmsetab <- aaply(S,1,sim_to_rmsetab,v=v,std=std,ret=ret)
    ## returns array { sim * AGQ * varcat } 
    ## meantab <- apply(rmsetab,c(2,3),mean,na.rm=TRUE)
    mediantab <- apply(rmsetab,c(2,3),median,na.rm=TRUE)
    if (!calc.range) {
        combtab <- mediantab
        r <- melt(mediantab)
        colnames(r)[3] <- "median"
    }  else {
        mintab <- apply(rmsetab,c(2,3),quantile,0.025,na.rm=TRUE)
        maxtab <- apply(rmsetab,c(2,3),quantile,0.975,na.rm=TRUE)
        combtab <- abind(median=mediantab,
                         min=mintab,max=maxtab,along=0)
        ## restore names(dimnames)
        names(dimnames(combtab)) <- c("type","AGQ","varcat")
        r <- dcast(melt(combtab),AGQ+varcat~type)
    }
    return(r)
    ## return(melt(rmsetab))
}


standardize_maxval <- function(d,v,
                               skipcols=c("(Intercept)",
                               "logLik","t.user","t.elapsed")) {

    d2 <- d[, -which(colnames(d) %in% skipcols),
            type="est"]
    mel_d <- melt(d2)
    mel2_d <- transform(mel_d,
                        varcat=ifelse(as.character(var) %in% v,
                        "RE","FE"))
    mel3_d <- ddply(mel2_d,"var",
                    transform, stdvalue=value/tail(value,1)-1)
    mel3_sum_d <- ddply(mel3_d,c("varcat","AGQ"),
                        summarise, rmse_stdvalue=sqrt(mean(stdvalue^2)))
    return(mel3_sum_d)
}

########### Graph 2############
##  center by 'correct' answer (max AGQ) and scale by max-AGQ SD

standardize_maxsd <- function(d,v,skipcols=c("(Intercept)",
                                  "logLik","t.user","t.elapsed")) {

    d2 <- d[, -which(colnames(d) %in% skipcols)]
    mel_d <- dcast(melt(d2),AGQ+var~type)
    mel2_d <- transform(mel_d,
                        varcat=ifelse(as.character(var) %in% v,
                        "RE","FE"))
    mel3_d <- ddply(mel2_d,"var",
                    transform,
                    stdvalue=(est-tail(est,1))/tail(se,1))
    mel3_sum_d <- ddply(mel3_d,c("varcat","AGQ"),
                        summarise, rmse_stdvalue_SE=sqrt(mean(stdvalue^2)))
    return(mel3_sum_d)
}

