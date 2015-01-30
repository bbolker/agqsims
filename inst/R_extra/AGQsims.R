library("agqsims")
f0 <- fit_gen(data=simfun_culc(n.ttt=2,beta=c(0,2)),
              formula= y~ttt,
              family="binomial",
              cluster="block",
              pkg="glmmPQL",
              maxAGQ=0,truevals=c(0,2,3.5))

AGQvec <-c(0)
nrepvec<- c(3,4,5,6,7,8,9,10,11,12,17,25,35,50,70,100,150,200)
thetavec <- 10^seq(-1,1,length=11)
nsims <- 100
dims <- c(nrep=length(nrepvec),
          theta=length(thetavec),
          sim=nsims,AGQ=length(AGQvec),
          var=dim(f0)[2],
          type=dim(f0)[3])
prod(dims[1:4])

sim_culcita_pql <- array(NA,dim=dims)
dimnames(sim_culcita_pql) <- list(
    nrep=nrepvec,theta=thetavec,
    sim=1:nsims,AGQ=AGQvec,var=dimnames(f0)[["var"]],
    type=dimnames(f0)[["type"]])
beta <- c(0,2)

for (i in seq_along(nrepvec) ) {
    nrep <- nrepvec[i]
    for (j in seq_along(thetavec)) {
        theta <- thetavec[j]
        for (k in 1:nsims) {
            cat(i,j,k,"\n")
            set.seed(1000+k)
            sim_culcita_pql[i,j,k,,,] <-
                fit_gen(data=simfun_culc(n.ttt=2,
                        n.rep=nrep,
                        beta=beta,
                        theta=theta),
                        nvar=4,
                        formula= y~ttt,
                        cluster="block",
                        family="binomial",
                        pkg="glmmPQL",
                        maxAGQ=0,
                        truevals=c(beta,theta))
        }
    }
}

f0 <- fit_gen(data=simfun_culc(n.ttt=2,beta=c(0,2)),
              formula= y~ttt+(1|block),
              family="binomial",
              maxAGQ=1,truevals=c(0,2,3.5))
AGQvec <-c(1,2,3,5,8,10,15,20,25,75,100)
nsims <- 2
nrepvec<- c(3,4,5,6,7,8,9,10,11,12,17,25,35,50,70,100,150,200)
thetavec <- 10^seq(-1,1,length=11)
dims <- c(nrep=length(nrepvec),
          theta=length(thetavec),
          sim=nsims,AGQ=length(AGQvec),
          var=dim(f0)[2],
          type=dim(f0)[3])
prod(dims[1:4])  ## total number of sims (735)

sim_culcita_lme4 <- array(NA,dim=dims)
dimnames(sim_culcita_lme4) <- list(
    nrep=nrepvec,theta=thetavec,
    sim=1:nsims,AGQ=AGQvec,var=dimnames(f0)[["var"]],
    type=dimnames(f0)[["type"]])
beta <- c(0,2)
for (i in seq_along(nrepvec) ) {
    nrep <- nrepvec[i]
    for (j in seq_along(thetavec)) {
        theta <- thetavec[j]
        for (k in 1:nsims) {
            cat(i,j,k,"\n")
            set.seed(1000+k)
            sim_culcita_lme4[i,j,k,,,] <-
                fit_gen(data=simfun_culc(n.ttt=2,
                        n.rep=nrep,
                        beta=beta,
                        theta=theta),
                        nvar=4,
                        formula= y~ttt+(1|block),
                        family="binomial",
                        AGQvec=AGQvec,
                        truevals=c(beta,theta))
        }
    }
}


