######### Generic data analysis  ###############

################## Toenail data ##############

## load data
data("toenail", package="DPpackage")

## fit for every package

## try to save some time by sampling more sparsely at higher AGQ
agqvec <- c(1:10,seq(11,30,by=2),seq(30,50,by=5),seq(60,100,by=10))
agqvec_glmmML <- agqvec[agqvec<75]

toenail_glmmML <- fit_gen(data=toenail,formula=outcome~treatment+visit,
                          cluster="ID",family="binomial",
                          pkg="glmmML",
                          agqvec=agqvec_glmmML)

## hmmm. slow ...
toenail_lme4 <- fit_gen(data=toenail,formula=outcome~treatment+visit+(1|ID),
                        family="binomial",
                        pkg="lme4",
                        agqvec=agqvec)

## input names of csv files and generate the same output
## toenail_SAS <- 


########### Culcita data #############

culcita_glmmML <- fit_gen(data=culcita_dat,formula=predation~ttt,
                          cluster="block",family="binomial",
                          pkg="glmmML",
                          agqvec=agqvec_glmmML)

culcita_lme4 <- fit_gen(data=culcita_dat, formula= predation~ttt+(1|block),
                        family="binomial",
                        pkg="lme4",
                        agqvec=agqvec)

########### cbpp data #############

cbpp$obs <- 1:nrow(cbpp) 
cbpp_glmmML <- fit_gen(data=cbpp,
                          formula=cbind(incidence,size-incidence)~period,
                          cluster="herd",family="binomial",
                          pkg="glmmML",
                          agqvec=agqvec_glmmML)


cbpp_lme4 <- fit_gen(data=cbpp, 
                      formula=cbind(incidence, size - incidence) ~ period +
                          (1 | herd),
                     family="binomial",
                     pkg="lme4",
                     agqvec=agqvec)


############ Contraception Data ##################

data("Contraception", package="mlmRev")
Contraception$ch <- factor(Contraception$livch != 0, labels = c("N","Y"))

contraception_glmmML <- fit_gen(data= Contraception,
                                formula= use~age+I(age^2)+urban+ch,
                                cluster="district", family= "binomial",
                                pkg="glmmML",
                                agqvec=agqvec_glmmML)

contraception_lme4 <- fit_gen(data=Contraception,
                               formula=use~age+I(age^2)+urban+ch+(1|district),
                               family= "binomial",
                               pkg="lme4",
                               agqvec=agqvec)

save(list=ls(pattern="_(lme4|glmmML)$"),
     file="AGQfits.RData")
