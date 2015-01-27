######### Generic data analysis  ###############

load("AGQfits.Rdata")
########### Graph 1 ############

## should restrict toenail_GLMM to AGQ<75;
##   AGQ=100 is *very* wonky and AGQ>75 is a bit wonky
graph1(toenail_glmmML)
graph1(toenail_lme4)

## check again: is glmmML doing adaptive GHQ or non-adaptive?

graph1(culcita_glmmML)
graph1(culcita_lme4)  

graph1(cbpp_glmmML)
graph1(cbpp_lme4)

graph1(contraception_glmmML)
graph1(contraception_lme4)

########### Graph 2 (Graph on same scale) ############
##  center and scale by 'correct' answer (max AGQ)

tmpf <- function(x,method) {
    data.frame(x,method=method)
}

mel_all_toenail <- rbind(tmpf(standardize_maxval(toenail_lme4,"REvar.ID"),"lme4"),
                         tmpf(standardize_maxval(toenail_glmmML[-(28:30),,],
                                     "REvar"),"glmmML"))
g_toenail <- ggplot(mel_all_toenail,
                    aes(x=AGQ,rmse_stdvalue,
                        color=varcat,linetype=method))+geom_line()  


##  center by 'correct' answer (max AGQ) and scale by max-AGQ SD

tmpf <- function(x,method) {
  data.frame(x,method=method)
}

mel_all_toenail_sd <- rbind(tmpf(standardize_maxsd(toenail_lme4,"REvar.ID"),"lme4"),
                         tmpf(standardize_maxsd(toenail_glmmML[-(28:30),,],
                                                 "REvar"),"glmmML"))

g_toenail_sd <- ggplot(mel_all_toenail_sd,
                    aes(x=AGQ,rmse_stdvalue,
                        color=varcat,linetype=method))+geom_line()  


g_toenail_sd+ggtitle("Standarized estimates center by 'correct' answer (max AGQ) and scale by max-AGQ SD")

### What should be the title of this graph?


##### center and scale by 'correct' answer (max AGQ)

mel_all_culcita <- rbind(tmpf(standardize_maxval(culcita_lme4,"REvar.block"),"lme4"),
                         tmpf(standardize_maxval(culcita_glmmML, "REvar"),"glmmML"))
g_culcita <- ggplot(mel_all_culcita,
                    aes(x=AGQ,rmse_stdvalue,color=varcat,linetype=method))+geom_line()  


mel_all_cbpp <- rbind(tmpf(standardize_maxval(cbpp_lme4,"REvar.herd"),"lme4"),
                         tmpf(standardize_maxval(cbpp_glmmML, "REvar"),"glmmML"))
g_cbpp =ggplot(mel_all_cbpp,aes(x=AGQ,rmse_stdvalue,color=varcat,linetype=method))+geom_line()  



mel_all_contraception <- rbind(tmpf(standardize_maxval(contraception_lme4,"REvar.district"),"lme4"),
                      tmpf(standardize_maxval(contraception_glmmML, "REvar"),"glmmML"))
g_contraception =ggplot(mel_all_contraception,aes(x=AGQ,rmse_stdvalue,color=varcat,linetype=method))+geom_line()  


g_toenail
g_culcita
g_cbpp
g_contraception

