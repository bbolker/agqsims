library("agqsims")
set.seed(101)
dd <- simfun_culc(N.binom=1)
f0 <- fit_gen(data=dd,formula= y~ttt+(1|block),
        family="binomial",
        maxAGQ=3)
f1 <- fit_gen(data=dd,formula= y~ttt,
        family="binomial",
        cluster="block",
        pkg="glmmPQL",
        maxAGQ=0)
f2 <- fit_gen(data=dd,formula= y~ttt,
        family="binomial",
        cluster="block",
        pkg="glmmML",
        maxAGQ=3)
