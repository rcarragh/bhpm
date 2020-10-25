RNGkind(sample.kind="Rounding")
library(bhpm)


set.seed(24290)
######################### First events, Severity 1, Model BB_dependent #####################
print("######################### First events, Severity 3, Model BB_dep (level 1) #####################")

data(demo.cluster.data)

gen.init <- bhpm.gen.initial.values(demo.cluster.data, level  = 1, hier = 2, model = "BB")
gen.init$gamma[1,]$value = -7

set.seed(20492)
raw = bhpm.pm(demo.cluster.data, hier = 2, nchains = 1, initial_values = gen.init)
conv = bhpm.convergence.diag(raw)
sink("conv.dat")
bhpm.print.convergence.summary(conv)
sink()
rm(conv)
gc()
summ = bhpm.summary.stats(raw)

sink("summary.dat")
bhpm.print.summary.stats(summ)
sink()


rm(summ)
gc()
ptheta = bhpm.ptheta(raw)
rm(raw)
gc()
write.table(ptheta, "ptheta.dat")
ptheta90 = ptheta[ptheta$ptheta > 0.90,] 
write.table(ptheta90, "ptheta90.dat")
#rm(conv)
rm(ptheta)
rm(ptheta90)
gc()
print("Finished.")

warnings()
