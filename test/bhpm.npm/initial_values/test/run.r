RNGkind(sample.kind="Rounding")
library(bhpm)
set.seed(25729)
######################### All events, Severity 1, Model 1a_independent #####################
print("######################### All events, Severity 1, Model 1a_independent #####################")
data(demo.cluster.data)
gen.init <- bhpm.gen.initial.values(demo.cluster.data, level  = 2, hier = 3)
gen.init$gamma[1,]$value = -7

raw = bhpm.npm(demo.cluster.data, lev = 2, initial_values = gen.init)
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

print("Removing objects...")
rm(summ)
gc()
ptheta = bhpm.ptheta(raw)
print("Removing objects...")
rm(raw)
gc()
print("Writing the ptheta probabilities....")
write.table(ptheta, "ptheta.dat")
ptheta95 = ptheta[ptheta$ptheta > 0.95,] 
write.table(ptheta95, "ptheta95.dat")
print("Removing objects...")
#rm(conv)
rm(ptheta)
rm(ptheta95)
gc()
print("Finished.")

warnings()
