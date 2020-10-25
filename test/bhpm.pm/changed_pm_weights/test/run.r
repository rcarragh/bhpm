RNGkind(sample.kind="Rounding")
library(bhpm)

set.seed(17902)
######################### All events, Severity 1, Model BB_dependent #####################
print("######################### All events, Severity 1, Model BB level 1 #####################")
data(demo.cluster.data)
pm = bhpm.pointmass.weights(demo.cluster.data)
pm = pm[pm$Outcome.Grp == "Group2" & pm$Outcome %in% c("Outcome5", "Outcome6", "Outcome7"),]
pm$weight_pm = 0.6

set.seed(17902)
raw = bhpm.pm(demo.cluster.data, nchains = 3, pm.weights = pm)
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
print("Removing objects...")
rm(raw)
gc()
write.table(ptheta, "ptheta.dat")
ptheta90 = ptheta[ptheta$ptheta > 0.90,] 
write.table(ptheta90, "ptheta90.dat")
print("Removing objects...")
#rm(conv)
rm(ptheta)
rm(ptheta90)
gc()
print("Finished.")

warnings()
