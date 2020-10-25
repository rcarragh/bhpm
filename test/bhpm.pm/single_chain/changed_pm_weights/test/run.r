RNGkind(sample.kind="Rounding")
library(bhpm)

set.seed(24002)
######################### All events, Severity 1, Model BB_dependent #####################
print("######################### All events, Severity 1, Model BB level 1 #####################")
data(demo.cluster.data)
pm = bhpm.pointmass.weights(demo.cluster.data)
pm = pm[pm$Outcome.Grp == "G00-G99" & pm$Outcome %in% c("G00-99_Outcome1", "G00-99_Outcome2", "G00-99_Outcome3"),]
pm$weight_pm = 0.6

set.seed(17902)
raw = bhpm.pm(demo.cluster.data, nchains = 1, pm.weights = pm)
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
