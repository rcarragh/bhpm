RNGkind(sample.kind="Rounding")
library(bhpm)

set.seed(18313)
######################### All events, Severity 1, Model BB_dependent #####################
print("######################### All events, Severity 1, Model BB level 1 #####################")
data(demo.cluster.data)

s.p = bhpm.sim.control.params(demo.cluster.data, "BB")
s.p1 = s.p[s.p$type == "SLICE" & s.p$Outcome.Grp == "I00-I99" & s.p$Outcome %in% c("I00-99_Outcome1", "I00-99_Outcome2"), ]
s.p1$value = 2
s.p1$control = 7
s.p2 = s.p[s.p$type == "MH" & s.p$variable == "theta" & s.p$Outcome.Grp == "G00-G99" & s.p$Outcome %in% c("G00-99_Outcome1", "G00-99_Outcome2", "G00-99_Outcome3"), ]
s.p2$value = 0.22
s.p = rbind(s.p1, s.p2)


set.seed(13518)
raw = bhpm.pm(demo.cluster.data, level = 2, nchains = 1, sim.params = s.p)
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
