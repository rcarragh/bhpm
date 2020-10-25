library(bhpm)

set.seed(1)

data(demo.cluster.data)

print(head(demo.cluster.data))

# 1. Fit the model:

mod.npm <- bhpm.npm(demo.cluster.data, nchains = 2)

# 2. Assess convergence
conv <- bhpm.convergence.diag(mod.npm)

# Printing a convergence summary will indicate if there are any obvious issues
# Any reported statistics that are greater than about 1.1 may indicate an issue.
bhpm.print.convergence.summary(conv)

print(max(conv$theta.conv.diag$stat))
# [1]  1.001697

# 3. If required calculate summary statistics (mean/median/hpi)
summ <- bhpm.summary.stats(mod.npm)
bhpm.print.summary.stats(summ)

# These may be accessed directly for model parameters:
print(head(summ$theta))
print(summ$theta[1,]$mean)
# [1] -0.2348521
hpi <- c(summ$theta[1,]$hpi_lower, summ$theta[1,]$hpi_upper)
print(hpi)
# [1] -0.6828268  0.1924472

# 4. Assuming the model have converged assess which AEs may be associated with treatment.
# The model paramter theta is used for this purpose.
theta.post.prob <- bhpm.ptheta(mod.npm)

# A large (posterior) probability that theta is > 0 is an indication that an adverse event is associated with treamtment.

print(theta.post.prob[ theta.post.prob$ptheta.pos > 0.95 | theta.post.prob$ptheta.neg > 0.95,])

#   Trt.Grp  Cluster Outcome.Grp  Outcome     ptheta ptheta.pos ptheta.zero   ptheta.neg
# 3        2 Cluster1      Group2 Outcome3 0.99973333 0.99973333           0 2.666667e-04
# 4        2 Cluster1      Group2 Outcome4 0.99995000 0.99995000           0 5.000000e-05
# 5        2 Cluster1      Group2 Outcome5 0.04705000 0.04705000           0 9.529500e-01
# 12       2 Cluster2      Group2 Outcome3 0.99993333 0.99993333           0 6.666667e-05
# 21       2 Cluster3      Group2 Outcome3 1.00000000 1.00000000           0 0.000000e+00
# 35       3 Cluster1      Group3 Outcome8 0.95175000 0.95175000           0 4.825000e-02
# 46       3 Cluster3      Group1 Outcome1 0.01590000 0.01590000           0 9.841000e-01
# 57       4 Cluster1      Group2 Outcome3 0.02540000 0.02540000           0 9.746000e-01
# 64       4 Cluster2      Group1 Outcome1 0.04168333 0.04168333           0 9.583167e-01

