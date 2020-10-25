library(bhpm)

set.seed(1)

data(demo.cluster.data)

print(head(demo.cluster.data))

# 1. Fit the model:

mod.pm <- bhpm.pm(demo.cluster.data, nchains = 2)

# 2. Assess convergence
conv <- bhpm.convergence.diag(mod.pm)

# Printing a convergence summary will indicate if there are any obvious issues
# Any reported statistics that are greater than about 1.1 may indicate an issue.
bhpm.print.convergence.summary(conv)

print(max(conv$theta.conv.diag$stat))
# [1] 1.093147

# 3. If required calculate summary statistics (mean/median/hpi)
summ <- bhpm.summary.stats(mod.pm)
bhpm.print.summary.stats(summ)

# These may be accessed directly for model parameters:
print(head(summ$theta))
print(summ$theta[1,]$mean)
# [1] -0.01232866
print(summ$theta[1,]$median)
# [1] 0

hpi <- c(summ$theta[1,]$hpi_lower, summ$theta[1,]$hpi_upper)
print(hpi)
# [1] -0.1629756169  0.0001113545

# 4. Assuming the model have converged assess which AEs may be associated with treatment.
# The model paramter theta is used for this purpose.
theta.post.prob <- bhpm.ptheta(mod.pm)

# A large (posterior) probability that theta is > 0 is an indication that an adverse event is associated with treamtment.

print(theta.post.prob[ theta.post.prob$ptheta > 0.80,])
print(theta.post.prob[ theta.post.prob$ptheta.pos > 0.80 | theta.post.prob$ptheta.neg > 0.80,])

   Trt.Grp  Cluster Outcome.Grp  Outcome    ptheta ptheta.pos ptheta.zero ptheta.neg
# 3        2 Cluster1      Group2 Outcome3 0.9998000  0.9998000   0.0002000          0
# 4        2 Cluster1      Group2 Outcome4 0.9993125  0.9993125   0.0006875          0
# 12       2 Cluster2      Group2 Outcome3 1.0000000  1.0000000   0.0000000          0
# 21       2 Cluster3      Group2 Outcome3 0.9999750  0.9999750   0.0000250          0
