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
# [1] 1.065125

# 3. If required calculate summary statistics (mean/median/hpi)
summ <- bhpm.summary.stats(mod.pm)
bhpm.print.summary.stats(summ)

# These may be accessed directly for model parameters:
print(head(summ$theta))
print(summ$theta[15,]$mean)
# [1] -0.04455222
print(summ$theta[1,]$median)
# [1] 0

hpi <- c(summ$theta[15,]$hpi_lower, summ$theta[15,]$hpi_upper)
print(hpi)
# [1] -0.4166467  0.0000000

# 4. Assuming the model has converged assess which AEs may be associated with treatment.
# The model paramter theta is used for this purpose.
theta.post.prob <- bhpm.ptheta(mod.pm)

# A large (posterior) probability that theta is > 0 is an indication that an adverse event is associated with treamtment.

print(theta.post.prob[ theta.post.prob$ptheta.pos > 0.80 | theta.post.prob$ptheta.neg > 0.80,])
#    Trt.Grp Cluster Outcome.Grp         Outcome    ptheta ptheta.pos ptheta.zero ptheta.neg
# 5        2  M/0-64     G00-G99 G00-99_Outcome1 0.9995125  0.9995125   0.0004375    5.0e-05
# 6        2  M/0-64     G00-G99 G00-99_Outcome2 0.9987875  0.9987875   0.0011875    2.5e-05
# 14       2 M/65-84     G00-G99 G00-99_Outcome1 0.9998875  0.9998875   0.0001125    0.0e+00
# 23       2   M/85+     G00-G99 G00-99_Outcome1 1.0000000  1.0000000   0.0000000    0.0e+00

