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
# [1] 1.000907

# 3. If required calculate summary statistics (mean/median/hpi)
summ <- bhpm.summary.stats(mod.npm)
bhpm.print.summary.stats(summ)

# These may be accessed directly for model parameters:
print(head(summ$theta))
print(summ$theta[1,]$mean)
# [1]  -0.08478039
hpi <- c(summ$theta[1,]$hpi_lower, summ$theta[1,]$hpi_upper)
print(hpi)
# [1] -0.4764439  0.3033665

# 4. Assuming the model has converged assess which AEs may be associated with treatment.
# The model paramter theta is used for this purpose.
theta.post.prob <- bhpm.ptheta(mod.npm)

# A large (posterior) probability that theta is > 0 is an indication that an adverse event is associated with treamtment.

print(theta.post.prob[ theta.post.prob$ptheta.pos > 0.95 | theta.post.prob$ptheta.neg > 0.95,])

#    Trt.Grp Cluster Outcome.Grp         Outcome     ptheta ptheta.pos ptheta.zero   ptheta.neg
# 5        2  M/0-64     G00-G99 G00-99_Outcome1 0.99966667 0.99966667           0 0.0003333333
# 6        2  M/0-64     G00-G99 G00-99_Outcome2 1.00000000 1.00000000           0 0.0000000000
# 7        2  M/0-64     G00-G99 G00-99_Outcome3 0.04766667 0.04766667           0 0.9523333333
# 14       2 M/65-84     G00-G99 G00-99_Outcome1 0.99986667 0.99986667           0 0.0001333333
# 23       2   M/85+     G00-G99 G00-99_Outcome1 1.00000000 1.00000000           0 0.0000000000
# 30       3  M/0-64       Bleed  Bleed_Outcome3 0.95205000 0.95205000           0 0.0479500000
# 53       3   M/85+     I00-I99 I00-99_Outcome1 0.01690000 0.01690000           0 0.9831000000
# 59       4  M/0-64     G00-G99 G00-99_Outcome1 0.02851667 0.02851667           0 0.9714833333
# 71       4 M/65-84     I00-I99 I00-99_Outcome1 0.04143333 0.04143333           0 0.9585666667
