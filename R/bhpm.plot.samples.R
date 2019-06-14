# bhpm.plot.posterior
# Model bhpm.BB
# R. Carragher
# Date: 29/06/2018

Id <- "$Id: bhpm.plot.samples.R,v 1.3 2019/04/18 15:01:21 clb13102 Exp clb13102 $"


bhpm.plot.samples <- function(samples, title) {
	if (is.null(samples)) {
		print("NULL sample");
		return(NULL)
	}

	chains = 0
	d = dim(samples)

	mcmc_obj <- list(NA)

	if (is.null(d)) {
		mcmc_obj[[1]] <- mcmc(samples)
	}
	else {

		l = length(d)

		if (l > 2) {
			print("Dimension length error");
			return(NULL)
		}

		chains = d[1]

		for (i in 1:chains) {
			mcmc_obj[[i]] <- mcmc(samples[i,])
		}
	}

	m = mcmc.list(mcmc_obj)
	plot(m, main = title)
}
