# bhpm.BB.ptheta
# Model bhpm.BB
# R. Carragher
# Date: 29/06/2018

Id <- "$Id: bhpm.ptheta.R,v 1.5 2019/04/28 13:51:47 clb13102 Exp clb13102 $"

bhpm.ptheta <- function(raw)
{
	if (is.null(raw)) {
		print("NULL raw data");
		return(NULL)
	}

	model = attr(raw, "model")

	if (is.null(model)) {
		print("Missing model attribute");
		return(NULL)
	}

	summ = bhpm.cluster.ptheta(raw)

	return(summ)
}
