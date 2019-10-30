#!/usr/bin/Rscript

# A simple script for qqplots using ggplot2
# contact: qcorbin@hsph.harvard.edu

qq <- function(pvals, facet = NULL, colour = NULL, group = NULL, nrow=NULL, thin_qt = 0.01, n_pts = 5000, ribbon = TRUE, confidence_level = 0.05, point.alpha = 1, point.size = 0.7, ribbon.alpha = 0.25, abline.colour = 'red', legend.title = NULL, theme.objects = NULL ){

	# REQUIRED ARGUMENTS
	# pvals = length(N) vector of p-values (or -log10(p-values))
	
	# GROUPING AND FACETING ARGUMENTS
	# facet = length(N) facetting variable 
	# colour = length(N) grouping and colouring variable 
	# nrow = number of rows in facet_grid; only used when `facet` != NULL
	
	# THINNING OPTIONS
	# thin_qt = thin points with p-values below thin_qt quantile (faster plotting)
	# n_pts = total number of points to show after thinning
	
	# PLOTTING OPTIONS
	# ribbon = show confidence band?
	# confidence_level = confidence level for confidence bands
	# point.alpha = alpha (transparency) of qq-plot points
	# point.size = size of qq-plot points
	# ribbon.alpha = alpha (transparency) of confidence bands
	# abline.colour = colour of null diagonal line
	# legend.title = title of legend (used only when `colour` != NULL)
	# theme.objects = extra ggplot theme objects, e.g., theme_minimal()
	

	require(data.table)
	require(ggplot2)

	`+.gg` <- `%+%`

	getDT <- function(p, group_label = '', conf_alpha = confidence_level ){
		n <- length(p)
		if( max(pvals) > 1 ){
			cat('\nAssuming input is -log10(p-value)\n')
			p <- 0.1^p
		}
		mlp <- (-1)*log10(sort(p))
		null_mlp <- (-1)*log10((1:n)/(n+1))
		null_min <- (-1)*log10( qbeta( conf_alpha/2, 1:n, n +1 - 1:n) )
		null_max <- (-1)*log10( qbeta( 1- conf_alpha/2, 1:n, n +1 - 1:n) )
		if( n*( 1 - thin_qt ) > 1.25*n_pts ){
			kp_0 <- ceiling(n*thin_qt)
			kp <- sort(c(
				1:(kp_0-1),
				sample( kp_0:n, n_pts, prob = 1 + mlp[kp_0:n] )
			))
		}else{
			kp <- 1:n
		}
		data.table('mlp'=mlp,'null_mlp'=null_mlp, 'null_min'=null_min, 'null_max'=null_max,'fc_group' = group_label)[kp,]
	}

	if( is.null(group) ) group <- 1

	if( !is.null(facet) ){
		dt <- data.table(pvals,facet,group)[,getDT(pvals,facet[1]),by=list(facet,group)]
		pl <- ggplot(dt, aes(y=mlp, ymin=null_min, ymax=null_max, x=null_mlp, colour=fc_group)) + facet_wrap(~fc_group,nrow=nrow)
	}else if( !is.null(colour)){
		dt <- data.table('pvals'=pvals,'colour'=colour,'group'=group)[,getDT(pvals,colour[1]),by=list(group,colour)]
		pl <- ggplot(dt, aes(y=mlp, ymin=null_min, ymax=null_max, x=null_mlp, colour=colour)) + guides(colour = guide_legend(title = legend.title))
	}else{
		dt <- data.table(pvals,group)[,getDT(pvals),by=list(group)]
		pl <- ggplot(dt, aes(y=mlp, x=null_mlp, ymin=null_min, ymax=null_max))
	}
	
	if( ribbon ) pl <- pl %+% geom_ribbon(colour=NA, alpha = ribbon.alpha)

	pl <- pl %+% geom_abline(slope = 1, intercept = 0, colour = abline.colour) %+% geom_point(size = point.size, alpha = point.alpha) %+% ylab(expression("Observed"~-log[10](p-value))) %+% xlab(expression("Expected"~-log[10](p-value)))
	
	if( !is.null(theme.objects) ) pl <- pl %+% theme.objects
	
	invisible(list('plot' = pl, 'data' = dt))
}
