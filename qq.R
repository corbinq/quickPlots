#!/usr/bin/Rscript

# A simple script for qqplots using ggplot2
# contact: qcorbin@hsph.harvard.edu

qq <- function(pvals, n_pvals = rep(1, length(pvals)), facet = NULL, colour = NULL, group = NULL, nrow=NULL, thin_qt = 0.01, n_digits = 2, n_pts = NA, ribbon = TRUE, confidence_level = 0.05, point.alpha = 1, point.size = 0.7, ribbon.alpha = 0.25, abline.colour = 'red', legend.title = NULL, theme.objects = NULL, print_plot = TRUE ){

	# REQUIRED ARGUMENTS
	# pvals = length(N) vector of p-values (or -log10(p-values))
	
	# Optional:
	# n_pvals = length(N) vector giving number of observations for each unique p-value
	
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
	# print_plot = show plot when function is called?
	

	require(data.table)
	require(ggplot2)

	`+.gg` <- `%+%`
	
	# fractional rounding to reduce number of unique points
	round_frac <- function(x, n, const = 10^n) round(x*const)/const

	getDT <- function(p, N, group_label = '', conf_alpha = confidence_level ){
		require(data.table)
		n <- length(p)
		if( max(pvals) > 1 ){
			cat('\nAssuming input is -log10(p-value)\n')
			p <- 0.1^p
		}
		ORD <- order(p)
		N <- N[ORD]
		Seq_N <- cumsum(N)
		N_t <- sum(N) + 1
		mlp <- (-1)*log10(p[ORD])
		null_mlp <- (-1)*log10(Seq_N/N_t)
		null_min <- (-1)*log10( qbeta( conf_alpha/2, Seq_N, N_t - Seq_N) )
		null_max <- (-1)*log10( qbeta( 1- conf_alpha/2, Seq_N, N_t - Seq_N) )
		if( !is.na(n_pts) & is.na(n_digits) ){
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
		}else if( !is.na(n_digits) ){
			rd <- function(x) round_frac(x, n_digits)
			unique(data.table('mlp'=rd(mlp),'null_mlp'=rd(null_mlp), 'null_min'=rd(null_min), 'null_max'=rd(null_max),'fc_group' = group_label))
		}
	}

	if( is.null(group) ) group <- 1

	if( !is.null(facet) ){
		dt <- data.table(pvals,n_pvals,facet,group)[,getDT(pvals,n_pvals,facet[1]),by=list(facet,group)]
		pl <- ggplot(dt, aes(y=mlp, ymin=null_min, ymax=null_max, x=null_mlp, colour=fc_group)) + facet_wrap(~fc_group,nrow=nrow)
	}else if( !is.null(colour)){
		dt <- data.table('pvals'=pvals,'n_pvals'=n_pvals,'colour'=colour,'group'=group)[,getDT(pvals,n_pvals,colour),by=list(group)]
		pl <- ggplot(dt, aes(y=mlp, ymin=null_min, ymax=null_max, x=null_mlp, colour=fc_group)) + guides(colour = guide_legend(title = legend.title))
	}else{
		dt <- data.table(pvals,n_pvals,group)[,getDT(pvals,n_pvals),by=list(group)]
		pl <- ggplot(dt, aes(y=mlp, x=null_mlp, ymin=null_min, ymax=null_max))
	}
	
	x_lims <- c(0, max(c(dt$null_mlp), na.rm = TRUE)*1.02)
	y_lims <- c(0, max(c(dt$mlp, dt$null_mlp), na.rm = TRUE)*1.02)
	
	if( ribbon ) pl <- pl %+% geom_ribbon(colour=NA, alpha = ribbon.alpha)

	pl <- pl %+% geom_abline(slope = 1, intercept = 0, colour = abline.colour) %+% geom_point(size = point.size, alpha = point.alpha) %+% ylab(expression("Observed"~-log[10]*'('*italic(p)*'-value)')) %+% xlab(expression("Expected"~-log[10]*'('*italic(p)*'-value)')) %+% coord_cartesian(ylim = y_lims, xlim = x_lims, expand = FALSE)
	
	if( !is.null(theme.objects) ) pl <- pl %+% theme.objects
	
	if( print_plot ) print(pl)
	
	invisible(list('plot' = pl, 'data' = dt))
}



