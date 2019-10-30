#!/usr/bin/Rscript

# Simple scripts for precision-recall and ROC curves using ggplot2
# contact: qcorbin@hsph.harvard.edu

library(ggplot2)
library(data.table)

## this just simplifies string concatenation
`+.default` <- .Primitive('+')
`+.character` <- paste0
`+.gg` <- `%+%`
`+` <- function (a, b) UseMethod('+')

getCurves <- function(x, y, group = NULL, extend = FALSE){
	
	require(data.table)
	
	if( is.null(group) ){
		dd <- data.table('x'=x,'uy'=y)[order(-x)]
		dd[,y:=mean(uy),by=x]

		sy <- sum(y)
		n <- nrow(dd)
		py <- sy/n

		dd <- dd[order(-x)]
		dd$Recall <- cumsum(dd$y)/sy

		dd$RE_SE <- sqrt( ( dd$Recall * (1-dd$Recall) + 1e-8) / ( n*py*(1-py) ) )

		dd$FPR <- cumsum(1-dd$y)/(n - sy)
		dd$Precision <- cumsum(dd$y)/seq(n)

		dd$PR_SE <- sqrt( py*(1-py) /seq(n))

		dd$FDR <- cumsum(1-dd$y)/seq(n)
		dd
	}else{
		ddat <- data.table('x'=x,'y'=y, 'z' = group)
		do.call(rbind, lapply(
		unique(group), function(gg){
			out <- with(
			subset(ddat, z == gg),
			getCurves(x,y,NULL,TRUE)
			)
			out$group <- gg
			out
		})
		)
	}
}


plotCurve <- function(eqn, data, line.size = 0.2, strat = NULL, strat.method = "max", facet = NULL, grouping = NULL, labs = NULL, cband = FALSE, palette = "Accent", wd_roc = 250, wd_pr = 35 , legend.position = 'bottom' , legend.width = 0.5, legend.space = 8, legend.size = 10, roc.xlim = NULL, stack.plots = FALSE , color.vals = NULL ){

	require(ggplot2); require(data.table); require(gridExtra); require(ggpubr); require(grid)

	cb_alpha <- ifelse(cband, 0.1, 0)

	vars <- as.character(eqn, "variables")[-1]
	outcome <- strsplit(gsub(" ", "", vars[1]), "\\+")[[1]]
	scores <- strsplit(gsub(" ", "", vars[2]), "\\+")[[1]]

	data <- as.data.table(data)
	
	eqt <- function(x,y=rep(1,length(x))) if(all(is.na(x[y>0])) | min(c(1,var(y,na.rm=TRUE),var(x,na.rm=TRUE)))==0 ){rep(as.numeric(NA),length(x))}else{ecdf(x[y==0])(x)}
	
	
	if( grepl('max', strat.method)  ){
		eqt <- function(x, y){
			if( all(is.na(x[y>0])) | min(c(1,var(y,na.rm=TRUE),var(x,na.rm=TRUE)))==0 ){
				rep(as.numeric(NA),length(x))
			}else{
				x <- x - ifelse(grepl('max.adj', strat.method), min(x, na.rm = TRUE)*(1-1/length(x)), 0)
				x/max(x, na.rm = TRUE)
			}
		}
	}

	eqt.std <- eqt
	
	if( !is.null(strat) ){
		strat <- paste0('`', strat, '`')
		strat <- paste(strat, sep = '', collapse = ',')
		eval(parse(text = gsub('%strat',strat,gsub('%y',outcome, "data[,MEAN_Y:=mean(`%y`),by=list(%strat)]"))));		data <- subset(data, MEAN_Y > 0 & MEAN_Y < 1)
		for(var in scores){ 
			eval(parse(text = 
				gsub('%y',outcome,gsub('%var',var,gsub('%strat',strat,"data[,`%var`:=eqt(`%var`,`%y`),by=list(%strat)]")))
			))
		}
		data2 <- data
		for(var in scores){
			eval(parse(text =
				gsub('%y',outcome,gsub('%var',var,gsub('%strat',strat,"data2[,`%var`:=eqt.std(`%var`,`%y`),by=list(%strat)]")))
			))
		}
		
		eval(parse(text = 
			gsub('%y',outcome, gsub('%strat',strat,
				paste0("top_rk <- data[,list(", paste("`",scores,"`=mean(`%y`[`",scores,"`>=max(`",scores,"`[`%y`<1], na.rm = TRUE)], na.rm = TRUE)",sep='',collapse=','), "),by=list(%strat)][,list(",paste("`",scores,"`=mean(`",scores,"`, na.rm = TRUE)",sep='',collapse=','),"),]")

			))
		))
		top_rk <-  data.table::melt(top_rk, measure.vars=scores, value.name = "top-ranked", variable.name="Predictor")
		
	}else{
		top_rk <- NULL
	}
	
	if( !is.null(facet) ){
		grouping <- facet
	}

	if(is.null(grouping)){
		groups <- NULL
	}else{
		groups <- data[[grouping]]
	}

	out <- as.data.table(do.call(rbind, 
		lapply(scores, function(sc){
			df <- getCurves(y = data[[outcome]], x = data[[sc]], group = groups)
			df$Predictor <- sc
			df$PR <- mean(df$y)
			df
		})
	))
	
	integr <- function(x, y, n = length(x)){
			y <- y[order(x)]; x <- sort(x)
			sum((y[-n] + y[-1])*(x[-1] - x[-n]))/2
	}

	AUC <- out[,list(
		"PRC" = integr(x=c(0,Recall),y=c(1,Precision)),
		"ROC" = integr(x=c(0,FPR,1),y=c(0,Recall,1)),
		"pROC10" = integr(x=c(0,FPR[FPR<=0.10]),y=c(0,Recall[FPR<=0.10])),
		"pROC05" = integr(x=c(0,FPR[FPR<=0.05]),y=c(0,Recall[FPR<=0.05])),
		"pROC01" = integr(x=c(0,FPR[FPR<=0.01]),y=c(0,Recall[FPR<=0.01]))
	),by=Predictor]
	if( !is.null(top_rk) ) AUC <- merge(AUC, top_rk, by = "Predictor")

	NUP <- length(unique(out$Predictor))
	
	roc.xmax <- 1
	if( !is.null(roc.xlim) ) roc.xmax <- max(roc.xlim)
	
	if(!is.null(labs)){
		out$Predictor <- factor(out$Predictor, levels = names(labs), labels = c(labs))
	}
	
	out[,`:=`(
		Precision=Precision,
		Recall=Recall,
		FPR=FPR
	),]		

	color_obj <- scale_colour_brewer(guide="none", palette = palette)
	color_legend_obj <- scale_colour_brewer(palette = palette)
	
	if( !is.null( color.vals ) ){
		color_obj <- scale_colour_manual(values = color.vals, guide = 'none')
		color_legend_obj <- scale_colour_manual(values = color.vals)			
	}
	
	coord_obj <- coord_cartesian(ylim=0:1, xlim=0:1, expand = FALSE, clip = 'off' )
	
	if(is.null(grouping) | !is.null(facet) ){
	
		PRC_OBJ <- ggplot(out, aes(y = Precision, x = Recall, colour = Predictor, ymin = Precision -1.96*PR_SE, ymax = Precision + 1.96*PR_SE, group = Predictor)) %+% geom_ribbon(fill='black',alpha=cb_alpha,colour=NA) %+% geom_line(size=line.size) %+% geom_abline(intercept=out$PR[1],slope=0,linetype=2) %+% theme_minimal() %+% color_obj %+% coord_obj %+% guides(colour=FALSE)

		ROC_OBJ <- ggplot(out, aes(y = Recall, x = FPR, colour = Predictor, ymin = Recall - 1.96*RE_SE, ymax = Recall + 1.96*RE_SE, group = Predictor)) %+% geom_ribbon(fill='black',alpha=cb_alpha,colour=NA) %+% geom_line(size=line.size) %+% geom_abline(intercept=0,slope=1,linetype=2) %+% theme_minimal() %+% color_obj %+% coord_obj %+% xlab("1 - Specificity") %+% ylab("Sensitivity") %+% guides(colour=FALSE)

		if( legend.position == 'bottom' ){
			LEGEND <- get_legend(ggplot(out, aes(y = Precision, x = Recall, colour = Predictor)) %+% color_legend_obj %+% geom_line(size=1.5) %+% theme_minimal() %+% guides(colour=guide_legend(title=NULL, nrow = ifelse(NUP > 3, 1, 1), byrow = (NUP > 4) )) %+% theme(legend.position="bottom")  )  
		}else{ 
			LEGEND <- get_legend(ggplot(out, aes(y = Precision, x = Recall, colour = Predictor)) %+% color_legend_obj %+% geom_line(size=1.5) %+% theme_minimal() %+% guides(colour=guide_legend(title=NULL )) %+% theme( legend.key = element_rect(size = legend.space, color = NA), legend.key.size = unit(legend.size, 'lines') )  ) 
		}
		
		if(!is.null(facet)){
			PRC_OBJ <- PRC_OBJ + facet_wrap(~group,nrow=1,scales='free') + scale_y_continuous(limits=0:1) + scale_x_continuous(limits=0:1)
			ROC_OBJ <- ROC_OBJ + facet_wrap(~group,nrow=1,scales='free') + scale_y_continuous(limits=0:1) + scale_x_continuous(limits=0:1) 

		}
		if( legend.position == 'bottom' ){
			PLOT <- grid.arrange(arrangeGrob(ROC_OBJ, PRC_OBJ, ncol = 2 - stack.plots), as_ggplot(LEGEND), ncol = 1, heights = c(8,1), widths = 16) 
		}else{
			PLOT <- grid.arrange(arrangeGrob(ROC_OBJ, PRC_OBJ, ncol = 2 - stack.plots), as_ggplot(LEGEND), ncol = 2, heights = 8, widths = c(16,8*legend.width))
		}
	}else{
		
		PRC_OBJ <- ggplot(out, aes(y = Precision, x = Recall, colour = Predictor)) %+% geom_line(alpha=0.25,aes(group=group)) %+% stat_summary(fun.y=mean,aes(x=I(round(Recall*wd_pr)/wd_pr)), geom='line', colour = 'black', size=line.size) %+% geom_abline(aes(intercept=PR),slope=0,linetype=2) %+% theme_minimal() %+% color_obj + coord_obj %+% facet_wrap(~Predictor,ncol=1) %+% guides(colour=FALSE)

		ROC_OBJ <- ggplot(out, aes(y = Recall, x = FPR, colour = Predictor)) %+% geom_line(alpha=0.25,aes(group=group)) %+% stat_summary(fun.y=mean,aes(x=I(round(FPR*wd_roc)/wd_roc)), geom='line', colour = 'black') %+% geom_abline(intercept=0,slope=1,linetype=2) %+% theme_minimal() %+% color_obj + coord_obj %+% xlab("1 - Specificity") %+% facet_wrap(~Predictor,ncol=1) %+% ylab("Sensitivity") %+% guides(colour=FALSE)

		if( legend.position == 'bottom' ){
			LEGEND <- get_legend(ggplot(out, aes(y = Precision, x = Recall, colour = Predictor)) %+% color_legend_obj %+% geom_line(size=1.5) %+% theme_minimal() %+% guides(colour=guide_legend(title=NULL, nrow = ifelse(NUP > 3, 1, 1), byrow = (NUP > 4) )) %+% theme(legend.position="bottom")  )  
		}else{ 
			LEGEND <- get_legend(ggplot(out, aes(y = Precision, x = Recall, colour = Predictor)) %+% color_legend_obj %+% geom_line(size=1.5) %+% theme_minimal() %+% guides(colour=guide_legend(title=NULL )) %+% theme( legend.key = element_rect(size = legend.space, color = NA), legend.key.size = unit(legend.size, 'lines') )  ) 
		}


		if( legend.position == 'bottom' ){
			PLOT <- grid.arrange(arrangeGrob(ROC_OBJ, PRC_OBJ, ncol = 2-stack.plots), as_ggplot(LEGEND), ncol = 1, heights = c(8,1), widths = 16) }else{ PLOT <- grid.arrange(ROC_OBJ, PRC_OBJ, as_ggplot(LEGEND), ncol = 3, heights = 8, widths = c(8,8,8*legend.width)) 
		}

	}

	PLOT

	invisible(list(
		data = out,
		AUC = AUC,
		plot = as_ggplot(PLOT)
	))
}


rprint <- function(D, dig = 3){
	print(as.data.table(lapply(as.data.frame(D), function(x){
		if(is.numeric(x)){
			if(!all(round(x)==x, na.rm = TRUE)){
				round(x,dig+1)
			}else{
				as.integer(x)
			}
		}else{
			x
		}
	}
	)), digits=dig)
}


getAUC <- function(eqn, data, line.size = 0.5, strat = NULL, facet = NULL, grouping = NULL, labs = NULL, cband = FALSE, palette = "Accent", wd_roc = 250, wd_pr = 35 , legend.position = 'bottom' , legend.width = 0.5, legend.space = 8, legend.size = 10, strat.method = ""){

	require(ggplot2); require(data.table); require(gridExtra); require(ggpubr); require(grid)

	cb_alpha <- ifelse(cband, 0.1, 0)

	vars <- as.character(eqn, "variables")[-1]
	outcome <- strsplit(gsub(" ", "", vars[1]), "\\+")[[1]]
	scores <- strsplit(gsub(" ", "", vars[2]), "\\+")[[1]]

	data <- as.data.table(data)
	
	eqt <- function(x,y=rep(1,length(x))) if(all(is.na(x[y>0])) | min(c(1,var(y,na.rm=TRUE),var(x,na.rm=TRUE)))==0 ){rep(as.numeric(NA),length(x))}else{ecdf(x)(x)}

	if( strat.method == "max" ){
		eqt <- function(x, y){
		if( all(is.na(x[y>0])) | min(c(1,var(y,na.rm=TRUE),var(x,na.rm=TRUE)))==0 ){
		rep(as.numeric(NA),length(x))
		}else{
		x <- x - min(x, na.rm = TRUE)
		x/max(x, na.rm = TRUE)
		}
		}
	}

	if( !is.null(strat) ){
		strat <- paste0('`', strat, '`')
		strat <- paste(strat, sep = '', collapse = ',')
		
		eval(parse(text = gsub('%strat',strat,gsub('%y',outcome, "data[,MEAN_Y:=mean(`%y`),by=list(%strat)]"))));		
		data <- subset(data, MEAN_Y > 0 & MEAN_Y < 1)

		for(var in scores){ 
			eval(parse(text = 
				gsub('%y',outcome,gsub('%var',var,gsub('%strat',strat,"data[,`%var`:=eqt(`%var`,`%y`),by=list(%strat)]")))
			))
		}
		
		eval(parse(text = 
			gsub('%y',outcome, gsub('%strat',strat,
				paste0("top_rk <- data[,list(", paste("`",scores,"`=mean(`",scores,"`[`%y`>0]>=max(`",scores,"`[`%y`<1], na.rm = TRUE), na.rm = TRUE)",sep='',collapse=','), "),by=list(%strat)][,list(",paste("`",scores,"`=mean(`",scores,"`, na.rm = TRUE)",sep='',collapse=','),"),]")
			))
		))
		top_rk <-  data.table::melt(top_rk, measure.vars=scores, value.name = "top-ranked", variable.name="Predictor")
		
	}else{
		top_rk <- NULL
	}
	
	if( !is.null(facet) ){
		grouping <- facet
	}

	if(is.null(grouping)){
		groups <- NULL
	}else{
		groups <- data[[grouping]]
	}


	out <- as.data.table(do.call(
	rbind, lapply(scores, function(sc){
			df <- getCurves(y = data[[outcome]], x = data[[sc]], group = groups)
			df$Predictor <- sc
			df$PR <- mean(df$y)
			df
	})
	))

	integr <- function(x, y, n = length(x)){
			y <- y[order(x)]; x <- sort(x)
			sum((y[-n] + y[-1])*(x[-1] - x[-n]))/2
	}

	AUC <- out[,list(
	"PRC" = integr(x=c(0,Recall),y=c(1,Precision)),
	"ROC" = integr(x=c(0,FPR,1),y=c(0,Recall,1)),
	"ROC10" = integr(x=c(0,FPR[FPR<0.10]),y=c(0,Recall[FPR<0.10])),
	"ROC05" = integr(x=c(0,FPR[FPR<0.05]),y=c(0,Recall[FPR<0.05])),
	"ROC01" = integr(x=c(0,FPR[FPR<0.01]),y=c(0,Recall[FPR<0.01]))
	),by=Predictor]
	if( !is.null(top_rk) ) AUC <- merge(AUC, top_rk, by = "Predictor")
	AUC
}
