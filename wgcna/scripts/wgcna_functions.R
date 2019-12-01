

swap_cols = function(oldCols, newCols, moduleColors, MEcols, dynamicCols){
  for (i in 1:nrow(swap)){
    oldCol=oldCols[i]
    newCol=newCols[i]
    oldME = paste('ME', oldCol, sep='')
    newME = paste('ME', newCol, sep='')
    
    #make old copies for reference
    oldMEcols = MEcols
    oldModuleColors = moduleColors
    

    #swap out any old for the new
    MEcols[oldMEcols==newME]<-oldME
    moduleColors[oldModuleColors==newCol]<-oldCol
    dynamicCols[dynamicCols==newCol]<-oldCol
    
    #swap out the targets
    MEcols[oldMEcols==oldME]<-newME
    moduleColors[oldModuleColors==oldCol]<-newCol
    dynamicCols[dynamicCols==oldCol]<-newCol
  }
  return(list(moduleColors, MEcols, dynamicCols))
}


ggVerboseScatterplot = function(x, y, sample = NULL, corFnc = "cor", corOptions = "use = 'p'", 
            main = "", xlab = NA, ylab = NA, cex = 1, cex.axis = 1.5, 
            cex.lab = 1.5, cex.main = 1.5, abline = FALSE, abline.color = 1, 
            abline.lty = 1, corLabel = corFnc, displayAsZero = 1e-05, 
            col = 1, bg = 0, pch = 1, lmFnc = lm, plotPriority = NULL, 
            ...) 
  {
    if (is.na(xlab)) 
      xlab = as.character(match.call(expand.dots = FALSE)$x)
    if (is.na(ylab)) 
      ylab = as.character(match.call(expand.dots = FALSE)$y)
    x = as.numeric(as.character(x))
    y = as.numeric(as.character(y))
    corExpr = parse(text = paste(corFnc, "(x, y ", prepComma(corOptions), 
                                 ")"))
    cor = signif(eval(corExpr), 2)
    if (abs(cor) < displayAsZero) 
      cor = 0
    corp = signif(corPvalueStudent(cor, sum(is.finite(x) & is.finite(y))), 
                  2)
    if (corp < 10^(-200)) 
      corp = "<1e-200"
    else corp = paste("=", corp, sep = "")
    if (!is.na(corLabel)) {
      mainX = paste(main, " ", corLabel, " = ", cor, sep = "")
    }
    else mainX = main
    if (length(col) < length(x)) 
      col = rep(col, ceiling(length(x)/length(col)))
    if (length(pch) < length(x)) 
      pch = rep(pch, ceiling(length(x)/length(pch)))
    if (length(cex) < length(x)) 
      cex = rep(cex, ceiling(length(x)/length(cex)))
    if (length(bg) < length(x)) 
      bg = rep(bg, ceiling(length(x)/length(bg)))
    if (is.null(plotPriority)) 
      plotPriority = rep(1, length(x))
    if (length(plotPriority) != length(x)) 
      stop("When given, length of 'plotPriority' must equal length of 'x'.")
    if (!is.null(sample)) {
      if (length(sample) == 1) {
        sample = sample(length(x), sample)
      }
      priority1 = plotPriority[sample]
      order1 = order(priority1, na.last = TRUE)
      
      
      # plot(x[sample][order1], y[sample][order1], main = mainX, 
      #      xlab = xlab, ylab = ylab, cex.axis = cex.axis, cex.lab = cex.lab, 
      #      cex.main = cex.main, col = col[sample][order1], bg = bg[sample][order1], 
      #      pch = pch[sample][order1], cex = cex[sample][order1], 
      #      ...)
      
      xvals=x[sample][order1]
      yvals=y[sample][order1]
      scat.df = data.frame(x=x,
                           y=y)
      plt=scat.df %>% 
        ggplot(aes(x=x,y=y)) +
        geom_point(color='black', fill=bg[order1], pch=21) +
        labs(x=xlab, y=ylab, subtitle=mainX)
      
      
    }
    else {
      order1 = order(plotPriority, na.last = TRUE)
      
      
      # plot(x[order1], y[order1], main = mainX, xlab = xlab, 
      #      ylab = ylab, cex.axis = cex.axis, cex.lab = cex.lab, 
      #      cex.main = cex.main, col = col[order1], bg = bg[order1], 
      #      pch = pch[order1], cex = cex[order1], ...)
      xvals=x[sample][order1]
      yvals=y[sample][order1]
      scat.df = data.frame(x=x,
                           y=y)
      plt=scat.df %>% 
        ggplot(aes(x=x,y=y)) +
        geom_point(color='black', fill=bg[order1], pch=21) +
        labs(x=xlab, y=ylab, subtitle=mainX)

    }
    if (abline) {
      lmFnc = match.fun(lmFnc)
      fit = lmFnc(y ~ x)
      abline(reg = fit, col = abline.color, lty = abline.lty)
    }
    invisible(sample)
    return(plt)
  }  
  