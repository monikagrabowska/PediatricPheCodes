plotManhattan <-
  function(d, suggestive.line=0.05, significant.line,
           annotate.level,
           y.axis.interval=5,
           ...) {
    if(sum(c("phenotype","p") %in% names(d))<2 ) stop("Data input must contain columns phenotype and p.")
    if(!length(d$OR)) stop("OR not provided in data input.")
    
    #Remove records with NA p values
    d=d[!is.na(d$p),]
    
    #Transform the significance thresholds
    if(missing(significant.line)) significant.line=suggestive.line/nrow(d)
    significant.line=-log10(significant.line)
    suggestive.line=-log10(suggestive.line)
    if(missing(annotate.level)) {
      if(is.na(significant.line)) {annotate.level=-log10(min(d$p,na.rm=T))}
      else{annotate.level=significant.line}
    }
    else {annotate.level=-log10(annotate.level)}
    
    #Nudge all of the p=0 results to the smallest double
    if(sum(d$p==0)>0)d[d$p==0,]$p=.Machine$double.xmin
    
    #Restrict to only those with appropriate p-values
    d=d[d$p>0 & d$p<=1,]
    
    #Create the - log p value for plotting. No errors given restriction to p>0
    d$value = -log10(d$p)
    
    max.y=max(ceiling(max(d$value)),4.4)
    
    y.axis.label=expression(-log[10](italic(p)))
    
    sizing=FALSE
    
    #Create OR direction
    d$direction = d$OR>=1
    
    plot=phenotypePlot(d,suggestive.line=suggestive.line,significant.line=significant.line,
                       annotate.level=annotate.level,
                       y.axis.interval=y.axis.interval,
                       ...)
    
  }
