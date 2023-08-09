phenotypePlot <-
  function(d, suggestive.line, significant.line,
           size.x.labels=9, size.y.labels=9,
           point.size = 3,
           annotate.phenotype=T,
           annotate.angle=0, annotate.size=3, annotate.level,
           annotate.list,
           lc.labels=F,
           y.axis.interval=y.axis.interval) {
    
    d=merge(d,PedPheCodes::annotate.phenotype.description,by.x="phenotype",by.y="phecode")
    
    d=d[!is.na(d$groupnum),]
    
    d$size=3
    
    #Remove lc.labels flag if no annotations
    if(!annotate.phenotype) lc.labels=F
    
    #Sort by the phenotype
    d=d[order(d$phenotype),]
    
    #Set the maximum x value to fit all phenotypes
    max.x = length(unique(d$phenotype))
    
    #Create the list of phenotypes, finding the best values for each phenotype
    phenotypes=aggregate(value ~ phenotype + groupnum, d,FUN=max)
    #Remove the least significant phenotypes; only has an effect if max.x was specified.
    phenotypes=phenotypes[order(phenotypes$value, decreasing=T),][1:min(nrow(phenotypes),max.x),]
    
    phenotypes=phenotypes[order(phenotypes$groupnum,phenotypes$phenotype),]
    
    phenotypes$seq = 1:nrow(phenotypes)
    
    #Limit to phenotype and seq, as they are the only relevant columns
    #Include value as min.value for annotation purposes
    phenotypes=phenotypes[,c("phenotype","seq","value")]
    names(phenotypes)[3]="min.value"
    
    #Add sequence information
    d=inner_join(phenotypes,d,by="phenotype")
    d=d[order(d$seq),]
    
    #Define the max y axis value if not provided
    max.y=ceiling(max(d$value))
    
    labels= summarize(group_by(d, groupnum), tick=mean(unique(seq)),label=as.character(group[1]))
    labels=labels[order(labels$tick),]
    
    color.palette = unique(d[order(d$seq),]$color)
    names(color.palette)=color.palette
    
    y.axis.label=expression(-log[10](italic(p)))
    x.axis.label="Phenotypes"
    
    #Generate the inital plot
    plot=ggplot(d,ylab=y.axis.label,xlab=x.axis.label) + theme(axis.title = element_text(face="bold",size=15))
    
    #Set the Y scale and labels
    plot=plot+scale_y_continuous(y.axis.label, limits=c(0,max.y), breaks=seq(0,max.y,y.axis.interval), expand=c(0,0))
    
    #Include lines for significance thresholds
    if (suggestive.line<=max.y && !missing(suggestive.line) && !is.na(suggestive.line)) {
      plot=plot+geom_hline(yintercept=suggestive.line,colour="blue", alpha=I(1/3),size=1)
    }
    if (significant.line<=max.y && !missing(significant.line) && !is.na(significant.line)) {
      plot=plot+geom_hline(yintercept=significant.line,colour="red",alpha=I(1/3),size=1)
    }
    
    plot=plot+aes(seq,value,size=size,colour=color)
    plot=plot+scale_size(range=c(point.size,point.size),guide="none")
    #Add points
    plot=plot+geom_point()
    
    #Color as defined
    plot = plot + scale_colour_manual(values= color.palette, guide="none")
    
    #Reduce label crowding
    labels$label[labels$label == "injuries & poisonings"] <- "injuries & \n poisonings"
    labels$label[labels$label == "immunologic & inflammatory disorders"] <- "immunologic & inflammatory"
    
    #Label the X axis with the groups
    plot=plot+scale_x_continuous(name=x.axis.label, limits=c(1,max.x), breaks=labels$tick, labels=labels$label, expand=c(0,0))
    
    #Set the default theme
    plot=plot+theme(
      panel.background=element_blank(),
      panel.grid.minor=element_blank(),
      axis.text.x=element_text(size=size.x.labels, colour="black", angle=-40, hjust=0, vjust=1),
      axis.text.y=element_text(size=size.y.labels, colour="black"),
      axis.line =element_line(colour="black"),
      axis.ticks=element_line(colour="black")
    )
    
    #Hide the legend by default
    plot = plot+theme(legend.position = "none")
    
    #Add OR information
    plot=plot+aes(shape = factor(direction), fill=color) +
      scale_shape_manual(values=c(25,24)) +
      scale_fill_manual(values= color.palette, guide="none")
    
    #If annotation present, start the definitions, otherwise skip it
    if (annotate.phenotype) 	{
      d$annotate=F
      #If provided with a list of phenotypes to annotate, select those.
      if(!missing(annotate.list)) d[d$phenotype %in% annotate.list, ]$annotate=T
      #Include those above the given threshold
      if((!missing(annotate.level) && !is.na(annotate.level)) && (sum(d$value>=annotate.level)>0)) d[d$value>=annotate.level, ]$annotate=T
      
      #Cap annotation length
      d$description = substr(d$description,1,100)
      #Add leading space
      d$description = paste0("  ",d$description)
      
      #If lower case labels are requested, lower case them.
      if(lc.labels) d$description=tolower(d$description)
      
      if(sum(d$annotate)==0) {
        message("No points met the annotation criteria.")
      }
      else {
        plot = plot + ggrepel::geom_text_repel(aes(label=description),colour="black",data=d[d$annotate,],size=annotate.size,angle=annotate.angle)
      }
    }
    plot
  }

