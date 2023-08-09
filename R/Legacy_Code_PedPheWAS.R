library(vroom)
library(tidyverse)
library(parallel)
options(ggrepel.max.overlaps = Inf)

# ideally would make these global variables?
vocabulary.map <- vroom::vroom("../Data files/mg_phecode_map_5.5.2023.csv",
                               .name = janitor::make_clean_names, delim = ",",
                               col_types = c(vocabulary_id = "c", code = "c", phecode = "c"))

rollup.map <- vroom::vroom("../Data files/mg_phecode_rollup_5.5.2023.csv",
                           .name = janitor::make_clean_names, delim = ",",
                           col_types = c(code = "c", phecode_unrolled = "c"))

exclusion.map <- vroom::vroom("../Data files/mg_phecode_exclude_5.5.2023.csv",
                               .name = janitor::make_clean_names, delim = ",",
                               col_types = c(code = "c", exclusion_criteria = "c"))

gender_restriction <- vroom::vroom("../Data files/mg_gender_restriction_5.5.2023.csv",
                              .name = janitor::make_clean_names, delim = ",",
                              col_types = c(phecode = "c", male_only = "logical", female_only = "logical"))

annotate.phenotype.description <- vroom::vroom("../Data files/mg_pheinfo_5.5.2023.csv",
                               .name = janitor::make_clean_names, delim = ",",
                               col_types = c(phecode = "c", description = "c", groupnum = "double", group = "c", color = "c"))

example_data_icd <- vroom::vroom("../Data files/person_icd.csv",
                             .name = janitor::make_clean_names, delim = ",",
                             col_types = c(id = "i", vocabulary_id = "c", code = "c", count = "i"))

example_data_gender <- vroom::vroom("../Data files/test_gender.csv",
                             .name = janitor::make_clean_names, delim = ",",
                             col_types = c(id = "i", sex = "c"))

example_data_genotypes <- vroom::vroom("../Data files/test_genotypes.csv",
                                    .name = janitor::make_clean_names, delim = ",",
                                    col_types = c(id = "i", rs_example = "i"))


aggregate.fun <- function(index) {
  if (is.numeric(index) == TRUE) {
    sum(index)
  }
  else {
    length(unique(index))
  }
}


restrictPedPhecodesBySex <- function(phenotypes,id.sex) {
  data=merge(phenotypes,id.sex,by=1,all.x=T)
  #Get the column of the sex
  g=dim(data)[2]
  #Get the restrictions found in the phenotypes data frame
  current_gender_restriction=gender_restriction[gender_restriction$phecode %in% names(phenotypes)[-1],]
  #Get male and female-only phenotypes
  male_only=current_gender_restriction[current_gender_restriction$male_only,"phecode"]
  female_only=current_gender_restriction[current_gender_restriction$female_only,"phecode"]
  #Set row column matches to NA where inds of a gender meet restricted phenotypes
  data[!is.na(data[,g])&data[,g]!="F",unlist(female_only)]=NA
  data[!is.na(data[,g])&data[,g]!="M",unlist(male_only)]=NA

  #Return everything, sans sex
  data[,-g]
}


createPedPhenotypes <-
  function(id.vocab.code.index,
           id.sex,
           min.code.count=2,
           full.population.ids=unique(id.vocab.code.index$id))
  {
    id.name=names(id.vocab.code.index)[1]

    #check to make sure numeric codes were not passed in
    if(!class(id.vocab.code.index[[3]]) %in% c("character","factor")) {stop("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
    names(id.vocab.code.index)=c("id","vocabulary_id","code","index")
    message("Mapping codes to phecodes...")
    phemapped=mapICDCodesToPedPhecodes(id.vocab.code.index) %>% transmute(id, code=phecode, index)

    message("Aggregating codes...")
    phecode=ungroup(summarize(group_by(phemapped,id,code),count=aggregate.fun(index)))
    phecode=phecode[phecode$count>0,]

    message("Mapping exclusions...")
    exclusions = merge(phecode %>% rename(exclusion_criteria=code), exclusion.map, by = "exclusion_criteria")
    exclusions = exclusions %>%  transmute(id, code, count=-1) %>% distinct()
    phecode=rbind(phecode,exclusions)

    #If there is request for a min code count, adjust counts to -1 if needed
    if(!is.na(min.code.count)&(max(!is.na(phecode$count)&phecode$count<min.code.count))) {
      phecode[!is.na(phecode$count)&phecode$count<min.code.count,]$count=-1
    }

    if(!is.na(min.code.count)) {
      message("Coalescing exclusions and min.code.count as applicable...")
      phecode=ungroup(summarize(group_by(phecode,id,code),count=max(count)))
    }

    message("Reshaping data...")
    phens=spread(phecode,code,count,fill=0)

    #Set exclusions to NA, preserving IDs just in case one is -1
    tmp_id=phens[,1]
    phens[phens==-1]=NA
    phens[,1]=tmp_id

    #Add in individuals present in input or the full population list, but without mapped phecodes
    missing_ids=setdiff(full.population.ids,unique(phens$id))
    if(length(missing_ids)>0) {
      empty_record=phens[1,-1]
      empty_record[]=0
      phens=rbind(phens,data.frame(id=missing_ids,empty_record,check.names=F))
    }

    #Change to logical if there is a min code count
    if(!is.na(min.code.count)) {phens[,-1]=phens[,-1]>0}

    phens=restrictPedPhecodesBySex(phens,id.sex)

    #Limit to full population ids
    phens = filter(phens, id %in% full.population.ids)

    #Rename the ID column to the input ID column name
    names(phens)[1]=id.name

    phens
  }


mapICDCodesToPedPhecodes <-
  function(input) {

    if(sum(names(input) %in% c("vocabulary_id","code"))!=2) {
      stop("Must supply a data frame with 'vocabulary_id' and 'code' columns")
    }
    if(!class(input[["code"]]) %in% c("character","factor")) {stop("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}

    #Perform the direct map
    withCallingHandlers(output <- merge(input,vocabulary.map,by=c("vocabulary_id","code")),
                        warning = function(w) { if (grepl("coercing into character vector", w$message)) {invokeRestart("muffleWarning")}})
    #Remove old columns
    output = output %>% select(-code,-vocabulary_id) %>% rename(code=phecode)

    #Make distinct
    output = distinct(output)

    #Perform the rollup
    withCallingHandlers(output <- merge(output ,rollup.map,by="code"),
                        warning = function(w) { if (grepl("coercing into character vector", w$message)) {invokeRestart("muffleWarning")}})
    output = output %>% select(-code) %>% rename(phecode=phecode_unrolled)

    #Make distinct (again)
    output = distinct(output)

    #Return the output
    output
  }


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


phenotypePlot <-
  function(d, suggestive.line, significant.line,
           size.x.labels=9, size.y.labels=9,
           point.size = 3,
           annotate.phenotype=T,
           annotate.angle=0, annotate.size=3, annotate.level,
           annotate.list,
           lc.labels=F,
           y.axis.interval=y.axis.interval) {

    d=merge(d,annotate.phenotype.description,by.x="phenotype",by.y="phecode")

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


phewas_ext <-
  function(phenotypes, genotypes, data, covariates=NA,cores=1, additive.genotypes=T,
           method="glm", strata=NA, factor.contrasts=contr.phewas,
           return.models=F, min.records=20, MASS.confint.level=NA) {
    #Require an input data frame
    if(missing(data)) {
      stop("A data frame must be supplied in 'data'")
    }
    #Rename outcomes and predictors parameters if used.
    if(missing(phenotypes)) {
      if(!missing(outcomes)) phenotypes=outcomes
      else stop("Either phenotypes or outcomes must be passed in.")
    }
    if(missing(genotypes)) {
      if(!missing(predictors)) genotypes=predictors
      else stop("Either genotypes or predictors must be passed in.")
    }
    #Convert covariates to a list if it is not one
    if(class(covariates)!="list") { covariates=list(covariates)}

    #Checks for each of the PheWAS methods
    if(method=="glm") {
      association_method=phe_as_ext
    } else if (method=="clogit") {
      association_method=phe_as_clogit
      #Check for a strata parameter
      if(is.na(strata)) { stop("clogit requires groups- please provide the column name in the 'strata' parameter.")}
    } else if (method=="lrt") {
      association_method=phe_as_lrt
      #Setup the genotypes as a single complete list if not done already
      if(class(genotypes)!="list") {genotypes=list(genotypes)}
    } else if (method == "logistf") {
      association_method=phe_as_logistf
    } else {
      stop("Method must be one of: 'glm', 'clogit', 'lrt', or 'logistf'.")
    }

    para=(cores>1)
    #Create the list of combinations to iterate over
    full_list=data.frame(t(expand.grid(phenotypes,genotypes,covariates,stringsAsFactors=F)),stringsAsFactors=F)

    #If parallel, run the parallel version.
    if(para) {
      #Check to make sure there is no existing phewas cluster.
      if(exists("phewas.cluster.handle")) {
        #If there is, kill it and remove it
        message("Old cluster detected (phewas.cluster.handle), removing...")
        try(stopCluster(phewas.cluster.handle), silent=T)
        rm(phewas.cluster.handle, envir=.GlobalEnv)
      }
      message("Starting cluster...")
      assign("phewas.cluster.handle", makeCluster(cores), envir = .GlobalEnv)
      message("Cluster created, finding associations...")
      clusterExport(phewas.cluster.handle,c("data"), envir=environment())
      clusterCall(phewas.cluster.handle,library,package="dplyr",character.only=T)
      #Loop across every phenotype- iterate in parallel
      result <-parLapplyLB(phewas.cluster.handle, full_list, association_method, additive.genotypes=additive.genotypes,
                           confint.level=MASS.confint.level, min.records=min.records,
                           return.models=return.models,factor.contrasts=factor.contrasts,strata=strata)
      #Once we have succeeded, stop the cluster and remove it.
      stopCluster(phewas.cluster.handle)
      rm(phewas.cluster.handle, envir=.GlobalEnv)
    } else {
      #Otherwise, just use lapply.
      message("Finding associations...")
      result=lapply(full_list,FUN=association_method, additive.genotypes=additive.genotypes,
                    confint.level=MASS.confint.level, my.data=data, min.records=min.records,
                    return.models=return.models,factor.contrasts=factor.contrasts,strata=strata)
    }

    if(return.models) {
      message("Collecting models...")
      models=lapply(result,function(x){attributes(x)$model})
      names(models)=sapply(models,function(x){paste0(as.character(terms(x))[c(2,1,3)],collapse=" ")})
    }

    message("Compiling results...")
    successful.phenotypes=na.omit(sapply(result,function(x){attributes(x)$successful.phenotype}))
    n.tests=length(successful.phenotypes)
    successful.phenotypes=unique(successful.phenotypes)
    successful.genotypes=unique(na.omit(sapply(result,function(x){attributes(x)$successful.genotype})))
    sig=bind_rows(result)

    #Report warning if any convergence errors
    if(max(grepl(pattern = "[Error: The model did not converge]", sig$note, fixed=TRUE))){
      warning("Not all models converged, check the notes column for details.")
    }

    message("Cleaning up...")

    if(return.models){sig=list(results=sig,models=models)}

    return(sig)
  }


phe_as_ext <-
  function(phe.gen, additive.genotypes=T,min.records=20,return.models=F,confint.level=NA, factor.contrasts=NA, my.data, ...) {
    if(!missing(my.data)) data=my.data
    #Retrieve the targets for this loop
    phe=phe.gen[[1]]
    gen=phe.gen[[2]]
    gens=gen
    cov=phe.gen[[3]]

    #Subset the data
    d=data %>% select(one_of(na.omit(unlist(c(phe,gen,cov)))))
    #Turn covariates into a string, if not NA
    if(!is.na(cov[1])) {covariates=paste(cov,collapse=",")}
    else {covariates=NA_character_} #Make sure it is a character NA for aggregation

    #Exclude the exclusions for the target phenotype
    d=d[!is.na(d[[phe]]),]
    n_no_snp=sapply(d %>% select(one_of(gen)),
                    FUN=function(x){sum(is.na(x))})
    #Exclude rows with missing data
    d=na.omit(d)
    n_total=nrow(d)

    n_cases=NA_integer_
    n_controls=NA_integer_
    allele_freq=NA_real_
    HWE_pval=NA_real_
    or=NA_real_
    se=NA_real_
    p=NA_real_
    beta=NA_real_
    type=NA_character_
    note=""
    model=NA
    formula.string=NA_character_
    expanded_formula=NA_character_
    gen_expansion=1:length(gen)

    #Drop columns with no variability
    drop.cols = names(d)[sapply(d, function(col) length(unique(col)))<=1]
    if(length(drop.cols>0)) {
      note=paste(note,"[Note: Column(s) dropped due to lack of variability: ",paste0(drop.cols,collapse=", "),"]")
      d=select(d, -one_of(drop.cols))
      #Remove dropped columns from covs- sticks around in the listed "covariates"
      cov=setdiff(cov,drop.cols)
    }

    if(n_total<min.records) {
      note=paste(note,"[Error: <",min.records," complete records]")
    } else if(sum(c(phe,gen) %in% names(d))!=length(c(phe,gen))) {
      note=paste(note,"[Error: non-varying phenotype or genotype]")
    } else {
      if(additive.genotypes) {
        snp.details=lapply(d %>% select(one_of(gen)),
                           FUN=function(x){
                             if(class(x) %in% c("numeric","integer")){
                               a.f=sum(x)/(2*n_total)
                               if(sum(!(na.omit(x) %in% 0:2))==0) {
                                 P=a.f
                                 Q=1-a.f
                                 AA=sum(x==2)
                                 xAA=P^2*n_total
                                 Aa=sum(x==1)
                                 xAa=2*P*Q*n_total
                                 aa=sum(x==0)
                                 xaa=Q^2*n_total
                                 hwe=pchisq((AA-xAA)^2/(xAA)+(Aa-xAa)^2/(xAa)+(aa-xaa)^2/(xaa),1)
                               } else {hwe=NA_real_}
                             } else {
                               a.f=NA_real_
                               hwe=NA_real_
                             }
                             c(a.f,hwe)
                           })


        allele_freq=sapply(snp.details,FUN=`[`,1)
        HWE_pval=sapply(snp.details,FUN=`[`,2)

        #Report a warning as needed.
        if(min(sapply(d %>% select(one_of(gen)),
                      FUN=function(x){class(x) %in% c("numeric","integer") & sum(!(na.omit(x) %in% 0:2))==0}))!=TRUE){
          note=paste(note,"[Warning: At least one genotype was not coded 0,1,2, but additive.genotypes was TRUE.]")}
      }
      #Alter factors to use special contrasts
      if(suppressWarnings(!is.na(factor.contrasts))) {
        d=data.frame(lapply(d,function(x){
          if("factor" %in% class(x)){
            x=droplevels(x)
            contrasts(x)=factor.contrasts(x)
          }
          x}),check.names=F)
      }
      #Create the formula:
      formula.string=paste0("`",phe,"` ~ `",paste(na.omit(c(gen,cov)),collapse = "` + `"),'`')
      my.formula = as.formula(formula.string)

      #Check if phenotype is logical (boolean)
      if(class(d[[phe]]) %in% c("logical")) {
        type = "logistic"
        #Create the logistic model
        n_cases=sum(d[[phe]])
        n_controls=n_total-n_cases
        if(n_cases<min.records|n_controls<min.records) {note=paste(note,"[Error: <",min.records," cases or controls]")}
        else {
          model = glm(my.formula, data=d, family=binomial)
          modsum= summary(model)
          #If the models did not converge, report NA values instead.
          if(model$converged) {
            #Find the observed genotype columns
            gen_expansion=attr(model.matrix(my.formula, data=d),"assign")
            gen_list=which(gen_expansion %in% 1:length(gen))
            gen_expansion=gen_expansion[gen_list]

            #Find the rows with results that gets merged across all loops
            gens=row.names(modsum$coef)[gen_list]
            or=exp(modsum$coef[gen_list,1])
            beta=modsum$coef[gen_list,1]
            se=modsum$coef[gen_list,2]
            p=modsum$coef[gen_list,4]
            expanded_formula=paste0(names(model$coefficients),collapse=" + ")
          } else {
            note=paste(note,"[Error: The model did not converge]")
          }
        }
      } else {
        type = "linear"
        if(n_total<min.records) {
          note=paste(note,"[Error: <",min.records," records with phenotype and genotype]")
        } else {
          model = glm(my.formula, data=d)

          modsum= summary(model)
          #If the models did not converge, report NA values instead.
          if(model$converged) {
            #Find the observed genotype columns
            gen_expansion=attr(model.matrix(my.formula, data=d),"assign")
            gen_list=which(gen_expansion %in% 1:length(gen))
            gen_expansion=gen_expansion[gen_list]

            #Find the rows with results that gets merged across all loops
            gen_list=grep(gen,row.names(modsum$coef))
            gens=row.names(modsum$coef)[gen_list]
            beta=modsum$coef[gen_list,1]
            se=modsum$coef[gen_list,2]
            p=modsum$coef[gen_list,4]
            expanded_formula=paste0(names(model$coefficients),collapse=" + ")
          } else {
            note=paste(note,"[Error: The model did not converge]")
          }
        }
      }
    }

    output=data.frame(phenotype=phe,snp=gens,
                      covariates=covariates,
                      beta=beta, SE=se,
                      OR=or,
                      p=p, type=type,
                      n_total=n_total, n_cases=n_cases, n_controls=n_controls,
                      HWE_p=HWE_pval[gen_expansion],allele_freq=allele_freq[gen_expansion],n_no_snp=n_no_snp[gen_expansion],
                      formula=formula.string,
                      expanded_formula = expanded_formula,
                      note=note, stringsAsFactors=F)

    #Add confidence intervals if requested.
    if(!is.na(confint.level)) {
      if(!is.na(model)[1]){
        suppressMessages(conf<-confint(model,c(1,gen_list),level=confint.level))
        lower=conf[-1,1]
        upper=conf[-1,2]
        if(type=="logistic") {
          lower=exp(lower)
          upper=exp(upper)
        }
      } else {
        lower=NA_real_
        upper=NA_real_
      }
      output$lower=lower
      output$upper=upper

      output=output[,c("phenotype","snp","beta","SE",
                       "lower","upper","OR","p","type",
                       "n_total","n_cases","n_controls",
                       "HWE_p","allele_freq","n_no_snp","formula","expanded_formula","note")]
    }

    #If the complete models were requested, add them as well.
    if(return.models) {attributes(output)$model=model}
    attributes(output)$successful.phenotype=ifelse(is.na(p),NA,phe)
    attributes(output)$successful.genotype=ifelse(is.na(p),NA,gen)
    #Return this to the loop to be merged.
    output
  }


phe_as_logistf <-
  function(phe.gen, additive.genotypes=T,min.records=20,return.models=F,confint.level=NA, factor.contrasts=NA, my.data, ...) {
    if(!missing(my.data)) data=my.data
    #Retrieve the targets for this loop
    phe=phe.gen[[1]]
    gen=phe.gen[[2]]
    gens=gen
    cov=phe.gen[[3]]

    #Subset the data
    d=data %>% select(one_of(na.omit(unlist(c(phe,gen,cov)))))
    #Turn covariates into a string, if not NA
    if(!is.na(cov[1])) {covariates=paste(cov,collapse=",")}
    else {covariates=NA_character_} #Make sure it is a character NA for aggregation

    #Set up confidence intervals
    return.confint=!is.na(confint.level)
    if(!return.confint) {confint.level=0.05}

    #Exclude the exclusions for the target phenotype
    d=d[!is.na(d[[phe]]),]
    n_no_snp=sapply(d %>% select(one_of(gen)),
                    FUN=function(x){sum(is.na(x))})
    #Exclude rows with missing data
    d=na.omit(d)
    n_total=nrow(d)
    n_cases=NA_integer_
    n_controls=NA_integer_
    allele_freq=NA_real_
    HWE_pval=NA_real_
    or=NA_real_
    se=NA_real_
    p=NA_real_
    beta=NA_real_
    type=NA_character_
    note=""
    model=NA
    formula.string=NA_character_
    expanded_formula=NA_character_
    gen_expansion=1:length(gen)

    #Drop columns with no variability
    drop.cols = names(d)[sapply(d, function(col) length(unique(col)))<=1]
    if(length(drop.cols>0)) {
      note=paste(note,"[Note: Column(s) dropped due to lack of variability: ",paste0(drop.cols,collapse=", "),"]")
      d=select(d, -one_of(drop.cols))
      #Remove dropped columns from covs- sticks around in the listed "covariates"
      cov=setdiff(cov,drop.cols)
    }

    if(n_total<min.records) {
      note=paste(note,"[Error: <",min.records," complete records]")
    } else if(sum(c(phe,gen) %in% names(d))!=length(c(phe,gen))) {
      note=paste(note,"[Error: non-varying phenotype or genotype]")
    } else {
      if(additive.genotypes) {
        snp.details=lapply(d %>% select(one_of(gen)),
                           FUN=function(x){
                             if(class(x) %in% c("numeric","integer")){
                               a.f=sum(x)/(2*n_total)
                               if(sum(!(na.omit(x) %in% 0:2))==0) {
                                 P=a.f
                                 Q=1-a.f
                                 AA=sum(x==2)
                                 xAA=P^2*n_total
                                 Aa=sum(x==1)
                                 xAa=2*P*Q*n_total
                                 aa=sum(x==0)
                                 xaa=Q^2*n_total
                                 hwe=pchisq((AA-xAA)^2/(xAA)+(Aa-xAa)^2/(xAa)+(aa-xaa)^2/(xaa),1)
                               } else {hwe=NA_real_}
                             } else {
                               a.f=NA_real_
                               hwe=NA_real_
                             }
                             c(a.f,hwe)
                           })


        allele_freq=sapply(snp.details,FUN=`[`,1)
        HWE_pval=sapply(snp.details,FUN=`[`,2)

        #Report a warning as needed.
        if(min(sapply(d %>% select(one_of(gen)),
                      FUN=function(x){class(x) %in% c("numeric","integer") & sum(!(na.omit(x) %in% 0:2))==0}))!=TRUE){
          note=paste(note,"[Warning: At least one genotype was not coded 0,1,2, but additive.genotypes was TRUE.]")}
      }
      #Alter factors to use special contrasts
      if(suppressWarnings(!is.na(factor.contrasts))) {
        d=data.frame(lapply(d,function(x){
          if("factor" %in% class(x)){
            x=droplevels(x)
            contrasts(x)=factor.contrasts(x)
          }
          x}),check.names=F)
      }
      #Create the formula:
      formula.string=paste0("`",phe,"` ~ `",paste(na.omit(c(gen,cov)),collapse = "` + `"),'`')
      my.formula = as.formula(formula.string)

      #Check if phenotype is logical (boolean) or 0/1
      if(class(d[[phe]]) %in% c("logical") | sum(d[[phe]] %in% c(0,1))==nrow(d)) {
        type = "firth logistic"
        #Create the logistic model
        n_cases=sum(d[[phe]])
        n_controls=n_total-n_cases
        if(n_cases<min.records|n_controls<min.records) {note=paste(note,"[Error: <",min.records," cases or controls]")}
        else {
          model = logistf::logistf(my.formula, data=d, dataout=F, alpha=confint.level)
          #If the models did not converge, report NA values instead.
          #if(model$converged) {
          #Find the observed genotype columns
          gen_expansion=attr(model.matrix(my.formula, data=d),"assign")
          gen_list=which(gen_expansion %in% 1:length(gen))
          gen_expansion=gen_expansion[gen_list]

          #Find the rows with results that gets merged across all loops
          gens=names(model$coef)[gen_list]
          or=exp(model$coef[gen_list])
          beta=model$coef[gen_list]
          se=(diag(model$var)^0.5)[gen_list]
          p=model$prob[gen_list]
          expanded_formula=paste0(names(model$coefficients),collapse=" + ")

          #Add confidence intervals if requested.
          if(return.confint) {
            lower=exp(model$ci.lower[gen_list])
            upper=exp(model$ci.upper[gen_list])
          } else {
            lower=NA_real_
            upper=NA_real_
          }

          #} else {
          #   note=paste(note,"[Error: The model did not converge]")
          # }
        }
      } else {
        paste(note,"[Error:  Firth's bias-reduced penalized-likelihood logistic regression requires logical or 1/0 input.]")
        type = "firth error"
      }
    }

    output=data.frame(phenotype=phe,snp=gens,
                      beta=beta, SE=se,
                      OR=or,
                      p=p, type=type,
                      n_total=n_total, n_cases=n_cases, n_controls=n_controls,
                      HWE_p=HWE_pval[gen_expansion],allele_freq=allele_freq[gen_expansion],n_no_snp=n_no_snp[gen_expansion],
                      formula=formula.string,
                      expanded_formula = expanded_formula,
                      note=note, stringsAsFactors=F)

    #Add confidence intervals if requested.
    if(return.confint) {
      output$lower=lower
      output$upper=upper

      output=output[,c("phenotype","snp","beta","SE",
                       "lower","upper","OR","p","type",
                       "n_total","n_cases","n_controls",
                       "HWE_p","allele_freq","n_no_snp","formula","expanded_formula","note")]
    }

    #If the complete models were requested, add them as well.
    if(return.models) {attributes(output)$model=model}
    attributes(output)$successful.phenotype=ifelse(is.na(p),NA,phe)
    attributes(output)$successful.genotype=ifelse(is.na(p),NA,gen)
    #Return this to the loop to be merged.
    output
  }


phe_as_lrt <-
  function(phe.gen, min.records=20,return.models=T, my.data, ...) {
    if(!missing(my.data)) data=my.data

    #Retrieve the targets for this loop
    phenotype=phe.gen[[1]]
    gen=phe.gen[[2]]
    gens=paste0(gen,collapse = ", ")
    cov=phe.gen[[3]]
    #Turn covariates into a string, if not NA
    if(!is.na(cov[1])) {covariates=paste(cov,collapse=",")}
    else {covariates=NA_character_} #Make sure it is a character NA for aggregation

    #Subset the data
    d=data[,na.omit(unlist(c(phenotype,gen,cov)))]
    d=na.omit(d)
    #Exclude rows with missing data
    n_total=nrow(d)
    n_cases=NA_integer_
    n_controls=NA_integer_
    type=NA_character_
    note=""
    lrt=NA
    p=NA

    #Drop columns with no variability
    drop.cols = names(d)[sapply(d, function(col) length(unique(col)))<=1]
    if(length(drop.cols>0)) {
      note=paste(note,"[Note: Column(s) dropped due to lack of variability: ",paste0(drop.cols,collapse=", "),"]")
      d=select(d, -one_of(drop.cols))
      #Remove dropped columns from covs- sticks around in the listed "covariates"
      cov=setdiff(cov,drop.cols)
    }

    if(n_total<min.records) {
      note=paste(note,"[Error: <",min.records," complete records]")
    } else if(!(phenotype %in% names(d))) {
      note=paste(note,"[Error: non-varying phenotype]")
    } else {

      #Check if phenotype is logical (boolean)
      if(class(d[,phenotype]) %in% c("logical")) {
        type = "logistic - LRT"
        #Create the logistic model
        n_cases=sum(d[,phenotype])
        n_controls=n_total-n_cases
        if(n_cases<min.records|n_controls<min.records) {note=paste(note,"[Error: <",min.records," cases or controls]")}
        else {
          model = glm(as.formula(paste("`",phenotype,"`"," ~ .", sep="", collapse="")), data=d, family=binomial)
          model.base = glm(as.formula(paste("`",phenotype,"`"," ~ .", sep="", collapse="")), data=d %>% select(-one_of(gen)), family=binomial)
          lrt=lrtest(model, model.base)
          p=lrt$Pr[2]
        }
      } else {
        type = "gaussian"
        if(n_total<min.records) {
          note=paste(note,"[Error: <",min.records," records with phenotype and genotype]")
        } else {
          model = lm(as.formula(paste(phenotype," ~ .", sep="", collapse="")), data=d)
          model.base = lm(as.formula(paste(phenotype," ~ .", sep="", collapse="")), data=d %>% select(-one_of(gen)))
          lrt=lrtest(model, model.base)
          p=lrt$Pr[2]
        }
      }
    }
    output=data.frame(phenotype=phenotype,genotype=gens,covariates=covariates,p=p,type=type,
                      n_total=n_total, n_cases=n_cases, n_controls=n_controls,
                      note=note, stringsAsFactors=F)



    #If the complete models were requested, add them as well.
    if(return.models) {attributes(output)$lrt=lrt}
    attributes(output)$successful.phenotype=ifelse(is.na(lrt),NA,phenotype)
    #Return this to the loop to be merged.
    output
  }


phe_as_clogit <-
  function(phe.gen, additive.genotypes=T,min.records=20,return.models=F,confint.level=NA, factor.contrasts=NA, strata, my.data, ...) {
    if(!missing(my.data)) data=my.data
    #Retrieve the targets for this loop
    phe=phe.gen[[1]]
    gen=phe.gen[[2]]
    gens=gen
    cov=phe.gen[[3]]
    #Turn covariates into a string, if not NA
    if(!is.na(cov[1])) {covariates=paste(cov,collapse=",")}
    else {covariates=NA_character_} #Make sure it is a character NA for aggregation

    #Subset the data
    d=data %>% select(one_of(na.omit(unlist(c(phe,gen,cov,strata)))))

    #Exclude the exclusions for the target phenotype
    d=d[!is.na(d[[phe]]),]
    n_no_snp=sapply(d %>% select(one_of(gen)),
                    FUN=function(x){sum(is.na(x))})
    #Exclude rows with missing data
    d=na.omit(d)
    #Exclude strata with only one predictor class (or value) represented
    strata.keep=d %>% group_by_at(.vars=strata) %>% summarize_at(.funs="n_distinct",na.rm=TRUE,.vars=gen) %>%
      filter_at(.vars=gen,.vars_predicate=all_vars(.>1)) %>% select(one_of(strata))
    d = inner_join(d,strata.keep,by=strata)

    n_total=nrow(d)
    n_cases=NA_integer_
    n_controls=NA_integer_
    n_strata=NA_integer_
    allele_freq=NA_real_
    HWE_pval=NA_real_
    or=NA_real_
    se=NA_real_
    p=NA_real_
    beta=NA_real_
    type=NA_character_
    note=""
    model=NA
    formula.string=NA_character_
    expanded_formula=NA_character_
    gen_expansion=1:length(gen)

    #Drop columns with no variability
    drop.cols = names(d)[sapply(d, function(col) length(unique(col)))<=1]
    if(length(drop.cols>0)) {
      note=paste(note,"[Note: Column(s) dropped due to lack of variability: ",paste0(drop.cols,collapse=", "),"]")
      d=select(d, -one_of(drop.cols))
      #Remove dropped columns from covs- sticks around in the listed "covariates"
      cov=setdiff(cov,drop.cols)
    }

    #Turn covariates into a string, if not NA
    if(!is.na(cov[1])) {covariates=paste(cov,collapse=",")}
    else {covariates=NA_character_} #Make sure it is a character NA for aggregation

    if(n_total<min.records) {
      note=paste(note,"[Error: <",min.records," complete records]")
    } else if(sum(c(phe,gen) %in% names(d))!=length(c(phe,gen))) {
      note=paste(note,"[Error: non-varying phenotype or genotype]")
    } else {
      if(additive.genotypes) {
        snp.details=lapply(d %>% select(one_of(gen)),
                           FUN=function(x){
                             if(class(x) %in% c("numeric","integer")){
                               a.f=sum(x)/(2*n_total)
                               if(sum(!(na.omit(x) %in% 0:2))==0) {
                                 P=a.f
                                 Q=1-a.f
                                 AA=sum(x==2)
                                 xAA=P^2*n_total
                                 Aa=sum(x==1)
                                 xAa=2*P*Q*n_total
                                 aa=sum(x==0)
                                 xaa=Q^2*n_total
                                 hwe=pchisq((AA-xAA)^2/(xAA)+(Aa-xAa)^2/(xAa)+(aa-xaa)^2/(xaa),1)
                               } else {hwe=NA_real_}
                             } else {
                               a.f=NA_real_
                               hwe=NA_real_
                             }
                             c(a.f,hwe)
                           })


        allele_freq=sapply(snp.details,FUN=`[`,1)
        HWE_pval=sapply(snp.details,FUN=`[`,2)

        #Report a warning as needed.
        if(min(sapply(d %>% select(one_of(gen)),
                      FUN=function(x){class(x) %in% c("numeric","integer") & sum(!(na.omit(x) %in% 0:2))==0}))!=TRUE){
          note=paste(note,"[Warning: At least one genotype was not coded 0,1,2, but additive.genotypes was TRUE.]")}
      }
      #Alter factors to use special contrasts
      if(suppressWarnings(!is.na(factor.contrasts))) {
        d=data.frame(lapply(d,function(x){
          if("factor" %in% class(x)){
            x=droplevels(x)
            contrasts(x)=factor.contrasts(x)
          }
          x}),check.names=F)
      }
      #Create the formula:
      formula.string=paste0("`",phe,"` ~ `",paste(c(gen,cov),collapse = "` + `"),'`'," + strata(`",strata,"`)")
      my.formula = as.formula(formula.string)

      #Check if phenotype is logical (boolean)
      if(class(d[[phe]]) %in% c("logical")) {
        type = "conditional logistic"
        #Create the logistic model
        n_cases=sum(d[[phe]])
        n_controls=n_total-n_cases
        n_strata=n_distinct(d[[strata]])
        if(n_cases<min.records|n_controls<min.records) {note=paste(note,"[Error: <",min.records," cases or controls]")}
        else {

          model = tryCatch(clogit(my.formula, data=d), warning = function(w) {w$message}, error = function(e) {e$message})
          #If the models did not converge, report NA values instead.
          if(class(model)[1]!="character") {
            #Find the observed genotype columns
            gen_expansion=attr(model.matrix(my.formula, data=d),"assign")
            gen_list=which(gen_expansion %in% 1:length(gen))
            gen_expansion=gen_expansion[gen_list]
            gen_list=gen_list-1
            #Create the model
            modsum= summary(model)
            #Find the rows with results that gets merged across all loops
            gens=row.names(modsum$coefficients)[gen_list]
            or=modsum$coefficients[gen_list,2]
            beta=modsum$coefficients[gen_list,1]
            se=modsum$coefficients[gen_list,3]
            p=modsum$coefficients[gen_list,5]
            expanded_formula=paste0(c(names(model$coefficients),paste0("strata(`","strata","`)")),collapse=" + ")
          } else {
            note=paste0(note,"[Error: Potential fitting problem (try changing covariates). clogit error: ",model,"]")
            model=NA
          }
        }
      } else {
        type = "linear"
        note=paste(note,"[Error: clogit requires a logical/Boolean outcome]")
      }
    }

    output=data.frame(phenotype=phe,snp=gens,
                      covariates=covariates,
                      beta=beta, SE=se,
                      OR=or,
                      p=p, type=type,
                      n_total=n_total, n_cases=n_cases, n_controls=n_controls,n_strata=n_strata,
                      HWE_p=HWE_pval[gen_expansion],allele_freq=allele_freq[gen_expansion],n_no_snp=n_no_snp[gen_expansion],
                      formula=formula.string,
                      expanded_formula = expanded_formula,
                      note=note, stringsAsFactors=F)

    #Add confidence intervals if requested.
    if(!is.na(confint.level)) {
      if(!is.na(model)[1]){
        modsum= summary(model,conf.int=confint.level)
        lower=modsum$conf.int[gen_list,3]
        upper=modsum$conf.int[gen_list,4]
      } else {
        lower=NA_real_
        upper=NA_real_
      }
      output$lower=lower
      output$upper=upper

      output=output[,c("phenotype","snp", "covariates","beta","SE",
                       "lower","upper","OR","p","type",
                       "n_total","n_cases","n_controls","n_strata",
                       "HWE_p","allele_freq","n_no_snp","formula","expanded_formula","note")]
    }

    #If the complete models were requested, add them as well.
    if(return.models) {attributes(output)$model=model}
    attributes(output)$successful.phenotype=ifelse(is.na(p),NA,phe)
    attributes(output)$successful.genotype=ifelse(is.na(p),NA,gen)
    #Return this to the loop to be merged.
    output
  }


contr.phewas=function(x){
  y=contrasts(x)
  colnames(y)=paste0('-',colnames(y))
  y
}


# simple example
ped_phenotypes <- createPedPhenotypes(example_data_icd, example_data_gender)

df_full <- inner_join(example_data_genotypes, ped_phenotypes, by = "id")

results <- phewas_ext(phenotypes=names(ped_phenotypes)[-1], genotypes=colnames(example_data_genotypes)[2], df_full, cores=4)

test_plot <- plotManhattan(results, annotate.phenotype = T, annotate.level = 0.005, y.axis.interval = 1)

test_plot
