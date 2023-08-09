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
    exclusions = merge(phecode %>% rename(exclusion_criteria=code), PedPheCodes::exclusion.map, by = "exclusion_criteria")
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
