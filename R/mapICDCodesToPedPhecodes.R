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