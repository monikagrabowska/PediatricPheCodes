restrictPedPhecodesBySex <- function(phenotypes,id.sex) {
  data=merge(phenotypes,id.sex,by=1,all.x=T)
  #Get the column of the sex
  g=dim(data)[2]
  #Get the restrictions found in the phenotypes data frame
  current_sex_restriction=PedPheCodes::sex.restriction[PedPheCodes::sex.restriction$phecode %in% names(phenotypes)[-1],]
  #Get male and female-only phenotypes
  male_only=current_sex_restriction[current_sex_restriction$male_only,"phecode"]
  female_only=current_sex_restriction[current_sex_restriction$female_only,"phecode"]
  #Set row column matches to NA where inds of a particular sex meet restricted phenotypes
  data[!is.na(data[,g])&data[,g]!="F",unlist(female_only)]=NA
  data[!is.na(data[,g])&data[,g]!="M",unlist(male_only)]=NA
  
  #Return everything, sans sex
  data[,-g]
}
