# read in the raw nielsen store data

weekstub = switch(yearnum, Year1 = "1114_1165", Year2 = "1166_1217", Year3 = "1218_1269",
                           Year4 = "1270_1321", Year5 = "1322_1373", Year6 = "1374_1426", Year7 = "1427_1478")

temppath = paste("C:/Users/matthew/Dropbox/iri data code/","r data/",sep="")

storechar = read.fwf(paste(rawdatapath,yearnum,"/External/",product,"/Delivery_Stores",sep=""),c(8,3,9,25,5,5,8),header=FALSE,skip=1,stringsAsFactors=FALSE)
colnames(storechar) = read.fwf(paste(rawdatapath,yearnum,"/External/",product,"/Delivery_Stores",sep=""),c(8,3,9,25,5,5,8),header=FALSE,nrow=1,stringsAsFactors=FALSE)
colnames(storechar) = gsub(" ","",colnames(storechar))
storechar$OU = gsub(" ","",storechar$OU)
storechar$Market_Name = gsub(" ","",storechar$Market_Name)
storechar$MskdName = gsub(" ","",storechar$MskdName)

storechar = storechar[order(storechar$IRI_KEY),]

#observe multiple observations in case of merger - very few obs are dropped

key.l = c(0,storechar$IRI_KEY[1:(nrow(storechar)-1)])

storechar = storechar[key.l != storechar$IRI_KEY,]

save(storechar,file=paste(outputpath,"storechar_",yearnum,"_",product,".RData",sep=""))

drug = try(read.table(paste(rawdatapath,yearnum,"/External/",product,"/",product,"_drug_",weekstub,sep=""),header=TRUE,stringsAsFactors=FALSE))
if(class(drug)=="try-error") {
  drug = try(read.table(paste(rawdatapath,yearnum,"/External/",product,"/",product,"_drug_",weekstub,".raw",sep=""),header=TRUE,stringsAsFactors=FALSE))
}

if(product == "saltsnck") {

  # sometimes this explodes so we should shrink the data first, then manipulate it

  storechar = storechar[storechar$Market_Name == "EAU CLAIRE" | storechar$Market_Name == "PITTSFIELD",]

  idgood = unique(storechar$IRI_KEY)

  drug = drug[drug$IRI_KEY %in% idgood,]

  nline = 10000000

  getout = FALSE

  groc = NULL
  indx=0

  while(!getout) {

    groc1 = read.table(paste(rawdatapath,yearnum,"/External/",product,"/",product,"_groc_",weekstub,sep=""),header=TRUE,stringsAsFactors=FALSE,skip=indx,nrows=nline)
    getout = nrow(groc1) != nline

    if(indx > 0) {
      colnames(groc1) = c("IRI_KEY","WEEK","SY","GE","VEND","ITEM","UNITS","DOLLARS","F","D","PR")
    }

    groc1 = groc1[groc1$IRI_KEY %in% idgood,]

    groc = rbind(groc,groc1)

    indx = indx + nline

  }

} else {
  groc = try(read.table(paste(rawdatapath,yearnum,"/External/",product,"/",product,"_groc_",weekstub,sep=""),header=TRUE,stringsAsFactors=FALSE))
  if(class(groc)=="try-error") {
    groc = try(read.table(paste(rawdatapath,yearnum,"/External/",product,"/",product,"_groc_",weekstub,".raw",sep=""),header=TRUE,stringsAsFactors=FALSE))
  }

}

stdata = rbind(drug,groc)

stdata = merge(stdata,storechar,by="IRI_KEY")

if(product != "saltsnck") {
  stdata = stdata[stdata$Market_Name == "EAU CLAIRE" | stdata$Market_Name == "PITTSFIELD",]
}

save(stdata,file=paste(outputpath,"stores_",yearnum,"_",product,".RData",sep=""))
