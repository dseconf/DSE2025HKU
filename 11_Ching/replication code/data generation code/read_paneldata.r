# read in the household panel data, year by year.

weekstub = switch(yearnum, Year1 = "1114_1165", Year2 = "1166_1217", Year3 = "1218_1269",
                           Year4 = "1270_1321", Year5 = "1322_1373", Year6 = "1374_1426", Year7 = "1427_1478")

ystub = substr(yearnum,5,5)

# Read in the product characteristics files

# note that these files were re-saved as csv, so I can read them into R

if(yearnum == "Year7") {
  brandinfo = read.csv(paste(rawdatapath,"parsed stub files 2007/prod_",product,".csv",sep=""),stringsAsFactors=FALSE)
} else {
  brandinfo = read.csv(paste(rawdatapath,"parsed stub files/prod_",product,".csv",sep=""),stringsAsFactors=FALSE)
}

brandinfo = brandinfo[,!(colnames(brandinfo) %in% c("L1", "Level", "GE", "VEND", "ITEM"))]
brandinfo$COLUPC = gsub("-","",brandinfo$UPC)

brandinfo$UPC <- NULL

save(brandinfo,file=paste(outputpath,"brands_",yearnum,"_",product,".RData",sep=""))

manent = read.csv(paste(rawdatapath,"demos trips external/manual store entry external.csv",sep=""),stringsAsFactors=FALSE)

static = read.csv(paste(rawdatapath,"demos trips external/static 1_7.csv",sep=""),stringsAsFactors=FALSE)

static = static[static$year == 1 & static$make_static == "yes",]

static$year <- NULL
static$make_static <- NULL

drugpanel = read.table(paste(rawdatapath,yearnum,"/External/",product,"/",product,"_PANEL_DR_",weekstub,".dat",sep=""),header=TRUE,stringsAsFactors=FALSE)

mapanel = read.table(paste(rawdatapath,yearnum,"/External/",product,"/",product,"_PANEL_MA_",weekstub,".dat",sep=""),header=TRUE,stringsAsFactors=FALSE)

grocpanel = read.table(paste(rawdatapath,yearnum,"/External/",product,"/",product,"_PANEL_GR_",weekstub,".dat",sep=""),header=TRUE,stringsAsFactors=FALSE)

hhpanel = rbind(drugpanel,grocpanel)
hhpanel = rbind(hhpanel,mapanel)

rm(drugpanel,mapanel,grocpanel)

# households that don't make the static will be dropped here
nobs1 = nrow(hhpanel)
nid1 = length(unique(hhpanel$PANID))
hhpanel = merge(hhpanel,static)
nobs2 = nrow(hhpanel)
nid2 = length(unique(hhpanel$PANID))

lupc = nchar(hhpanel$COLUPC)

#Add leading zeros into upc for merge

hhpanel$COLUPC[lupc==11] = paste("000",hhpanel$COLUPC[lupc==11],sep="")
hhpanel$COLUPC[lupc==12] = paste("0",substr(hhpanel$COLUPC[lupc==12],1,1),"0",substr(hhpanel$COLUPC[lupc==12],2,12),sep="")
hhpanel$COLUPC[lupc==13] = paste(substr(hhpanel$COLUPC[lupc==13],1,2),"0",substr(hhpanel$COLUPC[lupc==13],3,13),sep="")

hhpanel = merge(hhpanel,brandinfo,by.x="COLUPC",by.y="COLUPC")

load(paste(outputpath,"storechar_",yearnum,"_",product,".RData",sep=""))

load(paste(outputpath,"stores_",yearnum,"_",product,".RData",sep=""))

hhpanel = merge(hhpanel,storechar,all.x=TRUE)

hhpanel = merge(hhpanel,manent,by.x="IRI_KEY",by.y="IRI_Key",all.x=TRUE)

demodata = read.csv(paste(rawdatapath,"demos trips external/ads demo",ystub,".csv",sep=""),stringsAsFactors=FALSE)

oldnlist = c("Panelist.ID",demonames.orig)
newnlist = c("PANID",demovars)

i=0
for(n in oldnlist) {
  i=i+1
  colnames(demodata)[colnames(demodata)==n] = newnlist[i]
}

demodata = demodata[,c("PANID",demovars)]

hhpanel = merge(hhpanel,demodata)

save(hhpanel,file=paste(outputpath,"hhpanel_",product,"_",yearnum,".RData",sep=""))





