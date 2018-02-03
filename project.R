#############################################
# SPATIAL ANALYSIS OF EBOLA DATA REVISITED
# V.A. Suchar
# 04/28/2016
# Run through all the steps for CENTROID ROAD and SOCIAL weight matrices
############################################
library(GISTools)
library(rgdal)
library(spdep)
library(GWmodel)
#------------------------------------------------
# Get data
#-----------------------------------------------

# get incident rates
rates=read.csv("C:/Users/vasiles/ownCloud/Data Incubator Challenge/Diabetes INCIDENCE_2014.csv", 
               colClasses=c("NULL", NA, "NULL", "NULL", NA,"NULL","NULL","NULL","NULL","NULL", "NULL"),header=TRUE)
# get spatial data
#counties=readOGR(file.choose(), "cb_2016_us_county_20m") # in "C:\Users\vasiles\ownCloud\Data Incubator Challenge\us_counties
counties=readOGR(file.choose(), "UScounties") # in C:\Users\vasiles\ownCloud\Data Incubator Challenge\UScounties
proj4string(counties)

head(counties@data)
setdiff(counties@data$NAME,rates$County.1 )

dim(counties@data)

# add rates per 100 to the dataset
# reorder the rates by the order of counties in shp file
rates$FIPS=ifelse(rates$FIPS.Codes<10000, paste0(0, rates$FIPS.Codes), paste0(rates$FIPS.Codes))

counties@data=merge(counties@data, rates, by.x="FIPS", by.y="FIPS", all.x=TRUE)

#######################################################################
#######################################################################
# Make a map of the area
#######################################################################
#######################################################################
shade.rates <- auto.shading(c(min(counties@data$rate.per.1000),  max(counties@data$rate.per.1000)), cols = brewer.pal(5,"PRGn"))

choropleth(counties, counties@data$rate.per.1000, shading = shade.rates)
choro.legend(-110,70, shade.rates, title="2014 New cases/1000", cex=0.9)
title("Diabetes rates in 2014")

#######################################################################
# 1. GLOBAL MORAN'S I
#######################################################################
#create queen's case neighbors
counties.nb.queen=poly2nb(counties, queen=T)
summary(counties.nb.queen)

#turns nbobject into #weighted listwobject
counties.lw=nb2listw(counties.nb.queen, style="W", zero.policy=TRUE )


moran.I=moran.test(counties@data$rate.per.1000,listw=counties.lw,randomisation=T, zero.policy=TRUE, alternative="two.sided")
moran.I$estimate[1]
moran.I$p.value

#######################################################################
# 1. LOCAL MORAN'S I
#######################################################################

counties.li=localmoran(counties@data$rate.per.1000, counties.lw, zero.policy=TRUE, alternative="two.sided")

shade.q <- auto.shading(c(counties.li[,1], -counties.li[,1]), cols = brewer.pal(5,"PRGn"))

choropleth(counties, counties.li[,1], shading = shade.q)
choro.legend(-110,70, shade.q, title="Local Moran's I", cex=0.9)
title("Local Moran's I - Queen")


shade.pval <- shading(c(0.01, 0.05, 0.1), cols=rev(brewer.pal(4,"PuRd")))
#queen's map - holm
choropleth(counties, p.adjust(counties.li[,5], method='holm') ,shading = shade.pval)
choro.legend(-110,70, shade.pval, title="Local p-values", cex=0.8)
title(paste0("Local Moran's I p-values - Queen (Holm)", " ", "Week", " ",i-3))




