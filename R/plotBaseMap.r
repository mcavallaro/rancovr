source("utils.R")
library(sf)
library(spatstat)
 
A=st_read("shape_file_new/Areas.shp")

plot_map<-function(X.range, Y.range, main=NULL){
  plot(st_geometry(A), xlim=X.range, ylim = Y.range, axes=T, col='#e2e2e2', border='white', lwd=1, main=main)
  box(lwd=0.1)
}


PostcodeAreasBoundaries=st_read("shape_file_new/Areas.shp")
# Keep only Englandpostcode areas:
PostCodeAreas=as.character(PostcodeAreasBoundaries[[1]])
PostCodeAreasC = read.csv("Data/Book1.csv")
PostCodeAreasC = PostCodeAreasC[PostCodeAreasC$Country=='England', ]
id = sapply(PostCodeAreas, function(x){x %in% PostCodeAreasC[,1]})
PostcodeAreasBoundaries = st_geometry(PostcodeAreasBoundaries)[id]
rm(id)
RegionBoundaries = st_read("shape_files/RegionBoundaries.shp")
RegionBoundaries = st_transform(RegionBoundaries, crs=st_crs(PostcodeAreasBoundaries))
UK = st_read("Other maps/Map_UK.shp")
UK = st_transform(UK, crs=st_crs(PostcodeAreasBoundaries))
Wales. = st_union(UK[UK$NAME_1=='Wales',])




plotBaseMap<-function(main=NULL, add=T, Wales=TRUE, onlyregion=F,...){
  if (onlyregion){
    plot(st_geometry(UK), col=NA, border='white', lwd=0.8,  add=add,...)
  }else{
    plot(st_geometry(UK), col=NA, border='#e2e2e288', lwd=0.8,  add=add,...)    
  }
  if (Wales){
    plot(st_geometry(Wales.), axes=F, border='black', lwd=1, add=T)    
  }
  plot(st_geometry(RegionBoundaries), axes=F, border='black', lwd=1,  main=main, add=T)
  # box(lwd=0.1)
}


X.range = c(-6,2)
Y.range = c(50,56)
#EnglandBoundaries=st_union(st_geometry(RegionBoundaries))

is_inside<-function(i,j,x,y){
  latitude=y[i]
  longitude=x[j]
  stp=st_point(c(longitude, latitude))
  return(!is.empty(st_contains(st_geometry(EnglandBoundaries), stp)[[1]]))
}
