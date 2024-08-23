#*******************************************************************************
# Import and process point cloud data
# Sam Flake, 18 January 2023
#
# This script takes the raw point clouds from the Matter and Form 3D scanner,
# and calculates some variables like the length, width, height (curl), and the 
# volume of the convex hull enclosing the points, for use in analysis of 
# flammability of leaves.
#*******************************************************************************
library("geometry")
library("reshape2")
library("rgl")

# Import and process point clouds ----------------------------------------------

cloud_list <- dir("./LeafScans/xyz/")

hull_df <- data.frame(filename = cloud_list,
                      Species = character(length(cloud_list)),
                      IndCode = character(length(cloud_list)),
                      IndLeafCode = character(length(cloud_list)),
                      Dry = character(length(cloud_list)),
                      length = numeric(length(cloud_list)),
                      width = numeric(length(cloud_list)),
                      height = numeric(length(cloud_list)),
                      hull_volume = numeric(length(cloud_list)))
#751-754
system.time( #this takes 12 minutes or so
  for(i in 1:length(cloud_list)){
  cloud <- read.csv(paste0("./LeafScans/xyz/", cloud_list[i]), sep = " ", header = FALSE)
  pca <- princomp(cloud[, 1:3])
  hull_df[i, "length"] <- abs(min(pca$scores[, 1]) - max(pca$scores[, 1]))
  hull_df[i, "width"] <- abs(min(pca$scores[, 2]) - max(pca$scores[, 2]))
  hull_df[i, "height"] <- abs(min(pca$scores[, 3]) - max(pca$scores[, 3]))
  
  hull_df[i, "hull_volume"] <- convhulln(pca$scores, options = c("n", "FA"))$vol
  
  }
)

#to plot individual leaf
rgl::plot3d(cloud[, 1:3], angle = 45, type = "s", size = 0.8, lit = TRUE,
            xlim = c(-40, 40), ylim = c(0, 80), zlim = c(-60, 20))
plot3d(rgl::wire3d(rgl::mesh3d(hull$hull)))

#clean up some names
hull_df$Species <- toupper(substr(hull_df$filename, 1, 4))
hull_df[grep("casyf", hull_df$filename, ignore.case = TRUE), "Species"] <- "CASYf" 
hull_df[grep("casys", hull_df$filename, ignore.case = TRUE), "Species"] <- "CASYs"

hull_df$IndCode <- toupper(substr(hull_df$filename, 1, 5))
hull_df[grep("casyf", hull_df$filename, ignore.case = TRUE), "IndCode"] <- 
  paste0("CASYf", substr(hull_df[grep("casyf", hull_df$filename, ignore.case = TRUE), "filename"], 6, 6))  
hull_df[grep("casys", hull_df$filename, ignore.case = TRUE), "IndCode"] <- 
  paste0("CASYs", substr(hull_df[grep("casys", hull_df$filename, ignore.case = TRUE), "filename"], 6, 6))

hull_df$IndLeafCode <- toupper(substr(hull_df$filename, 1, 7))
hull_df[grep("casyf", hull_df$filename, ignore.case = TRUE), "IndLeafCode"] <- 
  paste0("CASYf", substr(hull_df[grep("casyf", hull_df$filename, ignore.case = TRUE), "filename"], 6, 8)) 
hull_df[grep("casys", hull_df$filename, ignore.case = TRUE), "IndLeafCode"] <- 
  paste0("CASYs", substr(hull_df[grep("casys", hull_df$filename, ignore.case = TRUE), "filename"], 6, 8)) 

hull_df$Dry <- "No"
hull_df[grep("dry", hull_df$filename), "Dry"] <- "Yes"

write.csv(hull_df, "./clean data/cloud data processed.csv")

#Restructure data --------------------------------------------------------------

hull_dry <- hull_df[hull_df$Dry == "Yes", ]
hull_fresh <- hull_df[hull_df$Dry == "No", ]

ind_leaves <- unique(c(hull_dry$IndLeafCode, hull_fresh$IndLeafCode))

hull_leaves <- data.frame(Species = character(length(ind_leaves)),
                          IndCode = character(length(ind_leaves)),
                          IndLeafCode = ind_leaves,
                          freshL3D = numeric(length(ind_leaves)),
                          freshW3D = numeric(length(ind_leaves)),
                          freshH3D = numeric(length(ind_leaves)),
                          freshHull = numeric(length(ind_leaves)),
                          dryL3D = numeric(length(ind_leaves)),
                          dryW3D = numeric(length(ind_leaves)),
                          dryH3D = numeric(length(ind_leaves)),
                          dryHull = numeric(length(ind_leaves)),
                          stringsAsFactors = FALSE)

for(i in 1:length(ind_leaves)){
  fresh_exist <- ind_leaves[i] %in% hull_fresh$IndLeafCode
  dry_exist <- ind_leaves[i] %in% hull_dry$IndLeafCode
  
  hull_leaves$Species[i] <- hull_df[hull_df$IndLeafCode == ind_leaves[i], "Species"][1]
  hull_leaves$IndCode[i] <- hull_df[hull_df$IndLeafCode == ind_leaves[i], "IndCode"][1]
  if(fresh_exist){
    hull_leaves$freshL3D[i] <- hull_fresh[hull_fresh$IndLeafCode == ind_leaves[i], "length"]
    hull_leaves$freshW3D[i] <- hull_fresh[hull_fresh$IndLeafCode == ind_leaves[i], "width"]
    hull_leaves$freshH3D[i] <- hull_fresh[hull_fresh$IndLeafCode == ind_leaves[i], "height"]
    hull_leaves$freshHull[i] <- hull_fresh[hull_fresh$IndLeafCode == ind_leaves[i], "hull_volume"]
    }
  if(dry_exist){
    hull_leaves$dryL3D[i] <- hull_dry[hull_dry$IndLeafCode == ind_leaves[i], "length"]
    hull_leaves$dryW3D[i] <- hull_dry[hull_dry$IndLeafCode == ind_leaves[i], "width"]
    hull_leaves$dryH3D[i] <- hull_dry[hull_dry$IndLeafCode == ind_leaves[i], "height"]
    hull_leaves$dryHull[i] <- hull_dry[hull_dry$IndLeafCode == ind_leaves[i], "hull_volume"]
  }
  
}

hull_leaves[, c(4:11)] <- apply(hull_leaves[, c(4:11)], c(1,2), 
                                FUN = function(x){ifelse(x == 0, NA, x)})

write.csv(hull_leaves, "./clean data/cloud_data_per_leaf.csv")
