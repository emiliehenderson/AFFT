### These function allow for creating overlapping and non-overlapping raster tiles.

### create_tiles_ext creates SpatVect objects of overlapping and nonoverlapping tile grids (returned as a list).
# This is handy for creating a grid that can be saved and reused. The overlapping grid can be used for creating
# overlapping tiles for processes that require neighborhood statistics. The non-overlapping grid can be used to trim
# overlapping tiles after processing. This is particularly handy when merging overlapping pixels does not make sense and
# trimming them off is a more valid approach. 

### create_raster_tiles generates tiles from in input grid. Such as one of the outputs from create_tiles_ext..

### trim_overlap trims overlapping area off of overlapping raster tiles. 

### commit_tile_index inserts raster tile info into our raster source table. 

### tile_index_OTF is similar to commit_tile_index but creates the index "On-The-Fly". No committing to a table. 

### check_overlay simply checks that tiled output overlays with a reference raster.

### Note: When these workflow were created the terra package had no option for created overlapping raster tiles. 
### A 'buffer' arg has since been added to terra::makeTiles() and could likely replace parts of this. Although, 
### for our purposes it's still convenient to have the overlapping and non-overlapping vector grids. Which will be used 
### in various parts of our data processing and modeling workflows. 


# Function to build extents (vector grids) for tiles
# Returns tile extents, both overlapping and nonoverlapping, and returns new overall raster ext. 
# This function was taken and altered from here:
# https://stackoverflow.com/questions/71709795/how-to-iterate-crop-raster-in-r-using-a-moving-window
create_tiles_ext <- function (x, n_row = 5000, n_col = 5000, overlap = 500) {
  
  require(terra)

  # Row and col index where the raster will be split
  rows <- seq(1, nrow(x), by = n_row)
  cols <- seq(1, ncol(x), by = n_col)
  
  # Index of cells that will be in the upper-left corner of the tiles
  cells <- cellFromRowColCombine(x, rows, cols)
  # Coordinates of the upper-left corner of the start cells
  xy_ul <- xyFromCell(x, cells)
  
  # Resolution
  rs <- res(x)
  
  ### Matrix of extents
  xy <- xyFromCell(x, cells)
  xy[,1] <- xy[,1] - rs[1]/2
  xy[,2] <- xy[,2] + rs[2]/2
  xy <- cbind(xy[,1], xy[,1] + n_col*rs[1], xy[,2] - n_row*rs[2], xy[,2])
  ###
  
  
  # Create matrix of extents (one ext per row)
  # xy <- cbind(
  #   xy_ul[,1],                   # x left
  #   xy_ul[,1] + n_col * rs[1],   # x right
  #   xy_ul[,2] - n_row * rs[2],   # y lower
  #   xy_ul[,2]                    # y upper
  # )
  # Matrix of extents with overlaps
  xy_ol <- cbind(
    xy[,1] - overlap * rs[1]/2,
    xy[,2] + overlap * rs[1]/2,
    xy[,3] - overlap * rs[2]/2,
    xy[,4] + overlap * rs[2]/2 
  )
  
  # Convert matrix to list of extents
  xy_lst    <- as.list(data.frame(t(xy)))
  xy_ol_lst <- as.list(data.frame(t(xy_ol)))
  
  # Expand raster to include external tiles overlaps
  ext(x) <- c(min(xy_ol[,1]), max(xy_ol[,2]), min(xy_ol[,3]), max(xy_ol[,4]))
  
  # Create IDs for each grid cell
  tileIDs <- data.frame(id = 1:length(xy_lst))
  tileIDs$id <- formatC(tileIDs$id, width = nchar(max(tileIDs$id)), flag = 0)

  # Cast non-overlapping grid to vect object with crs of x
  # # A lambda function with pipe forward. \() is equivalent to function(){}...Fancy!
  tiles <- lapply(xy_lst, \(i) ext(i) |> as.polygons()) |> 
    vect(crs = crs(x))
  # Add grid IDs
  values(tiles) <- tileIDs
  
  # Cast overlapping grid to vect object with crs of x
  OLtiles <- lapply(xy_ol_lst, \(i) ext(i) |> as.polygons()) |> 
    vect(crs = crs(x))
  values(OLtiles) <- tileIDs
  
  return(list(
    tiles    = tiles,
    OLtiles = OLtiles
    # new_ext   = ext(x)
  ))
  
} 


## Split raster into tiles
# r = input raster to be split into tiles
# grid = grid on which r should be tiled
# outpath = output location for tiles
# variable_name = name of raster variable. This will become the column name in the output table.
# prefix =  prefix for tile file names. Defaults to variable_name, but can be specified separately.

# Output:
# Raster tiles are written to disk. 
# a two-column dataframe is returned containing the grid ID from grid and the file path for each tile. 
create_raster_tiles <- function(r, grid, outpath, variable_name = "tile", prefix = NULL, dtype = NULL,overwrite = T,pfile = "progress.txt"){
  # browser()
  # if prefix is null set to variable_name 
  if(is.null(prefix)){
    prefix <- variable_name
  }
  # Get raster datatype
  if(is.null(dtype)){
    dtype <- datatype(r)[[1]]
  } 
  # Get the grid cell that actually overlap with the raster
  rext <- ext(r)
  rext <- vect(rext, crs = crs(r))
  rext <- st_as_sf(rext)
  # grid is transformed to the crs of the raster
  grid <- st_transform(grid, crs = st_crs(r))
  grid <- grid[rext,]
  
  
  # create empty dataframe for collecting grid ids and file paths.
  df <- data.frame("id" = NA, variable_name = NA)
  
  # extract each tile from the raster and write to disk
  # A version using apply was tested but turned out to be a little slower
  for(i in 1:nrow(grid)){
    cat('\r',"Tiling grid", grid[i,]$id, '', "\n")
    if(is.null(intersect(ext(r), ext(grid[i,])))){
      next
    }
    crp <- crop(r, ext(grid[i,]))
    crs(crp) <- crs(r)
    # Adding an extend function to handle cases where input data have smaller footprint than tile.
    crp <- extend(crp, ext(grid[i,]))
    cat("Writing tile as data type", dtype, '\n')
    names(crp) <- variable_name
    writeRaster(crp, filename = paste0(outpath, "/", prefix, "_", grid[i,]$id, ".tif"), datatype = dtype,overwrite = overwrite)
    
    write(paste(variable_name,". tile",grid[i,]$id),pfile,append = T)
    # get id and file path
    df[i,] <- data.frame(id = grid[i,]$id ,variable_name = paste0(outpath, "/", prefix, "_", grid[i,]$id, ".tif"))
  }
  # set field names and return tile index
  colnames(df) <- c("id", variable_name)
  return(df)
}



# trast <- test
# crt <- create_raster_tiles(trast, mgrid, outpath = "Y:/R1R2_VegMapping/Data/Raster/Source/elevation/testout",  
#                            prefix = "tiletest", return_index = T)


## Reduce overlapping tiles to non-overlapping extents.
# rlist = list of raster tiles (file names should contain grid ids)
# grid = non-overlapping grid vector file (defaults to master)
# outpath = output location
# prefix = character string to append to beginning of file names. 
trim_overlap <- function(rlist, grid = sf::st_read("Y:/MPSG_VegMapping/Data/Raster/spatial_database.gpkg", layer = "non-overlapping tiles"),
                                outpath, prefix, dtype = NULL){ # epsg = 5070
  library(stringr)
  library(dplyr)
  library(terra)
  library(sf)
  # Get raster datatype
  if(is.null(dtype)){
    dtype <- datatype(r)
  } 
  # ref <- paste0("EPSG:", epsg)
  lapply(1:length(rlist), function(i) {

    bn <- tools::file_path_sans_ext(basename(rlist[i]))
    # bn <- basename(rlist[i])
    # ndigits <- -str_count(bn, '\\d')
    ndigits <- -3
    getid <- str_sub(bn, start = ndigits)
    gselect <- filter(grid, id == getid)
    
    x <- rast(rlist[i])
    cat("Trimming overlap in tile", bn, "\n")
    # crs(x) <- crs(ref)
    names(x) <- prefix
    crop(x, ext(gselect), filename = paste0(outpath, "/", prefix, "_", gselect$id, ".tif"), datatype = dtype)
    return(unlist(paste0(outpath, "/", prefix, "_", gselect$id, ".tif")))
  })
}


## Create raster source table from scratch using a folder of tiles as template
initiate_rastersource_table <- function(tiles = "Y:/MPSG_VegMapping/Data/Raster/TemplateRaster/template_raster_tiles", 
                                        sourcedir = "Y:/MPSG_VegMapping/Data/Raster/Predictors",
                                        filename = "Y:/MPSG_VegMapping/Data/Raster/predictor_index.csv"){
  lf <- length(list.files(tiles))
  id <- 1:lf
  df <- data.frame(id)
  df$id <- sprintf("%03d", df$id)
  df$sourcedir <- sourcedir
  write.csv(df, filename, row.names = F, quote = F)
}


## Usage: for entering tile paths into raster source table
commit_rastersource <- function(tile_path, 
                                src_table = "Y:/MPSG_VegMapping/Data/Raster/predictor_index.csv",
                                replace_var = T,
                                save_to_file = F,
                                id_length = 3){
  
  # library(sf)
  library(stringr)
  library(dplyr)
  # The variable name in the raster source table column will match the respective directory name in the raster library
  variable_name <- basename(tile_path)
  # open raster source table
  src_t <- read.csv(src_table, colClasses = "character")
  # create dataframe of tile file paths sans extensions
  ls <- as.data.frame(tools::file_path_sans_ext(list.files(tile_path, pattern = ".tif$")))
  # set the variable name
  names(ls) <- variable_name
  
  # Extract the tile id from file names and add it as column
  ls$id <- str_sub(ls[,variable_name], start = -id_length)
  
  # Add file extension to file name
  ls[,variable_name] <- paste0(ls[,variable_name], ".tif")
  
  # Join the variable to the source table
  # If replace_var = T remove the existing version of the variable before joining
  if(replace_var){
    src_t <- src_t[, !names(src_t) %in% c(variable_name)]
  }
  
  j <- dplyr::left_join(src_t, ls, by = "id")
  
  if(save_to_file == T){
    # write backup copy
    # How many version have been written today?
    n_vers <- length(list.files("Y:/MPSG_VegMapping/Data/Raster/rastersource_backup", 
                                pattern = as.character(Sys.Date()), 
                                full.names = T))
    # write a copy to backup folder
    write.csv(src_t, paste0("Y:/MPSG_VegMapping/Data/Raster/rastersource_backup/rastersource_", Sys.Date(), "_", n_vers + 1, ".csv"),
              row.names = F, quote = F)
    
    # Replace Original
    write.csv(j, src_table, row.names = F, quote = F)
  }
  return(j)
}



tile_index_OTF <- function(tile_path, grid = sf::st_read("Y:/MPSG_VegMapping/Data/Raster/spatial_database.gpkg", layer = "non-overlapping tiles"),
                           variable_name = "variable"){
  
  # library(sf)
  library(sf)
  library(stringr)
  library(dplyr)
  # open precdictor index file
  # pindex <- read.csv(index_file, colClasses = "character")
  pindex <- st_drop_geometry(grid)
  # create dataframe of tile file paths sans extentions
  ls <- as.data.frame(tools::file_path_sans_ext(list.files(tile_path, pattern = ".tif$", full.names = T)))
  # set the variable name
  names(ls) <- variable_name
  bn <- basename(ls[1,])
  ndigits <- -str_count(bn, '\\d')
  ls$id <- str_sub(ls[,variable_name], start = ndigits)
  ls[,variable_name] <- paste0(ls[,variable_name], ".tif")
  # ls$id <- as.numeric(ls$id)
  j <- dplyr::left_join(pindex, ls, by = "id")
  return(j)
}


check_overlay <- function(x, template = "Y:/MPSG_VegMapping/Data/Raster/TemplateRaster/template_raster_tiles/template_raster_001.tif"){
  t <- rast(template)
  ch <- c(x, t)
  if(nlyr(ch) == 2){
    cat("Overlay check OK...")
  }
}
