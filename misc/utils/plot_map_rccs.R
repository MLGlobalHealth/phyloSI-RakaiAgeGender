library(rnaturalearth)
library(rnaturalearthdata)
library(osmdata)
require(ggpubr)
library(sf)
library(data.table)
library(ggspatial)

indir.deepsequence.data <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'map')

path.tests <- file.path(indir.deepsequence.data, 'RCCS_R15_R20', "all_participants_hivstatus_vl_220729.csv")
infile <- file.path(indir.deepsequence.data, 'RCCS_R15_R18', 'Rakai_community_geography_R15.rda')

# process dall
dall <- fread(path.tests)
dall <- dall[ROUND >= 15 & ROUND <= 18]

# rename variables according to Oli's old script + remove 1 unknown sex
setnames(dall,
         c('HIV_VL', 'COMM'),
         c('VL_COPIES', 'FC') )
dall[, HIV_AND_VL := ifelse( HIV_STATUS == 1 & !is.na(VL_COPIES), 1, 0)]
dall <- dall[! SEX=='']

DT <- copy(dall)
plotDT <- NA
plotcols <- 'PVLNSofHIV_MEAN'

if( ! is.data.table(plotDT))
{
  bool.plot.loc <- TRUE
}
if( is.data.table(plotDT) & !is.na(plotcols))
{
  bool.plot.loc <- FALSE
  stopifnot( all(plotcols %in% names(plotDT)))
  cols <- c('COMM_NUM', 'ROUND', 'SEX', plotcols)
  plotDT <- plotDT[, ..cols]
}

# theme
mytheme <- theme(
  panel.background = element_rect(fill = 'lightblue', color='black'),
  panel.grid.major = element_line(colour = "transparent")
)
noticks <- theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank()) 


palettes <- list(
  sex = wesanderson::wes_palette("Moonrise3")[c(2,1)],
  comm = wesanderson::wes_palette("Moonrise3")[c(1,3)],
  arvmed = c('#E03531', '#51B364'),
  cuarvmed = c('#E03531', '#51B364'),
  supp_status =  c('#8F190E', '#DBCB00', '#00DB25')
)


# Fix factors stored as integers (COMM_NUM)
# TODO: merge plotDT with rest
tmp <- new.env()
load(infile, envir=tmp)
ds <- as.data.table(tmp$comgps)
cols <- c('latitude', 'longitude')
ds[, (cols):=lapply(.SD, unlist), .SDcols=cols]
.f <- function(x) as.integer(as.character(x))
ds[, COMM_NUM := .f(COMM_NUM)]
tmp1 <- unique(DT[, .(COMM_NUM, FC)])
ds <- merge(ds, tmp1)
ds[is.na(longitude)]
setnames(ds, cols, paste0(toupper(cols), '_JITTER'))
ds_sf <- st_as_sf(ds, coords=c("LONGITUDE_JITTER", "LATITUDE_JITTER"), crs=4326)
ds_sf_t <- st_transform(ds_sf, crs=2163)



coord2 <- data.table(
  x=range(ds$LONGITUDE_JITTER)+c(-.1, .1), 
  y=range(ds$LATITUDE_JITTER) +c(-.1, .1)
)


# Latitude from -1.28538 to 3.66088 
# and longitude from 29.65 to 34.66659.

long_bounds=c(22, 43) # x
lati_bounds=c(-10, 10)  # y

# get uganda borders + lakes
# __________________________

uganda_borders <- c('uganda', 'kenya', "Tanzania", 'rwanda', 'democratic republic of the congo')


africa_sf <- ne_countries(scale=50, continent='africa', returnclass='sf')
uganda_sf <- ne_countries(scale=50, country=c('uganda', 'tanzania'), returnclass='sf')

# Define bounding box:
bb <- c(long_bounds[1], lati_bounds[1],
        long_bounds[2], lati_bounds[2])
lakes <- opq(bbox=bb, timeout=150) |>
  add_osm_feature(key='natural', value='water') |>
  osmdata_sf()

# Download data of interest
lakes10 <- ne_download(scale = 10, type = 'lakes', category = 'physical', returnclass='sf')
roads10 <- ne_download(scale = 10, type = 'roads', category = 'cultural', returnclass='sf')
rivers10 <- ne_download(scale=10, type='rivers_lake_centerlines', category='physical', returnclass='sf')

p1 <- ggplot() +
  geom_sf(data=africa_sf, fill='darkolivegreen4') + 
  geom_sf(data=lakes10,col= 'lightblue4', fill="lightblue") +
  geom_sf_text(data=africa_sf, aes(label=name), color = "black", fontface = "bold", check_overlap = FALSE) +
  geom_rect(data=coord2, aes(xmin = x[1], xmax = x[2], ymin = y[1], ymax = y[2]), color = "gold", fill = NA)  +
  coord_sf(xlim=long_bounds, ylim=lati_bounds) +
  labs(x='', y='') + 
  mytheme + 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()) + 
  annotation_scale(location = "tr", width_hint = 0.5) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.3, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)

p2 <- ggplot() +
  geom_sf(data=africa_sf, fill='darkolivegreen4', color = 'darkolivegreen4') +
  geom_sf(data=rivers10, color="lightblue", fill="lightblue") +
  geom_sf(data=roads10, color="grey80") +
  geom_sf(data=lakes10,  color="lightblue4", fill="lightblue") +
  geom_sf(data=ds_sf_t, aes(fill=FC, shape = FC),size=3) +
  coord_sf(xlim=coord2$x, ylim=coord2$y) +
  geom_rect(data=coord2, aes(xmin = x[1], xmax = x[2], ymin = y[1], ymax = y[2]), color = "gold", fill = NA)  +
  scale_fill_manual(values=c('fishing' = 'dodgerblue3', 'inland' = 'goldenrod3'), 
                     labels = c('Fishing communities', 'Inland communities'))  + 
  scale_shape_manual(values=c('fishing' = 24, 'inland' = 21), 
                    labels = c('Fishing communities', 'Inland communities')) +
  mytheme +
  theme(legend.position=c(0,0),
        legend.background = element_rect(fill = "white", color = "grey50"),
        legend.justification = c("left", "bottom"),
        panel.background = element_rect(fill = 'lightblue', color='lightblue4')) +
  labs(fill='Community type', shape='Community type')+ 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank())+
  annotation_scale(location = "tr", width_hint = 0.5) 

# CAN I PUT THE LEGEND ON NORTH-EAST CORNER to save space?
p <- gridExtra::grid.arrange(p1, p2, layout_matrix = rbind(c(NA,2), 
                                                           c(1, 2), 
                                                           c(NA, 2)), heights = c(0.09, 1.2, 0.03))

ggsave(p, file = file.path(outdir, 'map_RCCS_communities.pdf'), w = 9, h = 4.7)


