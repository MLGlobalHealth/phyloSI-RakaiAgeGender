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

long_bounds=c(25, 46) # x
lati_bounds=c(-12, 10)  # y

long_bounds=c(27, 40) # x
lati_bounds=c(-8, 6)  # y

# get uganda borders + lakes
# __________________________

uganda_borders <- c('uganda', 'kenya', "Tanzania", 'rwanda', 'democratic republic of the congo')


africa_sf <- ne_countries(scale=50, continent='africa', returnclass='sf')
uganda_sf <- ne_countries(scale=50, country=c('uganda', 'tanzania'), returnclass='sf')

# Define bounding box:
bb <- c(long_bounds[1], lati_bounds[1],
        long_bounds[2], lati_bounds[2])
# lakes <- opq(bbox=bb, timeout=150) |>
#   add_osm_feature(key='natural', value='water') |>
#   osmdata_sf()

# Download data of interest
lakes10 <- ne_download(scale = 10, type = 'lakes', category = 'physical', returnclass='sf')
roads10 <- ne_download(scale = 10, type = 'roads', category = 'cultural', returnclass='sf')
rivers10 <- ne_download(scale=10, type='rivers_lake_centerlines', category='physical', returnclass='sf')

p1 <- ggplot() +
  geom_sf(data=africa_sf, fill='palegoldenrod') + 
  geom_sf(data=lakes10,col= 'lightblue4', fill="lightblue") +
  geom_sf_text(data=africa_sf, aes(label=name), color = "black", fontface = "bold", check_overlap = FALSE) +
  geom_rect(data=coord2, aes(xmin = x[1], xmax = x[2], ymin = y[1], ymax = y[2]), color = "red", fill = NA, size = 1.2)  +
  coord_sf(xlim=long_bounds, ylim=lati_bounds) +
  labs(x='', y='') + 
  mytheme + 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()) 

 # m <- leaflet::leaflet() %>%
 #          setView(lng = mean(coord2$x), lat = mean(coord2$y), zoom = 12) %>%
 #         addTiles() %>%
 #         addProviderTiles(providers$CartoDB.Positron)
 # m
 
p2 <- ggplot() +
  geom_sf(data=africa_sf, fill='palegoldenrod', color = 'palegoldenrod') +
  geom_sf(data=rivers10, color="lightblue", fill="lightblue") +
  geom_sf(data=roads10, color="grey80") +
  geom_sf(data=lakes10,  color="lightblue4", fill="lightblue") +
  geom_sf(data=ds_sf_t, aes(shape = FC),size=3,fill='#876445') +
  coord_sf(xlim=coord2$x, ylim=coord2$y) +
  geom_rect(data=coord2, aes(xmin = x[1], xmax = x[2], ymin = y[1], ymax = y[2]), color = "red", fill = NA, size = 1.75)  +
  scale_shape_manual(values=c('fishing' = 24, 'inland' = 21), 
                    labels = c('Fishing communities', 'Inland communities')) +
  mytheme +
  theme(legend.position=c(0.6,0.85),
        legend.key.size = unit(0.6, 'pt'),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.justification = c("left", "bottom"),
        panel.background = element_rect(fill = 'lightblue', color='lightblue4')) +
  # labs(shape='')+ 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        legend.key = element_rect(fill = "white"),
        legend.title = element_blank())+
  annotation_scale(location = "br", width_hint = 0.5, pad_x = unit(0.1, 'cm'), pad_y = unit(0.05, 'cm')) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.3, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)


ggsave(p1, file = file.path(outdir, 'map_RCCS_communities_1.pdf'), w = 2.5, h = 2.5)
ggsave(p2, file = file.path(outdir, 'map_RCCS_communities_2.pdf'), w = 4, h = 4)

