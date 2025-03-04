
# load the appropriate R package libraries
library(rnaturalearth)
library(maps)
library(sp)
library(sf)
library(ggplot2)

#hudson
#hudson.n_UTM <- data.frame(x=c(357212.46,19314.01,658471.31,935303.78),y=c(4261346,4566676,5099985,4635884))
spRef_UTM_19 <- "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
spRef_DD <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"
hudson.n_UTM <- data.frame(x=c(357212.46,19314.01),y=c(4261346,4566676))
hudson.n_UTM.rl <- rbind(hudson.n_UTM, hudson.n_UTM[1,])
hudson.n_UTM.l <- Line(hudson.n_UTM.rl)
hudson.n_UTM.l <- Lines(list(hudson.n_UTM.l), ID = 'n')
hudson.n_UTM.spL <- SpatialLines(list(hudson.n_UTM.l))
proj4string(hudson.n_UTM.spL) <- CRS(spRef_UTM_19)
hudson.n <- spTransform(hudson.n_UTM.spL, CRS(spRef_DD))
hudson.n.sf <- st_as_sf(hudson.n,coords = c("x", "y"))
world <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

# contours <- sf::st_as_sf(sf::read_sf(("c:/work/shapefiles/Contours"), quiet=TRUE))
# contour_400m <- subset(contours, elev_m == -400)
# saveRDS(contour_400m, here::here("data","contour_400m.RDS"))
contour_400m <- readRDS(here::here("data","contour_400m.RDS"))
lat.range <- c(34.5,42)
lon.range <- c(-78,-70)
lat.range <- c(35,45)
lon.range <- c(-76,-66)

state_prov <- rnaturalearth::ne_states(c("united states of america", "canada"))

# ggmap <- ggplot(states) +
ggmap <- ggplot(state_prov) +
  geom_sf(alpha = 0.5) +
  # geom_sf(data = states, alpha = 0.5) +
  geom_sf(data = contour_400m, linewidth = 0.1) +
  geom_sf(data = hudson.n.sf, linewidth = 1, colour = "red") +
  xlim(lon.range) +
  ylim(lat.range) +
  theme_bw() +
  annotate("text", x = -70, y = 37, label = "Atlantic Ocean") +
  annotate("text", x = -68, y = 43, label = "Gulf\nof\nMaine") +
  theme(legend.title=element_blank(), axis.title = element_blank())

cairo_pdf(file.path("paper", "map.pdf"), width = 7, height = 7)
ggmap
dev.off()
