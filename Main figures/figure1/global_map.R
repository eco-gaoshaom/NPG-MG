library(ggplot2)
library(ggrepel)
library(maps)

site <- read.delim('capital_location.txt')

p <- ggplot() +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = rgb(220,220,220,max=255)) +  #世界地图模板
  theme_bw() +  
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(-150, -100, -50, 0, 50, 100, 150), expand = c(0, 0), 
                     labels = c('150°W', '100°W', '50°W', '0', '50°E', '100°E', '150°E')) +
  scale_y_continuous(breaks = c(-60, -30, 0, 30, 60), expand = c(0, 0), 
                     labels = c('60°S', '30°S', '0', '30°N', '60°N')) +
  labs(x = 'Longitude', y = 'Latitude', color = 'ocean and sea regions')

p
p1 <- p +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = rgb(220,220,220,max=255)) +
  geom_point(data = site, aes(x = long, y = lat, color = continent,shape= source), size = 2) +
  #geom_label_repel(data = site, aes(x = Longitude, y = Latitude, label = Station.identifier), size = 2, show.legend = FALSE) +  
  scale_color_manual(values = c('#F99233', '#3CA6CC', '#686687', '#01B250'))+
  scale_shape_manual(values = c(17,16))

p1
