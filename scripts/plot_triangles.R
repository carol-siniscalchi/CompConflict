
## modified from https://github.com/MaxTelford/MonoDeutData/blob/master/scripts/plot_triangles.R

library(ggtern)
library(plyr)

points<- data.frame(
  rbind(c( 1,1,0,0),
        c( 2,0.66666,0,0.33333),
        c( 3,0.66666,0.33333,0),
        c( 4,0,1,0),
        c( 5,0.33333,0.66666,0),
        c( 6,0,0.66666,0.33333),
        c( 7,0,0,1),
        c( 8,0,0.33333,0.66666),
        c( 9,0.33333,0,0.66666),
        c(10,0.16666,0.16666,0.66666),
        c(11,0.66666,0.16666,0.16666),
        c(12,0.16666,0.66666,0.16666)
  )
)

colnames(points) = c("IDPoint","L","T","R")

## obtained with the likelihood_transform.py script (https://github.com/MaxTelford/MonoDeutData/blob/master/scripts/likelihood_transform.py) 
triplets <- read.csv("~/LikelihoodMapping/triplets.tsv", sep = "\t")

polygon.labels<- data.frame(
  Label=c("T1",
          "T13",
          "T3",
          "T23",
          "T2",
          "T12",
          "T123"
  ))

#Assign each label an index
polygon.labels$IDLabel=1:nrow(polygon.labels)

#Create a map of polygons to points
polygons <- data.frame(
  rbind(c(1,1),c(1,2),c(1,3),
        c(2,3),c(2,11),c(2,12),c(2,5),
        c(3,5),c(3,4),c(3,6),
        c(4,6),c(4,12),c(4,10),c(4,8),
        c(5,8),c(5,7),c(5,9),
        c(6,9),c(6,10),c(6,11),c(6,2),
        c(7,11),c(7,12),c(7,10) )
)
#IMPORTANT FOR CORRECT ORDERING.
polygons$PointOrder <- 1:nrow(polygons)

#Rename the columns
colnames(polygons) = c("IDLabel","IDPoint","PointOrder")

#Merge the three sets together to create a master set.
df <- merge(polygons,points)
df <- merge(df,polygon.labels)
df <- df[order(df$PointOrder),]

#Determine the Labels Data library(plyr)
Labs = ddply(df,"Label",function(x){c(c(mean(x$T),mean(x$L),mean(x$R)))})
colnames(Labs) = c("Label","T","L","R")

pch = 21  

L_comp <- "#d40000"
R_comp <- "#aa4400"
Top_comp <- "#de8787"

LR <- "lightgrey"
LT <- "lightgrey"
TR <- "lightgrey"
LTR <- "grey"

my_named_vector_comp = c(
  "T1" = "black",
  "T2" = "black",
  "T3" = "black",
  "T12"  = "black",
  "T23" = "black",
  "T13" = "black",
  "T123" = "black")

ternary_comp <- ggtern(data=df, 
                  aes(L,T,R)) +
  geom_polygon(aes(fill=Label), color="black", alpha = 0.6) +
  geom_point(triplets, mapping = aes(triplets$T1, 
                                triplets$T2,
                                triplets$T3,
                                color=FinalT
             ),
             alpha=0.6,size = 0.2,
             show.legend = FALSE
  )+
  theme_nolabels() +
  theme_nogrid() +
  scale_fill_manual(values = c(L_comp, LR, LTR, LT, R_comp, TR, Top_comp)) +
  scale_color_manual(values = my_named_vector_comp, na.translate = FALSE)

ggsave(filename="triangle_comp.svg", plot=ternary_prot, width=10, height=8)

## topologies and labels added with Adobe Illustrator
