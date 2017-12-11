install.packages(vegan)


library(vegan)
library(ggplot2)


setwd("C:/Users/Zoe/Documents/Data/Lab\ Work/Dropbox/Labwork/Enzymes/Enzymes\ 2017")
map<- read.csv("Enzymes_MF.csv")

#check names
names(map)
str(map)
attach(map)
detach(map)

#read in normalized data

norm<- read.csv("2017Data_Normalized.csv")

names(norm)
str(norm)
attach(norm)
detach(norm) 

#make vector
sitevec<-norm$Sample
norm$Sample<-NULL
ord<-metaMDS(norm)
plot(ord)

#calculate ordination distances
pH = factor(pH)
Location = factor(Location)
moisture = factor(moisture)
orgc = factor(orgc)
ordvegdist1<- vegdist(norm, method = "bray")
adonis(ordvegdist1~ pH + moisture + orgc + Location, data = map, permutations = 1000)
ord.fit<-envfit(ord~., data=map, perm=1000, na.rm=TRUE)
ord.fit
ord.fit3<-envfit(ord~ pH + moisture + orgc, data=map, perm=1000, na.rm=TRUE)


#plotting ordination

sv1<-as.character(sitevec)
plot(ord)
fig1<-ordiplot(ord, type = "text")
fig1
NMDS1<-fig1$sites[,1]
NMDS1
NMDS2<-fig1$sites[,2]
NMDS2
plot(fig1$sites[,1],fig1$sites[,2],pch=sv1,cex=0.8,xlab='NMDS1',ylab='NMDS2')


#to plot elips of 95% confidence interval around site and treatment factors
with(map,ordiellipse(ord,Location,kind="se",conf=0.95))
#site_tr_vec3<-c(2,2,2,1,1,1,2,2,2,1,1,1)
#site_tr_vec2<-c(1,2,1,2)
with(enz_metadata,ordiellipse(ord,site_tr_vec,kind="se",conf=0.95,lty=1))
with(enz_metadata,ordiellipse(ord,Stand_Phase,kind="se",conf=0.95))
#to plot only variable with a certain p-value (p.max)
plot(ord.fit2, p.max=0.05, col="black", cex=0.8)
plot(ord.fit, p.max=0.05)


#make plot look good
grp.t <- data.scores[data.scores$Location == "T", ][chull(data.scores[data.scores$Location == 
                                                                   "T", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.m <- data.scores[data.scores$Location == "M", ][chull(data.scores[data.scores$Location == 
                                                                   "M", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.b <- data.scores[data.scores$Location == "B", ][chull(data.scores[data.scores$Location == 
                                                                   "B", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
hull.data <- rbind(grp.t, grp.m, grp.b)  #combine grp.a and grp.b
hull.data
data.scores<- as.data.frame(scores(ord))
data.scores
data.scores$site<- rownames(data.scores) 
data.scores$site2<- sitevec
data.scores$Location<- substr(data.scores$site2, start = 1, stop = 1)
a<- ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Location,group=Location),alpha=0.30) + # add the convex hulls
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=Location,colour=Location),size=4) + # add the point markers
  scale_colour_manual(values=c("T" = "blue", "M" = "green", "B" = "red")) +
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
a<- a + ggtitle("Differences in microbial enzymatic profiles with soil depth \n as a function of moisture, organic carbon, and pH" 
)
a<- a + theme(plot.title = element_text(20, hjust = 0.5))
