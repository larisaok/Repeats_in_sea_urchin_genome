library(data.table)
library(ggplot2)
library("RColorBrewer")
library(plyr)
library(reshape2)

setwd("~/Urchin_projects/Kimura_plot/Kimura_Spur/")
ce_rfam <- fread("Spur_5.0_genome.fna.align.landscape.Div.Rfam.tab")
ce_rname <- fread("Spur_5.0_genome.fna.align.landscape.Div.Rname.tab")
ce_rclass <- fread("Spur_5.0_genome.fna.align.landscape.Div.Rclass.tab")

names(ce_rclass) <- c("Class", c(1:(ncol(ce_rclass)-1)))
names(ce_rname) <- c("Name", "Class", "SubClass", c(1:(ncol(ce_rname)-3)))
names(ce_rfam) <- c("Class", "Family", c(1:(ncol(ce_rfam)-2)))

ce_rclass <- melt(ce_rclass, measure.vars=c(2:ncol(ce_rclass)), variable.name = "Divergence")
ce_rname <- melt(ce_rname, measure.vars=c(4:ncol(ce_rname)), variable.name = "Divergence")
ce_rfam <- melt(ce_rfam, measure.vars=c(3:ncol(ce_rfam)), variable.name = "Divergence")

ce_rclass$value <- ce_rclass$value/921855793*100
ce_rname$value <- ce_rname$value/921855793*100
ce_rfam$value <- ce_rfam$value/921855793*100

setDT(ce_rclass)
setDT(ce_rname)
setDT(ce_rfam)
ce_rclass <- ce_rclass[value > 0]
ce_rname <- ce_rname[value > 0]
ce_rfam <- ce_rfam[value > 0]

colourCount = length(unique(ce_rclass$Class))
palette_rclass <- colorRampPalette(brewer.pal(12, "Paired"))(colourCount)

rclass_plot <- ggplot(data=ce_rclass, aes(x=Divergence, y=value, fill=Class))+
  geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(name="DNA Families", values = palette_rclass) + 
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion (%)") +
  theme_bw()
rclass_plot
ggsave("rep_class.png", plot = rclass_plot, device = "png", dpi = 700) 

colourCount = length(unique(ce_rname$Name))
palette_rname <- colorRampPalette(brewer.pal(12, "Paired"))(colourCount)

rname_plot <- ggplot(data=ce_rname, aes(x=Divergence, y=value, fill=Name))+
  geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(name="DNA Families", values = palette_rname) + 
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion (%)") +
  theme_bw() +
  theme(legend.position = "none")
rname_plot
ggsave("rep_name.png", plot = rname_plot, device = "png", dpi = 700) 

colourCount = length(unique(ce_rname$SubClass))
palette_rSubClass <- colorRampPalette(brewer.pal(12, "Paired"))(colourCount)

rSubClass_plot <- ggplot(data=ce_rname, aes(x=Divergence, y=value, fill=SubClass))+
  geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(name="DNA Families", values = palette_rSubClass) + 
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion (%)") +
  theme_bw()
rSubClass_plot
ggsave("rep_SubClass.png", plot = rSubClass_plot, device = "png", dpi = 700) 


colourCount = length(unique(ce_rfam$Family))
palette_rfam <- colorRampPalette(brewer.pal(12, "Paired"))(colourCount)

rfam_plot <- ggplot(data=ce_rfam, aes(x=Divergence, y=value, fill=Family))+
  geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(name="DNA Families", values = palette_rfam) + 
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion (%)") +
  theme_bw()
rfam_plot
ggsave("rep_Families.png", plot = rfam_plot, device = "png", dpi = 700) 
