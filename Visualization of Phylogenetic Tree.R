library(data.table)
library(tidyverse)
library(ggthemes)
library(treeio)
library(ape)
library(ggtree) 
library(ggtreeExtra) 
library(ggnewscale)


tree <- read.newick("Plant-1.INF.treefile", node.label = "support")

digital <- fread("Plant-1.INF.digital.txt", sep = " ", header = F) %>% as.data.frame()
dim(digital)
names(digital) <- c("leaf", seq(1,1030))

digital <- digital %>%
                    mutate_if(is.numeric,as.character) %>% 
                    pivot_longer(cols = -1, names_to = "pos", values_to = "type") %>%
                    mutate(type=case_when(type=="A" ~ "A-T", type=="B" ~ "A-C", type=="C" ~ "A-G",
                                          type=="D" ~ "T-A", type=="E" ~ "T-C", type=="F" ~ "T-G",
                                          type=="G" ~ "C-A", type=="H" ~ "C-T", type=="I" ~ "C-G",
                                          type=="J" ~ "G-A", type=="K" ~ "G-T", type=="L" ~ "G-C",
                                          type=="d" ~ "Indel", type=="0" ~ "unchange"))

fulllength_append <- data.frame(leaf=unique(digital$leaf), pos="1", type="FL")


mutation <- rbind(digital, fulllength_append) %>% filter(type!="unchange") %>% 
                  mutate(pos=as.numeric(pos)) %>% 
                  mutate(length=ifelse(type=="FL",1030,2)) %>% 
                  mutate(sample=str_split_fixed(leaf,"_",2)[,1]) %>% 
                  mutate(shoot=str_split_fixed(sample,"[.]",4)[,2]) %>% 
                  arrange(factor(leaf, levels = tree@phylo$tip.label), pos) 


ggtree(tree, layout = "roundrect", alpha=0) %<+% mutation[c("leaf","shoot")] +
      geom_tree(aes(color=shoot), layout = "roundrect") +
      geom_nodelab(aes(label=paste(support, node, sep = "|")), size=1) +
      geom_tiplab(aes(label=sample, color=shoot), size=1) +
      scale_color_manual(values = shoot_color) +

      new_scale_color() +
      geom_facet(panel = "mutation", data = mutation, geom = geom_segment,
                  aes(x=pos, xend=pos+length, y=y, yend=y, color=type)) +
      scale_color_manual(values = seg_color) +
      scale_x_continuous(expand=c(0,0)) +
   theme_tree2()



