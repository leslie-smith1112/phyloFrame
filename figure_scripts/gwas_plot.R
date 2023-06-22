### - FIGURE 1B - ### 


library(ggplot2)
library(dplyr)
time <- as.numeric(rep(seq(1,7),each=7))  # x Axis
value <- runif(49, 10, 100)               # y Axis
group <- rep(LETTERS[1:7],times=7)        # group, one shape per group
value <- round(value)

data <- data.frame(time, value, group)
ggplot(data, aes(x=time, y=value, fill=group)) + geom_area()
# data downloaded from https://www.ebi.ac.uk/gwas/docs/file-downloads # 
data <- read.table("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/remake_gwas_ancestry_info2.tsv", sep='\t', header=T)
# ancestries will be European, Asian, African, America, Oceanian 
colnames(data) <- c("year", "ancestry", "value")

table(data$ancestry)
data$ancestry <- factor(data$ancestry, levels =c("European","Asian","African","America",
                                                 "Oceanian") )
# GWAS data plot Figure 1A 
ggplot(data, aes(x=year, y=value, fill=ancestry, color = ancestry)) + geom_area(alpha=0.9, linewidth=.25, colour="white") + theme_minimal() +
  scale_fill_manual(name = "", values = c("European" = "#CA4136", "Asian" = "#496849", "African" = "#1F619E", "America" = "#FFA500", 
                                          "Oceanian" = "#BF812D")) +
  scale_color_manual(name = "", values = c("European" = "#CA4136", "Asian" = "#496849", "African" = "#1F619E", "America" = "#FFA500", 
                                          "Oceanian" = "#BF812D")) 

total.pop <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/total_world_pop_2023.xlsx")
#%% downloaded from: https://population.un.org/wpp/Download/Standard/MostUsed/ ## 

total.pop$simple <- "temp"
#change country names to be uniform 
total.pop$simple[total.pop$Ancestry == "Countries with Access to the Sea: Africa"] <- "Africa"
total.pop$simple[total.pop$Ancestry == "Countries with Access to the Sea: Asia"] <- "Asia"
total.pop$simple[total.pop$Ancestry == "Countries with Access to the Sea: Europe"] <- "Europe"
total.pop$simple[total.pop$Ancestry == "Countries with Access to the Sea: Latin America and the Caribbean"] <- "America"

total.pop$simple[total.pop$Ancestry == "Countries with Access to the Sea: Northern America"] <- "America"
total.pop$simple[total.pop$Ancestry == "Countries with Access to the Sea: Oceania"] <- "Oceania"
total.pop$simple[total.pop$Ancestry == "LLDC: Africa"] <- "Africa"
total.pop$simple[total.pop$Ancestry == "LLDC: Asia"] <- "Asia"
total.pop$simple[total.pop$Ancestry == "LLDC: Europe"] <- "Europe"
total.pop$simple[total.pop$Ancestry == "LLDC: Latin America"] <- "America"
total.pop$simple[total.pop$Ancestry == "LLDC: Africa"] <- "Africa"
total.pop <- total.pop[!total.pop$simple == "temp",]
#simple = the simplified country name (Removed the LLDC and or the Countries with Access to the Sea)
new <- total.pop  %>% group_by(Year, simple) %>% summarise(n = sum(totalpopulation)) %>% mutate(percentage = n / sum(n))
# get current world population proportions
new.y <- new[new$Year == 2023,]
new.y$simple <- factor(new.y$simple, levels = c("Europe", "Asia", "Africa", "America", "Oceania"))

ggplot(new.y, aes(fill=simple, y=n, x=Year, color = simple)) + 
  geom_bar(position="fill", stat="identity") + theme_minimal() +
  scale_fill_manual(name = "", values = c("Europe" = "#CA4136", "Asia" = "#496849", "Africa" = "#1F619E", "America" = "#FFA500",
                                          "Oceania" = "#BF812D")) +
  scale_color_manual(name = "", values = c("Europe" = "#CA4136", "Asia" = "#496849", "Africa" = "#1F619E", "America" = "#FFA500",
                                           "Oceania" = "#BF812D")) + scale_x_continuous(breaks = 2023)




