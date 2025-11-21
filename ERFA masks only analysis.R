# clear the environment: 
rm(list = ls())
# clear all the plots: 
dev.off()
# clear the console: ctrl+L

#### Packages for this script ####
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggsignif")
install.packages("ggpubr")
install.packages("ggbeeswarm")

library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(ggbeeswarm)
#### end ####

# Set location of data and where plots will be saved to:
mainDir <- "C:/Users/jonatmt/OneDrive - Universitetet i Oslo/Desktop/Amalie masks"
setwd(mainDir)
# Create the path where the data is stored
DataPath <- "C:/Users/jonatmt/OneDrive - Universitetet i Oslo/Desktop/Amalie masks"

# Set the pixel size and focal adhesion sizes
pixel_length_x <- 0.1083
pixel_length_y <-0.1083
pixel_area <- pixel_length_x*pixel_length_y
MinFAsize <- 0.15
MaxFAsize <- 10

#### Read in the data and keep focal adhesion data ####
MeanPixelLength <- (pixel_length_x+pixel_length_y)/2
## Read in all the csv files and store in a list of dataframes
AllRawData <- list()
# Find all csv files within the DataPath directory
Files <- list.files(DataPath, pattern = "*.csv")
# Generate a list of Filepaths to be read in with the directory followed by the
# csv filename
FilePaths <- list()
for (i in 1:length(Files)){
  FilePaths[[i]] <- paste(DataPath,Files[i], sep ="/")
}
for (i in 1:length(FilePaths)){
  # Read the data from each csv file
  df <- read.csv(FilePaths[[i]])
  # Change the data to be in ??m and ??m^2 instead of pixels
  df$area <- df$area*pixel_area
  df$area_bbox <- df$area_bbox*pixel_area
  df$area_convex <- df$area_convex*pixel_area
  df$area_filled <- df$area_filled*pixel_area
  df$equivalent_diameter_area <- df$equivalent_diameter_area*pixel_area
  df$ERThresh_intersection_area <- df$ERThresh_intersection_area*pixel_area
  df$axis_major_length <- df$axis_major_length*MeanPixelLength
  df$axis_minor_length <- df$axis_minor_length*MeanPixelLength
  df$feret_diameter_max <- df$feret_diameter_max*MeanPixelLength
  df$perimeter <- df$perimeter*MeanPixelLength
  df$perimeter_crofton <- df$perimeter*MeanPixelLength
  # Assign data to the raw data list
  AllRawData[[i]] <- df
  # Name each dataframe in the list after the csv file it came from
  names(AllRawData)[i] <- Files[i]
  rm(df)
}
rm(FilePaths, Files, i)

FARawData <- list()
for (s in 1:length(AllRawData)) {
  df <- AllRawData[[s]]
  # Select only FAs that are between the minimum and maxium defined size
  subset_df <- df[df$area > MinFAsize & df$area < MaxFAsize, ]
  # Count how many FAs are in the cell
  subset_df$NumberOfAdhesions <- nrow(subset_df)
  # Calculate the sum area of focal adhesions
  subset_df$TotalFAarea <- sum(subset_df$area)
  # Add a column containing the name of the cell
  subset_df <- mutate(subset_df, Name = names(AllRawData)[s])
  # Add a column containing the surface condition the cells were grown on
  subset_df <- mutate(subset_df, Condition = strsplit(names(AllRawData)[s], "_")[[1]][1])
  # Add a column containing the date of the experiment
  subset_df <- mutate(subset_df, Date = strsplit(names(AllRawData)[s], "_")[[1]][2])
  FARawData[[s]] <- subset_df
  names(FARawData)[s] <- names(AllRawData)[s]
  rm(df, subset_df)
}

AllFAData <- FARawData[[1]]
for (s in 2:length(FARawData)) {
  Nextdf <- FARawData[[s]]
  AllFAData <- bind_rows(AllFAData, Nextdf)
  rm(Nextdf)
}

MeanCellData <- data.frame()
# Calculate column means for each dataframe
for (s in 1:length(AllRawData)) {
  df <- AllRawData[[s]]
  # Keep only FAs within size limit
  subset_df <- df[df$area > MinFAsize & df$area < MaxFAsize, ]
  # Count how many FAs are in the cell
  subset_df$NumberOfAdhesions <- nrow(subset_df)
  # Number of FAs which also intersect with ER
  subset_df$NumFAintersectER <- nrow(subset_df[subset_df$ERThresh_num_intersection_areas >= 1, ])
  # Calculate the sum area of focal adhesions
  subset_df$TotalFAarea <- sum(subset_df$area)
  # Convert subset_df to a data frame
  subset_df <- as.data.frame(subset_df)
  # Calculate column means using summarize
  mean_df <- subset_df %>% summarize(across(everything(), \(x) mean(x, na.rm = TRUE)))
  mean_df$PercERFAintersection <- (mean_df$NumFAintersectER/mean_df$NumberOfAdhesions)*100
  mean_df <- mutate(mean_df, Name = names(AllRawData)[s])
  mean_df <- mutate(mean_df, Condition = strsplit(names(AllRawData)[s], "_")[[1]][1])
  mean_df <- mutate(mean_df, Date = strsplit(names(AllRawData)[s], "_")[[1]][2])
  MeanCellData <- bind_rows(MeanCellData, mean_df)
  rm(subset_df, mean_df, df)
}

rm(s)
#### end ####

#### Plotting functions ####
## Define a basic theme for the plots
theme_jmt <- function(){theme(
  text = element_text(family = "sans", size = 25), # set default font to be Arial
  plot.title = element_text(hjust=0.5, face = 'plain'), # align title to centre and make bold
  panel.grid.major = element_line(colour = "grey80", size = 0.1), # strip major gridlines
  panel.grid.minor = element_blank(), # strip minor gridlines
  panel.background = element_blank(), # strip background
  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),  # tilt the x axis
  axis.line = element_line(colour = "black", size = 1.5), # add axis lines in black
  axis.ticks = element_line(colour = "black", size = 1.5),
  axis.ticks.length=unit(2,"mm")) # add tick marks in black
}

## Where required, labels are the Title, X axis and Y axis labels defined in a character string (in order)
## e.g. Labels <- c(title = "Normalised to Background value", x = "Distance (um)", y = "Normalised fluorescence")
## these are then read into the function.

## Define a function to plot the relative intensities of the bands
Boxviolin <- function(Data, Labels){
  ggplot(Data, aes(Condition,area,col=Condition,shape=Condition)) + theme_jmt() +
    geom_violin() + geom_boxplot(width=0.1, show.legend = FALSE,colour = "Grey60",alpha = 0.5, outlier.shape = NA) + geom_jitter() +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

Crossbar <- function(Data, y_var, Labels){
  ggplot(Data,aes(Condition,{{y_var}}, col=Condition,shape=Condition)) + theme_jmt() +
    geom_jitter(position = position_jitterdodge()) +
    stat_summary(fun="mean", geom="crossbar", width=0.5, color = "black") + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

Superplot_statistics <- function(Plot, Data, y_var, Yaxis_Max, Yaxis_Min, Yaxis_step){
  ReplicateAverages <- Data %>% group_by(Condition, Date) %>% summarise_all(mean)
  unique_comparators <- unique(ReplicateAverages$Condition)
  
  NumberOfComparisons <- (0.5 * length(unique_comparators) - 0.5) * length(unique_comparators)
  Yaxis_length <- Yaxis_Max - Yaxis_Min
  TopY <- Yaxis_Min + (Yaxis_length * 0.9)
  y_pos_values <- rev(seq((TopY - (Yaxis_length * 0.1 * (NumberOfComparisons - 1))), TopY, by = (Yaxis_step * 0.5)))
  CompNum <- 1
  
  # Nested loop for geom_signif
  for (loc1 in 1:(length(unique_comparators) - 1)) {
    for (loc2 in (loc1 + 1):length(unique_comparators)) {
      comparison_loc1 <- unique_comparators[loc1]
      comparison_loc2 <- unique_comparators[loc2]
      
      # Define aes mapping for geom_signif
      aes_mapping <- aes(x = Condition, y = {{y_var}}, color = factor(Date))
      
      geom_signif_call <- substitute(
        geom_signif(data = ReplicateAverages, comparisons = list(c(loc1, loc2)),
                    map_signif_level = TRUE, test = "wilcox.test",
                    color = "black", y_position = y_pos, na.rm = TRUE),
        list(loc1 = comparison_loc1,
             loc2 = comparison_loc2,
             y_pos = y_pos_values[CompNum])
      )
      # Add geom_signif with labels
      Plot <- Plot + eval(geom_signif_call)
      # Add a blank annotation to preserve the original labels
      Plot <- Plot + annotate("text", x = 1, y = 1, label = "")
      CompNum <- CompNum + 1
    }
    
  }
  #Plot <- Plot + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  print(Plot)
}

Superplot_beeswarm <- function(Data, y_var, Labels){
  ReplicateAverages <- Data %>% group_by(Condition, Date) %>% summarise_each(list(mean))
  
  ggplot(Data, aes(x=Condition,y={{y_var}},color=factor(Date))) + 
    geom_beeswarm(alpha = 0.2, cex=0.5, dodge.width = 0.5, size = 1) + scale_colour_brewer(palette = "Set2") + 
    stat_summary(fun="mean", geom="crossbar", width=0.5, color = "black") + 
    geom_beeswarm(data=ReplicateAverages, size=5, alpha =0.9) +
    geom_beeswarm(data=ReplicateAverages, size=5, color="black", shape = 1) +
    #stat_compare_means(data=ReplicateAverages, comparisons = list(c("ESS1.5 polyK", "ESS15 polyK")), method="t.test", paired=TRUE) + 
    theme(legend.position="none") + theme_jmt() +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]])
}

Superplot_CellMeanStDev <- function(Data, y_var, Labels){
  # Subset the data based on the y_var
  CellMeanStDev <- Data %>%
    select(Condition, Date, Name, {{y_var}})
  
  CellMeanStDev <- CellMeanStDev %>%
    group_by(Name, Condition) %>%
    summarize(
      Date = first(Date),
      mean_var = mean({{y_var}}),
      stdev_var = sd({{y_var}}),
      sem_var = stdev_var / sqrt(n())
    )
  ReplicateAverages <- CellMeanStDev %>% group_by(Condition, Date) %>% summarise_each(list(mean))
  
  ggplot(CellMeanStDev, aes(x = Condition, y = mean_var)) +
    geom_point(aes(color = Date), alpha = 0.5, position = position_dodge2(width = 0.5, padding = 0.1), size = 3) +
    geom_errorbar(aes(ymin = mean_var - stdev_var, ymax = mean_var + stdev_var, color = Date), alpha = 0.5, width = 0.5, position = position_dodge2(width = 0.5, padding = 0.1)) +
    stat_summary(fun="mean", geom="crossbar", width=0.5, color = "black") + 
    geom_point(data=ReplicateAverages, aes(color = factor(Date)), position = position_dodge2(width = 0.25, padding = 0.1), alpha = 0.9, size=5) +
    geom_point(data=ReplicateAverages, size=5, color="black", position = position_dodge2(width = 0.25, padding = 0.1), shape = 1) +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) +
    theme_jmt() + theme(legend.position="none") + scale_colour_brewer(palette = "Set2")
}

Superplot_CellMeanStDev_NoStat <- function(Data, y_var, Labels, Yaxis_Max, Yaxis_Min, Yaxis_step){
  # Subset the data based on the y_var
  CellMeanStDev <- Data %>%
    select(Condition, Date, Name, {{y_var}})
  
  CellMeanStDev <- CellMeanStDev %>%
    group_by(Name, Condition) %>%
    summarize(
      Date = first(Date),
      mean_var = mean({{y_var}}),
      stdev_var = sd({{y_var}}),
      sem_var = stdev_var / sqrt(n())
    )
  
  ReplicateAverages <- CellMeanStDev %>% group_by(Condition, Date) %>% summarise_each(list(mean))
  
  Plot <- ggplot(CellMeanStDev, aes(x = Condition, y = mean_var)) +
    geom_point(aes(color = Date), alpha = 0.5, position = position_dodge2(width = 0.5, padding = 0.1), size = 3) +
    geom_errorbar(aes(ymin = mean_var - stdev_var, ymax = mean_var + stdev_var, color = Date), alpha = 0.5, width = 0.5, position = position_dodge2(width = 0.5, padding = 0.1)) +
    stat_summary(fun="mean", geom="crossbar", width=0.5, color = "black") + 
    geom_point(data=ReplicateAverages, aes(color = factor(Date)), position = position_dodge2(width = 0.25, padding = 0.1), alpha = 0.9, size=5) +
    geom_point(data=ReplicateAverages, size=5, color="black", position = position_dodge2(width = 0.25, padding = 0.1), shape = 1) +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) +
    theme_jmt() + theme(legend.position="none") + scale_colour_brewer(palette = "Set2") +
    scale_y_continuous(limits = c(Yaxis_Min, Yaxis_Max), breaks = seq(Yaxis_Min, Yaxis_Max, Yaxis_step))
  
  print(Plot)
}

Superplot_CellMeanStDev_Stat <- function(Data, y_var, Labels, Yaxis_Max, Yaxis_Min, Yaxis_step){
  # Subset the data based on the y_var
  CellMeanStDev <- Data %>%
    select(Condition, Date, Name, {{y_var}})
  
  CellMeanStDev <- CellMeanStDev %>%
    group_by(Name, Condition) %>%
    summarize(
      Date = first(Date),
      mean_var = mean({{y_var}}),
      stdev_var = sd({{y_var}}),
      sem_var = stdev_var / sqrt(n())
    )
  
  ReplicateAverages <- CellMeanStDev %>% group_by(Condition, Date) %>% summarise_each(list(mean))
  
  Plot <- ggplot(CellMeanStDev, aes(x = Condition, y = mean_var)) +
    geom_point(aes(color = Date), alpha = 0.5, position = position_dodge2(width = 0.5, padding = 0.1), size = 3) +
    geom_errorbar(aes(ymin = mean_var - stdev_var, ymax = mean_var + stdev_var, color = Date), alpha = 0.5, width = 0.5, position = position_dodge2(width = 0.5, padding = 0.1)) +
    stat_summary(fun="mean", geom="crossbar", width=0.5, color = "black") + 
    geom_point(data=ReplicateAverages, aes(color = factor(Date)), position = position_dodge2(width = 0.25, padding = 0.1), alpha = 0.9, size=5) +
    geom_point(data=ReplicateAverages, size=5, color="black", position = position_dodge2(width = 0.25, padding = 0.1), shape = 1) +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) +
    theme_jmt() + theme(legend.position="none") + scale_colour_brewer(palette = "Set2") +
    scale_y_continuous(limits = c(Yaxis_Min, Yaxis_Max), breaks = seq(Yaxis_Min, Yaxis_Max, Yaxis_step))
  
  unique_comparators <- unique(ReplicateAverages$Condition)
  
  NumberOfComparisons <- (0.5 * length(unique_comparators) - 0.5) * length(unique_comparators)
  Yaxis_length <- Yaxis_Max - Yaxis_Min
  TopY <- Yaxis_Min + (Yaxis_length * 0.9)
  y_pos_values <- rev(seq((TopY - (Yaxis_length * 0.1 * (NumberOfComparisons - 1))), TopY, by = (Yaxis_step * 0.5)))
  CompNum <- 1
  
  # Nested loop for geom_signif
  for (loc1 in 1:(length(unique_comparators) - 1)) {
    for (loc2 in (loc1 + 1):length(unique_comparators)) {
      comparison_loc1 <- unique_comparators[loc1]
      comparison_loc2 <- unique_comparators[loc2]
      
      # Define aes mapping for geom_signif
      aes_mapping <- aes(x = Condition, y = {{y_var}}, color = factor(Date))
      
      geom_signif_call <- substitute(
        geom_signif(data = ReplicateAverages, comparisons = list(c(loc1, loc2)),
                    map_signif_level = TRUE, test = "t.test",
                    color = "black", y_position = y_pos, na.rm = TRUE),
        list(loc1 = comparison_loc1,
             loc2 = comparison_loc2,
             y_pos = y_pos_values[CompNum])
      )
      # Add geom_signif with labels
      Plot <- Plot + eval(geom_signif_call)
      # Add a blank annotation to preserve the original labels
      Plot <- Plot + annotate("text", x = 1, y = 1, label = "")
      CompNum <- CompNum + 1
    }
  }
  print(Plot)
}

#### end ####

Labs <- c(title = "Focal Adhesion areas", x = "", y = expression(paste("Focal Adhesion area (", mu, "m"^2, ")")))
CellFAareaSuperplot <- Superplot_CellMeanStDev(AllFAData, area, Labs)
ggsave(paste("CellMeanStDev of AllFocalAdhesions.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_Stat(AllFAData, area, Labs, 8, 0, 1)
ggsave(paste("CellMeanStDev_Stat of AllFocalAdhesions.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Area of all Focal Adhesions", x = "", y = expression(paste("Focal adhesion area (", mu, "m"^2, ")")))
AllFaDataPlot <- Superplot_beeswarm(AllFAData, area, Labs)
Superplot_statistics(AllFaDataPlot, AllFAData, area, 8, 0, 1)
ggsave(paste("Beeswarm of AllFocalAdhesions.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "All Focal Adhesions", x = "", y = expression(paste("Focal adhesion area (", mu, "m"^2, ")")))
Crossbar(AllFAData, area, Labs) + theme(legend.position = "none") +
  geom_signif(comparisons = list(c(unique(AllFAData$Condition)[1],unique(AllFAData$Condition)[2])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 8) +
  geom_signif(comparisons = list(c(unique(AllFAData$Condition)[2],unique(AllFAData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 9) +
  geom_signif(comparisons = list(c(unique(AllFAData$Condition)[1],unique(AllFAData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 10) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0,11))
ggsave(paste("Crossbar of AllFocalAdhesions.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "FA-ER area intersection (%)")
Crossbar(AllFAData, ERThresh_intersection_percentage, Labs) + theme(legend.position = "none") +
  geom_signif(comparisons = list(c(unique(AllFAData$Condition)[1],unique(AllFAData$Condition)[2])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 100) +
  geom_signif(comparisons = list(c(unique(AllFAData$Condition)[2],unique(AllFAData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 107) +
  geom_signif(comparisons = list(c(unique(AllFAData$Condition)[1],unique(AllFAData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 114) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0,120))
ggsave(paste("Crossbar of ERThresh intersection percentage.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Superplot_CellMeanStDev_Stat(AllFAData, ERThresh_intersection_percentage, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev_Stat of FA-ER area intersection for AllFocalAdhesions.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Overlap of ER \nwith FAs", x = "", y = "FA-ER area intersection (%)")
Crossbar(MeanCellData, ERThresh_intersection_percentage, Labs) + theme(legend.position = "none") +
  geom_signif(comparisons = list(c(unique(AllFAData$Condition)[1],unique(AllFAData$Condition)[2])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 100) +
  geom_signif(comparisons = list(c(unique(AllFAData$Condition)[2],unique(AllFAData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 107) +
  geom_signif(comparisons = list(c(unique(AllFAData$Condition)[1],unique(AllFAData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 114) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0,120))
ggsave(paste("Crossbar of ERThresh intersection percentage for MeanCellData.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Superplot_CellMeanStDev_Stat(MeanCellData, ERThresh_intersection_percentage, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev_Stat of FA-ER area intersection for MeanCellData.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Focal adhesion eccentricity", x = "", y = "FA eccentricity")
Crossbar(AllFAData, eccentricity, Labs) + theme(legend.position = "none") +
  geom_signif(comparisons = list(c(unique(AllFAData$Condition)[1],unique(AllFAData$Condition)[2])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 1.00) +
  geom_signif(comparisons = list(c(unique(AllFAData$Condition)[2],unique(AllFAData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 1.07) +
  geom_signif(comparisons = list(c(unique(AllFAData$Condition)[1],unique(AllFAData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 1.14) +
  scale_y_continuous(breaks = seq(0, 1.00, by = 0.2), limits = c(0,1.20))
ggsave(paste("Focal adhesion eccentricity.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Mean Focal Adhesion \narea in each cell", x = "", y = expression(paste("Focal adhesion area (", mu, "m"^2, ")")))
Crossbar(MeanCellData, area, Labs) + theme(legend.position = "none") +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[1],unique(MeanCellData$Condition)[2])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 1.5) +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[2],unique(MeanCellData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 2) +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[1],unique(MeanCellData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 2.5) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), limits = c(0,3))
ggsave(paste("Crossbar of Mean Focal Adhesion per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Mean Focal Adhesion eccentricity in each cell", x = "", y = "FA eccentricity")
Crossbar(MeanCellData, eccentricity, Labs) + theme(legend.position = "none") +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[1],unique(MeanCellData$Condition)[2])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 1.00) +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[2],unique(MeanCellData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 1.07) +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[1],unique(MeanCellData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 1.14) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1.2))
ggsave(paste("Crossbar of Mean Focal Adhesion eccentricity per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Cell areas", x = "", y = expression(paste("Cell area (", mu, "m"^2, ")")))
Crossbar(MeanCellData, Cell.area, Labs) + theme(legend.position = "none") +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[1],unique(MeanCellData$Condition)[2])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 3500) +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[2],unique(MeanCellData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 3250) +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[1],unique(MeanCellData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 3750) +
  scale_y_continuous(breaks = seq(0, 4000, by = 500), limits = c(0,4000))
ggsave(paste("Cell areas.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Superplot_CellMeanStDev_Stat(MeanCellData, Cell.area, Labs, 4000, 0, 500)
ggsave(paste("CellMeanStDev_Stat of Crossbar of Cell areas.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Total Focal Adhesion \narea in each cell", x = "", y = expression(paste("Total FA area (", mu, "m"^2, ")")))
Crossbar(MeanCellData, TotalFAarea, Labs) + theme(legend.position = "none") +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[1],unique(MeanCellData$Condition)[2])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 15) +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[2],unique(MeanCellData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position =12) +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[1],unique(MeanCellData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 18) +
  scale_y_continuous(breaks = seq(0, 25, by = 5), limits = c(0,25))
ggsave(paste("Crossbar of Total Focal Adhesion per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Superplot_CellMeanStDev_Stat(MeanCellData, TotalFAarea, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev_Stat of Crossbar of Total Focal Adhesion per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Number of Focal \nAdhesions", x = "", y = expression(paste("Number of FAs per cell")))
Crossbar(MeanCellData, NumberOfAdhesions, Labs) + theme(legend.position = "none") +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[1],unique(MeanCellData$Condition)[2])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 75) +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[2],unique(MeanCellData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 85) +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[1],unique(MeanCellData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 95) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0,100))
ggsave(paste("Number of Focal Adhesion per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Superplot_CellMeanStDev_Stat(MeanCellData, NumberOfAdhesions, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev_Stat of Number of Focal Adhesion per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Proportion of FAs \nintersecting with ER", x = "", y = "FAs occupied by ER (%)")
Crossbar(MeanCellData, PercERFAintersection, Labs) + theme(legend.position = "none") +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[1],unique(MeanCellData$Condition)[2])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 100) +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[2],unique(MeanCellData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 107) +
  geom_signif(comparisons = list(c(unique(MeanCellData$Condition)[1],unique(MeanCellData$Condition)[3])), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 114) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0,120))
ggsave(paste("Proportion of FAs intersecting with ER.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")



Labs <- c(title = "Overlap of ER \nwith FAs", x = "", y = "FA-ER area \nintersection (%)")
Superplot_CellMeanStDev_Stat(AllFAData, ERThresh_intersection_percentage, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev_Stat of FA-ER area intersection for AllFocalAdhesions.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(AllFAData, ERThresh_intersection_percentage, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev of FA-ER area intersection for AllFocalAdhesions.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Cell areas", x = "", y = expression(paste("Cell area (", mu, "m"^2, ")")))
Superplot_CellMeanStDev_Stat(MeanCellData, Cell.area, Labs, 4000, 0, 500)
ggsave(paste("CellMeanStDev_Stat of Crossbar of Cell areas.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, Cell.area, Labs, 4000, 0, 500)
ggsave(paste("CellMeanStDev of Crossbar of Cell areas.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Total Focal Adhesion \narea in each cell", x = "", y = expression(paste("Total FA area (", mu, "m"^2, ")")))
Superplot_CellMeanStDev_Stat(MeanCellData, TotalFAarea, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev_Stat of Crossbar of Total Focal Adhesion per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, TotalFAarea, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev of Crossbar of Total Focal Adhesion per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Proportion of FAs \nintersecting with ER", x = "", y = "FAs occupied by ER (%)")
Superplot_CellMeanStDev_Stat(MeanCellData, PercERFAintersection, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev_Stat of Proportion of FAs intersecting with ER.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, PercERFAintersection, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev of Proportion of FAs intersecting with ER.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Number of FAs \nper 100 µm^2", x = "", y = "FA number (per 100 µm^2)")
Superplot_CellMeanStDev_Stat(MeanCellData, FAnum_per_100um, Labs, 10, 0, 2)
ggsave(paste("CellMeanStDev_Stat of FA number per 100um.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, FAnum_per_100um, Labs, 10, 0, 2)
ggsave(paste("CellMeanStDev of FA number per 100um.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Area of FAs \nper 100 µm^2", x = "", y = "FA area (µm^2)")
Superplot_CellMeanStDev_Stat(MeanCellData, TotalFAarea_prop_100um, Labs, 10, 0, 2)
ggsave(paste("CellMeanStDev_Stat of FA area per 100um.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, TotalFAarea_prop_100um, Labs, 10, 0, 2)
ggsave(paste("CellMeanStDev of FA area per 100um.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Number of FAs \nintersecting ER", x = "", y = "FA number")
Superplot_CellMeanStDev_Stat(MeanCellData, NumFAintersectER, Labs, 50, 0, 10)
ggsave(paste("CellMeanStDev_Stat of FA number intersecting ER.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, NumFAintersectER, Labs, 50, 0, 10)
ggsave(paste("CellMeanStDev of FA number intersecting ER.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Number of FAs \nper cell", x = "", y = "FA number")
Superplot_CellMeanStDev_Stat(MeanCellData, NumberOfAdhesions, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev_Stat of FA number per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, NumberOfAdhesions, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev of FA number per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Area of FAs", x = "", y = "FA area (µm^2)")
Superplot_CellMeanStDev_Stat(MeanCellData, area, Labs, 2, 0, 0.4)
ggsave(paste("CellMeanStDev_Stat of FA area.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, area, Labs, 2, 0, 0.4)
ggsave(paste("CellMeanStDev of FA area.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")