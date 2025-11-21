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
install.packages("gridExtra")
install.packages("grid")
install.packages("lmerTest")

library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(gridExtra)
library(grid)
library(lmerTest)
#### end ####

# Set location of data and where plots will be saved to:
mainDir <- "C:/Users/jonatmt/OneDrive - Universitetet i Oslo/Desktop/Amalie masks"
#mainDir <- "/Users/amalie/Done"
setwd(mainDir)
# Create the path where the data is stored
DataPath <- "C:/Users/jonatmt/OneDrive - Universitetet i Oslo/Desktop/Amalie masks/csv files"
#DataPath <- "/Users/amalie/Done/csv files"

# Set the pixel size and focal adhesion sizes
pixel_length_x <- 0.10300000916700082
pixel_length_y <-0.10300000916700082
pixel_area <- pixel_length_x*pixel_length_y
MinFAsize <- 0.15
MaxFAsize <- 15

# Set the order of conditions along the x axis for plots
x_axis_order = c("Scr-control", "control", "WT-BFP", "L24Q-BFP")

#### Read in the data and keep focal adhesion data ####
MeanPixelLength <- (pixel_length_x+pixel_length_y)/2
Cell_Areas <- read.csv(paste(mainDir,"Cell_Area_results.csv", sep ="/"))
#Cell_Areas$Name <- sub("^[^-]*-", "", Cell_Areas$Name)
Cell_Areas$Filename <- sub("\\.tif$", "_NoERsheets.csv", Cell_Areas$Filename)
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
  
  # Get the matching name (remove extension from filename)
  #file_base_name <- tools::file_path_sans_ext(Files[i])
  
  # Find the matching Cell_Area value
  cell_area_value <- Cell_Areas$Area[Cell_Areas$Filename == Files[i]]
  
  # If found, add it as a new column to the df
  if (length(cell_area_value) == 1) {
    df$Cell_Area <- cell_area_value
  } else {
    df$Cell_Area <- NA  # Or handle mismatches however you prefer
    warning(paste("No matching Cell_Area found for:", Files[i]))
  }
  # Add focal adhesion size proportional to the cell area
  df$FApropCA <- df$area/df$Cell_Area
  
  # Assign data to the raw data list
  AllRawData[[i]] <- df
  # Name each dataframe in the list after the csv file it came from
  names(AllRawData)[i] <- Files[i]
  rm(df)
}
rm(FilePaths, Files, i, cell_area_value, Cell_Areas)

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
  subset_df <- mutate(subset_df, Condition = strsplit(names(AllRawData)[s], "_")[[1]][4])
  # Add a column containing the date of the experiment
  subset_df <- mutate(subset_df, Date = strsplit(names(AllRawData)[s], "_")[[1]][1])
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
AllFAData$Condition <- factor(AllFAData$Condition, levels = x_axis_order)

MeanCellData <- data.frame()
# Calculate column means for each dataframe
for (s in 1:length(AllRawData)) {
  df <- AllRawData[[s]]
  # Keep only FAs within size limit
  subset_df <- df[df$area > MinFAsize & df$area < MaxFAsize, ]
  # Count how many FAs are in the cell
  subset_df$NumberOfAdhesions <- nrow(subset_df)
  # Number of focal adhesions per 100 um
  subset_df$FAnum_per_100um <- (subset_df$NumberOfAdhesions / subset_df$Cell_Area) * 100
  # Number of FAs which also intersect with ER
  subset_df$NumFAintersectER <- nrow(subset_df[subset_df$ERThresh_num_intersection_areas >= 1, ])
  # Number of FAs which also intersect with ER per 100 um
  subset_df$NumFAintersectER_per_100um <- (subset_df$NumFAintersectER / subset_df$Cell_Area) * 100
  # Calculate the sum area of focal adhesions
  subset_df$TotalFAarea <- sum(subset_df$area)
  # Total FA area per 100 um^2
  subset_df$TotalFAarea_prop_100um <- (subset_df$TotalFAarea / subset_df$Cell_Area) * 100
  # Convert subset_df to a data frame
  subset_df <- as.data.frame(subset_df)
  # Calculate column means using summarize
  mean_df <- subset_df %>% summarize(across(everything(), \(x) mean(x, na.rm = TRUE)))
  mean_df$PercERFAintersection <- (mean_df$NumFAintersectER/mean_df$NumberOfAdhesions)*100
  mean_df <- mutate(mean_df, Name = names(AllRawData)[s])
  mean_df <- mutate(mean_df, Condition = strsplit(names(AllRawData)[s], "_")[[1]][4])
  mean_df <- mutate(mean_df, Date = strsplit(names(AllRawData)[s], "_")[[1]][1])
  MeanCellData <- bind_rows(MeanCellData, mean_df)
  rm(subset_df, mean_df, df)
}
MeanCellData$Condition <- factor(MeanCellData$Condition, levels = x_axis_order)
rm(s)

summary_df <- MeanCellData %>%
  group_by(Condition, Date) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
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
    #geom_errorbar(aes(ymin = mean_var - stdev_var, ymax = mean_var + stdev_var, color = Date), alpha = 0.5, width = 0.5, position = position_dodge2(width = 0.5, padding = 0.1)) +
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
    #geom_errorbar(aes(ymin = mean_var - stdev_var, ymax = mean_var + stdev_var, color = Date), alpha = 0.5, width = 0.5, position = position_dodge2(width = 0.5, padding = 0.1)) +
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


Superplot_CellMeanStDev_ANOVA <- function(Data, y_var, Labels, Yaxis_Max, Yaxis_Min, Yaxis_step){
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
    #geom_errorbar(aes(ymin = mean_var - stdev_var, ymax = mean_var + stdev_var, color = Date), alpha = 0.5, width = 0.5, position = position_dodge2(width = 0.5, padding = 0.1)) +
    stat_summary(fun="mean", geom="crossbar", width=0.5, color = "black") + 
    geom_point(data=ReplicateAverages, aes(color = factor(Date)), position = position_dodge2(width = 0.25, padding = 0.1), alpha = 0.9, size=5) +
    geom_point(data=ReplicateAverages, size=5, color="black", position = position_dodge2(width = 0.25, padding = 0.1), shape = 1) +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) +
    theme_jmt() + theme(legend.position="none") + scale_colour_brewer(palette = "Set2") +
    scale_y_continuous(limits = c(Yaxis_Min, Yaxis_Max), breaks = seq(Yaxis_Min, Yaxis_Max, Yaxis_step))
  
  # Run one-way ANOVA
  fit <- aov(mean_var ~ Condition, data = ReplicateAverages)
  print(summary(fit))   # ANOVA table
  print(shapiro.test(residuals(fit)))
  print(leveneTest(mean_var ~ Condition, data = ReplicateAverages))
  
  # Tukey’s HSD posthoc
  tukey_res <- TukeyHSD(fit)
  print(tukey_res)
  
  print(Plot)
}
#### end ####


# Define the dates you want to keep
keep_values <- c("240321", "240322", "240419", "240426", "250123")  # Replace these with your actual date values
# Subset the MeanCellData to keep only the desired dates
select_dates_df <- subset(MeanCellData, Date %in% keep_values)
# Check if the data was correctly filtered (optional)
#print(select_dates_df)


## Compare the areas of the cells
Labs <- c(title = "Mean cell area", x = "", y = expression(paste("Cell area (", mu, "m"^2, ")")))
Superplot_CellMeanStDev_Stat(MeanCellData, Cell_Area, Labs, 4000, 0, 500)
ggsave(paste("CellMeanStDev_Stat of Crossbar of Cell areas.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_ANOVA(MeanCellData, Cell_Area, Labs, 4000, 0, 500)
ggsave(paste("CellMeanStDev of Crossbar of Cell areas.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

## Compare the numbers and areas of the focal adhesions
Labs <- c(title = "Mean number of FAs \nper cell", x = "", y = "FA number")
Superplot_CellMeanStDev_Stat(MeanCellData, NumberOfAdhesions, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev_Stat of FA number per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, NumberOfAdhesions, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev of FA number per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Mean number of FAs \nper 100 µm^2", x = "", y = "FA number (per 100 µm^2)")
Superplot_CellMeanStDev_Stat(MeanCellData, FAnum_per_100um, Labs, 10, 0, 2)
ggsave(paste("CellMeanStDev_Stat of FA number per 100um.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, FAnum_per_100um, Labs, 10, 0, 2)
ggsave(paste("CellMeanStDev of FA number per 100um.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Mean Area of FAs", x = "", y = "FA area (µm^2)")
Superplot_CellMeanStDev_Stat(MeanCellData, area, Labs, 2, 0, 0.4)
ggsave(paste("CellMeanStDev_Stat of FA area.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, area, Labs, 2, 0, 0.4)
ggsave(paste("CellMeanStDev of FA area.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Mean Area of FAs \nper 100 µm^2", x = "", y = "FA area (µm^2)")
Superplot_CellMeanStDev_Stat(MeanCellData, TotalFAarea_prop_100um, Labs, 10, 0, 2)
ggsave(paste("CellMeanStDev_Stat of FA area per 100um.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, TotalFAarea_prop_100um, Labs, 10, 0, 2)
ggsave(paste("CellMeanStDev of FA area per 100um.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

Labs <- c(title = "Total Focal Adhesion \narea in each cell", x = "", y = expression(paste("Total FA area (", mu, "m"^2, ")")))
Superplot_CellMeanStDev_Stat(MeanCellData, TotalFAarea, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev_Stat of Crossbar of Total Focal Adhesion per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, TotalFAarea, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev of Crossbar of Total Focal Adhesion per cell.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")



## Number of FAs that have any kind of contact with the ER
Labs <- c(title = "Mean number of FAs \nintersecting ER", x = "", y = "FA number")
Superplot_CellMeanStDev_Stat(MeanCellData, NumFAintersectER, Labs, 50, 0, 10)
ggsave(paste("CellMeanStDev_Stat of FA number intersecting ER.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, NumFAintersectER, Labs, 50, 0, 10)
ggsave(paste("CellMeanStDev of FA number intersecting ER.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

## The same as previous plot, but prortional to the total number of FAs in the cell
Labs <- c(title = "Proportion of FAs \nintersecting with ER", x = "", y = "FAs occupied by ER (%)")
Superplot_CellMeanStDev_Stat(MeanCellData, PercERFAintersection, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev_Stat of Proportion of FAs intersecting with ER.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(MeanCellData, PercERFAintersection, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev of Proportion of FAs intersecting with ER.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

## The area of FAs which also have ER, as a percentage of the area of the FA
Labs <- c(title = "Area Overlap of ER \nwith FAs", x = "", y = "FA-ER area \nintersection (%)")
Superplot_CellMeanStDev_Stat(AllFAData, ERThresh_intersection_percentage, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev_Stat of FA-ER area intersection for AllFocalAdhesions.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")
Superplot_CellMeanStDev_NoStat(AllFAData, ERThresh_intersection_percentage, Labs, 100, 0, 20)
ggsave(paste("CellMeanStDev of FA-ER area intersection for AllFocalAdhesions.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 15, units = "cm")

#### Statistics Tests and Plots ####
subDir <- "Stats"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
rm(subDir)

## The first thing we look for is normality of the data. Note that anything with
## a sample size greater than 30 can be considered normal because of the central
## limit theorem.

## We can do density plots of the band intensities to make sure the data is
## normally distributed

## A quantile-quantile plot can also be used to check for normal distribution
## the data should align with the 45 degree line if it is normal.

Summary_plots_by_condition_date <- function(df, MainTitlePrefix) {
  # Ensure Date is treated as a factor
  df$Date <- as.factor(df$Date)
  
  # Identify variables to plot
  variableCols <- names(df)[sapply(df, is.numeric) & !(names(df) %in% c("Date", "Condition", "Name"))]
  
  # Unique values for layout
  Dates <- unique(df$Date)
  Conditions <- unique(df$Condition)
  n_rows <- length(Dates)
  n_cols <- length(Conditions)
  
  # Plot types and file suffixes
  plot_types <- c("QQ", "Density", "Boxplot")
  
  for (plot_type in plot_types) {
    pdf_filename <- paste0(MainTitlePrefix, "_", plot_type, ".pdf")
    pdf(pdf_filename, width = n_cols * 4, height = n_rows * 3)
    
    for (var in variableCols) {
      Plots <- list()
      
      for (d in seq_along(Dates)) {
        for (c in seq_along(Conditions)) {
          subset_df <- df %>%
            filter(Date == Dates[d], Condition == Conditions[c])
          
          plot_title <- paste("Date:", Dates[d], "| Condition:", Conditions[c])
          
          if (nrow(subset_df) >= 3) {
            if (plot_type == "QQ") {
              p <- ggqqplot(subset_df, x = var) +
                ggtitle(plot_title) +
                theme_minimal()
            } else if (plot_type == "Density") {
              p <- ggplot(subset_df, aes_string(x = var)) +
                geom_density(fill = "skyblue", alpha = 0.5) +
                theme_minimal() +
                ggtitle(plot_title)
            } else if (plot_type == "Boxplot") {
              p <- ggplot(subset_df, aes_string(x = "Condition", y = var)) +
                geom_boxplot(fill = "tomato", alpha = 0.6) +
                theme_minimal() +
                ggtitle(plot_title)
            }
          } else {
            p <- ggplot() + 
              theme_void() + 
              ggtitle(paste(plot_title, "\n(Insufficient Data)"))
          }
          
          Plots[[length(Plots)+1]] <- p
        }
      }
      
      # Combine and title each page
      grid_title <- textGrob(var, gp = gpar(fontsize = 16, fontface = "bold"))
      arranged <- arrangeGrob(grobs = Plots, ncol = n_cols, top = grid_title)
      grid.newpage() 
      grid.draw(arranged)
    }
    
    dev.off()
  }
}
Summary_plots_by_condition_date(MeanCellData, "SummaryPlots")

## Perform the Shapiro-Wilk test of normality on samples. In this test the null
## hypothesis is that the data is normal and therefore if the p value <0.05 the 
## data should be considered to NOT be normal. The test works best on small sample
## sizes and should be considered alongside the graphs.

# This function will create a csv file with some summary statistics for a list
# including the Shapiro value. It also creates a list of these summaries for
# viewing in R.

DataStatSumTest <- function(DataFrame, Filename) {
  # Split the dataframe by unique Date values
  DataList <- split(DataFrame, DataFrame$Date)
  
  # Set columns to summarise (excluding metadata)
  Columns <- setdiff(names(DataFrame), c("Name", "Condition", "Date"))
  
  OutputList <- list()
  
  for (p in 1:length(DataList)) {
    df <- DataList[[p]]
    dfStatistics <- data.frame()
    
    for (col in Columns) {
      for (g in unique(df$Condition)) {
        Testdf <- dplyr::filter(df, Condition == g)
        values <- as.numeric(Testdf[[col]])
        
        MeanValue <- tryCatch(mean(values), error = function(e) NaN)
        MedianValue <- tryCatch(median(values), error = function(e) NaN)
        VarValue <- tryCatch(var(values), error = function(e) NaN)
        StandardDeviation <- tryCatch(sd(values), error = function(e) NaN)
        StandardError <- tryCatch(StandardDeviation / sqrt(length(values)), error = function(e) NaN)
        SWpvalue <- tryCatch(shapiro.test(values)$p.value, error = function(e) NaN)
        # A One sample t test compares the mean to the mu value. Probably not needed.
        #OneSamplep <- tryCatch(t.test(values, mu = 1)$p.value, error = function(e) NaN)
        
        StatsResults <- data.frame(
          Condition = g,
          IntensityColumn = col,
          Mean = MeanValue,
          Median = MedianValue,
          Variance = VarValue,
          Standard_Deviation = StandardDeviation,
          Standard_Error = StandardError,
          Shapiro_Wilk_p_value = SWpvalue
          #`One sample T-test p value` = OneSamplep
        )
        
        dfStatistics <- rbind(dfStatistics, StatsResults)
      }
    }
    
    OutputList[[p]] <- dfStatistics
    names(OutputList)[p] <- names(DataList)[p]
  }
  
  SaveResult <- do.call(rbind, OutputList)
  write.csv(SaveResult, file = paste(Filename, ".csv", sep = ""), row.names = FALSE)
  return(OutputList)
}

SumMeanCellData <- DataStatSumTest(MeanCellData,"Summary of MeanCellData")

## Finally we want to generate a set of p-values for the data which can be checked for
## significance and compared to the summaries generated above to know which is most
## relevant. For normal data with equal variance it is the t test, unequal variance use
## the Welch's t test. If the data is nonparametric (a.k.a. not normal), then it is
## the Mann-Whitney test.

# Effect size functions
cohen_d <- function(x, y) {
  pooled_sd <- sqrt(((length(x) - 1) * sd(x)^2 + (length(y) - 1) * sd(y)^2) / (length(x) + length(y) - 2))
  return(abs(mean(x) - mean(y)) / pooled_sd)
}

cliffs_delta <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  rank_diff <- sum(outer(x, y, FUN = ">")) - sum(outer(x, y, FUN = "<"))
  return(rank_diff / (n1 * n2))
}

ConditionStatsTests <- function(DataFrame, Filename) {
  Columns <- setdiff(names(DataFrame), c("Name", "Condition", "Date"))
  dfStatistics <- data.frame()
  
  AllConditions <- unique(DataFrame$Condition)
  if (length(AllConditions) < 2) {
    message("Not enough conditions to compare.")
    return(NULL)
  }
  
  conditionPairs <- combn(AllConditions, 2, simplify = FALSE)
  
  for (Col in Columns) {
    # Step 1: Summarise per-date means
    Summarised <- DataFrame |>
      dplyr::group_by(Date, Condition) |>
      dplyr::summarise(MeanValue = mean(.data[[Col]], na.rm = TRUE), .groups = "drop")
    
    for (pair in conditionPairs) {
      cond1 <- pair[1]
      cond2 <- pair[2]
      
      # Filter only the two conditions
      SubData <- Summarised[Summarised$Condition %in% c(cond1, cond2), ]
      
      # Extract per-date means
      values1 <- SubData$MeanValue[SubData$Condition == cond1]
      values2 <- SubData$MeanValue[SubData$Condition == cond2]
      
      n1 <- length(values1)
      n2 <- length(values2)
      nDates <- length(unique(SubData$Date))
      
      if (n1 < 2 | n2 < 2 | nDates < 2) next  # Need at least 2 values per group and at least 2 dates
      
      # Classic tests on per-date means
      StudentT <- tryCatch(t.test(values1, values2, var.equal = TRUE)$p.value, error = function(e) NA)
      WelchT <- tryCatch(t.test(values1, values2, var.equal = FALSE)$p.value, error = function(e) NA)
      MannWhitney <- tryCatch(wilcox.test(values1, values2)$p.value, error = function(e) NA)
      cohen_d_val <- tryCatch(cohen_d(values1, values2), error = function(e) NA)
      cliffs_delta_val <- tryCatch(cliffs_delta(values1, values2), error = function(e) NA)
      
      # Mixed effects model: Use original (non-averaged) data for this
      LMM_Data <- DataFrame[DataFrame$Condition %in% c(cond1, cond2), ]
      LMM <- tryCatch(
        lmerTest::lmer(as.formula(paste(Col, "~ Condition + (1 | Date)")), data = LMM_Data),
        error = function(e) NA
      )
      LMM_p <- tryCatch(coef(summary(LMM))[2, "Pr(>|t|)"], error = function(e) NA)
      
      # Store results
      SummaryVec <- data.frame(
        ValueColumn = Col,
        Condition1 = cond1,
        Condition2 = cond2,
        SampleSize1 = n1,
        SampleSize2 = n2,
        NumDates = nDates,
        StudentT = StudentT,
        WelchT = WelchT,
        MannWhitney = MannWhitney,
        LMM_p = LMM_p,
        Cohen_d = cohen_d_val,
        Cliff_delta = cliffs_delta_val
      )
      
      dfStatistics <- rbind(dfStatistics, SummaryVec)
    }
  }
  
  write.csv(dfStatistics, file = paste0(Filename, ".csv"), row.names = FALSE)
  return(dfStatistics)
}

StatsMeanCellData <- ConditionStatsTests(MeanCellData,"SigStats of MeanCellData")

setwd(mainDir)
#### end ####