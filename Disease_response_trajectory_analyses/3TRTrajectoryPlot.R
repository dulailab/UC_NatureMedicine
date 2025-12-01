## ---------------------------
##
## Script name: Trajectory Visualization
## Purpose of script: Generate a visualization of trajectory data using ggplot2.
## Author: Dr. Scott Jelinsky
## Date Created: 2025-01-05
## Version: 1.1 with error checking
## Email: Scott.Jelinsky@Pfizer.com
##
## ---------------------------
##
## Notes:
##   - Applies LOESS smoothing to highlight trends.
##   - Annotates plot with relevant details for each subject.
##   - Saves the final plot as a PDF.
##
## ---------------------------
##
## Function: TrajectoryVisualization
##
## Purpose: Generate a ggplot visualization of Trajectory Data.
## Inputs:
##   - data: A data frame with the following columns:
##       * Subject: Identifier for each subject
##       * Cluster: Cluster group
##       * CDC: Additional data column (not used in visualization)
##       * Time: Time in weeks (numeric)
##       * Value: Measured value
##    - Example Input File 
##        A tibble: 6 × 5
##       Subject  Cluster CDC    Time Value
##       <chr>    <fct>   <chr> <dbl> <dbl>
##     1 33221470 A       DC        0     5
##     2 33221470 A       DC        2     1
##     3 33221470 A       DC        4     1
##     4 33221470 A       DC        8     0
##     5 33221470 A       DC       12     0
##     6 33221470 A       DC       16     0
##   - output_file: File path for saving the plot as a PDF
##
## Outputs:
##   - A ggplot object representing the visualization
##   - Saves the plot to the specified file path
## Parameters:
##   data        - A data frame containing the trajectory data with columns: Subject, Cluster, CDC, Time, Value.
##   output_file - A string specifying the file path to save the plot.
##   Title       - (Optional) A string for the plot title. Default is "Vibrato Study - Brepocitinib, Ritlecitinib".
##   add_label   - (Optional) A boolean indicating whether to add cluster labels to the plot. Default is TRUE.
##   savePDF     - (Optional) A boolean indicating whether to save the plot as a PDF. Default is TRUE.
##
## ---------------------------

TrajectoryVisualization <- function(data, output_file, 
                                    Title = "Vibrato Study - Brepocitinib, Ritlecitinib",
                                    StudyName=StudyName,
                                    add_label = TRUE,
                                    Include_Subtitles=FALSE,
                                    savePDF = TRUE) {
  # Load necessary libraries
  libraries_needed <- c("dplyr", "ggplot2", "ggrepel", "Cairo")
  missing_libs <- libraries_needed[!libraries_needed %in% installed.packages()[,"Package"]]
  if(length(missing_libs)) {
    stop(paste("The following required packages are missing:", paste(missing_libs, collapse = ", ")))
  }
  
  # Load libraries
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(Cairo)
  })
  
  # Error Checking
  
  ## Check if 'data' is a data frame
  if(!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  
  ## Required columns
  required_cols <- c("Subject", "Cluster", "Time", "Value", "Trajectory")
  missing_cols <- setdiff(required_cols, colnames(data))
  if(length(missing_cols) > 0) {
    stop(paste("The following required columns are missing from 'data':", paste(missing_cols, collapse = ", ")))
  }
  
  ## Check data types
  if(!is.numeric(data$Time)) {
    stop("'Time' column must be numeric.")
  }
  if(!is.numeric(data$Value)) {
    stop("'Value' column must be numeric.")
  }
  
  ## Check 'output_file' is a non-empty string
  if(!is.character(output_file) || length(output_file) != 1 || nchar(output_file) == 0) {
    stop("'output_file' must be a non-empty string specifying the file path.")
  }
  
  ## Check 'Title' is a string
  if(!is.character(Title) || length(Title) != 1) {
    stop("'Title' must be a single string.")
  }
  
  ## Check 'add_label' is boolean
  if(!is.logical(add_label) || length(add_label) != 1) {
    stop("'add_label' must be a single boolean value (TRUE or FALSE).")
  }
  
  ## Check 'savePDF' is boolean
  if(!is.logical(savePDF) || length(savePDF) != 1) {
    stop("'savePDF' must be a single boolean value (TRUE or FALSE).")
  }
  
  # Error Handling: Check if the specified "Trajectory" exists and  contains only expected values (A-E)
  valid_values <- c("T1", "T2", "T3", "T4")
  unique_values <- unique(data$Trajectory)
  unexpected_values <- setdiff(unique_values, valid_values)
  
  if(length(unexpected_values) > 0) {
    warning(paste("The following unexpected values were found in the column Trajectory",
                  "and will be set to NA:", paste(unexpected_values, collapse = ", ")))
  }
  
  # Ensure 'Time' is numeric for accurate plotting
  data$Time <- as.numeric(data$Time)
  
  # Define specific colors for all categories
  custom_colors <- c(
    "T2" = "#FF9999",        # Light Red
    "lightyellow" = "#FFFF99",       # Light Yellow
    "T1" = "#FF66CC",       # Hot Pink
    "T4" = "#CC9933",         # Brownish
    "Elderberry" = "#9933FF",   # Purple
    "T3" = "orange",         
    "Lightpurple" = "#9966FF",        # Light Purple
    "mintgreen" = "#99FFCC"      # Mint Green
  )
  
  # Create the initial ggplot object
  plot <- data %>%
    ggplot(aes(
      x = Time,                  # Map 'Time' to the x-axis
      y = Value,                 # Map 'Value' to the y-axis
      #color = Subject,           # Color lines by 'Subject'
      #group = Subject            # Group data points by 'Subject' for plotting
    )) +
    # Add LOESS smoothing for each cluster to identify trends
    geom_smooth(
      aes(group = Cluster, color=Trajectory),  # Group and color by 'Cluster'
      linewidth=1.5,
      method = "loess",                        # Use LOESS smoothing method
      se = TRUE,
      
      #color = "red"
      # Display standard error bands
    ) +
    scale_color_manual(values = custom_colors)+
    # Apply a clean and classic theme to the plot
    theme_classic() +
    # Add descriptive labels and titles
    labs(
      y = "PRO2",                                    # Y-axis label
      x = "Time (Weeks)",                               # X-axis label
      title = Title,                                    # Plot title
      subtitle = if_else(Include_Subtitles, 
                         "KML Profile", 
                         ""),                        # Plot subtitle
      caption =  if_else(Include_Subtitles, 
                         paste0("Source Data: ", StudyName),
                         "")# Add a caption for data source
    ) +
    # Customize theme elements for better readability
    theme(
      legend.position = "none",         # Hide the legend for clarity
      text = element_text(size = 12),   # Set base text size
      plot.title = element_text(size = 16, face = "bold",
                                vjust = -5.5),  # Title styling
      axis.title.x = element_text(size = 14),                # X-axis title styling
      axis.title.y = element_text(size = 14),                # Y-axis title styling
      axis.text.x = element_text(size = 12),                 # X-axis text size
      axis.text.y = element_text(size = 12),                 # Y-axis text size
      legend.title = element_blank(),                        # Remove legend title
      legend.text = element_text(size = 12),                  # Legend text size
      
    )
  
  # Prepare data for labeling clusters at the maximum time point
  labelData <- data %>% ungroup%>%
    select(Cluster, Time, Value, Trajectory) %>%      # Select relevant columns
    filter(complete.cases(.))             # Remove rows with missing values
  
  # Calculate predicted values at the maximum time point for each cluster
  labelInfo <- labelData %>%
    group_by(Cluster) %>%
    reframe(
      predAtMax = predict(loess(Value ~ Time, span = 0.8),
                          newdata = data.frame(Time = max(Time, na.rm = TRUE))),
      MaxTime = max(Time, na.rm = TRUE), 
      Trajectory=Trajectory
    ) %>%
    # ungroup() %>%
    mutate(
      Label = Trajectory,   # Assign cluster names as labels
      Subject = Cluster  # Assign cluster names for coloring
    ) %>%
    unique()
  
  # Add labels to the plot if 'add_label' is TRUE
  if(add_label) {
    plot <- plot + 
      geom_text_repel(
        data = labelInfo,
        aes(x = MaxTime, y = predAtMax, label = Label, color = Label),
        nudge_x = 5  # Slightly shift labels to the right for clarity
      )
  }
  #plot+ scale_color_manual(values = custom_colors)
  
  # Save the plot as a PDF if 'savePDF' is TRUE
  if (savePDF) {
    tryCatch({
      ggsave(output_file, plot = plot, width = 5, height = 5, units = "in")
      message(paste("Plot successfully saved to", output_file))
    }, error = function(e) {
      warning(paste("Failed to save plot to", output_file, ":", e$message))
    })
  }
  
  return(plot)  # Return the ggplot object for further manipulation or display
}

# Example usage:
# Generate example table with 100 rows
example_table <- data.frame(
  Subject = rep(3322147:3322157, each = 9),
  Cluster = sample(LETTERS[1:5], 99, replace = TRUE),
  CDC = sample(c("CDC", "DC"), 99, replace = TRUE),
  Time = rep(c(0, 2, 4, 8, 12, 16, 20, 24, 32, 36), 10)[1:99],
  Value = sample(0:5, 99, replace = TRUE), 
  Trajectory = sample(c("T1", "T2", "T3"), 99, replace=TRUE)
)

# Generate plot with error checking
plot <- TrajectoryVisualization(data=example_table, 
                                output_file= "Example.pdf", 
                                Title = "Brepocitinib, Ritlecitinib",
                                StudyName="Vibrato",
                                add_label = TRUE, 
                                Include_Subtitles=FALSE,
                                savePDF = FALSE)
print(plot)

ada <- read.csv("C:/Users/ibi9245/OneDrive - Takeda/Projects/VARSITY/adam/ada.csv")
vedo <- read.csv("C:/Users/ibi9245/OneDrive - Takeda/Projects/VARSITY/adam/vedo.csv")
gemini <- read.csv("C:/Users/ibi9245/OneDrive - Takeda/Projects/VARSITY/adam/gemini.csv")

plot <- TrajectoryVisualization(data=ada, 
                                output_file= "ada.pdf", 
                                Title = "VARSITY - Adalimumab",
                                StudyName="VARSITY",
                                add_label = TRUE, 
                                Include_Subtitles=FALSE,
                                savePDF = TRUE)
plot <- TrajectoryVisualization(data=vedo, 
                                output_file= "vedo.pdf", 
                                Title = "VARSITY - Vedolizumab",
                                StudyName="VARSITY",
                                add_label = TRUE, 
                                Include_Subtitles=FALSE,
                                savePDF = TRUE)
plot <- TrajectoryVisualization(data=gemini, 
                                output_file= "gemini.pdf", 
                                Title = "GEMINI - Vedolizumab",
                                StudyName="GEMINI",
                                add_label = TRUE, 
                                Include_Subtitles=FALSE,
                                savePDF = TRUE)
