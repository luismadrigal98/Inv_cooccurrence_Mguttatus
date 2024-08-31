##' Inversion dictionary creation
##'
##' This script reads a detailed inversion regions file, processes the data to 
##' select unique inversion IDs and their corresponding chromosomes, and writes 
##' the processed data to a new file.
##'
##' Input: 
##' - C:\\Users\\usuario\\OneDrive - University of Kansas\\Projects\\Inversion paper\\summary.new.inversion.csv"
##'
##' Output:
##' - "C:\\Users\\usuario\\OneDrive - University of Kansas\\Projects\\Co-occurence of inversions\\anchorwave.inversions.txt"
##'
##' Steps:
##' 1. Read the input file containing detailed inversion regions.
##' 2. Select relevant columns (INV_ID, chrom, alt_lines) and ensure uniqueness.
##' 3. Arrange the data by INV_ID and chrom.
##' 4. Write the processed data to the output file.
##'
##' Author: Luis Javier Madrigal Roca
##' Date: 08/31/2024
##' ____________________________________________________________________________
##' >>>>>>>>>>>>>>>>>>>>>>>>>>>>> Start <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

library(dplyr)
library(gtools)

Data <-  read.csv("C:\\Users\\usuario\\OneDrive - University of Kansas\\Projects\\Inversion paper\\summary.new.inversion.csv",
                    header = TRUE)

Data <- Data |> select(INV_ID, Chr, Line) |>
  distinct(.keep_all = TRUE)

line_condenser <- function(df, reference, field_to_condense) 
{
  ##' This function will take all the observations with the same reference value,
  ##' and will parse the field_to_condense values into a single string.
  ##' 
  ##' @param df A data frame with the data to be condensed.
  ##' @param reference The column name that will be used to group the data.
  ##' @param field_to_condense The column name that will be condensed.
  ##' __________________________________________________________________________
  
  df <- df |> arrange(!!sym(reference), !!sym(field_to_condense))
  
  df <- df |> group_by(!!sym(reference), Chr) |> 
    summarise(!!sym(field_to_condense) := paste(!!sym(field_to_condense), 
                                                collapse = ";"), 
              .groups = 'drop')
  
  return(df)
}

Data <- line_condenser(Data, "INV_ID", "Line")

Data <- Data |> arrange(INV_ID, Chr)

Data <- Data[mixedorder(Data$INV_ID),]

write.table(Data, "C:\\Users\\usuario\\OneDrive - University of Kansas\\Projects\\Co-occurence of inversions\\list.of.inversions.txt", 
      sep = "\t", row.names = FALSE, quote = FALSE)

##' >>>>>>>>>>>>>>>>>>>>>>>>>>> End of script <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<