##***************************
## Assignment 01 - Group 10 
## 2026-11-02
##
##***************************
## _ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())

rm(list = ls())

setwd("./r/")

##Import TE annotation datasets -----

## Import T. patagoniensis dataset
df_patagoniensis <- read_csv("../data/GCA_001447655.1.csv")
df_patagoniensis 


## Import T. spiralis dataset
df_spiralis = read_csv("../data/trichinella_spiralis_te.csv")
df_spiralis 

## Clean dataset columns -----

##__James: First Edit----

#I combined the two data frames at the start, before further modifications
#This reduces the lines of code by roughly half, since two separate data frames are no longer being modifed separately before they are joined

#I also deleted the code that used str_replace_all to remove brackets around the numbers in the qry_left, in_repeat_end, and no_bp_in_complement columns before joining the two data frames together because these columns are not used for downstream analyses anyways, so it was easier to just filter out these columns by selecting for the columns of interest which was matching_class and qry_width before joining the data frames
df_combined = bind_rows(
    df_spiralis |>
      select(matching_class, qry_width) |>
      mutate(species = "T. spiralis"),
    df_patagoniensis |>
      select(matching_class, qry_width) |>
      mutate(species = "T. patagoniensis"))

View(df_combined) #End of my first edit

## Define non-TE classes to exclude -----

ls_non_te = c("simple_repeat", "Simple_repeat",
              "Low_complexity", "Satellite",
              "rRNA", "snRNA", "ARTEFACT",
              "Unknown")

##__James: Second Edit----

#I created a list of TE classes that are found in both species, and filtered matching class by this list

#Justification: I noticed that the two species do not contain exactly the same TE classes as each other
#For example, in the bar plot made later on, you can see that some TE classes are missing a bar completely in either species
#I thought that by only using TE classes in common between the two species, it would eliminate the potential confounding effects of different TE classes between species on how the relationship between TE annotation count and total annotated sequence length differs between the two species
#Therefore, I used the following code to create a list of TE classes in common between the two species, and then further down I filter the matching class by this list

common_class = df_spiralis |>
  filter(matching_class %notin% ls_non_te) |>
  filter(matching_class %in% df_patagoniensis$matching_class) |>
  select(matching_class) |>
  unique() |>
  pull() |>
  as.list()

## Summarize TE annotations -----

df_summary = df_combined |>
  filter(matching_class %in% common_class) |> #End of my second edit
  group_by(species, matching_class) |>
  summarise(
    annotation_count = n() ,
    total_length = sum(qry_width, na.rm = TRUE)) |>
  arrange(desc(annotation_count))

View(df_summary)

## Rank TE classes by annotation count -----

df_summary |>
  group_by(species) |>
  arrange(desc(annotation_count), .by_group = TRUE)


## Rank TE classes by total annotated length -----

df_summary |>
  group_by(species) |>
  arrange(desc(total_length), .by_group = TRUE)


# Top 3 TE classes by total length for each species
df_summary |> 
  group_by(species) |> 
  arrange(desc(total_length), .by_group = TRUE) |> 
  slice_head(n = 3)

## Plot distribution of TE annotation lengths -----

df_combined |>
  filter(matching_class %in% common_class) |>
  ggplot(aes(x = qry_width)) +
  geom_histogram() +
  facet_grid(~ species) +
  labs(
    title = "Distribution of TE Annotation Lengths",
    x = "Annotated TE length (bp)",
    y = "Frequency"
  )

## Plot annotation counts against total annotated lengths -----

df_summary |> 
  ggplot(aes(x = annotation_count, y = total_length)) +
  geom_point() +
  geom_smooth() +
  facet_grid(~ species) +
  labs(title = "Comparison of TE Abundance & Length",
       x = "annotated count",
       y = "annotated sequence length(bp)")


## Plot average TE length by class -----

##__James: Third Edit----

#I modified the bar plot to show the annotation count via a colour gradient for each TE, in addition to already showing the average annotated length

#Justification: This is the only plot of the three to use average TE length as a variable. Since total TE length is expected to increase with annotation count, I was interested whether this trend would still hold with average TE length
#This provides further insight into the research question by potentially illuminating whether an observed relationship between TE annotation count and total annotated length holds even when the total length is averaged out by the number of contributing TEs
  
df_summary |>
  ggplot(aes(x = matching_class,
             y = total_length / annotation_count,
             fill = annotation_count)) +
  geom_col() +
  coord_flip() +
  facet_grid(~ species) +
  labs(
    title = "Average TE Length by Class",
    x = "TE class",
    y = "Average annotated length (bp)",
    fill = "Annotation Count"
  )

