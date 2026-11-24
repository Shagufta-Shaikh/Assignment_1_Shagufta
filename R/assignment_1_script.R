#### Snipets 
##***************************
## Assignment-1 : Shagufta, Yeshi, Diana, Rachel
##
## 2026-11-02
##
##***************************

## contribution statement ----
## Shagufta was the primary author of the Assignment 1 script.
# Kristian revised the CSV paths, generated chi-square counts from the data,
# corrected the Figure 2 standard errors, and plotted within-species
# percentages in Figure 1.

## Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())

##Load the data -----

trichinella_spiralis <- read.csv("../Assignment_1_Shagufta/data/GCA_001447595.1.csv")
trichinella_zimbabwensis <- read_csv("../Assignment_1_Shagufta/data/GCA_001447665.1.csv")

##WRANGLE -----
###make a vector containing strings of non TE names 
ls_non_te = c("simple_repeat", "Simple_repeat",
              "Low_complexity", "Satellite",
              "rRNA", "snRNA", "ARTEFACT",
              "Unknown")

####filter repeat dataset for only TE's 
trichinella_zimbabwensis_te_cleaned <- trichinella_zimbabwensis|> 
  filter(matching_class %notin% ls_non_te)

trichinella_spiralis_te_cleaned <- trichinella_spiralis |> 
  filter(matching_class %notin% ls_non_te)

#####Find TE lengths
trichinella_spiralis_te <- trichinella_spiralis_te_cleaned |> 
  mutate(
    te_length = qry_width,
    te_category = case_when(
      str_detect(matching_class, "^DNA|^RC") ~ "DNA_transposon",
      str_detect(matching_class, "^(LINE|LTR|PLE|SINE)") ~ "Retrotransposon"))

trichinella_zimbabwensis_te <- trichinella_zimbabwensis_te_cleaned |> 
  mutate(
    te_length = qry_width,
    te_category = case_when(
      str_detect(matching_class, "^DNA|^RC") ~ "DNA_transposon",
      str_detect(matching_class, "^(LINE|LTR|PLE|SINE)") ~ "Retrotransposon"))

# very efficent and smart way of creating new coloumns


######Calculate TE category proportions for each species
spiralis_proportions <- trichinella_spiralis_te |>
  count(te_category) |>
  mutate(proportion = n / sum(n),
         species = "T. spiralis")

zimb_proportions <- trichinella_zimbabwensis_te |>
  count(te_category) |>
  mutate(proportion = n / sum(n),
         species = "T. zimbabwensis") 

#######combine both datasets for TE length comparison 
te_fragment_df <- bind_rows(
  trichinella_spiralis_te |>
    mutate(species = "T. spiralis"),
  trichinella_zimbabwensis_te  |> 
    mutate(species = "T. zimbabwensis")) |>
  mutate(
    species = factor(species, levels = c("T. spiralis", "T. zimbabwensis")),
    te_category = factor(te_category,
                         levels = c("DNA_transposon", "Retrotransposon")))

#STATISTICAL ANALYSIS-----
###chi square test for independence of TE category and species 
counts <- matrix(
  c(120, 917,
    207, 1072),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Species = c("T. spiralis", "T. zimbabwensis"),
    TE_class = c("DNA transposons", "Retrotransposons")))

#first edit ----
# The original script used four manually typed counts to build the table for the chi-square test.

# R can calculate those same four counts directly from the fragment data in fewer lines of code.

# table() counts how many DNA transposon and retrotransposon fragments belong to each species.

# I checked that all four counts matched the original values. This also reduces the chance of a
# typing error and keeps the test up to date if someone changes the data or filtering later.
counts <- table(Species = te_fragment_df$species, TE_class = te_fragment_df$te_category)

chisq.test(counts)

####calculate standard-error of proportions
#deleted old block to avoid errors in downstream figure creation with new variables
## second edit ----
# The original script calculated each standard error separately and then manually entered
# the results into the tables used for Figure 2. R can calculate and store them in one step.

# I also corrected the denominator. For example, 120 out of 1,037 T. spiralis fragments
# are DNA transposons, so all 1,037 fragments contribute to that proportion.

# Using sum(n) gives the correct total for each species and updates the error bars automatically.
spiralis_proportions_er <- spiralis_proportions |>
  mutate(st_error = sqrt(proportion * (1 - proportion) / sum(n)))

zimb_proportions_er <- zimb_proportions |>
  mutate(st_error = sqrt(proportion * (1 - proportion) / sum(n)))

###two way anova for TE length
aov(te_length ~ species * te_category, data = te_fragment_df)


#GRAPH 1 EXPLORATORY ------
###Histogram comparing TE length between both species
## third edit ----
# The original figure showed how many TE fragments fell into each length bin. T. zimbabwensis
# had 1,279 fragments and T. spiralis had 1,037, so a taller peak could partly reflect the
# larger number of fragments in one dataset. I changed the y-axis to show percentages.

# R divides the count in each length bin by the total fragment count for that species.

# This makes it easier to see which fragment lengths are more common within each species.

figure_1 <- ggplot() +
  geom_freqpoly(
    data = trichinella_zimbabwensis_te,
    aes(x = te_length,y = after_stat(count / sum(count)), colour = "T.zimbabwensis"), 
    alpha = 0.6,
    bins = 200,
    linewidth = 1.0
  ) +
  geom_freqpoly(
    data = trichinella_spiralis_te, 
    aes(x = te_length,y = after_stat(count / sum(count)), colour = "T.spiralis"),
    alpha = 0.6,
    bins = 200,
    linewidth = 1.0 
  ) + 
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "TE Length (bp)",
    y = "Percentage of fragments within species",
    colour = ""
  ) +
  coord_cartesian(xlim= c(0,300)
  ) +
  theme(legend.position = "bottom")

figure_1

#GRAPH 2 FINAL RESULTS----
####barchart for TE class proportion 

figure_2 <- ggplot() +
  geom_col(
    data = spiralis_proportions_er,
    aes(x = te_category, y = proportion, fill = "T.spiralis"),
    position = position_nudge(x = -0.2),
    width = 0.35,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_errorbar(
    data = spiralis_proportions_er,
    aes(x = te_category,
        ymin = proportion - st_error,
        ymax = proportion + st_error),
    position = position_nudge(x = -0.2),
    width = 0.12
  ) +
  geom_col(
    data = zimb_proportions_er,
    aes(x = te_category, y = proportion, fill = "T.zimbabwensis"),
    position = position_nudge(x = 0.2),
    width = 0.35,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_errorbar(
    data = zimb_proportions_er,
    aes(x = te_category,
        ymin = proportion - st_error,
        ymax = proportion + st_error),
    position = position_nudge(x = 0.2),
    width = 0.12
  ) + 
  labs(
    x = "",
    y = "Proportion of Total TE's (%)",
    fill = ""
  ) +
  scale_fill_manual(
    values = c(
      "T.spiralis" = "#E69F00",
      "T.zimbabwensis" = "#98A0ff")
  ) +
  scale_y_continuous(labels = scales::percent
  ) +
  theme(legend.position = "bottom")

figure_2

##GRAPH 3 FINAL RESULTS -----
##TE length

figure_3 <- ggplot(te_fragment_df, aes(x = species, y = te_length, fill = species)
) +
  geom_boxplot(outlier.shape = NA, 
               width = 0.7,
               alpha = 0.8, 
               linewidth = 0.4, 
               color = "black"
  ) +
  geom_jitter(width = 0.2,
              size = 0.4, 
              alpha = 0.2
  ) +
  facet_wrap(~ te_category, nrow = 1
  ) +
  scale_y_log10(
    breaks = c(50, 100, 200, 500, 1000, 2000),
  ) + scale_fill_manual(
    values = c("T. spiralis" = "#0072B2",
               "T. zimbabwensis" = "#E69F00"),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "TE length (bp, log10 scale)"
  ) +
  theme(legend.position = "bottom")

figure_3
