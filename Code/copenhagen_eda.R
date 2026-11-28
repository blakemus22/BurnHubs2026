library(readr)
library(dplyr)
library(data.table)
library(tidyr)
library(doParallel)
library(igraph)
library(grid)
library(scales)
library(fitdistrplus)
library(colorspace)
library(clue)
library(fields)
library(ggplot2)

edges <- read_csv("~/Fall 2026/BurnHubs2026/Data/edges.csv")
nodes <- read_csv("~/Fall 2026/BurnHubs2026/Data/nodes.csv")

# bad names, I rename them
edges <- edges %>%
  rename(i = `# source`) %>%
  rename(j = target)

# just checking out the data here
range(edges$rssi)
range(edges$timestamp)


edges <- edges %>%
  rename(i = `# source`) %>%
  rename(j = target)


edges2 <- edges %>%
  filter(i != j) %>%
  mutate(a = pmin(i, j), b = pmax(i, j), ts = as.numeric(timestamp)) %>%
  distinct(a, b, timestamp)


edges2$hour = edges2$timestamp/3600

ggplot(edges2, aes(x = hour)) +
  geom_histogram(bins = 672, fill = "#4a235a", color = "#4a235a") +
  scale_x_continuous(
    breaks = seq(0, 672, by = 24),
    labels = function(x) x / 24,
    limits = c(0,672)
  ) +
  labs(x = "Day", y = "Number of contacts") +
  theme_minimal()

contacts_day4 <- contacts_day3 %>%
  filter(duration < 30000)

dim(contacts_day3[contacts_day3$duration >= 30000,])
p <- ggplot(contacts_day4, aes(x = duration)) +
  scale_x_continuous(
    breaks = seq(0, 30000, by = 3000),
    labels = function(x) x/60
  ) +
  geom_histogram(bins = 30, fill = "#4a235a", color = "#4a235a") +
  labs(x = "Duration(minutes)", y = "Number of Contacts") + 
  theme_minimal()
p

bins <- ggplot_build(p)$data[[1]]

ggplot(contacts_day4, aes(x = log10(duration/60))) +
  geom_histogram(bins = 30, fill = "#4a235a", color = "#4a235a") +
  labs(x = "Log(base 10) of Duration", y = "Number of Contacts") + 
  theme_minimal()


mean(contacts_day3$duration/60)
