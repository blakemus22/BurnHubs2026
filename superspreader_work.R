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

edges <- read_csv("~/Fall 2026/BurnHubs2026/edges.csv")
nodes <- read_csv("~/Fall 2026/BurnHubs2026/nodes.csv")
source("Fall 2026/temporal_networks/code/temporal_functions.R")

## I want to try and replicate some of the plots made about super spreading in Pung et al.

# These plots illustrate how often represent 80% of the contacts in a network

#### CODE FROM RETENTION_INDEX.R ####
step_len <- 600
res <- 300
edges <- edges %>%
  rename(i = `# source`) %>%
  rename(j = target) %>% 
  dplyr::filter(rssi >= median(rssi))

## I also want to structure the data following the methods of Pung et al.
# set a time step based on the median contact delay.
contacts <- edges %>%
  filter(i != j) %>%
  mutate(a = pmin(i, j), b = pmax(i, j), ts = as.numeric(timestamp)) %>%
  distinct(a, b, ts) %>%
  arrange(a, b, ts) %>%
  group_by(a, b) %>%
  mutate(new_contact = is.na(lag(ts)) | (ts - lag(ts)) > res,
         contact_id  = cumsum(new_contact)) %>%
  group_by(a, b, contact_id) %>%
  summarise(time_start = min(ts), time_end = max(ts) + res, .groups = "drop") %>%
  mutate(duration = time_end - time_start) %>%
  arrange(a, b, time_start) %>%
  group_by(a, b) %>%
  mutate(inter_duration = lead(time_start) - time_end) %>%
  ungroup()

quantile(contacts$duration,       probs = c(0.25, 0.5, 0.75))
quantile(contacts$inter_duration, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
# the median delay between two contacts is 600 seconds
# by the logic of Pung et al., we should set the time period to be 600 seconds. 

# filter out all the contacts under 600 seconds
contacts2 <- contacts %>%
  dplyr::filter(duration >= median(duration)) %>%
  dplyr::select(a, b, time_start, time_end, duration)

# separate into time steps: each row corresponds to a pair that had contact during a time step
contacts_steps <- contacts2 %>%
  mutate(step_first = floor((time_start) / step_len) + 1,
         step_last  = floor((time_end - 1) / step_len) + 1,   # time_end is exclusive
         n_steps    = step_last - step_first + 1) %>%
  uncount(n_steps, .id = "k") %>%                 # one row per step spanned
  mutate(step       = step_first + k - 1,
         step_begin = (step - 1) * step_len,
         step_stop  = step_begin + step_len) %>%
  dplyr::select(i = a, j = b, step)
# here, we process the data to calculate r_bar:
# make contacts symmetric:
edges_sym <- bind_rows(
  contacts_steps,
  contacts_steps %>% rename(i = j, j = i)    # swap the two columns
)

# now, we need to calculate the mean number of retained contacts per time step
node_list <- sort(unique(edges_sym$i))


n_steps <- 4032

# filter to one edge per pair per step - a = min of pair, b = max of pair
edges_step <- edges_sym %>%
  filter(i != j) %>%
  mutate(a = pmin(i, j),
         b = pmax(i, j)) %>%
  distinct(a, b, step)


#### New code ####
deg_long <- bind_rows(
  edges_step %>% dplyr::select(node = a, step),
  edges_step %>% dplyr::select(node = b, step)
) %>%
  count(node, step, name = "k0")

# data frame of nodes that are responsible for the top 80% of connections at each step
p80 <- deg_long %>%
  group_by(step) %>%
  mutate(total = sum(k0)) %>%
  filter(total > 1) %>%                       
  arrange(desc(k0), .by_group = TRUE) %>%
  mutate(cum_before = lag(cumsum(k0), default = 0)) %>%
  filter(cum_before < 0.8 * total) %>%           # keep a node if we haven't reached 80% before it
  ungroup()

# find the proportion of steps that each node is a super spreader
node_top80 <- p80 %>%
  group_by(node) %>%
  summarise(prop = n()/4032) %>%
  ungroup()

prop_thresh <- function(x){
  # finds the proportion of nodes that are superspreaders in more than x steps
  return(sum(node_top80$prop > x/200)/680)
}
prop_thresh_vec <- unlist(lapply(0:200, prop_thresh))

plot_df <- data.frame(threshold = (0:200) / 200,
                      prop      = prop_thresh_vec)

plot_df <- plot_df[plot_df$prop > 0, ]


## lets make some reference curves too
        
# rewired every step: exact binomial
thr_fine <- seq(0, 1, length.out = 201)
prop_df <- data.frame(
  threshold = thr_fine,
  prop_low     = 1 - pbinom(ceiling(thr_fine * 100) - 1, 100, .2),
  prop_high    = 1 - pbinom(ceiling(thr_fine * 4032) - 1, 4032, .2),
  prop_real    = prop_thresh_vec
)

prop_df$prop_real[prop_df$prop_real == 0] <- NA
prop_df$prop_low[prop_df$prop_low == 0] <- NA
prop_df$prop_high[prop_df$prop_high == 0] <- NA


ggplot(prop_df) +
  geom_line(aes(x = prop_real, y = threshold), colour = "purple", linewidth = 1.5, alpha = 0.8) +
 ylim(0,1) + 
  xlim(0,1) +
  labs(title = "",
       x = "Proportion of Nodes",
       y = "Proportion of Steps") +
  theme_minimal(base_size = 14) + 
  theme(aspect.ratio = 1)


