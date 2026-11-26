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

unique_nodes <- contacts_steps %>%
  group_by(step) %>%
  summarise(tot_nodes = n_distinct(c(i,j))) %>%
  ungroup()
mean(unique_nodes$tot_nodes)
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

N_total <- n_distinct(c(edges_step$a, edges_step$b))

# offset the step by 1 to see what edges are present at the next step
next_edges <- edges_step %>%
  mutate(step = step - 1,       
         retained = 1)
# rejoin with previous data frame
# if a contact is present in consecutive steps, it will have retained = 1
edges_step <- edges_step %>%
  left_join(next_edges, by = c("a", "b", "step")) %>%
  mutate(retained = ifelse(is.na(retained), 0L, retained))

long <- bind_rows(
  edges_step %>% dplyr::select(node = a, step, retained),
  edges_step %>% dplyr::select(node = b, step, retained)
)

# remove the last step, then calculate the number of connections + retained per node per step
node_step <- long %>%
  filter(step < n_steps) %>%
  group_by(node, step) %>%
  summarise(k0 = n(),
            r = sum(retained),
            .groups = "drop")

# mean retained contacts per step, among ALL nodes(even those with no contacts)
# r_stat is the mean degree of the network at each step
mean_r_by_step <- node_step %>%
  group_by(step) %>%
  summarise(r_bar_temp = sum(r)/N_total,
            r_bar_stat = sum(k0)/N_total,
            n_active = n(),
            .groups = "drop")

# lastly, we need to calculate r_dynam
# r_dynam: retained contacts if the network got randomly rewired every step

# degree of every node at every step
deg <- bind_rows(
  edges_step %>% dplyr::select(node = a, step),
  edges_step %>% dplyr::select(node = b, step)
) %>%
  count(node, step, name = "k")

# shift back one step so the t+1 values line up with step t
deg_next <- deg %>% mutate(step = step - 1L)
E_next   <- edges_step %>% count(step, name = "E1") %>% mutate(step = step - 1L)


rand_next <- edges_step %>%
  # set the degree at the NEXT step
  left_join(deg_next %>% rename(k1_a = k), by = c("a" = "node", "step")) %>%
  left_join(deg_next %>% rename(k1_b = k), by = c("b" = "node", "step")) %>%
  # set the total degree for each step
  left_join(E_next, by = "step") %>%
  # replace NA values, p = probability of edge between a and b given degree distributions
  mutate(k1_a = coalesce(as.numeric(k1_a), 0),
         k1_b = coalesce(as.numeric(k1_b), 0),
         E1   = coalesce(as.numeric(E1), 0),
         p    = ifelse(E1 > 0, pmin(1, k1_a * k1_b / (2 * E1)), 0)) %>%
  group_by(step) %>%
  # sum over the existing edges at each step
  summarise(exp_retained_edges = sum(p), .groups = "drop") %>%
  # find the expected retained edges per node
  mutate(mean_r_rand_next = 2 * exp_retained_edges / N_total)

rand_next <- data.frame(step = 1:(n_steps - 1)) %>%
  left_join(rand_next, by = "step") %>%
  mutate(mean_r_rand_next = coalesce(mean_r_rand_next, 0))

mean_r_by_step$r_bar_dynam <- rand_next$mean_r_rand

# now we can finally compute r_bar!
mean_r_by_step$r_bar <- (mean_r_by_step$r_bar_temp - mean_r_by_step$r_bar_dynam)/(mean_r_by_step$r_bar_stat - mean_r_by_step$r_bar_dynam)


ggplot(mean_r_by_step, aes(x = r_bar)) +
  geom_density(fill = "purple", colour = "purple4", alpha = 0.5, linewidth = 0.8) +
  geom_vline(aes(xintercept = median(r_bar, na.rm = TRUE)),
             linetype = "dashed", linewidth = 0.8) +
  xlim(0,1) +
  labs(x = expression(bar(r)),
       y = "Density") +
  theme_minimal(base_size = 14) +
  theme(axis.text = element_text(size = 18),
        axis.text.y = element_blank())

    
