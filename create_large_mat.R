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

edges <- read_csv("~/Fall 2026/BurnHubs2026/edges.csv")
nodes <- read_csv("~/Fall 2026/BurnHubs2026/nodes.csv")
tz = "UTC"

# bad names, I rename them
edges <- edges %>%
  rename(i = `# source`) %>%
  rename(j = target)

# just checking out the data here
range(edges$rssi)
range(edges$timestamp)


contacts_day <- edges %>%
  filter(i != j) %>%
  mutate(a   = pmin(i, j),
         b   = pmax(i, j),
         ts  = as.numeric(timestamp),
         day = as.Date(as.POSIXct(ts, origin = "1970-01-01", tz = tz))) %>%
  distinct(a, b, ts, day) %>%
  arrange(a, b, ts) %>%
  group_by(a, b) %>%
  mutate(
    # consecutive reading, but on a different day: this is the cut we are adding
    bridged     = !is.na(lag(ts)) & (ts - lag(ts)) <= res & day != lag(day),
    new_contact = is.na(lag(ts)) | (ts - lag(ts)) > res | day != lag(day),
    contact_id  = cumsum(new_contact)
  ) %>%
  group_by(a, b, contact_id, day) %>%
  summarise(time_start = min(ts),
            time_end   = max(ts) + res,
            bridged    = first(bridged),    # TRUE if this piece is the 2nd half of a split contact
            .groups = "drop") %>%
  mutate(duration = time_end - time_start) %>%
  arrange(a, b, time_start) %>%
  group_by(a, b) %>%
  mutate(inter_duration = ifelse(lead(bridged, default = FALSE),
                                 NA,                                  # artificial gap at the split
                                 lead(time_start) - time_end)) %>%
  ungroup()

## Cleaning the data to have a better stucture
# Following the methods of pung et al, we filter 
# to only include the contacts greater than 600 seconds

# Then, we set a time step of one day
# Each row corresponds to an individual i's contacts on one day t.
# Each column corresponds to indiviual j.
# the value of an entry corresponds to how many minutes 
# inviduals i and j were connected on day t

# filter out the negligible contacts and created a "day" variable
contacts_day2 <- contacts_day %>%
  dplyr::filter(duration >= median(duration)) %>%
  mutate(day = floor(time_start/(60*60*24)) + 1,
         duration = duration /60)
  
# the weights represent how many minutes the pair was in contact
contacts_weight <- contacts_day2 %>%
  group_by(a, b, day) %>%
  summarize(count = sum(duration), .groups = "drop")

# make contacts symmetric:
contacts_weight_sym <- bind_rows(
  contacts_weight,
  contacts_weight %>% rename(a = b, b = a)    # swap the two columns
)
contacts_weight_sym$day <- as.integer(contacts_weight_sym$day)

# grid with all possible i, j pairs on all possible days
ids  <- sort(unique(c(contacts_weight_sym$a, contacts_weight_sym$b)))
days <- 1:28
grid <- CJ(a = ids, b = ids, day = days)

# now we join the grid up with the count data
# this step can be edited to produce an unweighted graph instead of weighted graph
contacts_weight_sym2 <- left_join(grid, contacts_weight_sym, by = c("a", "b", "day")) %>%
  mutate(count = replace_na(count, 0))


# expand this into a matrix. Each entry represents the number of contacts
# that a pair (i,j) had during that day. 
large_contact_mat <- contacts_weight_sym2 %>%
  pivot_wider(names_from = "b",
              names_prefix = "j_",
              values_from = count)
write.csv(large_contact_mat, "C:/Users/blake/Documents/Fall 2026/BurnHubs2026/copenhagen_contact_mat.csv", row.names = TRUE)
