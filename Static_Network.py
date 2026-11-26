
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# INPUTS

FILENAME = "copenhagen_contact_mat.csv"
p = 1
GAMMA = 0.2
INITIAL_INFECTED = 5
chosen_day = 1
NUM_SIMS = 100

all_histories = []
weekend_histories = []
weekday_histories = []
weekend_infected = []
weekday_infected= []
total_infected = []

def load_contact_network(filename):

    df = pd.read_csv(filename)

    person_ids = np.sort(df["i"].unique())
    days = np.sort(df["day"].unique())

    contact_cols = [f"j_{i}" for i in person_ids]

    N = len(person_ids)
    contacts = np.zeros((len(days), N, N), dtype=float)

    for k, day in enumerate(days):
        daily = (df.loc[df["day"] == day].set_index("i").reindex(person_ids))

        C = daily[contact_cols].to_numpy(dtype=float)
        contacts[k] = C / 288

        np.fill_diagonal(contacts[k], 0) # No self loops

    return contacts


def outbreak_sim_static(contacts,p,gamma,initial_infected,chosen_day):

    rng = np.random.default_rng()

    D, N, N2 = contacts.shape
    
    C = contacts[chosen_day]
    
    max_days = D

    state = np.zeros(N, dtype=int) # 0 is sus, 1 is infected, 2 is recovered

    infection_day = np.full(N, np.nan)
    recovery_time = np.full(N, np.nan)

    initial_ids = rng.choice(N, size=initial_infected,replace=False)

    state[initial_ids] = 1
    infection_day[initial_ids] = 0

    recovery_time[initial_ids] = 1/gamma #rng.exponential(scale = 1/gamma, size=initial_infected)

    records = []

    for day in range(max_days + 1):

        records.append({ "day": day,"S": np.sum(state == 0),"I": np.sum(state == 1),"R": np.sum(state == 2)})

        if day == max_days:
            break

        infected = (state == 1)
        susceptible = (state == 0)

        force_infection = p * C[:, infected].sum(axis=1)
        p_infection = 1 - np.exp(-force_infection)

        new_infections = (susceptible & (rng.random(N) < p_infection))
        state[new_infections] = 1
        infection_day[new_infections] = day + 1
        
        new_recoveries = (infected & (recovery_time <= day + 1))
        state[new_recoveries] = 2

        n_new = np.count_nonzero(new_infections)

        if n_new > 0:
            recovery_time[new_infections] = (day + 1 + 1/gamma)

    history = pd.DataFrame(records)

    return history


contacts = load_contact_network(FILENAME)

contacts = contacts[1:]

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
bins = np.linspace(0, 450, 31)
axes = axes.flatten()

for day in range(28):

    total_infected = []

    for sim in range(NUM_SIMS):

        history = outbreak_sim_static(
            contacts, p, GAMMA, INITIAL_INFECTED, day
        )

        infected = history["I"].iloc[-1] + history["R"].iloc[-1]

        total_infected.append(infected)

        if (day % 7 == 0) or (day % 7 == 6):
            weekend_infected.append(infected)
            weekend_histories.append(history)
        else:
            weekday_infected.append(infected)
            weekday_histories.append(history)

axes[0].hist(weekend_infected, bins=bins, edgecolor="black")
axes[0].set_title("Weekend")
axes[0].set_xlabel("Total Infected")
axes[0].set_ylabel("Frequency")

axes[1].hist(weekday_infected, bins=bins, edgecolor="black")
axes[1].set_title("Weekday")
axes[1].set_xlabel("Total Infected")
axes[1].set_ylabel("Frequency")

for ax in axes:
    ax.set_xlim(0, 500)
    ax.set_ylim(0,600)
    ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.show()
