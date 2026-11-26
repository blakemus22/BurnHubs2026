import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

#INPUTS

FILENAME = Path(__file__).resolve().parent / "copenhagen_contact_mat.csv"
p = 1 # Probability of transmission: Beta = p*C
GAMMA = 0.2
INITIAL_INFECTED = 5

NUM_SIMS = 100

all_histories = []
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


def outbreak_sim(contacts,p,gamma,initial_infected):

    rng = np.random.default_rng()

    D, N, N2 = contacts.shape

    max_days = D

    state = np.zeros(N, dtype=int) # 0 is sus, 1 is infected, 2 is recovered

    infection_day = np.full(N, np.nan)
    recovery_time = np.full(N, np.nan)

    initial_ids = rng.choice(N, size=initial_infected)

    state[initial_ids] = 1
    infection_day[initial_ids] = 0

    recovery_time[initial_ids] = 1/gamma #rng.exponential(scale = 1/gamma, size=initial_infected)

    records = []

    r_eff = np.zeros(max_days + 1)

    for day in range(max_days + 1):

        records.append({ "day": day,"S": np.sum(state == 0),"I": np.sum(state == 1),"R": np.sum(state == 2)})

        if day == max_days:
            break

        C = contacts[day]

        infected = (state == 1)
        susceptible = (state == 0)

        beta = p * C[:, infected].sum(axis=1)
        # p_infection = 1 - np.exp(-force_infection)

        new_infections = (susceptible & (rng.random(N) < beta))
        state[new_infections] = 1
        infection_day[new_infections] = day + 1

        r_eff[day] = 1/gamma * new_infections.sum() / np.count_nonzero(infected) if np.count_nonzero(infected) > 0 else 0


        new_recoveries = (infected & (recovery_time <= day + 1))
        state[new_recoveries] = 2

        n_new = np.count_nonzero(new_infections)

        if n_new > 0:
            recovery_time[new_infections] = (day + 1/gamma + 1)

    epidemic_size = np.sum(state == 2) + np.sum(state == 1)


    history = pd.DataFrame(records)

    return history, epidemic_size, r_eff

contacts = load_contact_network(FILENAME)

contacts = contacts[1:]

history, epidemic_size, r_eff = outbreak_sim(contacts, p, GAMMA, INITIAL_INFECTED)

for sim in range(NUM_SIMS):
    history, epidemic_size, r_eff = outbreak_sim(contacts, p, GAMMA, INITIAL_INFECTED)
    infected = history["I"].iloc[-1] + history["R"].iloc[-1]
    history["R_eff"] = r_eff

    total_infected.append(infected)
    all_histories.append(history)

total_infected = np.array(total_infected)


plt.figure(figsize=(10,6))

plt.hist(total_infected, bins=40, density=True, alpha=0.7, color="blue")
plt.title("Distribution of Epidemic Sizes")
plt.xlabel("Final Epidemic Size")
plt.ylabel("Density")
plt.grid(axis="y", alpha=0.3)
plt.show()


plt.figure(figsize=(10,6))
plt.title("R_eff distribution over time")
for sim in range(NUM_SIMS):
    plt.plot(all_histories[sim][["R_eff"]], alpha=0.1, color="blue")
plt.xlabel("Day")
plt.ylabel("R Effective")
plt.grid(axis="y", alpha=0.3)
plt.show()


print(f"Final epidemic size: {epidemic_size} out of {contacts.shape[1]} individuals")
print(f"R_eff: {r_eff.tolist()}")