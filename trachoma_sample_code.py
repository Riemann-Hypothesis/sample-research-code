# A small sample of the ~1000 lines of code used in my trachoma research project.

import csv
import numpy as np
import math
import networkx as nx
from multiprocessing import Pool

class IterateObject:
    def __init__(self, state_vector, treatment_frequency):
        self.state_vector = state_vector
        self.treatment_frequency = treatment_frequency

class Village:
    def __init__(self, id, visits):
        self.id = id
        self.visits = visits
        self.first_agent_info = {}
        # List of unique agents: index of agents in the list is their index in the corresponding correlation matrix
        self.unique_agents = []

    def numAgents(self):
        return len(self.unique_agents)

    # In this version, I'm making the process of generating an adjacency matrix for a population a function with input of VID (as an integer)
    def matrices(self, base_transmission_rate, base_recovery_rate, h, ti, o):
        numAgents = self.numAgents()
        beta_initial = base_transmission_rate
        gamma_initial = base_recovery_rate
        # Note: this line is no longer needed because num_unique_agents is now global
        # num_unique_agents = len(village.unique_agents)

        # Initialize the matrix; it will be NxN where N is the number of unique agents in the village
        # Initializing an NxN matrix where all entries are 0.
        #matrix = np.zeros((num_unique_agents, num_unique_agents))

        # Initialize infection matrix.
        infection_matrix = np.zeros((numAgents, numAgents))

        # Initialize recovery matrix.
        recovery_matrix = np.zeros((numAgents, numAgents))

        # Note: matrix[i,j] represents the probability that j infects i
        for potential_carrier in self.unique_agents:
            j = self.unique_agents.index(potential_carrier)
            for potential_victim in self.unique_agents:
                i = self.unique_agents.index(potential_victim)

                # Check to see if the agents are in the same household. Since household ID is simply the first 7 digits of the unique agent ID, we check if the two agents have the same first 7 digits.
                #if potential_carrier[:7] == potential_victim[:7]:
                #    same_household = 1
                #else:
                #    same_household = 0
                # If we are updating a diagonal entry in the matrix, fill it with the base recovery rate, gamma_initial.
                if potential_carrier == potential_victim:
                    recovery_matrix[i,j] = gamma_initial

                # If we are updating any other entry in the matrix, fill it with the base transmission rate, beta_initial, and modify it with the parameters
                else:
                    # When same_household = 1, h comes into effect
                    # matrix[i,j] = beta_initial + (h * same_household)
                    folded_visits = []
                    same_household = self.first_agent_info[potential_victim].hid == self.first_agent_info[potential_carrier].hid
                    infection_matrix[i,j] = beta_initial
                    if same_household:
                        infection_matrix[i,j] *= h

                    if self.first_agent_info[potential_carrier].ti:
                        infection_matrix[i,j] *= ti

                    if self.first_agent_info[potential_carrier].ocular:
                        infection_matrix[i,j] *= o

        # save the matrices to files
        #np.save("matrix_saves/vid_{}_infection_matrix_bin_household_ti_ocular.npy".format(self.id), infection_matrix)
        #np.savetxt("matrix_saves/vid_{}_infection_matrix_household_ti_ocular.txt".format(self.id), infection_matrix)

        #np.save("matrix_saves/vid_{}_recovery_matrix_bin_household_ti_ocular.npy".format(self.id), recovery_matrix)
        #np.savetxt("matrix_saves/vid_{}_recovery_matrix_household_ti_ocular.txt".format(self.id), recovery_matrix)

        return infection_matrix, recovery_matrix

    def makeRandomInitialVector(self, initial_infected):
        random_initial_vector = []

        # length of the vector must be the length of the village that it corresponds to; that is, it must be num_unique_agents long

        # initial_infected of the population starts infected
        initial_infected_count = round(initial_infected * len(self.unique_agents))
        initial_infected_list = [1 for i in range(initial_infected_count)]

        # 1 - initial_infected of the population starts susceptible
        initial_susceptible_count = len(self.unique_agents) - initial_infected_count
        initial_susceptible_list = [0 for i in range(initial_susceptible_count)]

        # Now, we combine these to create our full initial state vector.
        random_initial_vector = initial_infected_list + initial_susceptible_list

        # Now, we shuffle the list to randomize the initial state vector.
        random.shuffle(random_initial_vector)

        return random_initial_vector


class Agent:
    def __init__(self, id, infected, age, ocular, tf, ti, gender, nasal, flies):
        self.id = id
        # household is contained in the ID (first 7 numbers)
        self.hid = id[:7]

        self.infected = infected

        # adding in other parameters [NOTE THAT THEY DO NOT WORK YET]
        # important parameters
        self.age = age
        self.ocular = ocular
        self.tf = tf
        self.ti = ti

        # less important parameters
        self.gender = gender
        self.nasal = nasal
        self.flies = flies

# Pull and sort the data in the spreadsheet
with open ("trachomadata.csv") as csvfile:
    trachomadata = list(csv.reader(csvfile, delimiter=","))

    # Find the index of the column with the village IDs (VID)
    vid_idx = trachomadata[0].index("VID")

    # Find the index of the column with the unique agent identifiers (CID)
    cid_idx = trachomadata[0].index("CID")

    # Find the index of the column with the statuses of the agents at each timestep
    infected_idx = trachomadata[0].index("IndividualPCR")

    # Find the index of the column with the visit #s
    visit_idx = trachomadata[0].index("\ufeffvisit")

    # Find the index of the column with the ages
    age_idx = trachomadata[0].index("age")

    # Find the index of the column tracking ocular discharge presence
    ocular_idx = trachomadata[0].index("ocular")

    # Find the index of the column with TF (to stay consistent with my regressions, I will track TF by the data for the right eye: RTF)
    tf_idx = trachomadata[0].index("RTF")

    # Find the index of the column with TI (to stay consistent with my regressions, I will track TI by the data for the right eye: RTI)
    ti_idx = trachomadata[0].index("RTI")

    # Find the index of the column with the ages
    gender_idx = trachomadata[0].index("gender")

    # Find the index of the column tracking ocular discharge presence
    nasal_idx = trachomadata[0].index("nasal")

    # Find the index of the column tracking flies on face
    flies_idx = trachomadata[0].index("flies")

    # Loop through the spreadsheet (every row except the 0th row, which contains the labels)
    for row in trachomadata[1:]:
        # The village ID of a data entry (row in the spreadsheet) is the number in the VID column
        current_vid = int(row[vid_idx])

        # If the current village ID isn't in the dictionary of villages, add it to the dictionary with the key being that village ID
        # For any number of villages, the code will construct a master dictionary that sorts data based on village ID
        # Good compatibility
        if current_vid not in villages.keys():
            villages[current_vid] = Village(current_vid, {})
        current_village = villages[current_vid]
        current_visit_num = int(row[visit_idx])
        if current_visit_num not in current_village.visits.keys():
            current_village.visits[current_visit_num] = {}

        # Note: cid is now a string instead of an integer
        cid = str(row[cid_idx])

        # Limit the number of agents in village 587 to 100 for the sake of time
        if current_vid == 587:
            if len(current_village.unique_agents) <= 100:
                if cid not in current_village.unique_agents:
                    current_village.unique_agents.append(cid)
                current_village.visits[current_visit_num][cid] = Agent(cid, True if row[infected_idx] == "P" else False, row[age_idx], True if row[ocular_idx] == "1" else False, row[tf_idx], True if row[ti_idx] == "1" else False, row[gender_idx], row[nasal_idx], row[flies_idx])
                for visit in current_village.visits.values():
                    for agent in visit.values():
                        if agent.id not in current_village.first_agent_info:
                            current_village.first_agent_info[agent.id] = agent
            else:
                pass

        else:
            if cid not in current_village.unique_agents:
                current_village.unique_agents.append(cid)
            current_village.visits[current_visit_num][cid] = Agent(cid, True if row[infected_idx] == "P" else False, row[age_idx], True if row[ocular_idx] == "1" else False, row[tf_idx], True if row[ti_idx] == "1" else False, row[gender_idx], row[nasal_idx], row[flies_idx])
            for visit in current_village.visits.values():
                for agent in visit.values():
                    if agent.id not in current_village.first_agent_info:
                        current_village.first_agent_info[agent.id] = agent

def getVillageData(village_of_interest):
    with open(str(village_of_interest) + "_agent_states.csv", "w") as file:
        agent_states = {}
        for agent in villages[village_of_interest].unique_agents:
            # make all states question marks
            agent_states[agent] = [str(agent)] + ["?" for _ in villages[village_of_interest].visits]
        timestep = 1
        # overwrite the question marks (all data that we actually have will overwrite question marks, while the data that we don't have will remain question marks)
        for visit in villages[village_of_interest].visits.values():
            '''
            ### get the agents that weren't filled in
            ### which agents at this visit were not filled in?
            # below is the list of agents that were not recorded at this visit
            missing_agents = list(set(villages[village_of_interest].unique_agents) - set(visit.keys()))

            # prevalence among observed agents, as a fraction
            observed_prevalence = sum([1 for i in ]visit.values()) / len(visit.values())

            # make the list of infected status
            # prevalence: percent of population infected
            inf = round(observed_prevalence * len(missing_agents))
            inf_list = [1 for i in range(inf)]

            sus = len(missing_agents) - inf
            sus_list = [0 for i in range(sus)]

            # Now, we combine these to create our list of infected/susceptible, with one entry per agent.
            sim_states = inf_list + sus_list

            # Step 2: Shuffle this list.
            random.shuffle(sim_states)
            '''
            for agent in visit.values():
                agent_states[agent.id][timestep] = str(agent.infected)
            timestep += 1

        writer = csv.writer(file)
        print(list(agent_states.values()))
        writer.writerows(list(agent_states.values()))
