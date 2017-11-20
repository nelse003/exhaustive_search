import giant.xray.edstats as ed
import matplotlib.pyplot as plt

# Looping over multiple files

# Running Edstats
edstats, summary = ed.score_file_with_edstats("refine_1.mtz","refine_1.pdb")

#edstats.scores.to_csv("edstat_test.csv")

RSR = edstats.scores.loc['Ra']

# Splitting RSR score into required chain
RSR_chain = edstats.scores.loc['Ra',(slice(None),'A',slice(None),slice(None))]
ordered_chain  = RSR_chain.sort_index(level=2)

# Plotting of RSR vs residue for one chain
fig = plt.figure()
plt.plot(ordered_chain.index.get_level_values(2).values,ordered_chain.values)
plt.ylabel("RSR")
plt.xlabel("Residue")
plt.savefig("RSR_test")
plt.close()

# RSR of multiple xtals