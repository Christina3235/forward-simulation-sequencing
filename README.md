# simulations and real data statistics
#%%
#trying to take samples in different time periods (in generations)-MANUALLY
import msprime
import tskit
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import SVG
rr=1e-6
demography = msprime.Demography()
demography.add_population(name="A", initial_size=2000)
demography.add_population_parameters_change(time=60, population="A", initial_size=100)
demography.add_population_parameters_change(time=100, population="A", initial_size=2000)
ts = msprime.sim_ancestry(samples= [
     msprime.SampleSet(5,time=90),
     msprime.SampleSet(5,time=80),
     msprime.SampleSet(5,time=70)], demography=demography, recombination_rate= rr, sequence_length= 1000, random_seed=123)
mts = msprime.sim_mutations(ts, rate=1e-5, random_seed=5678)
mts.draw_svg(y_axis=True, size=(1000, 1000))
#%%

#%%
#Sampling in every 10 generations from a list of samplesets
import msprime
import tskit
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import SVG
rr=1e-6
Samples=[]
for i in range (60,100,10):
    Samples.append(msprime.SampleSet(5,time=i))
    
demography = msprime.Demography()
demography.add_population(name="A", initial_size=2000)
demography.add_population_parameters_change(time=60, population="A", initial_size=100)
demography.add_population_parameters_change(time=100, population="A", initial_size=2000)
ts = msprime.sim_ancestry(samples= Samples, demography=demography, recombination_rate= rr, sequence_length= 1000, random_seed=123)
mts = msprime.sim_mutations(ts, rate=1e-5, random_seed=5678)
mts.draw_svg(y_axis=True, size=(2000, 2000))
#%%


#%%
#Exponetnial growth
#sampling every 10 generations for 40 generations
import msprime
import tskit
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import SVG
rr=1e-6
Samples=[]
for i in range (60,200,10):
    Samples.append(msprime.SampleSet(5,time=i))
demography = msprime.Demography()
demography.add_population(name="A", initial_size=100)
demography.add_population_parameters_change(time=100, growth_rate=0.01, population="A")
ts = msprime.sim_ancestry(samples=Samples, demography=demography, recombination_rate= rr, sequence_length= 10000, random_seed=123)
mts = msprime.sim_mutations(ts, rate=1e-6, random_seed=5678)
mts.draw_svg(y_axis=True, size=(6000, 6000))
#%%

#%%
#small initial size, exponential growth, computing statistics
#something problemaic with the growth rate
import msprime
import tskit
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import SVG
rr=1e-6
Samples=[]
for i in range (60,200,10):
    Samples.append(msprime.SampleSet(5,time=i,population='A'))
demography = msprime.Demography()
demography.add_population(name="A", initial_size=3000)
demography.add_population_parameters_change(time=100, growth_rate=1, population="A")
ts = msprime.sim_ancestry(samples=Samples, demography=demography, recombination_rate= rr, sequence_length= 10000, random_seed=123)
mts = msprime.sim_mutations(ts, rate=1e-6, random_seed=5678)
mts.draw_svg(y_axis=True, size=(6000, 6000))
print(Samples)
#%%

#%%
print(Samples)

#%%
#Exponetnial growth
#sampling every 10 generations for 40 generations
import msprime
import tskit
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import SVG
rr=1e-6
Samples=[]
for i in range (60,200,10):
    Samples.append(msprime.SampleSet(5,time=i))
demography = msprime.Demography()
demography.add_population(name="A", initial_size=100)
demography.add_population_parameters_change(time=100, growth_rate=0.01, population="A")
ts = msprime.sim_ancestry(samples=Samples, demography=demography, recombination_rate= rr, sequence_length= 10000, random_seed=123)
mts = msprime.sim_mutations(ts, rate=1e-6, random_seed=5678)
mts.draw_svg(y_axis=True, size=(6000, 6000))
#%%

#%%
print(mts.tables.mutations)
print(mts.tables.nodes)
print(mts)
print(mts.pairwise_diversity())
print(mts.diversity())
print(mts.Tajimas_D())
SVG(mts.first().draw_svg())
tree=mts.first()
print(tree.newick(precision=3))
#%%

#%%
#haplotype computation

import numpy as np
import msprime

tables = mts.dump_tables()

missing_sites = np.setdiff1d(np.arange(tables.sequence_length), tables.sites.position)
for pos in missing_sites:
    tables.sites.add_row(position=pos, ancestral_state="A")  # Add sites at every pos
tables.sort()
missing_ts = tables.tree_sequence()
SVG(missing_ts.draw_svg())
#%%

#%%
for i, h in enumerate(missing_ts.haplotypes()):
    print(f"Sample {i}: {h}")
#%%

#%%
import numpy as np
np.set_printoptions(linewidth=200)
print("Genotypes")
for v in mts.variants():
    print(f"Site {v.site.id}: {v.genotypes}")
    if v.site.id>=4:
        print("...")
        break
#%%

#%%
#use natural logarithm fubction to compute the rate ofthe exponential growth
import math 
import msprime
import tskit
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import SVG
rr=1e-6
a=(math.log((30000/5)/100))
Samples=[]
for i in range (60,200,10):
    Samples.append(msprime.SampleSet(5,time=i))
demography = msprime.Demography()
demography.add_population(name="A", initial_size=30000)
demography.add_population_parameters_change(time=100, growth_rate=a, population="A")
ts = msprime.sim_ancestry(samples=Samples, demography=demography, recombination_rate= rr, sequence_length= 10000, random_seed=123)
mts = msprime.sim_mutations(ts, rate=1e-6, random_seed=5678)
mts.draw_svg(y_axis=True, size=(6000, 6000))
#%%

#%%
#I tried to solve the equation to find the growth rate if we started with a population of 5 individuals
#and ended up with a population of 30000 individuals
import math 
import msprime
import tskit
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import SVG
rr=1e-6

Samples=[]
for i in range (60,200,10):
    Samples.append(msprime.SampleSet(5,time=i))
demography = msprime.Demography()
demography.add_population(name="A", initial_size=30000)
demography.add_population_parameters_change(time=100, growth_rate=0.25, population="A")
ts = msprime.sim_ancestry(samples=Samples, demography=demography, recombination_rate= rr, sequence_length= 10000, random_seed=123)
mts = msprime.sim_mutations(ts, rate=1e-6, random_seed=5678)
mts.draw_svg(y_axis=True, size=(6000, 6000))

#%%

#%%
print(mts.pairwise_diversity())
print(mts.diversity())
print(mts.Tajimas_D())
#%%

#%%
import numpy as np
import matplotlib



mutat_color_dict = {}
cmap = matplotlib.cm.get_cmap("viridis")
for i in range(30000):
    rgba = cmap(i / 30000)
    mutat_color_dict[i] = np.array(rgba) * 255

edge_colors = {}
string = ".tree .edge {stroke-width:4px}"
pop = mts.tables.nodes.population
for node in mts.nodes():
    node = node.id
    edge_colors[node] = mutat_color_dict[pop[node]]
    string += (
        ".node.n"
        + str(node)
        + "> .sym {fill: rgb"
        + str(tuple(mutat_color_dict[pop[node]][0:3]))
        + "}"
    )
    string += (
        ".node.n"
        + str(node)
        + " > .edge {stroke: rgb"
        + str(tuple(mutat_color_dict[pop[node]][0:3]))
        + "}"
    )
for mutat in mts.migrations():
    edge_colors[mutat.node] = mutat_color_dict[mutat.source]
    string += (
        ".node.n"
        + str(mutat.node)
        + "> .sym {fill: rgb"
        + str(tuple(mutat_color_dict[pop[mutat.node]][0:3]))
        + "}"
    )
    string += (
        ".node.n"
        + str(m.node)
        + " > .edge {stroke: rgb"
        + str(tuple(mutat_color_dict[mutat.dest][0:3]))
        + "}"
    )
node_labels = {node.id: "" for node in mts.nodes()}
SVG(
    mts.first().draw_svg(
        size=(500, 400), style=string, node_labels=node_labels, symbol_size=10
    )
)

#%%

#%%
#trying to make a graph with Tajima's D values from samples of each time
import math 
import msprime
import tskit
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import SVG
rr=1e-6

gen_start=100
gen_end=200
size_start=30000
size_end=5
growth_rate=-math.log((size_end/size_start), math.e)/(gen_end-gen_start)
growth_rate

Samples=[]
for i in range (60,200,10):
    Samples.append(msprime.SampleSet(5,time=i))
demography = msprime.Demography()
demography.add_population(name="A", initial_size=30000)
demography.add_population_parameters_change(time=100, growth_rate=growth_rate, population="A")
demography.add_population_parameters_change(time=200,growth_rate=0, population='A')
ts = msprime.sim_ancestry(samples=Samples, demography=demography, recombination_rate= rr, sequence_length= 10000, random_seed=123)
mts = msprime.sim_mutations(ts, rate=1e-6, random_seed=5678)
mts.draw_svg(y_axis=True, size=(6000, 6000))

#If it is a list of lists of samples we return an array for each window in the output, 
#which contains the value of the statistic separately for each of sample_sets in the order they are given.



print(mts.samples(time = 150))


print(mts.Tajimas_D(mts.samples(time = 60)))
print(mts.pairwise_diversity(mts.samples(time = 60)))


print(mts.tables)

#%%

#%%
print(Samples)
mts.genotype_matrix()
#%%
mts.genotype_matrix()[0:, mts.samples(time=60)]
#%%

#%%
tajimas=[]
for i in range (60,190,10):
    tajimas.append(mts.Tajimas_D(mts.samples(time=i)).flatten()[0])
print(tajimas)
#%%

#%%
import matplotlib.pyplot as plt
x=range(60,190,10)
y=tajimas

mts.Tajimas_D(mts.samples(time=i)).flatten()[0]
plt.xlabel("times of sampling in generations")
plt.ylabel("Tajima's D")
plt.scatter(list(x),y)
print(plt.show())
#%%

#%%
#remove empty values from the list tajimas
while ("" in tajimas):
    tajimas.remove("")
print(tajimas)
#%%

#%% 
pairwise_div=[]
for i in range (60,190,10):
    pairwise_div.append(mts.pairwise_diversity(mts.samples(time=i)))
print(pairwise_div)
#%%

#%%
while ("null" in pairwise_div):
    pairwise_div.remove("null")
print(pairwise_div)
#%%

#%%
import matplotlib.pyplot as plt
x=range(60,190,10)
y=pairwise_div
plt.xlabel("times of sampling in generations")
plt.ylabel("Pairwise diversity")


plt.scatter(list(x),y)
print(plt.show())
#%%
#Tajima's D
#Tajima's D is a population genetic test statistic 
#created by and named after the Japanese researcher Fumio Tajima.[1] 
#Tajima's D is computed as the difference between two measures of genetic diversity: 
#the mean number of pairwise differences and the number of segregating sites, 
#each scaled so that they are expected to be the same in a neutrally evolving population of constant size.


#Allele frequency spectrum
#Στην πληθυσμιακή γενετική, το φάσμα συχνότητας αλληλόμορφων, 
#που μερικές φορές ονομάζεται φάσμα συχνότητας θέσης, 
#είναι η κατανομή των συχνοτήτων αλληλόμορφων ενός δεδομένου συνόλου τόπων σε έναν πληθυσμό ή δείγμα















