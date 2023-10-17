import numpy as np
from cpolymer.polymer import Polymer
from cpolymer.lsimu import LSimu
from cpolymer.constrain import Box, Sphere, Point

from cpolymer.halley.constrain import Spherical, Nowhere
from cpolymer.halley.vectors import V

len_chrom = [46, 162, 63, 306, 115, 54, 218, 112, 87, 149, 133, 365, 184, 156, 218, 189]
dist_centro = [30, 47, 22, 89, 30, 29, 99, 21, 71, 87, 88, 30, 53, 125, 65, 111]

Radius = 16.6
Mt = 0.3 * 16.6
Nchromosomes = 16

nucleus = Sphere(position=[0, 0, 0], radius=Radius)

Sim = LSimu()
dnuc = 3
bead_type = 1  # All the beads are going to be the same type
liaison = {"1-1": [1, 1], "1-2": [1, 2], "1-3": [1, 3], "1-4": [(dnuc + 1) / 2., 4], "1-5": [0, 5],
           "2-2": [1, 6], "2-3": [1, 7], "2-4": [(dnuc + 1) / 2., 8], "2-5": [0, 9],
           "3-3": [1, 10], "3-4": [(dnuc + 1) / 2., 11], "3-5": [Mt, 12],
           "4-4": [dnuc, 13], "4-5": [0, 14],
           "5-5": [0, 15]}

for X in range(Nchromosomes):
    # We need to define geometrical constrain:
    # The centromere must be at a distance mt of the spb positioned at (-Radius,0,0)
    # We use the module halley for that
    S_spb = Spherical(V(-Radius, 0, 0), radius=Mt)  # we define a sphere centered on the spb
    S_centered = Spherical(V(0, 0, 0), radius=Radius - Mt * 0.9)  # a sphere centered that intersect the spb sphere

    circle = S_spb * S_centered
    centromere = circle.get_random()

    # We must then construct a sphere centered on the centromere with a radius sqrt(Nbead)
    # and look at its intersection with the nucleus
    Nucleus = Spherical(V(0, 0, 0), radius=Radius * 0.95)  # a sphere centered that intersect the spb sphere
    d1 = dist_centro[X]
    Telo1_possible = Spherical(centromere, radius=np.sqrt(d1)) * Nucleus
    telo1 = Telo1_possible.get_random()

    d2 = len_chrom[X] - dist_centro[X]
    Telo2_possible = Spherical(centromere, radius=np.sqrt(d2)) * Nucleus
    telo2 = Telo1_possible.get_random()

    if X != 11:
        Sim.add(Polymer(N=len_chrom[X], type_bead=[2] + [1] * (d1 - 2) + [3] + [1] * (d2 - 1) + [2],
                        liaison=liaison,
                        gconstrain=[nucleus],
                        lconstrain=[Point(index=0, position=telo1._v),
                                    Point(index=d1, position=centromere._v),
                                    Point(index=len_chrom[X] - 1, position=telo2._v)]))
    else:
        # This chromosome is different because it has a nucleole
        Sim.add(Polymer(N=len_chrom[X], type_bead=[2] + [1] * (d1 - 2) + [3] + [1] * (90 - d1) + [4] * 150 + \
                                                  [1] * (len_chrom[X] - 150 - 90 - 1) + [2],
                        liaison=liaison,
                        gconstrain=[nucleus],
                        lconstrain=[Point(index=0, position=telo1._v),
                                    Point(index=d1, position=centromere._v),
                                    Point(index=90 + 75, position=(0.66 * Radius, 0, 0)),
                                    # We add a new constrain: the center
                                    # of the nucleole must be at 2/3 of
                                    # the radius at the
                                    # opposite of the spb:
                                    Point(index=len_chrom[X] - 1, position=telo2._v)]))

# Then Add the spb


Sim.add(Polymer(N=1, type_bead=5, liaison=liaison))
Sim.molecules[-1].coords = np.array([[-Radius, 0, 0]])

for i, c in enumerate(dist_centro, 1):
    if c is not None:
        Sim.add_extra_bond(mol1=[len(Sim.molecules), 1], mol2=[i, c], typeb=liaison["3-5"][1])

simsoft = LSimu()
print
liaison
for k, (bond_size, bond_type) in liaison.items():
    simsoft.add_bond(typeb="harmonic", idbond=bond_type, K=350, R0=bond_size)

    if bond_size != 0:
        Sim.add_bond(typeb="fene", idbond=bond_type, K=30. / (bond_size * bond_size),
                     R0=1.5 * bond_size, epsilon=1, sigma=bond_size)
    else:
        Sim.add_bond(typeb="harmonic", idbond=bond_type, K=0, R0=0)

    bead1, bead2 = map(int, k.split("-"))
    if bead1 == 5 or bead2 == 5:
        simsoft.add_pair(typep="soft", idpair1=bead1, idpair2=bead2, A=0, rc=bond_size)

        Sim.add_pair(typep="lj/cut", idpair1=bead1, idpair2=bead2, epsilon=0, sigma=bond_size, cutoff1=bond_size * 1.15)
    else:
        simsoft.add_pair(typep="soft", idpair1=bead1, idpair2=bead2, A=5, rc=bond_size)

        Sim.add_pair(typep="lj/cut", idpair1=bead1, idpair2=bead2, epsilon=1, sigma=bond_size, cutoff1=bond_size * 1.15)

Sim.add_box(Box([-1.1 * Radius, -1.1 * Radius, -1.1 * Radius], [1.1 * Radius, 1.1 * Radius, 1.1 * Radius]))

#Now we have to create a template for the simulation.
Template="""
################################
#Template
#Must contain the variables
#  typecell 
#  outtraj
#  outfile
#  radius
#  interaction
#  run_length
#  samplingrate
#  particle


# VARIABLES
variable tcell index $typecell    # Define the variable from the tempalte
variable fname index ${tcell}conf2.txt    # configuration initiale


# Initialization
#correspond to x=y=z=1
lattice fcc 4
units		lj
boundary	f f f
atom_style	molecular
log 		log.txt
read_data	${fname}



neighbor 2.0 multi


include $softinteractions



group particle type 1 2 3 4
group normal type 1 3
group telo type 2
group ribo type 4

compute hic particle  pair/local dist
compute hicp particle property/local patom1 patom2



dump  init all dcd $samplingrate $outtraj.${tcell}.comp.dcd



###########################################################
#Definiton of nucleus and its interaction
#the telomere part is added when the nuceus has the right size

variable rad equal $radius

region mySphere sphere 0.0 0.0 0.0 v_rad side in

fix wall1 normal wall/region mySphere lj126 1 0.5 0.56 
fix wall2 ribo wall/region mySphere lj126 1 1.62 1.81 
fix wall telo wall/region mySphere  lj93  4 2 6


#####################################################
# Equilibration (Langevin dynamics )

velocity 	particle create 1.0 1231
fix		1 particle nve/limit 0.0005
fix		lang particle langevin 1.0 1.0 1.0 904297
run 30000

unfix 1
fix		1 particle nve/limit 0.005
run 30000

unfix 1
fix		1 particle nve/limit 0.05
run 30000

include $interaction
thermo_style	custom step temp 
thermo          10000
fix		1 particle nve/limit 0.05
timestep	0.005 
run		$run_length


write_data $outfile
"""
with open("tscript","w") as f:
    f.writelines(Template)

# Then let's generate all the files needed to run the simulation:
REP = "./"
cell = "nucleus_yeast"
Sim.generate_xyz(REP + "/%sconf2.txt" % cell, Mass={"1": 1, "2": 1, "3": 1, "4": 1, "5": 1})
Sim.generate_pdb(REP + "/%snoyau2.pdb" % cell,
                 shift=1)  # Not necellary to run the simulation but usefull for the analysis
Sim.generate_interactions(REP + "/interactions", info_bond=["special_bonds fene"])
simsoft.generate_interactions(REP + "/softinteractions")

Sim.generate_script(REP + "/nucleus_init.txt", template_name="./tscript", outfile="final.xyz",
                    outtraj="dump_init", samplingrate=10000, run_length=1000000,
                    interaction="interactions",
                    softinteractions="softinteractions", typecell=cell,
                    radius=Radius)

r = Sim.run("nucleus_init.txt")
print(r)
