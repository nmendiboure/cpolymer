import numpy as np
from cpolymer.polymer import Polymer
from cpolymer.lsimu import LSimu
from cpolymer.constrain import Box, Sphere

diameter = 0.030
Radius = 5 / diameter
Nchromosomes = 20
Schromosome = int(3000000 / 5000 / 20)
nucleus = Sphere(position=[0, 0, 0], radius=Radius)

Sim = LSimu()

bead_type = 1  # ALl the beads are going to be the same type
for X in range(Nchromosomes):
    bead_size = 1
    bond_type = 1
    Sim.add(Polymer(N=Schromosome, type_bead=1,
                    liaison={"{0}-{0}".format(bead_type): [bead_size, bond_type]},
                    gconstrain=[nucleus]))

Sim.add_bond(typeb="harmonic", idbond=bond_type, K=350, R0=bead_size)

Sim.add_pair(typep="lj/cut", idpair1=bead_type, idpair2=bead_type, epsilon=1, sigma=bead_size, cutoff1=1.15)

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


include $interaction





group particle type 1 

compute hic particle  pair/local dist
compute hicp particle property/local patom1 patom2



dump  init all dcd $samplingrate $outtraj.${tcell}.comp.dcd



###########################################################
#Definiton of nucleus and its interaction
#the telomere part is added when the nuceus has the right size

variable rad equal $radius

region mySphere sphere 0.0 0.0 0.0 v_rad side in

fix wall1 particle wall/region mySphere lj126 1 1 1.12 


#####################################################
# Equilibration (Langevin dynamics )

velocity 	particle create 1.0 1231
fix		1 particle nve/limit 0.0005
fix		lang particle langevin 1.0 1.0 1.0 904297
run 20000

unfix 1

thermo_style	custom step temp 
thermo          10000
fix		1 particle nve/limit 0.05
timestep	0.005 
run		$run_length


write_data $outfile
"""
with open("tscript","w") as f:
    f.writelines(Template)


REP = "./"
cell = "simple_yeast"
Sim.generate_xyz(REP + "/%sconf2.txt" % cell, Mass={str(bead_type): 1})
Sim.generate_pdb(REP + "/%snoyau2.pdb" % cell,
                 shift=1)  # Not necellary to run the simulation but usefull for the analysis
Sim.generate_interactions(REP + "/interactions")
Sim.generate_script(REP + "/nucleus_init.txt", template_name="./tscript", outfile="final.xyz",
                    outtraj="dump_init", samplingrate=1000, run_length=20000,
                    interaction="interactions", typecell=cell,
                    radius=Radius)

r = Sim.run("nucleus_init.txt")