# coding=utf-8
"""
Adds hydrogens to 2HYY-noH.pdb
"""

from protons import app
from simtk.openmm import openmm
from simtk.unit import *
from protons.app.integrators import GBAOABIntegrator

# Load relevant template definitions for modeller, forcefield and topology

app.Modeller.loadHydrogenDefinitions('imatinib-hydrogens.xml')
forcefield = app.ForceField('amber10-constph.xml', 'gaff.xml', 'imatinib.xml', 'tip3p.xml', 'ions_tip3p.xml')

pdb = app.PDBFile('2HYY-noH.pdb')

modeller = app.Modeller(pdb.topology, pdb.positions)

# The pdb contains solvent but not the right ions.
# This would mean we need to equilibrate again even if we just add new ions.
# In this case its easiest to just delete and re-add the solvent with the right amount of ions

modeller.deleteWater()
ions = [ atom for atom in modeller.topology.atoms() if atom.element.symbol in ['Cl', 'Na'] ]
modeller.delete(ions)


modeller.addHydrogens(forcefield=forcefield)
modeller.addSolvent(forcefield, model='tip3p', padding=1.0*nanometers, positiveIon='Na+', negativeIon='Cl-', ionicStrength=120*millimolar, neutralize=True)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * nanometers,
                                 constraints=app.HBonds, rigidWater=True,
                                 ewaldErrorTolerance=0.0005)
system.addForce(openmm.MonteCarloBarostat(1.0 * atmosphere, 300.0 * kelvin))
simulation = app.Simulation(modeller.topology, system, GBAOABIntegrator() )
simulation.context.setPositions(modeller.positions)
rename_res = {"HIS": "HIP", "ASP":"AS4", "GLU":"GL4"}
for resi in modeller.topology.residues():
    if resi.name in rename_res:
        resi.name = rename_res[resi.name]

app.PDBFile.writeFile(modeller.topology, simulation.context.getState(getPositions=True).getPositions(),
                       open('2HYY-H.pdb', 'w'))

