# coding=utf-8
"""
Adds hydrogens to 2HYY-noH.pdb
"""
from protons import app
from simtk.openmm import  openmm
from simtk.unit import *
from protons.app.integrators import GBAOABIntegrator

# Load relevant template definitions for modeller, forcefield and topology
app.Modeller.loadHydrogenDefinitions('imatinib-hydrogens.xml')
forcefield = app.ForceField('amber10-constph.xml','gaff.xml', 'imatinib.xml', 'tip3p.xml', 'ions_tip3p.xml')

pdb = app.PDBFile('2HYY-noH.pdb')

modeller = app.Modeller(pdb.topology, pdb.positions)

# The pdb contains solvent but not the right ions.
# This would mean we need to equilibrate again even if we just add new ions.
# In this case its easiest to just delete and re-add the solvent with the right amount of ions

residues = [ resi for resi in modeller.topology.residues() if not resi.name in ['STI'] ]
modeller.delete(residues)


modeller.addHydrogens(forcefield=forcefield)
modeller.addSolvent(forcefield, model='tip3p', padding=1.0*nanometers, neutralize=False)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * nanometers,
                                 constraints=app.HBonds, rigidWater=True,
                                 ewaldErrorTolerance=0.0005)
system.addForce(openmm.MonteCarloBarostat(1.0 * atmosphere, 300.0 * kelvin))
simulation = app.Simulation(modeller.topology, system, GBAOABIntegrator())
simulation.context.setPositions(modeller.positions)
# simulation.minimizeEnergy()

app.PDBFile.writeFile(modeller.topology, simulation.context.getState(getPositions=True).getPositions(),
                       open('imatinib-solvated-for-calibration.pdb', 'w'))

