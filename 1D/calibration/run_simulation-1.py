from __future__ import print_function
from protons import app
from simtk import unit, openmm as mm
from protons.app import ConstantPHCalibration, ForceFieldProtonDrive, NCMCProtonDrive
from protons.app import MetadataReporter, TitrationReporter, NCMCReporter, SAMSReporter
import shutil
import signal
import json
from protons.app.logger import log, logging
import numpy as np
from openmmtools.integrators import LangevinIntegrator, ExternalPerturbationLangevinIntegrator
log.setLevel(logging.DEBUG)
import os
import netCDF4
import saltswap
# Define what to do on timeout
class TimeOutError(RuntimeError):
    """This error is raised when an operation is taking longer than expected."""
    pass

def timeout_handler(signum, frame):
   log.warn("Script is running out of time. Attempting to exit cleanly.")
   raise TimeOutError("Running out of time, shutting down!")

# Register the timeout handling
signal.signal(signal.SIGALRM, timeout_handler)
# Input files
input_pdb_file = "imatinib-solvated-for-calibration.pdb"
ligand_xml = "imatinib.xml"

# Load the PDB file and the forcefield files
pdb_object = app.PDBFile(input_pdb_file)
forcefield = app.ForceField('amber10-constph.xml', 'gaff.xml', ligand_xml, 'tip3p.xml', 'ions_tip3p.xml')

# Prepare the Simulation
topology = pdb_object.topology
positions = pdb_object.positions

# Quick fix for histidines in topology
for residue in topology.residues():
    if residue.name == 'HIS':
        residue.name = 'HIP'
# Naming the output files
basename = "Imatinib-solvent-1"
name_netcdf = '%s.nc' % basename
dcd_output_name = '%s.dcd' %basename
weights_txt_name = '%s-weights.txt' % basename
output_context_xml = "resume-%s-state.xml" % basename
output_drive_xml = "resume-%s-drive.xml" % basename
output_calibration_json = "resume-%s-calibration.json" % basename
# Integrator options
timestep = 2.0 * unit.femtosecond
constraint_tolerance = 1.e-7
collision_rate = 1.0 / unit.picosecond
number_R_steps = 1

class ExternalGBAOABIntegrator(ExternalPerturbationLangevinIntegrator):
    """
    Implementation of the gBAOAB integrator which tracks external protocol work.

    Parameters
    ----------
        number_R: int, default: 1
            The number of sequential R steps.  For instance V R R O R R V has number_R = 2
        temperature : simtk.unit.Quantity compatible with kelvin, default: 298*unit.kelvin
           The temperature.
        collision_rate : simtk.unit.Quantity compatible with 1/picoseconds, default: 1.0/unit.picoseconds
           The collision rate.
        timestep : simtk.unit.Quantity compatible with femtoseconds, default: 1.0*unit.femtoseconds
           The integration timestep.


    """

    def __init__(self, number_R_steps=1, temperature=298.0 * unit.kelvin,
                 collision_rate=1.0 / unit.picoseconds,
                 timestep=1.0 * unit.femtoseconds,
                 constraint_tolerance=1e-7
                 ):
        Rstep = " R" * number_R_steps

        super(ExternalGBAOABIntegrator, self).__init__(splitting="V{0} O{0} V".format(Rstep),
                                                       temperature=temperature,
                                                       collision_rate=collision_rate,
                                                       timestep=timestep,
                                                       constraint_tolerance=constraint_tolerance,
                                                       measure_shadow_work=False,
                                                       measure_heat=False,
                                                       )
# Steps of MD before starting the main loop
num_thermalization_steps = 10000
# Steps of MD in between MC moves
steps_between_updates = 10000
ncmc_steps_per_trial = 10000  # 20 ps / 20 fs
prop_steps_per_trial = 1
total_iterations = 1000
modulo_ligand_update = 1 # Update ligand every n iterations

pre_run_minimization_tolerance = 1e-5 * unit.kilojoule / unit.mole
minimization_max_iterations = 0
# SAMS settings
beta_burnin = 0.5
flatness_criterion = 0.2
# Script specific settings
script_timeout = 428400 # 119 hours

# Platform Options
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed', 'DeterministicForces': 'true', 'CudaDeviceIndex': os.environ['CUDA_VISIBLE_DEVICES']}
# System Configuration
nonbondedMethod = app.PME
constraints = app.HBonds
rigidWater = True
ewaldErrorTolerance = 1.e-5
barostatInterval =  25
switching_distance = 0.85 * unit.nanometers
nonbondedCutoff = 1.0 * unit.nanometers
pressure = 1.0 * unit.atmosphere
temperature = 300.0 * unit.kelvin

system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, constraints=constraints,
                                 rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance, nonbondedCutoff=nonbondedCutoff)

#
for force in system.getForces():
    if isinstance(force, mm.NonbondedForce):
        force.setUseSwitchingFunction(True)

        force.setSwitchingDistance(switching_distance)

# NPT simulation
system.addForce(
    mm.MonteCarloBarostat(
        pressure,
        temperature,
        barostatInterval))


integrator = ExternalGBAOABIntegrator(number_R_steps=number_R_steps, temperature=temperature, collision_rate=collision_rate, timestep=timestep, constraint_tolerance=constraint_tolerance)
ncmc_propagation_integrator = ExternalGBAOABIntegrator(number_R_steps=number_R_steps, temperature=temperature, collision_rate=collision_rate, timestep=timestep, constraint_tolerance=constraint_tolerance)

# Define a compound integrator
compound_integrator = mm.CompoundIntegrator()
compound_integrator.addIntegrator(integrator)
compound_integrator.addIntegrator(ncmc_propagation_integrator)
compound_integrator.setCurrentIntegrator(0)

driver = ForceFieldProtonDrive(temperature, topology, system, forcefield, ['amber10-constph.xml', ligand_xml], pressure=pressure,
                                       perturbations_per_trial=ncmc_steps_per_trial, propagations_per_step=prop_steps_per_trial)

# Initializing the weights using values from previous simulations
g_initial = {'TYR': [0.0, 126.7],
                 'AS4': [0.0, -63.2, -65.1, -63.1, -69.5],
                 'GL4': [0.0, -33.8, -39.7, -36.1, -38.5],
                 'HIP': [0.0, 27.5, 29.6],
                 'CYS': [0.0, 154.4],
                 'LYS': [0.0, -6.8]}
driver.import_gk_values(g_initial)
driver.adjust_to_ph(7.4)
# Assumes ligand is always the last titration group
ligand_titration_group_index = len(driver.titrationGroups) - 1

# Define residue pools
pools = {'ligand' : [ligand_titration_group_index]}
driver.define_pools(pools)
# Create SAMS sampler
simulation = app.ConstantPHCalibration(topology, system, compound_integrator, driver, group_index=ligand_titration_group_index, platform=platform, platformProperties=properties)
simulation.context.setPositions(positions)
simulation.minimizeEnergy(tolerance=pre_run_minimization_tolerance, maxIterations=minimization_max_iterations)
simulation.step(num_thermalization_steps)

dcdreporter = app.DCDReporter(dcd_output_name, int(steps_between_updates/30))
ncfile = netCDF4.Dataset(name_netcdf, "w")
metdatarep = MetadataReporter(ncfile, shared=True)
ncmcrep = NCMCReporter(ncfile,1,shared=True)
titrep = TitrationReporter(ncfile,1,shared=True)
simulation.reporters.append(dcdreporter)
simulation.update_reporters.append(metdatarep)
simulation.update_reporters.append(ncmcrep)
simulation.update_reporters.append(titrep)
samsrep = SAMSReporter(ncfile,1,shared=True)
simulation.calibration_reporters.append(samsrep)

# Raises an exception if the simulation runs out of time, so that the script can be killed cleanly from within python
signal.alarm(script_timeout)

try:
    for i in range(total_iterations):
        log.info("Iteration %i", i)
        if i == 5:
            log.info("Simulation seems to be working. Suppressing debugging info.")
            log.setLevel(logging.INFO)
        simulation.step(steps_between_updates)                
        simulation.update(1, pool='ligand')
        simulation.adapt()
    # Reset timer
    signal.alarm(0)

except TimeOutError:
    log.warn("Simulation ran out of time, saving current results.")

finally:
    # export the context
    serialize_state(simulation.context,output_context_xml)
    # export the driver
    serialize_drive(simulation.drive, output_drive_xml)
    # export the calibration status
    serialize_sams_status(simulation, output_calibration_json)

    ncfile.close()

# End of script
