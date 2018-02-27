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
import sys
import netCDF4
from saltswap.wrappers import Salinator
import json


# Define what to do on timeout
class TimeOutError(RuntimeError):
    """This error is raised when an operation is taking longer than expected."""
    pass


def timeout_handler(signum, frame):
    """Handle a timeout."""
    log.warn("Script is running out of time. Attempting to exit cleanly.")
    raise TimeOutError("Running out of time, shutting down!")

def serialize_state(context, outputfile):
    """
    Serialize the simulation state to xml.
    """
    xmls = mm.XmlSerializer
    statexml = xmls.serialize(context.getState(getPositions=True, getVelocities=True))
    with open(outputfile,'w') as statefile:
        statefile.write(statexml)

def serialize_drive(drive, outputfile):
    """
    Serialize the drive residues to xml.
    """
    drivexml = drive.serialize_titration_groups()
    with open(outputfile, 'wb') as drivefile:
        drivefile.write(drivexml)

def serialize_sams_status(calibration, outputfile):
    """
    Serialize the state of the SAMS calibration as json
    """
    samsProperties = calibration.export_samsProperties()
    samsjson = json.dumps(samsProperties, sort_keys=True, indent=4, separators=(',', ': '))
    with open(outputfile, 'w') as samsfile:
        samsfile.write(samsjson)

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


def main(jsonfile):
    """Main simulation loop."""

    settings = json.load(open(jsonfile))
    prms = settings["parameters"]
    # Register the timeout handling
    signal.signal(signal.SIGALRM, timeout_handler)
    # Input files
    inp = settings["input"]
    idir = inp["dir"].format(**prms)
    input_pdbx_file = os.path.join(idir, inp["calibration_system"].format(**prms))
    custom_xml_provided = False
    if "ffxml" in inp:
        custom_xml_provided = True

    if custom_xml_provided:
        custom_xml = os.path.join(idir, inp["ffxml"].format(**prms))
        custom_xml = os.path.abspath(custom_xml)
        forcefield = app.ForceField('amber10-constph.xml', 'gaff.xml', custom_xml, 'tip3p.xml', 'ions_tip3p.xml')
    else:
        forcefield = app.ForceField('amber10-constph.xml', 'gaff.xml', 'tip3p.xml', 'ions_tip3p.xml')

    # Load the PDBxfile and the forcefield files
    pdb_object = app.PDBxFile(input_pdbx_file)

    # Prepare the Simulation
    topology = pdb_object.topology
    positions = pdb_object.positions

    # Quick fix for histidines in topology
    # Openmm relabels them HIS, which leads to them not being detected as
    # titratable. Renaming them fixes this.
    for residue in topology.residues():
        if residue.name == 'HIS':
            residue.name = 'HIP'

    # Naming the output files
    out = settings["output"]
    odir = out["dir"].format(**prms)

    if not os.path.isdir(odir):
        os.makedirs(odir)
    lastdir = os.getcwd()
    os.chdir(odir)

    name_netcdf = out["netcdf"].format(**prms)
    dcd_output_name = out["dcd"].format(**prms)

    # Files for resuming simulation
    resumes = out["resume_files"]

    output_context_xml = resumes["state"].format(**prms)
    output_drive_xml = resumes["drive"].format(**prms)
    output_calibration_json = resumes["calibration"].format(**prms)

    # Integrator options
    integrator_opts = prms["integrator"]
    timestep = integrator_opts["timestep_fs"] * unit.femtosecond
    constraint_tolerance = integrator_opts["constraint_tolerance"]
    collision_rate = integrator_opts["collision_rate_per_ps"] / unit.picosecond
    number_R_steps = 1

    # Steps of MD before starting the main loop
    num_thermalization_steps = int(prms["num_thermalization_steps"])
    # Steps of MD in between MC moves
    steps_between_updates = int(prms["steps_between_updates"])

    ncmc = prms["ncmc"]
    counterion_method = ncmc["counterion_method"].lower()
    if counterion_method not in ["chen-roux", "chenroux", "background"]:
        raise ValueError("Invalid ncmc counterion method, {}. Please pick Chen-Roux or background.".format(counterion_method))
    ncmc_steps_per_trial = int(ncmc["steps_per_trial"])
    prop_steps_per_trial = int(ncmc["propagations_per_step"])
    total_iterations = int(prms["total_attempts"])

    modulo_ligand_update = 1 # Update ligand every n iterations

    # settings form minimization
    minimization = prms["minimization"]
    pre_run_minimization_tolerance = float(minimization["tolerance_kjmol"]) * unit.kilojoule / unit.mole
    minimization_max_iterations = int(minimization["max_iterations"])

    # SAMS settings
    sams = prms["SAMS"]
    beta_burnin = sams["beta_sams"]
    flatness_criterion = sams["flatness_criterion"]

    # Script specific settings
    script_timeout = 428400 # 119 hours

    # Platform Options

    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed', 'DeterministicForces': 'true', 'CudaDeviceIndex': os.environ['CUDA_VISIBLE_DEVICES']}
    # System Configuration
    sysprops = prms["system"]
    nonbondedMethod = app.PME
    constraints = app.HBonds
    rigidWater = True
    ewaldErrorTolerance = float(sysprops["ewald_error_tolerance"])
    barostatInterval =  int(sysprops["barostat_interval"])
    switching_distance = float(sysprops["switching_distance_nm"]) * unit.nanometers
    nonbondedCutoff = float(sysprops["nonbonded_cutoff_nm"]) * unit.nanometers
    pressure = float(sysprops["pressure_atm"]) * unit.atmosphere
    temperature = float(sysprops["temperature_k"]) * unit.kelvin

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

    if custom_xml_provided:
        driver = ForceFieldProtonDrive(temperature, topology, system, forcefield, ['amber10-constph.xml', custom_xml], pressure=pressure,
                                           perturbations_per_trial=ncmc_steps_per_trial, propagations_per_step=prop_steps_per_trial)
    else:
        driver = ForceFieldProtonDrive(temperature, topology, system, forcefield, ['amber10-constph.xml'], pressure=pressure,
                                           perturbations_per_trial=ncmc_steps_per_trial, propagations_per_step=prop_steps_per_trial)
    
    # Assumes ligand is always the last titration group
    ligand_titration_group_index = len(driver.titrationGroups) - 1

    # Define residue pools
    pools = {'ligand' : [ligand_titration_group_index]}
    driver.define_pools(pools)
    # Create SAMS sampler
    simulation = app.ConstantPHCalibration(topology, system, compound_integrator, driver, group_index=ligand_titration_group_index, platform=platform, platformProperties=properties, samsProperties=sams)
    simulation.context.setPositions(positions)
    salinator = Salinator(context=simulation.context,
                            system=system,
                            topology=topology,
                            ncmc_integrator=compound_integrator.getIntegrator(1),
                            salt_concentration=0.150 * unit.molar,
                            pressure=pressure,
                            temperature=temperature)
    salinator.neutralize()
    salinator.initialize_concentration()
    swapper = salinator.swapper
    # if chen-roux required, attach swapper. Else, use neutralizing background charge
    if counterion_method in ["chenroux", "chen-roux"]:
        simulation.drive.attach_swapper(swapper)

    simulation.minimizeEnergy(tolerance=pre_run_minimization_tolerance, maxIterations=minimization_max_iterations)
    simulation.step(num_thermalization_steps)

    dcdreporter = app.DCDReporter(dcd_output_name, int(steps_between_updates/10))
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
        os.chdir(lastdir)

    # End of script


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Please provide a single json file as input.")
    else:
        # Provide the json file to main function
        main(sys.argv[1])
