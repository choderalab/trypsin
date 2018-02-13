# coding=utf-8
"""
This script parametrizes 1D
"""

from sys import argv
from protons.app.ligands import *
from protons.app import logger
from protons.app.logger import log
from lxml import etree
import json
import sys
import os

log.setLevel(logger.logging.DEBUG)


def main(args):
    """
    Run the program
    Parameters
    ----------
    args - cmd line arguments

    """
    import os

    if len(args) != 2:
        print("Please provide a single json input file.")
        sys.exit(1)

    with open(args[1].strip(), 'r') as settingsfile:
        settings = json.load(settingsfile)

    # Retrieve parameter fields
    prms = settings["parameters"]
    
    pH = float(prms["pH"])
    max_penalty = float(prms["max_penalty"])
    tautomerize = bool(prms["tautomerize"])
    resname = prms["pdb_resname"]

    #retrieve input fields
    inp = settings["input"]
    idir = inp["dir"].format(**prms)
    iepik = inp["epik"].format(**prms)
    icalib = inp["calibration"].format(**prms)
    iepik_path = os.path.abspath(os.path.join(idir, iepik))
    ical_path = os.path.abspath(os.path.join(idir, icalib))
    
    if not os.path.isfile(iepik_path):
        raise IOError("Could not find epik input at {}.".format(**locals()))

    # retrieve output field
    out = settings["output"]
    odir = out["dir"].format(**prms)
    offxml = out["ffxml"].format(**prms)
    oepik = out["epik"].format(**prms)
    ohxml = out["hxml"].format(**prms)
    ocalib = out["calibration_system"].format(**prms)
     
    if not os.path.isdir(odir):
        os.makedirs(odir)

    # Debugging/intermediate files
    debug = out["debug"]
    dhydrogen_fix = debug["fixed_hydrogens"].format(**prms)
    dext_res = debug["extracted_residue"].format(**prms)
    dkeep_files = bool(debug["keep_files"])

    lastdir=os.getcwd()    
    
    # Begin processing

    # run epik
    generate_epik_states(iepik_path, oepik, pH=pH, max_penalty=max_penalty, workdir=odir)
    os.chdir(odir)
    # process into mol2
    epik_results_to_mol2(oepik,dhydrogen_fix)
    # retrieve states
    isomer_info = retrieve_epik_info(oepik)

    # parametrize
    generate_protons_ffxml(dhydrogen_fix, isomer_info, offxml, pH, resname=resname)
    # create hydrogens
    create_hydrogen_definitions(offxml, ohxml)

    # set up calibration system
    extract_residue(ical_path,dext_res,resname=resname)
    
    # prepare solvated system
    prepare_calibration_system(dext_res, ocalib, offxml, ohxml)

    if not dkeep_files:
        os.remove(dext_res)
        os.remove(dhydrogen_fix)

    os.chdir(lastdir)

if __name__ == "__main__":
    main(argv)
