#! /usr/bin/env python
"""This script generates input JSON files, and a submission script for calibration runs"""

import os
import json

# The simulation will need a json input file with settings/file paths et cetera
amino_acid_template = "amino_acid_calibration.json"
ligand_template = "ligand_calibration.json"

# The simulation needs to be submitted
submit_template = "submit_calibration.sh"
# Read it as as string for formatting later
with open(submit_template) as submit_file:
    submit_str = submit_file.read()

# This script will submit everything after you're done. Use responsibly
# No template needed. Essentially its a list of bsub {scriptname} commands
master_submit = "master_submit.sh"
master_str = "# These jobs will be submitted \n"

# This is the script that will read the json file
# NOTE this does not get edited, it just reads the json.
run_script = "run_calibration.py"
runscript_path = os.path.abspath(run_script)
# Amino acids for which mmCIF input structures have been prepared
amino_acids = ["asp", "his", "glu", " lys", "tyr", "cys"]

# Small molecules
ligs = ["1D"]


# For every amino amino_acid, ligand
for system in amino_acids + ligs:
    # For every charge method (background/chenroux)
    for method in ("chenroux", "background"):
        jobname = "{}-{}-calibration".format(system, method)
        outdirname = "{}-results".format(jobname)
        # Make a new directory, ignore if exists already
        if not os.path.exists(outdirname):
            os.makedirs(outdirname)

        # Take a copy of the template json and fill it in
        if system in ligs:
            templatename = ligand_template
        else:
            templatename = amino_acid_template
        jsontemplate = json.load(open(templatename))
        jsontemplate["output"]["dir"] = outdirname
        jsontemplate["parameters"]["basename"] = system
        jsontemplate["parameters"]["methodname"] = method
        jsontemplate["parameters"]["ncmc"]["counterion_method"] = method
        os.chdir(outdirname)
        with open("settings.json", "w") as outfile:
            json.dump(jsontemplate, outfile, indent=4, sort_keys=True)
        setting_path = "{}/settings.json".format(outdirname)
        # Take a copy of the run submit template, and fill it in.
        with open("submit.sh", "w") as submitfile:
            submitfile.write(submit_str.format(jobname, runscript_path, setting_path))
        os.chdir("..")
        # add a line to the master submit script
        master_str += "bsub < {}/submit.sh\n".format(outdirname)

with open(master_submit, "w") as masterfile:
    masterfile.write(master_str)
