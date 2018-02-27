from protons.app.ligands import extract_residue, prepare_calibration_system
from protons import app
import os 


for res, rsn in [('asp',"AS4"), ('glu',"GL4"), ('his',"HIS"), ('lys', 'LYS'), ('tyr', 'TYR'), ('cys', 'CYS')]:
    vacname = "{}.pdb".format(res)
    prepare_calibration_system(vacname, "{}-solvated.cif".format(res),delete_old_H=False)
