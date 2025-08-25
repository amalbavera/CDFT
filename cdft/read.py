#############################################################################################################################

##    ## ##      ##  ######  ##     ## ######## ##     ## 
###   ## ##  ##  ## ##    ## ##     ## ##       ###   ### 
####  ## ##  ##  ## ##       ##     ## ##       #### #### 
## ## ## ##  ##  ## ##       ######### ######   ## ### ## 
##  #### ##  ##  ## ##       ##     ## ##       ##     ## 
##   ### ##  ##  ## ##    ## ##     ## ##       ##     ## 
##    ##  ###  ###   ######  ##     ## ######## ##     ##

#############################################################################################################################

import os

import numpy as np

from pathlib import Path as path

#############################################################################################################################
### Output serach strings for NWCHEM files

nwchem_total_energy_str:str   = "Total DFT energy ="

nwchem_orb_str:str            = "DFT Final Molecular Orbital Analysis"
nwchem_up_orb_str:str         = "DFT Final Alpha Molecular Orbital Analysis"
nwchem_dw_orb_str:str         = "DFT Final Beta Molecular Orbital Analysis"
nwchem_occ_str:str            = "Occ="

nwchem_geom_str:str           = "XYZ format geometry"

au_to_ev = 27.211399

#############################################################################################################################

class nwchem:

    __slots__ = ("atoms", "energy", "elements", "homo", "lumo", "convert")
#
### Initialize
#
    def __init__(self, outfile:str|path, convert:float=au_to_ev) -> None:
        
        self.atoms    = None
        self.homo     = None
        self.lumo     = None
        self.energy   = None
        self.elements = None
        self.convert  = convert
        
        if isinstance(outfile,str): outfile = path(outfile)
        
        fileio(outfile, option="r")
        
        self._get_output_data(outfile)
#
### Extract all data from output file
#
    def _get_output_data(self, data) -> None:
        
        skipped_lines           = 0

        alpha_orbital_energy    = []
        beta_orbital_energy     = []
        alpha_orbital_occup     = []
        beta_orbital_occup      = []

        elements                = []
        #
        ### Runtime variables
        #
        found_spin_up           = False
        found_spin_dw           = False

        found_geometry          = False
        skip_geom_lines         = 3

        with open(data, "r") as out:
            for line in out:
                line = line.strip()
                #
                ### Elements
                #
                if nwchem_geom_str in line:
                    skipped_lines  = 0
                    found_geometry = True
                    continue

                if   (found_geometry) and (skipped_lines < skip_geom_lines):
                    skipped_lines += 1
                    continue

                if   (found_geometry) and (not line):
                    skipped_lines  = 0
                    found_geometry = False
                    continue
                    
                elif (found_geometry) and (skipped_lines == skip_geom_lines):
                    elements += [line.split()[0]]
                #
                ### Total Energy
                #
                if nwchem_total_energy_str in line: self.energy = float(line.split()[4])*self.convert
                #
                ### Orbital Energies and Occupations
                #
                if (nwchem_up_orb_str in line) or (nwchem_orb_str in line):
                    found_spin_up = True
                    found_spin_dw = False
                    continue

                if nwchem_dw_orb_str in line:
                    found_spin_up = False
                    found_spin_dw = True
                    continue

                if nwchem_occ_str in line:
                    data_line       = line.split()
                    negative_energy = len(data_line) < 5

                    occupation_data, energy_data = data_line[2:] if negative_energy else data_line[2::2]

                    energy     = float(energy_data.partition("=")[2].replace("D","e")) if negative_energy else float(energy_data.replace("D","e"))
                    occupation = float(occupation_data.partition("=")[2].replace("D","e"))

                    if   found_spin_up:
                        alpha_orbital_energy += [energy]
                        alpha_orbital_occup  += [occupation]
                        continue
                        
                    elif found_spin_dw:
                        beta_orbital_energy  += [energy]
                        beta_orbital_occup   += [occupation]
                        continue

        alpha_orbital_energy = np.asarray(alpha_orbital_energy, dtype=float)*self.convert
        beta_orbital_energy  = np.asarray(beta_orbital_energy,  dtype=float)*self.convert
        
        alpha_orbital_occup  = np.asarray(alpha_orbital_occup,  dtype=float)
        beta_orbital_occup   = np.asarray(beta_orbital_occup,   dtype=float)

        if beta_orbital_occup.size:
            self.homo = max(alpha_orbital_energy[alpha_orbital_occup==1.0][-1], beta_orbital_energy[beta_orbital_occup==1.0][-1])
            self.lumo = min(alpha_orbital_energy[alpha_orbital_occup==0.0][ 0], beta_orbital_energy[beta_orbital_occup==0.0][ 0])

        else:
            self.homo = alpha_orbital_energy[alpha_orbital_occup==2.0][-1]
            self.lumo = alpha_orbital_energy[alpha_orbital_occup==0.0][ 0]
    
        self.elements = np.asarray(elements,   dtype=str)
        self.atoms    = np.ones_like(elements, dtype=int)

#############################################################################################################################

##     ##  #######  ########  ##     ## ##       ########  ######
###   ### ##     ## ##     ## ##     ## ##       ##       ##    ##
#### #### ##     ## ##     ## ##     ## ##       ##       ##
## ### ## ##     ## ##     ## ##     ## ##       ######    ######
##     ## ##     ## ##     ## ##     ## ##       ##             ##
##     ## ##     ## ##     ## ##     ## ##       ##       ##    ##
##     ##  #######  ########   #######  ######## ########  ######

#############################################################################################################################
### IO operations

def fileio(filename=None, option="read"):

    if option in ("w", "write"):

        try:
            return True if os.path.getsize(filename) == 0 else False

        except OSError:
            return True

    elif option in ("r","s","read","silent"):

        try:
            if os.path.getsize(filename) == 0:
                if option in ("r","read"):
                    raise RuntimeError(f"{filename} is empty")
                return False
            else:
                return True

        except OSError:
            if option in ("r","read"):
                raise FileNotFoundError(f"{filename} not found")
            return False

#############################################################################################################################
