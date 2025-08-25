#############################################################################################################################

##       #### ########  ########     ###    ########  #### ########  ######
##        ##  ##     ## ##     ##   ## ##   ##     ##  ##  ##       ##    ##
##        ##  ##     ## ##     ##  ##   ##  ##     ##  ##  ##       ##
##        ##  ########  ########  ##     ## ########   ##  ######    ######
##        ##  ##     ## ##   ##   ######### ##   ##    ##  ##             ##
##        ##  ##     ## ##    ##  ##     ## ##    ##   ##  ##       ##    ##
######## #### ########  ##     ## ##     ## ##     ## #### ########  ######

#############################################################################################################################

#
### Local libraries
#
from .model import two_parabola as tpm
from .model import generalized  as gqm

#############################################################################################################################

   ##   ########  ##     ## 
 ####   ##     ## ###   ### 
   ##   ##     ## #### #### 
   ##   ########  ## ### ## 
   ##   ##        ##     ## 
   ##   ##        ##     ## 
 ###### ##        ##     ##

#############################################################################################################################

def one_parabola(cation=None, neutral=None, anion=None, ref_cation=None, ref_neutral=None, ref_anion=None):
    #
    ### Reference species (B)
    #
    ionization_potential_ref = ref_cation.energy  - ref_neutral.energy
    electron_affinity_ref    = ref_neutral.energy - ref_anion.energy
    
    hardness_ref             = ionization_potential_ref - electron_affinity_ref
    #
    ### Species of interest (A)
    #
    ionization_potential = cation.energy  - neutral.energy
    electron_affinity    = neutral.energy - anion.energy
    
    hardness             = ionization_potential - electron_affinity        
    #
    ### Charge transfer
    #
    reciprocal_hardness          = 1.0/(hardness + hardness_ref)
    
    charge_transfer_electrophile = -0.5*(electron_affinity_ref - ionization_potential)*reciprocal_hardness
    charge_transfer_nucleophile  =  0.5*(electron_affinity - ionization_potential_ref)*reciprocal_hardness
    
    return charge_transfer_nucleophile, charge_transfer_electrophile

#############################################################################################################################

 #######  ########  ##     ## 
##     ## ##     ## ###   ### 
       ## ##     ## #### #### 
 #######  ########  ## ### ## 
##        ##        ##     ## 
##        ##        ##     ## 
######### ##        ##     ##

#############################################################################################################################

def two_parabola(cation=None, neutral=None, anion=None, ref_cation=None, ref_neutral=None, ref_anion=None):
    #
    ### Reference species (B)
    #    
    potential_minus_ref, potential_plus_ref, hardness_ref = tpm(cation=ref_cation, neutral=ref_neutral, anion=ref_anion)
    #
    ### Species of interest (A)
    #   
    potential_minus, potential_plus, hardness = tpm(cation=cation, neutral=neutral, anion=anion)
    #
    ### Charge transfer
    #
    reciprocal_hardness          = 1.0/(hardness + hardness_ref)
    
    charge_transfer_electrophile = (potential_minus_ref - potential_plus)*reciprocal_hardness
    charge_transfer_nucleophile  = (potential_plus_ref - potential_minus)*reciprocal_hardness

    return charge_transfer_nucleophile, charge_transfer_electrophile

#############################################################################################################################

 ######    #######  ##     ## 
##    ##  ##     ## ###   ### 
##        ##     ## #### #### 
##   #### ##     ## ## ### ## 
##    ##  ##  ## ## ##     ## 
##    ##  ##    ##  ##     ## 
 ######    ##### ## ##     ##

#############################################################################################################################

def generalized(cation=None, neutral=None, anion=None, ref_cation=None, ref_neutral=None, ref_anion=None):
    #
    ### Reference species (B)
    #    
    potential_minus_ref, potential_plus_ref, hardness_minus_ref, hardness_plus_ref = gqm(cation=ref_cation, neutral=ref_neutral, anion=ref_anion)
    #
    ### Species of interest (A)
    #    
    potential_minus, potential_plus, hardness_minus, hardness_plus = gqm(cation=cation, neutral=neutral, anion=anion)
    #
    ### Charge transfer
    #
    reciprocal_hardness_plus     = 1.0/(hardness_plus + hardness_minus_ref)
    reciprocal_hardness_minus    = 1.0/(hardness_minus + hardness_plus_ref)
    
    charge_transfer_electrophile = (potential_minus_ref - potential_plus)*reciprocal_hardness_plus
    charge_transfer_nucleophile  = (potential_plus_ref - potential_minus)*reciprocal_hardness_minus

    return charge_transfer_nucleophile, charge_transfer_electrophile

#############################################################################################################################
