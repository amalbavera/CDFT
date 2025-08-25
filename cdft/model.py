#############################################################################################################################

##       #### ########  ########     ###    ########  #### ########  ######
##        ##  ##     ## ##     ##   ## ##   ##     ##  ##  ##       ##    ##
##        ##  ##     ## ##     ##  ##   ##  ##     ##  ##  ##       ##
##        ##  ########  ########  ##     ## ########   ##  ######    ######
##        ##  ##     ## ##   ##   ######### ##   ##    ##  ##             ##
##        ##  ##     ## ##    ##  ##     ## ##    ##   ##  ##       ##    ##
######## #### ########  ##     ## ##     ## ##     ## #### ########  ######

#############################################################################################################################

#############################################################################################################################

   ##   ########  ##     ## 
 ####   ##     ## ###   ### 
   ##   ##     ## #### #### 
   ##   ########  ## ### ## 
   ##   ##        ##     ## 
   ##   ##        ##     ## 
 ###### ##        ##     ##

#############################################################################################################################

def one_parabola(cation=None, neutral=None, anion=None):

    ionization_potential = cation.energy - neutral.energy
    electron_affinity    = neutral.energy - anion.energy
    
    chemical_potential   = -0.5*(ionization_potential + electron_affinity)
    chemical_hardness    = ionization_potential - electron_affinity
    
    return chemical_potential, chemical_hardness

#############################################################################################################################

 #######  ########  ##     ## 
##     ## ##     ## ###   ### 
       ## ##     ## #### #### 
 #######  ########  ## ### ## 
##        ##        ##     ## 
##        ##        ##     ## 
######### ##        ##     ##

#############################################################################################################################

def two_parabola(cation=None, neutral=None, anion=None,):

    ionization_potential     = cation.energy - neutral.energy
    electron_affinity        = neutral.energy - anion.energy
    
    chemical_potential_plus  = -0.25*(ionization_potential + 3.0*electron_affinity)
    chemical_potential_minus = -0.25*(3.0*ionization_potential + electron_affinity)
    
    chemical_hardness        = 0.5*(ionization_potential - electron_affinity)

    return chemical_potential_minus, chemical_potential_plus, chemical_hardness

#############################################################################################################################

 ######    #######  ##     ## 
##    ##  ##     ## ###   ### 
##        ##     ## #### #### 
##   #### ##     ## ## ### ## 
##    ##  ##  ## ## ##     ## 
##    ##  ##    ##  ##     ## 
 ######    ##### ## ##     ##

#############################################################################################################################

def generalized(cation=None, neutral=None, anion=None):

    ionization_potential     = cation.energy - neutral.energy
    electron_affinity        = neutral.energy - anion.energy
    
    chemical_potential_plus  = neutral.lumo
    chemical_potential_minus = neutral.homo
    
    chemical_hardness_plus   = -2.0*(electron_affinity    + neutral.lumo)
    chemical_hardness_minus  =  2.0*(ionization_potential + neutral.homo)

    return chemical_potential_minus, chemical_potential_plus, chemical_hardness_minus, chemical_hardness_plus

#############################################################################################################################
