import numpy as np
from numpy import log as ln
import matplotlib.pyplot as plt
from scipy.optimize import root
from lmfit import Model
from lmfit.models import LinearModel

class FloryHuggins():

    def __init__(self,N1=1,N2=1,rho=1400):
        """Simple Flory Huggins Model

        Keyword arguments:
        N1 -- number of lattice sites occupied by one polymer i.e number of amino acids (default 1)
        N2 -- number of lattice sites occupied by solvent (default 1)
        rho -- density of protein in mg/mL (default 1400)

        """
        self.N1 = N1
        self.N2 = N2
        #self.mw_aa = mw_aa # average Mw of amino acid (110 g/mol)
        self.rho = rho # density of protein (g/cm^3)
        # calculate protein molecular weight
        #if mw_protein is None:
        #    self.mw_protein = self.N1 * self.mw_aa
        #else:
        #    self.mw_protein = mw_protein

    def fmix(self,phi1,chi):
        """Flory Huggins free energy of mixing"""
        phi2 = 1.-phi1
        return phi1/self.N1*ln(phi1) + phi2/self.N2*ln(phi2) + chi*phi1*phi2

    def dfmix_by_dphi(self,phi1,chi):
        """Derivative with respect to phi of Flory Huggins free energy of mixing"""
        phi2 = 1.-phi1
        return (1.+ln(phi1))/self.N1 - (1.+ln(phi2))/self.N2 + chi*(1.-2.*phi1)

    def _find_common_tangent(self,p,chi):
        """functions for finding common tangent of free energy curves"""
        phi1,phi2 = p
        equal_derivs = self.dfmix_by_dphi(phi1,chi) - self.dfmix_by_dphi(phi2,chi) # = 0
        common_slope = (self.fmix(phi1,chi) - self.fmix(phi2,chi)) / (phi1 - phi2) - self.dfmix_by_dphi(phi1,chi)# = 0 
        return [equal_derivs,common_slope]

    def find_common_tangent(self,p,chi):
        """Find the common tangent line between two free energy minima
        
        Keyword Arguments:
        p -- start_params as list, [phi1,phi2], these have to be carefully selected otherwise algorithm can fail
        chi -- chi parameters
        """
        phi1,phi2 = p
        return root(self._find_common_tangent,[phi1,phi2],args=(chi),method="lm")

    def phi_protein(self,c):
        """convert protein concentration in mg/mL into volume fraction
        1 mgmL is eq to 1 kg/m^3
        density is in SI units of kg/m^3  
        """
        return c / self.rho

    def phi_to_conc(self,phi):
        """convert phi value to protein concentration in mg/mL."""
        return phi * self.rho 

    def chi(self,phi1,phi2):
        """calculate chi from ph1 and phi2
        this is a rearrangement of the equation dfmix_by_dphi(phi1,chi)=dfmix_by_dphi(phi2,chi)
        for chi where phi1 and phi2 are the volume fractions of protein in the dilute and condensed
        phases, respectively.
        """
        return (1./self.N1*ln(phi2/phi1) + 1./self.N2*ln((1-phi1)/(1-phi2)))/(2.*phi2-2.*phi1)

    def chiError(self,phi1,dphi1,phi2,dphi2):
        # propagate errors in calculation of chi parameter
        # needs checking
        expr1 = ((-1./(phi1*self.N1)-1./((1-phi1)*self.N2))*(2*phi2-2*phi1)+2*((1./self.N1)*ln(phi2/phi1)+(1./self.N2)*ln((1-phi1)/(1-phi2))))/(2*phi2-2*phi1)**2.
        expr2 = ((1./(phi2*self.N1)+1./((1-phi2)*self.N2))*(2*phi2-2*phi1)+2*((1./self.N1)*ln(phi2/phi1)+(1./self.N2)*ln((1-phi1)/(1-phi2))))/(2*phi2-2*phi1)**2.
        err = np.sqrt(expr1**2.*dphi1**2. + expr2**2.*dphi2**2.)
        return err
        
    def __str__(self):
        s = (self.N1, self.N2, self.rho)
        return ("Flory Huggins model\n N1=%d\n N2=%d\n rho=%.1f\n"%s)

#############################################
# Functions for fitting chi from phi values #
#############################################
def f_to_temp(phi12,params,chi):
    # converts phi1 and phi2 values to temperature
    # chi is an instance of flory-huggins model.chi
    phi1, phi2 = phi12
    A,B = params
    return B/chi(phi1,phi2)-A

def residual(p,phi12,chi):
    # calculate difference between temperatures and model
    return T - f_to_temp(phi12,p,chi)

def calcT(x,dS,dH):
    # rearrangement of chi to give temperature
    return dH/(x-dS)

def calcChi(x,dS,dH):
    # returns Chi value given x=T and dS and dH
    return 1./x * dH + dS


##########################################
# Function to generate coexistence curve #
##########################################

def generate_coexistence_curve(FHmodel,temps,chis,p0,threshold=1e-6,plot=True):
    """ Generate coexistence curve (phase diagram)

    Keyword arguments:
    FHmodel -- initialised FloryHuggins class object
    p0 -- starting phi1,phi2 values e.g. [0.001,0.7]
    threshold -- if difference between phi1 and phi2 is less than this values are discarded (default 1e-6)
    plot -- make matplotlib plot of phase diagram (default True)
    
    Returns:
    four numpy arrays containing phi1s,phi2s,temperature and chi_values
    for which the difference between phi1 and phi2 was more than the threshold
    """
    phi1s = []
    phi2s = []
    temperatures = []
    chi_values = []

    for t,chi in zip(temps,chis):
        # iterate through temp,and min_chi,chi,max_chi to find phi values
        # finds common tangent using p0 and chi
        solve = FHmodel.find_common_tangent(p0,chi)
        print(solve.x)

        if abs(solve.x[0]-solve.x[1]) < threshold:
            # checks whether you are near the critical point
            pass
        else:
            phi1s.append(solve.x[0])
            phi2s.append(solve.x[1])
            temperatures.append(t)
            chi_values.append(chi)
    if plot:
        phis = np.concatenate([phi1s,phi2s[::-1]])
        Ts = np.concatenate([temperatures,temperatures[::-1]])
        plt.plot(phis,Ts,"k-")
            
    phi1s = np.array(phi1s)
    phi2s = np.array(phi2s)
    temperatures = np.array(temperatures)
    chi_values = np.array(chi_values)

    return phi1s, phi2s, temperatures, chi_values
