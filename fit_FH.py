import os
import yaml
import numpy as np
from lmfit import Model
from lmfit.models import LinearModel
import matplotlib.pyplot as plt

from flory_huggins import FloryHuggins, generate_coexistence_curve, calcT, calcChi

configs = yaml.load(open("fit_FH.yml",'r'))
for i in configs.keys():
    FHmodel = configs[i]["model"]
    p0 = configs[i]["p0"]
    lowerT, upperT = configs[i]["temp_range"] 
    temps = np.arange(lowerT,upperT,0.2)+273.15 # temperature range to calulate phase diagram
    outpath = configs[i]["outpath"]
    if os.path.exists(outpath):
        pass
    else:
        os.system("mkdir %s"%outpath)

    print("Fitting...\n %s for "%FHmodel)
    for fname in configs[i]["data"]:
        print(fname)
        # Reading experimental data
        fname = fname
        fname_ext = os.path.splitext(fname)[0]

        data = np.genfromtxt(fname,comments="#")
        T = data[:,0]+273.15
        phi1 = FHmodel.phi_protein(data[:,1])
        phi1_err = FHmodel.phi_protein(data[:,2])
        phi2 = FHmodel.phi_protein(data[:,3])
        phi2_err = FHmodel.phi_protein(data[:,4])

        # Fit chis for dH and dS from phi1 and phi2 values
        x0 = [-1.,100.]
        #mod = Model(calcT)
        mod = Model(calcChi)
        chis = FHmodel.chi(phi1,phi2)
        chi_errs = FHmodel.chiError(phi1,phi1_err,phi2,phi2_err)
        #out = mod.fit(T,x=chis,A=-1,B=10,kws=(FHmodel.chi))
        out = mod.fit(chis,x=T,dS=-1.,dH=10,kws=(FHmodel.chi))
        print(out.fit_report())

        # extract fit params
        dH = out.params["dH"].value
        std_dH = out.params["dH"].stderr
        dS = out.params["dS"].value
        std_dS = out.params["dS"].stderr

        # plot dH and dS fits
        plt.plot(1./T,mod.func(T,dS=dS,dH=dH),"k--",
                label=r"$\Delta H = %.3f \pm %.3f, \Delta S = %.3f \pm %.3f$"%(dH,std_dH,dS,std_dS))
        #plt.errorbar(1./T,chis,yerr=chi_errs,fmt="ro")
        plt.errorbar(1./T,chis,fmt="ro")
        plt.ylabel(r"$\chi$")
        plt.xlabel(r"1/T")
        plt.legend(loc=0)
        plt.savefig(os.path.join(outpath,"Chi_vs_T_%s.pdf"%fname_ext))
        plt.show()

        # plot expt. phi values for coexistence curve 
        plt.errorbar(phi1,T,fmt="ro",xerr=phi1_err,label=fname)
        plt.errorbar(phi2,T,fmt="ro",xerr=phi2_err) 
        plt.ylim(273.15,373.15)
        plt.ylabel(r"T (K)")
        plt.xlabel(r"$\phi$")
    #    plt.show()
        
        # calc chis for range of temps for construction of coexistence curve
        chis = dH/temps+dS
        
        phi1s, phi2s, temperatures, chi_values = generate_coexistence_curve(FHmodel,temps,chis,p0,threshold=1e-5,plot=False)
        plt.plot(phi1s,temperatures,"k-",phi2s,temperatures,"k-")

        # save data
        np.savetxt(os.path.join(outpath,"coexistence_curves_%s.txt"%fname_ext),
                np.column_stack([temperatures,
                              chi_values,
                              phi1s,
                              phi2s,
                              FHmodel.phi_to_conc(phi1s),
                              FHmodel.phi_to_conc(phi2s)]),
                fmt="%.1f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
                header="T\tChi\tphi1\tphi2\tphi1_mgml\tphi2_mgml\n \
                        dH = %.3f +/- %.3f, dS = %.3f +/- %.3f\n \
                        fitted from %s"%(dH,std_dH,dS,std_dS,fname))
        #plt.legend()
        plt.savefig(os.path.join(outpath,"phase_diagram_%s.pdf"%(fname_ext)))
        plt.show()
