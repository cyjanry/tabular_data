#!/usr/bin/env python
#######################################################################
#    This code is used to select the right resolution of the tabuler  #
# which is used in look up table real gas properties.                 #
#                                                                     #
#    Author :   Jianhui Qi        j.qi@uq.edu.au                      #
#    Advisor:   Ingo H. J. Jahn   i.jahn@uq.edu.au                    #
#    Date   :   01-06-2017                                            #
#    version:   1.                                                    #
#######################################################################

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt 
from pyRefpropMania import *

########################################
# Define the limitation for plotting:
Tmin        =   600.0                    #[K] The lower limit of the calculation temperature
Tmax        =   1000.0                   #[K] The upper limit of the calculation temperature
            
pmin        =   8.0e6                    #[kPa] The lower limit of the pressure need to be compared
pmax        =   22.0e6                   #[kPa] The upper limit of the pressure need to be compared


nTh         =   100.                      # How many data points for temperature in high resolution grid
nPh         =   100.                      # How many data points for pressure in high resolution grid

nTcmax      =   7.                       # How many data points for temperature  in coarse resolution grid
nPcmax      =   8.                      # How many data points for pressure in coarse resolution grid

tol         =   1.e-7                    # Usding the tolerance to adjust the convergence. Usually 1.e-7 is fine.
########################################

specie     =  'D'    # which property you want to compare with?

########################################
# Actaviate the fluid parameter, with REFEPROP:

fluidType  =   'CO2'                     # Fluid name
myFluid    =   refpropFluid(fluidType)   # Actaviate the REFPROP data base
R          =   8.3144598                 # [J/(K mol)] --> SI unit [kg m2 s-2 K-1 mol-1]
mW         =   44.01e-3                  # [kg/mol] molecure weight 
########################################



# calculation reference high resolution grid
T_high   = np.arange(Tmin, Tmax + (Tmax - Tmin)/(nTh) , (Tmax - Tmin)/(nTh-1.))
p_high   = np.arange(pmin, pmax + (pmax - pmin)/(nPh) , (pmax - pmin)/(nPh-1.))
Ph, Th = np.meshgrid(p_high,T_high)   # T row and p colume



def interpolation_2D(nTc,nPc):

    #print 'starting to interpolation progress'
    T_coarse = np.arange(Tmin, Tmax + (Tmax - Tmin)/(nTc) , (Tmax - Tmin)/(nTc-1.))
    p_coarse = np.arange(pmin, pmax + (pmax - pmin)/(nPc) , (pmax - pmin)/(nPc-1.))
    Pc, Tc = np.meshgrid(p_coarse,T_coarse)
    #print "---------------\n",Pc,"\n------------------------\n",Tc

    Zh           = []
    Zc           = []
    Z            = []
    rhoh,rhoc    = [],[]
    muh,muc      = [],[]
    kappah,kappac= [],[]
    Cph,Cpc      = [],[]
    Cvh,Cvc      = [],[]
    hh,hc        = [],[]
    #cpMcv_list   = []

    # calculation to grep the Z axis data.
    #
    #         j(P)
    #      ___________
    #      | . . . . 
    # i(T) | . . . .
    #      | . . . .
    #      | . . . . 
    #      |
    #
    #
    #
    # for Z list 0
    for i in range(len(Ph)):
        rho_temp = []
        mu_temp  = []
        kappa_temp = []
        Cp_temp  = []
        Cv_temp  = []
        h_temp   = []
        for j in range(len(Ph[i])):
            if specie == 'D':  rho_temp.append(   myFluid.getProps('D','T',Th[i][j]-273.15,'P',Ph[i][j]/1000.)       )
            if specie == 'V':  mu_temp.append(    myFluid.getProps('V','T',Th[i][j]-273.15,'P',Ph[i][j]/1000.)/1.e6  )  #[Pa*s] Or [kg/(m s)]
            if specie == 'K':  kappa_temp.append( myFluid.getProps('K','T',Th[i][j]-273.15,'P',Ph[i][j]/1000.)       )
            if specie == 'C':  Cp_temp.append(    myFluid.getProps('C','T',Th[i][j]-273.15,'P',Ph[i][j]/1000.)*1000. )  #[J/(kg K)]
            if specie == 'O':  Cv_temp.append(    myFluid.getProps('O','T',Th[i][j]-273.15,'P',Ph[i][j]/1000.)*1000. )  #[J/(kg K)]
            if specie == 'H':  h_temp.append(     myFluid.getProps('H','T',Th[i][j]-273.15,'P',Ph[i][j]/1000.)*1000. )  #[J/(kg K)]
        if   specie == 'D':rhoh.append(rho_temp)
        elif specie == 'V':muh.append(mu_temp)
        elif specie == 'K':kappah.append(kappa_temp)
        elif specie == 'C':Cph.append(Cp_temp)
        elif specie == 'O':Cvh.append(Cv_temp)
        elif specie == 'H':hh.append(h_temp)
    if   specie == 'D':     Zh = rhoh
    elif specie == 'V':     Zh = muh
    elif specie == 'K':     Zh = kappah
    elif specie == 'C':     Zh = Cph
    elif specie == 'O':     Zh = Cvh
    elif specie == 'H':     Zh = hh

    for i in range(len(Pc)):
        rho_temp = []
        mu_temp  = []
        kappa_temp = []
        Cp_temp  = []
        Cv_temp  = []
        h_temp   = []
        for j in range(len(Pc[i])):
            if specie == 'D':  rho_temp.append(   myFluid.getProps('D','T',Tc[i][j]-273.15,'P',Pc[i][j]/1000.)       )
            if specie == 'V':  mu_temp.append(    myFluid.getProps('V','T',Tc[i][j]-273.15,'P',Pc[i][j]/1000.)/1.e6  )  #[Pa*s] Or [kg/(m s)]
            if specie == 'K':  kappa_temp.append( myFluid.getProps('K','T',Tc[i][j]-273.15,'P',Pc[i][j]/1000.)       )
            if specie == 'C':  Cp_temp.append(    myFluid.getProps('C','T',Tc[i][j]-273.15,'P',Pc[i][j]/1000.)*1000. )  #[J/(kg K)]
            if specie == 'O':  Cv_temp.append(    myFluid.getProps('O','T',Tc[i][j]-273.15,'P',Pc[i][j]/1000.)*1000. )  #[J/(kg K)]
            if specie == 'H':  h_temp.append(     myFluid.getProps('H','T',Tc[i][j]-273.15,'P',Pc[i][j]/1000.)*1000. )  #[J/(kg K)]
        if   specie == 'D':rhoc.append(rho_temp)
        elif specie == 'V':muc.append(mu_temp)
        elif specie == 'K':kappac.append(kappa_temp)
        elif specie == 'C':Cpc.append(Cp_temp)
        elif specie == 'O':Cvc.append(Cv_temp)
        elif specie == 'H':hc.append(h_temp)
    if   specie == 'D':     Zc = rhoc
    elif specie == 'V':     Zc = muc
    elif specie == 'K':     Zc = kappac
    elif specie == 'C':     Zc = Cpc
    elif specie == 'O':     Zc = Cvc
    elif specie == 'H':     Zc = hc


    # Staring to interpolate:
    Zi = [] # initialize an empty value list 
    #        bj__________
    #      a|  fx |
    #      i|fy   |
    #       |_____|
    #       |     A
    #       |
    #       |
    for a in range(int(nTh)):

        Ztemp_list = []
        for b in range(int(nPh)):
            #print "a=",a,"b=", b
            #print "--------A----------"

            for i in range(int(nTc)):
                #print "i=",i
                if abs(T_high[a] - T_coarse[i])<tol: # check if T_high[a] is equal to the lower limit value of the T_coarse
                    fy = 0.
                    break
                elif abs(T_high[a] - T_coarse[i+1])<tol: # check if T_high[a] is equal to the upper limit value of the T_coarse
                    #i = i+1 # make sure that the T_coarse is picking the upper limit value
                    fy = 1.
                    break
                elif  (T_high[a]-T_coarse[i]) * (T_high[a]-T_coarse[i+1]) < 0.:
                    fy = (T_high[a] - T_coarse[i])/(T_coarse[i+1] - T_coarse[i])
                    break
                else:
                    pass
            #print "T_high[",a,"]=",T_high[a]," and, T_coarese[",i,"]=",T_coarse[i]," and T_coarse[",i+1,"]=",T_coarse[i+1]
            #print "--------B----------"
            for j in range(int(nPc)):
                #print "j=",j

                #print "p_high[",b,"]=",p_high[b]," and, p_coarese[",j,"]=",p_coarse[j], ' and ',p_high[b] - p_coarse[j]#," and p_coarse[",j+1,"]=",p_coarse[j+1]
                if abs(p_high[b] - p_coarse[j])<tol:
                    #print p_high[b],p_coarse[j]

                    fx = 0.
                    break
                elif abs(p_high[b] - p_coarse[j+1])<tol:
                    #print p_high[b]
                    fx = 1.
                    #j = j+1 # make sure that the p_coarse is picking up the upper limit value.
                            # Because p_high[b] is equal to p_coarse[j+1], the j should be add 1 to j+1
                    break
                elif (p_high[b]-p_coarse[j]) * (p_high[b]-p_coarse[j+1]) < 0.:

                    fx = (p_high[b] - p_coarse[j])/(p_coarse[j+1] - p_coarse[j])
                    break
                else:
                    pass
            #print "fx = ",fx,"fy = ", fy, " and i=",i," j=",j , " and T_high[a]=",T_high[a]," and P_high[b]=",p_high[b]
            #print "a=",a," and b=",b         
            A11 = Zc[i][j]    #;  print "A11 = ", A11
            A12 = Zc[i][j+1]  #;  print "A12 = ", A12
            A21 = Zc[i+1][j]  #;  print "A21 = ", A21
            A22 = Zc[i+1][j+1]#;  print "A22 = ", A22
            temp0 = A11 + (A12 - A11)*fx
            temp1 = A21 + (A22 - A21)*fx
            Ztemp = temp0 + (temp1-temp0)*fy

            Ztemp_list.append(Ztemp)
        Zi.append(Ztemp_list)
    #print "End the interpolation progress"
    return Zh,Zc,Zi

	#print Zc
	#print Th, Tc,"\n-------\n",Ph,Pc
	#print Zh,Zc
	#print "-------\n",Zi[0],"\n-----------------------\n",Zh[0]

	######################################################################################
	# checking if the code is working via plot the difference between the first line data:
	# During calculation, you can bulk comment them all.
	'''
	temp = []
	for m in range(len(Zi[0])):
	    temp.append(Zi[0][m] - Zh[0][m])
	print temp

	#And then plot the error:
	fig = plt.figure()
	ax  = fig.add_subplot(111)
	x   = np.arange(0,len(temp),1)
	ax.plot(x, temp)
	plt.show()
	'''
	######################################################################################

	# calculation the percent of error for the interpolation methods


def err_calculation(Zi,Zh):

	err = []
	for m in range(len(Zi)):
	    err_temp = []
	    for n in range(len(Zi[m])):
	        err_temp.append( abs((Zi[m][n] - Zh[m][n])/Zh[m][n])*100. )
	    err.append(err_temp)

	maxmium_error = np.amax(err)
	print maxmium_error
	return maxmium_error


def plot_err_surface(X,Y,err_list):

	fig = plt.figure()
	ax = fig.add_subplot(111, projection = '3d')
	# Plot the surface.
	#surf = ax.plot_surface(Pc,Tc, Zc, #cmap=cm.coolwarm,
	#                       linewidth=0, antialiased=False)
	#ax.scatter(Ph,Th,Zh,c ='b',alpha = 0.5)
	#ax.scatter(Pc,Tc,Zc,c='r',marker = '^',alpha = 1.0)
	ax.plot_surface(X,Y, err_list,linewidth=0, antialiased=False)
	ax.scatter(X,Y,err_list,c='r')
	
	


	ax.set_xlabel('Number of Temperature arguments')
	ax.set_ylabel('Numebr of Pressure arguments')
	ax.set_zlabel('Error/[%]')
	plt.show()



def main():

    Tc_list  = np.arange(2,nTcmax+1,1)
    Pc_list  = np.arange(2,nPcmax+1,1)
    PCL,TCL  = np.meshgrid(Pc_list,Tc_list)
    #print PCL,"\n",TCL
    Zh_list  = []
    Zi_list  = []
    err_list = []
    print "############################################################"
    print "    This code will help you to determin the resolution "
    print "    you need to perform the real gas properties with"
    print "    an maximum error below 1%"
    print "    The reference resolution is:\n"
    print "    Temperature number: ",nTh," and Pressure number:",nPh,"\n"
    print "    Then the figure will show " 
    print "############################################################\n"


    for i in range(len(Tc_list)):
    	err_temp =[]
    	for j in range(len(Pc_list)):
    		print "\n---------------------------\n"
    		print "Stating the interpolation for nTc=",Tc_list[i]," and nPc=",Pc_list[j],"..."
    		Zh,Zc,Zi = interpolation_2D(Tc_list[i],Pc_list[j])

    		max_err = err_calculation(Zi,Zh)
#    		Zh_list.append(Zh)
#    		Zi_list.append(Zi)
    		err_temp.append(max_err)
    #print TCL,PCL,err_list
        err_list.append(err_temp)
    plot_err_surface(TCL,PCL,err_list)
    return 0

###
if __name__ == "__main__":
    #userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    #uoDict = dict(userOptions[0])
    
    #if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        #printUsage()
        #sys.exit(1)
    

    # execute the code
    try:
        main()
        print "\n \n"
        print "SUCESS."
        print "\n \n"
    except:
    	print "This run has gone bad."
    	sys.exit(1)
    	

    #except MyError as e:
    #    print "This run has gone bad."
    #    print e.value
    #    sys.exit(1)
