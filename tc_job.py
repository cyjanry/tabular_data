########################################
# Define the limitation for plotting:
Tmin        =   500.0                    #[K] The lower limit of the calculation temperature
Tmax        =   1500.0                   #[K] The upper limit of the calculation temperature
            
pmin        =   5.0e6                    #[kPa] The lower limit of the pressure need to be compared
pmax        =   32.0e6                   #[kPa] The upper limit of the pressure need to be compared


nTh         =   150.                      # How many data points for temperature in high resolution grid
nPh         =   100.                      # How many data points for pressure in high resolution grid


nTcmin      =   10.
nTcmax      =   20.

nPcmin      =   4.                       # How many data points for temperature  in coarse resolution grid
nPcmax      =   10.                      # How many data points for pressure in coarse resolution grid

percent_tol =   0.5                        # [%]
tol         =   1.e-7                    # Usding the tolerance to adjust the convergence. Usually 1.e-7 is fine.
########################################



########################################
# Actaviate the fluid parameter, with REFEPROP:
specie     =   'D'    # which property you want to compare with?
fluidType  =   'CO2'                     # Fluid name
########################################