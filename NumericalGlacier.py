import numpy as np
import matplotlib.pyplot as plt
from math import exp # exponential function
matplotlib.use('qt5agg')

# Flow law and physical constants
nf = 3 # Power in Glens flow law
g = 9.82 # Gravity
rho = 918 # Ice density
isotemp = 270  # Isotermal
activation_energy = 60000 # Activation energy for flow law (J/mol)
gas_constant = 8.3143  # Gas constant (J/(mol*K))
# Flow law rate factor
rate_factor = 1.35*10**(-5)*exp(-activation_energy/gas_constant/isotemp)
flow_law_coefficient  = 2*rate_factor*(rho*g)**nf/(nf+2)

rho = []
rho_water = 997 # Density of water
rho_ice = 917 # Density of ice

# Set up model domain
num_gridpoints = 1001 # Number of gridpoints
grid_spacing = 1000 # Grid spacing
# Grid along x-axis
x = (np.arange(0,num_gridpoints) * grid_spacing).reshape(-1,1)
# Staggered grid along x-axis
staggered_x_grid = (np.arange(0,num_gridpoints-1) * grid_spacing + grid_spacing/2).reshape(-1,1)

# Time evolution parameters
dt = 0.01 # Time step
num_timesteps = 50000 # Number of timsteps
time = 0 # Start time
start_timesteps = 0 # Number of timestep at start (for restarting a run)

# Surface mass balance
surface_mass_balance = np.zeros((num_gridpoints,1))
surface_mass_balance[0:]=0.5

# Bottom topography
# Flat or gaussian shaped
bottom_topography  = np.zeros((num_gridpoints,1))
b0 = 300
slope = 0.005
lambd = 300
x_shift = 130000
sigma = 15000
for i in range(0, num_gridpoints):
    bottom_topography [i] = b0-slope*x[i]+lambd*np.exp(-((x[i]-x_shift)/(sigma))**2)

# Initial profile
# No ice sheet:
ice_thickness = np.zeros((num_gridpoints,num_timesteps+1))
# Vialov-profilen:
# hdiv=((0.3/(2*rate_factor))^(1/nf)*(nf+2)^(1/nf)*2/(rho*g))^(1/(2+2/nf))*500000^0.5;
# ice_thickness(1:501,1)=(1-(x(1:501)/500000).^(4/3)).^(3/8)*hdiv;
# hdiv=((0.3/(2*rate_factor))^(1/nf)*(nf+2)^(1/nf)*2/(rho*g))^(1/(2+2/nf))*200000^0.5;
# ice_thickness(501:701,1)=(1-((x(501:701)-x(501))/200000).^(4/3)).^(3/8)*hdiv;
# ice_thickness(301:501,1)=ice_thickness(701:-1:501);

surface_elevation = ice_thickness + bottom_topography     # Surface elevation

# Calculate the total volume
total_volume = 0
total_mass_balance = 0
for i in range(0, num_gridpoints-1):
    total_mass_balance = total_mass_balance + 0.5*(surface_mass_balance[i]+surface_mass_balance[i+1])*grid_spacing
    total_volume = total_volume + 0.5*(ice_thickness[i]+ice_thickness[i+1])*grid_spacing

# Visualization setup
plt.ion()
fig = plt.figure(1)
ax = fig.gca()
ax.plot(x,surface_elevation[:,0])
ax.plot(x,bottom_topography )
plt.axis([0, 3000, 0, 400000])

# Arrays for calculations
hsg = np.zeros(num_gridpoints-1)
dsdx = np.zeros(num_gridpoints-1)
diff = np.zeros(num_gridpoints-1)
flux = np.zeros(num_gridpoints-1)
timestep = np.zeros((num_timesteps,1))
difftest = []

# Time iteration
for j in range(start_timesteps, num_timesteps+start_timesteps):
    if (j % 500) == 0: 
        plt.plot(x[:,0], surface_elevation[:,j], color="blue")
        plt.plot(x[:,0], x[:,0]*0, color="black")
        plt.axis([0, 500000, -4750, 3300])
        plt.draw()
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.gcf().show()
        
    # Explicit iteration scheme
    # thickness, slopes og diffusion (staggered grid)
    for i in range(num_gridpoints-1):
        hsg[i] = (ice_thickness[i+1,j] + ice_thickness[i,j]) / 2
        dsdx[i] = (surface_elevation[i+1,j] - surface_elevation[i,j]) / grid_spacing
        diff[i] = flow_law_coefficient * hsg[i]**5 * dsdx[i]**2
        flux[i] = -diff[i]*dsdx[i]

    # Calculate maximum stable time step
    difftest = np.max(diff)
    if difftest > 0:
        dt = grid_spacing*grid_spacing/difftest/3
        dt = np.min([dt, 1.0])
        dt = np.max([dt, 0.0001])
    else:
        dt = 0.0001
        
    # Update time
    time = time + dt
    timestep[j] = dt
    # Matrix coefficients for explicit time iteration
    c = np.zeros((num_gridpoints,num_gridpoints))
    # Boundary conditions
    diffid = flow_law_coefficient /2*ice_thickness[0,j] * ice_thickness[1,j]**4*((surface_elevation[2,j]-surface_elevation[0,j])/(2*grid_spacing))**2
    c[0,0] = -diffid
    c[0,2] = diffid
  
    for i in range(1, num_gridpoints-1):
        c[i,i-1] = diff[i-1]
        c[i,i] = -diff[i-1]-diff[i]
        c[i,i+1] = diff[i]
    c[num_gridpoints-1, num_gridpoints-1] = 1

    # Calculate new surface
    ice_thickness[:,j+1] = (ice_thickness[:,j].reshape(-1,1) + dt*(surface_mass_balance + 1/(grid_spacing**2)* np.matmul(c, ice_thickness[:,j]).reshape(-1,1) )).reshape(1,-1) 

    # Ensure ice thickness doesn't become negative
    for i in range(num_gridpoints):
        if ice_thickness[i,j+1] < 0:
            ice_thickness[i,j+1]=0
            
    # Check for water-ice interface
    for i in range(num_gridpoints-1):
        if -rho_water * bottom_topography [i,0] >= rho_ice * ice_thickness[i,j]:
            ice_thickness[i+1,j+1]=0
               
    # Update surface
    surface_elevation[:,j+1] = bottom_topography [:,0] + ice_thickness[:,j+1]
    
    
# End of time iteration
# Final plot
plt.figure(2)
plt.plot(x, surface_mass_balance)
plt.show()
