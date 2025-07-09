#This script will make the data files required for the ternary potential

import numpy as np
import sys,os
import subprocess
from pylab import *
from scipy.optimize import minimize
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
import shutil
import configparser

#----------------------------------------------------------------------------------------------------------------------------------
# User inputs
#----------------------------------------------------------------------------------------------------------------------------------

run_label='test_1' # A label to distinguish these runs in data_log and the saved lowest files

# System properties
NPOSTX  = 1  #Number of patches in x direction
NPOSTY  = 1  #Number of patches in y direction    ##############################################
WIDTHX  = 1000  #!X-width of patch
WIDTHY  = 10  #!y-width of patch                  ##############################################
HEIGHT1 = 10 #!Porous suface depth
GRIDX=200 #!System size in x-direction [L.U.] Set GRIDX=1 for a 2D simulation
GRIDY=20*NPOSTY #!System size in y-direction [L.U.]
GRIDZ=100   #!System size in z-direction [L.U.]

#Fluid properties
N_PHASE=4 #!Number of phases to simulate
alpha=1.0 #Interface width [L.U.]

#Interfacial tensions between the fluids
g12=1.33
g13=1.9
g14=1.0
g23=0.5
g24=0.8
g34=0.4 

LAMBDA=5.0E1 #Boyer's stabilisation parameter (50 recommended, OR 0 if only partial wetting occurs)

#Fluid constraints
const_label_1=1 #Constrain of phase 1: 0=no constraint, 1=volume constraint, 2=pressure constraint
const_value_1=200000 #Value we wish to constrain phase 1 to 
const_strength_1=0.001 #Magnitude of contraint strength (only used for volume constraints)

const_label_2=2 #Constrain of phase 2: 0=no constraint, 1=volume constraint, 2=pressure constraint
const_value_2=-0.03 #Value we wish to constrain phase 2 to
const_strength_2=0.00 #Magnitude of contraint strength (only used for volume constraints)

const_label_3=2 #Constrain of phase 3: 0=no constraint, 1=volume constraint, 2=pressure constraint
const_value_3=-0.03 #Value we wish to constrain phase 3 to
const_strength_3=0.0 #Magnitude of contraint strength (only used for volume constraints)

const_label_4=0 #Constraint of phase 4: 0=no constraint, 1=volume constraint, 2=pressure constraint
const_value_4=0 #Value we wish to constrain phase 3 to
const_strength_4=0.00 #Magnitude of contraint strength (only used for volume constraints)

#----------------------------------------------------------------------------------------------------------------------------------
# Additional inputs
#----------------------------------------------------------------------------------------------------------------------------------
startstat=True # startstat=True if we have not previously made the LIS

PS=[1,1,1,1] #0=fix this phase so it is not changed during minimisation, 1=vary this phase

dirname='.' #Name of directory to make the data folder in
#----------------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------------------------
# Parameter conversions
#----------------------------------------------------------------------------------------------------------------------------------

# A and B are parameters in Panter's complete wetting model. Currently, Boyer's stabilisation parameter is prefered, so A and B are
# set to zero. 
A=0.0
B=0.0

g13=g13+A
g14=g14+A/2.0
g24=g24+A/2.0
g34=g34+A/2.0

g34=g34+B
g12=g12+B/2.0
g13=g13+B/2.0
g14=g14+B/2.0

#----------------------------------------------------------------------------------------------------------------------------------
# Update the input file name
input_file = f"input.txt"

# Optional: Print or use the input file as needed
print(f"Input file: {input_file}")

config = configparser.ConfigParser()
config.read(input_file)

# Read system properties
NPOSTX = int(config['System properties']['NPOSTX'])
NPOSTY = int(config['System properties']['NPOSTY'])
WIDTHX = int(config['System properties']['WIDTHX'])
WIDTHY = int(config['System properties']['WIDTHY'])
HEIGHT1 = int(config['System properties']['HEIGHT1'])
GRIDX = int(config['System properties']['GRIDX'])
GRIDY = int(config['System properties']['GRIDY'])
GRIDZ = int(config['System properties']['GRIDZ'])

# Read fluid properties
N_PHASE = int(config['Fluid properties']['N_PHASE'])
alpha = float(config['Fluid properties']['alpha'])

# Read interfacial tensions
g12 = float(config['Interfacial tensions']['g12'])
g13 = float(config['Interfacial tensions']['g13'])
g14 = float(config['Interfacial tensions']['g14'])
g23 = float(config['Interfacial tensions']['g23'])
g24 = float(config['Interfacial tensions']['g24'])
g34 = float(config['Interfacial tensions']['g34'])

# Read constraints
const_label_1 = int(config['Constraints']['const_label_1'])
const_value_1 = float(config['Constraints']['const_value_1'])
const_strength_1 = float(config['Constraints']['const_strength_1'])

const_label_2 = int(config['Constraints']['const_label_2'])
const_value_2 = float(config['Constraints']['const_value_2'])
const_strength_2 = float(config['Constraints']['const_strength_2'])

const_label_3 = int(config['Constraints']['const_label_3'])
const_value_3 = float(config['Constraints']['const_value_3'])
const_strength_3 = float(config['Constraints']['const_strength_3'])

const_label_4 = int(config['Constraints']['const_label_4'])
const_value_4 = float(config['Constraints']['const_value_4'])
const_strength_4 = float(config['Constraints']['const_strength_4'])

gamma=np.array([[0.0, g12, g13, g14],[g12, 0.0, g23, g24],[g13, g23, 0.0, g34],[g14, g24, g34, 0.0]])

#----------------------------------------------------------------------------------------------------------------------------------
# Functions
#----------------------------------------------------------------------------------------------------------------------------------

# Make the input parameters
def make_data():
	
	#Write the data files

	#Check if 'data' folder exists
	if not os.path.exists(dirname+'/data'):
		os.mkdir(dirname+'/data')
		
	#Write kappa
	with open(dirname+'/data'+'/gamma.in','w') as f:
		for i in range(0,N_PHASE):
			for j in range(0,N_PHASE):
				f.write("%s\n"%gamma[i,j])
		
	#Write A	
	with open(dirname+'/data'+'/A.in','w') as f:
		f.write("%s\n"%A)
		f.write("%s\n"%B)
		
	#Write wetting parameters
	g_1=0
	g_2=0
	g_3=0
	g_4=0
	with open(dirname+'/data'+'/wetenergy.in','w') as f:
		f.write("%s\n"%g_1)
		f.write("%s\n"%g_2)
		f.write("%s\n"%g_3)
		f.write("%s\n"%g_4)
		
	#Write constraints
	with open(dirname+'/data'+'/constraints.in','w') as f:
		f.write("%s %s %s\n"%(const_label_1,const_value_1,const_strength_1))
		f.write("%s %s %s\n"%(const_label_2,const_value_2,const_strength_2))
		f.write("%s %s %s\n"%(const_label_3,const_value_3,const_strength_3))
		f.write("%s %s %s\n"%(const_label_4,const_value_4,const_strength_4))
		
	#Write phase switches
	with open(dirname+'/data'+'/phaseswitch.in','w') as f:
		for i in range(0,N_PHASE):
			f.write("%s\n"%PS[i])

	return

# Write the data.in file
def write_data():

	#Some variables that are not used currently, but still need to be specified
	TOPEXX  = 0  #!x-extension of reentrant cap
	TOPEXY  = 0  #!y-extension of reentrant cap
	HEIGHT2 = 0  #!Cap thickness
	LIPX    = 0  #!Doubly reentrant lip x-thickness
	LIPY    = 0  #!Doubly reentrant lip y-thickness
	LIPZ    = 0  #!Doubly reentrant lip depth

	with open("data.in",'w') as f:
		f.write("%s\n"%NPOSTX)
		f.write("%s\n"%NPOSTY)
		f.write("%s\n"%WIDTHX)
		f.write("%s\n"%WIDTHY)
		f.write("%s\n"%TOPEXX)
		f.write("%s\n"%TOPEXY)
		f.write("%s\n"%HEIGHT1)
		f.write("%s\n"%HEIGHT2)
		f.write("%s\n"%LIPX)
		f.write("%s\n"%LIPY)
		f.write("%s\n"%LIPZ)
		f.write("%s\n"%GRIDX)
		f.write("%s\n"%GRIDY)
		f.write("%s\n"%GRIDZ)

		f.write("%s\n"%N_PHASE)
		f.write("%s\n"%LAMBDA)
	return
	

		
# Make the inital coords
def make_coords(startstat):
	coords=np.zeros(GRIDX*GRIDY*GRIDZ*(N_PHASE-1),dtype=float)
		
	if startstat==True:
		
		#Make the 2-liquid infused surface
		PERIOD_X=int(GRIDX/(NPOSTX))
		PERIOD_Y=int(GRIDY/(NPOSTY))	
	
		for p1 in range(0,NPOSTX):
			for p2 in range(0,NPOSTY):
				XSTART = 0
				XSTOP = int(XSTART+WIDTHX/2)
				YSTART = int(p2*PERIOD_Y + (PERIOD_Y-WIDTHY)/2)
				YSTOP=YSTART+WIDTHY
				
				for j1 in range(0,PERIOD_X):
					for j2 in range(p2*PERIOD_Y,(p2+1)*PERIOD_Y):				
						for j3 in range(0,HEIGHT1):
							cur=j1*GRIDY*GRIDZ+j2*GRIDZ+j3       
							if ( (j1<XSTOP) and (j2>=YSTART) and (j2<YSTOP)):
								coords[3*cur:3*cur+3]=np.array([0.0,1.0,0.0])
							else:
								coords[3*cur:3*cur+3]=np.array([0.0,0.0,1.0])											

	else:
		
		#Import the minimised coordinates of just the surface
		fid=open("lowests_grid",'r')
		coords=np.loadtxt(fid,skiprows=0)
		fid.close

		#Initialize the droplet as a half-cylinder based on the constant volume.
		angle = pi/2
		V0=const_value_1/GRIDY
		Radius=np.sqrt(V0/(angle-sin(angle)*cos(angle)))
		xc=GRIDX/2
		zc=HEIGHT1

		#Save the phase distribution before initializing the droplet.
		with open("coords_before",'w') as f:
			s1=0
			for j1 in range(0,GRIDX):
				for j2 in range(0,GRIDY):
					for j3 in range(0,GRIDZ):
						for j4 in range(0,N_PHASE-1):				
							f.write("%s\n"%coords[s1])
							s1=s1+1

		#Initialize the droplet
		for j1 in range(0,GRIDX):
			for j2 in range(0,GRIDY):
				for j3 in range(0,GRIDZ):
					cur=j1*GRIDY*GRIDZ+j2*GRIDZ+j3
					dis=(j1-xc)*(j1-xc)+(j3-zc)*(j3-zc)
					dis=np.sqrt(dis)
					if(dis<=Radius) and (j3>HEIGHT1): 
						coords[3*cur]=1.0
						coords[3*cur+1]=0.0
						coords[3*cur+2]=0.0

		#Save the phase distribution after initializing the droplet
		with open("coords_after",'w') as f:
			s1=0
			for j1 in range(0,GRIDX):
				for j2 in range(0,GRIDY):
					for j3 in range(0,GRIDZ):
						for j4 in range(0,N_PHASE-1):				
							f.write("%s\n"%coords[s1])
							s1=s1+1	

	#Write coords to file					
	with open("coords",'w') as f:
		s1=0
		for j1 in range(0,GRIDX):
			for j2 in range(0,GRIDY):
				for j3 in range(0,GRIDZ):
					for j4 in range(0,N_PHASE-1):				
						f.write("%s\n"%coords[s1])
						s1=s1+1
	return 

#Extract the energy of the current system
def get_energy():
	with open('lowest','r') as f:
		lines=f.readlines()
		line=lines[2]
		words=line.split()
		energy=words[2]
		energy = float(energy)
	return energy

# Strip the header lines from the output file 'lowest'
def process_lowest():
	fid=open("lowest",'r')
	lowests=np.loadtxt(fid,skiprows=4)
	fid.close
	np.savetxt('lowests',lowests)	
	return	
#----------------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------------------------
# Compile minimisation code
#----------------------------------------------------------------------------------------------------------------------------------
subprocess.run(["../../GMIN","-n"])
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
# Main
#----------------------------------------------------------------------------------------------------------------------------------

#Write header line of data file
f_data=open('data_log','a')
f_data.write(("%s "*15 + "\n") %("const_value_1", "energy", "g12", "g13", "g14", "g23", "g24", "g34", "GRIDX", "GRIDY", "GRIDZ", "WIDTHX", "WIDTHY", "const_value_2", "const_value_3"))
f_data.close()

# Make the surface if it has not been made in a previous run
if startstat==True:
	#Save original const labels
	cl1=const_label_1
	cl2=const_label_2
	cl3=const_label_3
	
	const_label_1=0
	const_label_2=0
	const_label_3=0
	
	#Make the LIS surface
	make_data()
	write_data()
	print('Making initial LIS')
	make_coords(startstat)
	subprocess.run(["./gmin"])
	process_lowest()
	subprocess.run(["cp","lowests","lowests_grid"])

	#Restore values
	const_label_1=cl1
	const_label_2=cl2
	const_label_3=cl3
	

# Now, make add the drop and minimise

#Update and write the data files
make_data()
write_data()	

#Initialise the system
make_coords(False)

#Run the minimisation
subprocess.run(["./gmin"])
energy=get_energy()

#Save the minimised state
lowestname='lowest_'+run_label
subprocess.run(["cp","lowest",str(lowestname)])

#Save the key simulation parameters
f_data = open('data_log', 'a')
f_data.write(("%s " * 15 + "\n") % (const_value_1, energy, g12, g13, g14, g23, g24, g34, GRIDX, GRIDY, GRIDZ, WIDTHX, WIDTHY, const_value_2, const_value_3))
f_data.close()







