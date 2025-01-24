#This script is designed such that the magnetofriction can ITERATIVELY pick the best omega to match the lower boundary nicely at all times.

#This is going to be bloody difficult...

#PLAN:

#Set up MF parameters as normal in here. Use a best guess for the ideal vaule of omega
#Initially that will be a bit stilly

#Establish that reading in existing steps is fine -- snapshots will need to align with mags, unfortunately.

#At every mag step, run with an omega for ONE step. Make logs of what happens.
#Determine ideal Omega based on this, and stop when necessary.

#EXTRAS relative to the original: mag_start, mag_end? As in, code runs from the start to the end magnetogram.

#Can remove snapshots after the fact if there are too many.



import os
import shutil
import numpy as np
import sys
from numpy import random
import time
import matplotlib.pyplot as plt

from init import compute_initial_condition
from write_electric import compute_electrics
from scipy.io import netcdf_file

from helicity_tools import compute_helicity, compare_fields

if len(sys.argv) > 1:
    run = int(sys.argv[1])
else:
    run = 0

if len(sys.argv) > 2:
    nprocs = int(sys.argv[2])
else:
    nprocs = 1

try:
    if os.uname()[1] == 'brillouin.dur.ac.uk':
        hflag = 0
        winflag = 0
    else:
        hflag = 1
        winflag = 0
except:
    hflag = 0
    winflag = 1

start = 0
#DYNAMIC SYSTEM PARAMETERS
#-------------------------------------
voutfact = 5.0
shearfact = 1.0#3.7e-5   #factor by which to change the imported 'speed'
eta0 = 0.0

tmax = 750.0
tstart = 0.0

nx = 64
ny = 64
nz = 64

ndiags = 750
nmags = max(500, 500*tmax/250.0) #Number of magnetograms used.
nplots = nmags


nu0 = 10.0#np.geomspace(1.0,50.0,10)[run]
eta = 5e-4*nu0

x0 = -130.0; x1 = 130.0
y0 = -130.0; y1 = 130.0
z0 = 0.0; z1 = 130.0

init_number = run
#omega = np.geomspace(1e-4,1e-1,10)[run]
#omega = 1e-2[1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,1e-1][run]
omega = 1e-2

backfield_angle = 0.1#Angle of background field in degrees.
#Variables for the pressure term
decay_type = 0  #Decay types -- 0 for none, 1 for exponential, 2/3 for tanh. Same as the 2D cases.

if decay_type == 0: #No pressure
    zstar = 0.0; a = 0.0; b = 0.0; deltaz = 0.0

if decay_type == 1: #exponential decay
    zstar = run#np.linspace(0.0,0.3,10)[run//50]*z1
    if zstar > 0:
        b = zstar/np.log(2)
        a = 0.5
        deltaz = 0.0*z1
    else:
        decay_type = 0
        zstar = 0.0; a = 0.0; b = 0.0; deltaz = 0.0

if decay_type == 2: #smooth tanh
    a = 0.25; b = 1.0
    zstar = np.linspace(0.0,0.3,10)[run]*z1
    deltaz = 0.1*z1

if decay_type == 3: #sharp tanh
    a = 0.25; b = 1.0
    zstar = np.linspace(0.0,0.3,7)[run]*z1
    deltaz = 0.02*z1

#SAVE VARIABLES TO BE READ BY FORTRAN
#-------------------------------------

variables = np.zeros((40))
#Basics
variables[0] = run
variables[1] = nx
variables[2] = ny
variables[3] = nz
variables[4] = tmax
#Outputs
variables[5] = nplots
variables[6] = ndiags
#Parameters
variables[7] = voutfact
variables[8] = shearfact
variables[9] = eta
variables[10] = nu0
variables[11] = eta0
#Grid things
variables[12] = x0
variables[13] = x1
variables[14] = y0
variables[15] = y1
variables[16] = z0
variables[17] = z1

#Pressure things
variables[18] = a
variables[19] = b
variables[20] = deltaz
variables[21] = zstar

variables[22] = hflag

variables[23] = decay_type

#Background field angle (degrees from vertical)
variables[24] = backfield_angle

#Number of imported magnetograms
variables[25] = nmags
variables[26] = tstart

variables[27] = init_number #Code to give the magnetograms and electric fields so they don't need to be done every time'

variables[28] = omega   #Just for plotting and things. Not used in the Fortran

#SOME FOLDER ADMIN
#-------------------------------------

if not hflag:
    data_directory = '/extra/tmp/trcn27/mf3d/%03d/' % run
else:
    data_directory = '/nobackup/trcn27/mf3d0/%03d/' % run

if winflag:
    data_directory = '/Data/'
    
if not os.path.isdir('./inits'):
    os.mkdir('./inits')

if not os.path.isdir('./parameters'):
    os.mkdir('./parameters')

if not os.path.isdir('./diagnostics'):
    os.mkdir('./diagnostics')

if not os.path.isdir(data_directory[:-4]):
    os.mkdir(data_directory[:-4])

if os.path.isdir(data_directory) and start == 0:
    for i in range(1000):
        if os.path.isfile('%s%04d.nc' % (data_directory, i)):
            os.remove('%s%04d.nc' % (data_directory, i))
else:
    os.mkdir(data_directory)

if not os.path.isdir('./mf_mags/'):
    os.mkdir('./mf_mags/')

#Make magnetogram filenames if necessary
if os.path.isdir('./mf_mags/%03d/' % run) and start == 0:
    for i in range(1000):
        if os.path.isfile('./mf_mags/%03d/%04d.nc' % (run, i)):
            os.remove('./mf_mags/%03d/%04d.nc' % (run, i))
else:
    os.mkdir('./mf_mags/%03d/' % run)


#Create initial condition using new init.py (potential field with arbitrary lower boundary and domain dimensions)
#-------------------------------------

class Grid():
    def __init__(self):
        self.x0 = x0; self.x1 = x1
        self.y0 = y0; self.y1 = y1
        self.z0 = z0; self.z1 = z1
        self.nx = nx ; self.ny = ny; self.nz = nz
        self.dx = (self.x1 - self.x0)/nx
        self.dy = (self.y1 - self.y0)/ny
        self.dz = (self.z1 - self.z0)/nz
        self.xs = np.linspace(self.x0,self.x1,self.nx+1)
        self.ys = np.linspace(self.y0,self.y1,self.ny+1)
        self.zs = np.linspace(self.z0,self.z1,self.nz+1)

mag_max = 500
#All the above needs to be set up once -- the following needs to happen in a complex loop
grid= Grid()
mag_import = int((nmags-1)*tstart/tmax)
print('Initial mag to build field from = ', mag_import)

#Find lbound from magnetogram file
fname = './magnetograms/%04d.nc' % (mag_import)

try:
    data = netcdf_file(fname, 'r', mmap=False)
    print('File', fname, 'found')

except:
    print('File', fname, 'not found')
    raise Exception('Initial boundary data not found')

bz = np.swapaxes(data.variables['bz'][:],0,1)

#Compute initial condition at the start -- but not later on as this will be obtained from a snapshot

print('Calculating initial condition...')

init = compute_initial_condition(grid, bz, run, background_strength = 0.0, background_angle = backfield_angle, boundary_error_limit = 1e-6, init_filename = './inits/init%03d.nc' % init_number)

print('Calculating initial boundary conditions...')

#compute_electrics(run, init_number, omega = omega, start = 0, end = 500)

hrefs = np.load('./hdata/h_ref.npy')
ts    = np.load('./hdata/ts.npy')

halls = []; omegas = []; tplots = []

if hflag < 0.5:
    os.system('make')

def find_minimum(xs, ys):
    #Fits a quadratic function to these three data points and finds the find_minimum
    #Hopefully should be second-order convergence, but we'll see

    print('Finding minimum', xs, ys)
    a = 0
    b = 0
    for i in range(3):
        #Do 'a' first
        den = 1
        for j in range(3):
            if j != i:
                den *= (xs[i] - xs[j])
        a += ys[i]/den

        #Then 'b'
        num = 0
        for j in range(3):
            if j != i:
                num -= xs[j]
        b += ys[i]*num/den

    c = ys[0] - a*xs[0]**2 - b*xs[0]
    x = np.linspace(0,np.max(xs)*1.1,100)
    y = a*x**2 + b*x + c

    plt.plot(x, y)
    plt.scatter(xs, ys)

    print('a', a, 'b', b)

    minx = -b/(2*a)
    plt.scatter(minx, a*minx**2 + b*minx + c)
    plt.close()

    if a > 0:
        return minx

    else:
        return 1e6

nmags_per_run = 5

find_intercept = False

for mag_start in range(start, 500, nmags_per_run):

    mag_end = mag_start + nmags_per_run

    xs = []; ys = []  #For the function interpolation

    go  = True
    stopnext = False
    while go:
        #Find the ideal omega

        if stopnext:
            go = False

        maxomega = 1.0
        minomega = 1e-4

        omega = max(minomega, omega)
        omega = min(maxomega, omega)

        print('Running step', mag_start, 'omega = ', omega)

        compute_electrics(run, init_number, omega = omega, start = mag_start, end = mag_end)

        variables[29] = mag_start
        variables[30] = mag_end

        np.savetxt('parameters/variables%03d.txt' % run, variables)   #variables numbered based on run number (up to 1000)

        #print('Using output directory "%s"' % (data_directory))
        if nprocs <= 4:
            os.system('/usr/lib64/openmpi/bin/mpiexec ffpe-summary=none -np %d ./bin/mf3d %d' % (nprocs, run))
        else:
            os.system('/usr/lib64/openmpi/bin/mpiexec -np %d --oversubscribe ./bin/mf3d %d' % (nprocs, run))

        #Check helicity against reference

        #Want to add the option to use other measures here, so probably shouldn't call it helicity

        if False:
            check = compute_helicity(run, mag_end)
            target = hrefs[mag_end]

        else:
            check = compare_fields(run, mag_end)
            target = 0

        xs.append(omega); ys.append(check - target)

        #THIS ASSUMES A ZERO INTERCEPT RATHER THAN MINIMISING ANYTHING

        if find_intercept:
            if abs(check - target) < 1e0:   #Close enough
                go = False
                print('GOOD ENOUGH',  omega, check-target, xs, ys)

            else:   #not close enough, make a better estimate
                print('NOT GOOD ENOUGH', omega, check-target, xs, ys)
                if len(xs) == 1:
                    if check > target and omega > minomega:
                        omega = omega/1.1
                    elif check < target and omega < maxomega:
                        omega = omega*1.1
                    else:
                        go = False
                else:
                    #Check for anomalies...
                    if (ys[-1] - ys[-2]) / (xs[-1] - xs[-2]) < 0:
                        go = False

                    else:
                        nr_target = xs[-1] - ys[-1]*((xs[-1] - xs[-2])/(ys[-1] - ys[-2]))
                        if nr_target > maxomega:
                            omega = maxomega
                            stopnext = True
                            print('OMEGA BIG', nr_target)
                        elif nr_target < minomega:
                            omega = minomega
                            stopnext = True
                            print('OMEGA SMALL', nr_target)

                        else:
                            print('TARGET',nr_target)
                            omega = nr_target
        else:   #Find the minimum. This is a bit more intense and will always take a little while
            fact = 1.1
            #Estalishing when minimum is found is tricky.
            if len(xs) == 1:  #Only one point -- make a guess
                print('Increasing Omega')
                omega = omega*fact
            elif len(xs) == 2:  #Two points -- make a slightly more informed guess
                if ys[-1] > ys[-2]:
                    print('Further from minimum, trying the other way')
                    omega = omega/(fact**2)
                else:
                    print('Correct drection, extending')

                    omega = omega*fact

            elif len(xs) > 2:  #Three points -- can make a pretty good guess!
                print('Interpolating...')

                min_target = find_minimum(xs[-3:], ys[-3:])

                if min_target > 1e5:
                    print('Mininum not nearby, increasing')
                    omega = omega*fact
                elif min_target < 0:
                    print('Mininum too low, guessing')
                    omega = omega/fact
                else:
                    print('Minimum found, using that')
                    omega = min_target

                print(xs, ys, min_target, omega)

                dydx = abs(ys[-1] - ys[-2])/abs(xs[-1] - xs[-2])
                print('dydx', dydx)
                if dydx < 1e-5:
                    go = False


    halls.append(check)
    omegas.append(omega)
    tplots.append(mag_end*0.5)

    print(mag_end, ts[mag_start + 1], hrefs[mag_start + 1], check)

    np.save('./hdata/halls.npy', halls)
    np.save('./hdata/omegas.npy', omegas)
    np.save('./hdata/tplots.npy', tplots)

#Finish off by running the remaining

mag_start = 499
mag_end = min(mag_start, int(tmax*2))

compute_electrics(run, init_number, omega = omega, start = mag_start, end = mag_end)

variables[29] = mag_start
variables[30] = mag_end

np.savetxt('parameters/variables%03d.txt' % run, variables)   #variables numbered based on run number (up to 1000)

#print('Using output directory "%s"' % (data_directory))
if nprocs <= 4:
    os.system('/usr/lib64/openmpi/bin/mpiexec ffpe-summary=none -np %d ./bin/mf3d %d' % (nprocs, run))
else:
    os.system('/usr/lib64/openmpi/bin/mpiexec -np %d --oversubscribe ./bin/mf3d %d' % (nprocs, run))

#RUN CODE
#-------------------------------------

#os.system('python diagnostics.py %d' % run)

#os.system('rm parameters/variables%03d.txt' % run)
