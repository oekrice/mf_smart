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
from scipy.optimize import curve_fit

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

#DYNAMIC SYSTEM PARAMETERS
#-------------------------------------
voutfact = 5.0
shearfact = 1.0#3.7e-5   #factor by which to change the imported 'speed'
eta0 = 0.0

tmax = 750.0
tstart = 0.0

start = 0

nx = 128
ny = 128
nz = 128

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
elif not os.path.isdir(data_directory):
    os.mkdir(data_directory)

if not os.path.isdir('./mf_mags/'):
    os.mkdir('./mf_mags/')

#Make magnetogram filenames if necessary
if os.path.isdir('./mf_mags/%03d/' % run) and start == 0:
    for i in range(1000):
        if os.path.isfile('./mf_mags/%03d/%04d.nc' % (run, i)):
            os.remove('./mf_mags/%03d/%04d.nc' % (run, i))
elif not os.path.isdir('./mf_mags/%03d/' % run):
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

nmags_per_run = 1

find_intercept = False

test_range_init = 10   #This should be changed.

test_range = test_range_init

fig_num = 0

omega_init = 1e-2

for mag_start in range(start, 500, nmags_per_run):

    mag_end = mag_start + nmags_per_run

    go  = True

    stopnext = False

    ntests = 5
    stepcount = 0


    min_omega = 1e-4
    max_omega = 1.0

    xs = []; ys = []  #For the function interpolation

    #Put in a two-step process -- initial 5-point quadratic thingummy over the whole range and go from there

    omegas_range = np.geomspace(1e-4,2.5e-2,5) #Or something like that


    for test_omega in omegas_range:
        print('Running step', mag_start, 'omega = ', test_omega)

        compute_electrics(run, init_number, omega = test_omega, start = mag_start, end = mag_end)

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

        check = compare_fields(run, mag_end)
        target = 0

        #Use some 'random' samples and narrow in on it from there.
        #Let's be reasonably smart.

        #1. Establish a range within which lies a minimum
        #2. Narrow that down somehow and do a poly within it. Hmmm...

        #Perhaps. Hard though.
        #Could fit a poly to it and find the minimum of that? Would be smooth at least

        xs.append(test_omega); ys.append(check - target)

    print('Initial step', xs, ys)

    #Fit a curve to this -- if quite small then the minimum might not be so
    def func(x, a, b, c): #Parameters to be checked
        lx = x#np.log(x)   #I think this is probably more reasonable
        return a + b*lx + c*lx**2

    popt, pcov = curve_fit(func, xs[:], ys[:])

    xs_all = np.linspace(min(xs[:]), max(xs[:]), 500)
    ys_all = func(xs_all, *popt)

    #Check for any variation at all
    if np.max(ys_all) - np.min(ys_all) < 1e-3:
        omega_init = 1e-2

    #Check if minimum is at either end, if so just use that
    elif np.abs(ys_all[0] - np.min(ys_all)) < 1e-3:
        omega_init = xs_all[0]

    elif np.abs(ys_all[-1] - np.min(ys_all)) < 1e-3:
        omega_init = xs_all[-1]

    else: #Is somewhere in the middle, so do a nice interpolation
        print('Minimum found somewhere...')
        min_index = np.where(ys_all == np.min(ys_all))[0][0]
        omega_init = xs_all[min_index]

        #Create range based on this
        omegas_range = np.geomspace(omega_init/2.0,omega_init*2.0,5)

        for test_omega in omegas_range:
            print('Running step', mag_start, 'omega = ', test_omega)

            compute_electrics(run, init_number, omega = test_omega, start = mag_start, end = mag_end)

            variables[29] = mag_start
            variables[30] = mag_end

            np.savetxt('parameters/variables%03d.txt' % run, variables)   #variables numbered based on run number (up to 1000)

            if nprocs <= 4:
                os.system('/usr/lib64/openmpi/bin/mpiexec ffpe-summary=none -np %d ./bin/mf3d %d' % (nprocs, run))
            else:
                os.system('/usr/lib64/openmpi/bin/mpiexec -np %d --oversubscribe ./bin/mf3d %d' % (nprocs, run))

            check = compare_fields(run, mag_end)
            target = 0

            xs.append(test_omega); ys.append(check - target)

        popt, pcov = curve_fit(func, xs[:], ys[:])

        xs_all = np.linspace(min(xs[:]), max(xs[:]), 500)
        ys_all = func(xs_all, *popt)

        min_index = np.where(ys_all == np.min(ys_all))[0][0]
        omega_init = xs_all[min_index]  #Maybe that will suffice...


        xs_all = np.linspace(min(xs[:]), max(xs[:]), 100)
        ys_all = func(xs_all, *popt)
        plt.scatter(xs[:], ys[:])
        plt.scatter(xs[:], ys[:], c = 'red')

        plt.plot(xs_all, ys_all)
        plt.title(mag_start)
        #plt.xscale('log')
        plt.savefig('testfigs/fig%d.png' % fig_num)
        plt.close()
        fig_num += 1

    print('Optimum value found:', omega_init)
    #Run final step with this optimum value

    compute_electrics(run, init_number, omega = omega_init, start = mag_start, end = mag_end)

    variables[29] = mag_start
    variables[30] = mag_end

    np.savetxt('parameters/variables%03d.txt' % run, variables)   #variables numbered based on run number (up to 1000)

    if nprocs <= 4:
        os.system('/usr/lib64/openmpi/bin/mpiexec ffpe-summary=none -np %d ./bin/mf3d %d' % (nprocs, run))
    else:
        os.system('/usr/lib64/openmpi/bin/mpiexec -np %d --oversubscribe ./bin/mf3d %d' % (nprocs, run))

    hel = compute_helicity(run, mag_end)

    halls.append(hel)
    omegas.append(omega_init)
    tplots.append(mag_end*0.5)

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
