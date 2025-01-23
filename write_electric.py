#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27

Generates .nc files of the electric field, to be read-in to the Fortran code.
Needs to deal with different resolutions as the inputs are all the same (192 I think)
Read in the 'magnetograms' from the .nc files
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib
from scipy.io import netcdf_file
from scipy.interpolate import RegularGridInterpolator
from scipy.fft import fft, ifft2, fft2, ifft
import os
from scipy.ndimage import gaussian_filter

class Grid():
    """In the interest of doing it properly, put grid parameters in here"""
    def __init__(self, run):
        
        paras = np.loadtxt('parameters/variables%03d.txt' % run)

        import_resolution = 128

        #Define the grid onto which the electric field should be outputted
        self.nx = int(paras[1])
        self.ny = int(paras[2])

        self.x0 = paras[12]; self.x1 = paras[13]
        self.y0 = paras[14]; self.y1 = paras[15]
                
        self.xs = np.linspace(self.x0, self.x1, self.nx+1)
        self.ys = np.linspace(self.y0, self.y1, self.ny+1)
        
        self.dx = self.xs[1] - self.xs[0]
        self.dy = self.ys[1] - self.ys[0]
        
        self.xc = np.linspace(self.x0 - self.dx/2, self.x1 + self.dx/2, self.nx+2)
        self.yc = np.linspace(self.y0 - self.dy/2, self.y1 + self.dy/2, self.ny+2)

        self.dx = self.xs[1] - self.xs[0]
        self.dy = self.ys[1] - self.ys[0]
        
        self.xs_import = np.linspace(self.x0, self.x1, import_resolution+1)
        self.ys_import = np.linspace(self.y0, self.y1, import_resolution+1)

        
        self.dx_import = self.xs[1] - self.xs[0]
        self.dy_import = self.ys[1] - self.ys[0]
        
        self.xc_import = np.linspace(self.x0 - self.dx/2, self.x1 + self.dx/2, import_resolution+2)
        self.yc_import = np.linspace(self.y0 - self.dy/2, self.y1 + self.dy/2, import_resolution+2)

        
    def lap_xribs(self, xribs):
        """Calculates laplacian of a quantity saved on the x ribs, using staggered as appropriate"""
        lap_x = np.zeros((self.nx+2, self.ny+1))
        #x direction
        lap_x[1:-1,1:-1] += (xribs[2:,1:-1] -  2*xribs[1:-1,1:-1] + xribs[:-2,1:-1])/self.dx**2
        lap_x[1:-1,1:-1] += (xribs[1:-1,2:] -  2*xribs[1:-1,1:-1] + xribs[1:-1,:-2])/self.dy**2
        return lap_x

    def lap_yribs(self, yribs):
        """Calculates laplacian of a quantity saved on the y ribs, using staggered as appropriate"""
        lap_y = np.zeros((self.nx+1, self.ny+2))
        #x direction
        lap_y[1:-1,1:-1] += (yribs[2:,1:-1] -  2*yribs[1:-1,1:-1] + yribs[:-2,1:-1])/self.dx**2
        lap_y[1:-1,1:-1] += (yribs[1:-1,2:] -  2*yribs[1:-1,1:-1] + yribs[1:-1,:-2])/self.dy**2
        return lap_y


    def curl_inplane(self, C):
        """Calculates the in-plane curl of the quantity C(x,y)"""
        """Outputs on the respective ribs, hopefully"""
        
        curl_x = np.zeros((self.nx+2, self.ny+1))
        curl_y = np.zeros((self.nx+1, self.ny+2))

        curl_x = (C[:,1:] - C[:,:-1])/self.dy
        curl_y = -(C[1:,:] - C[:-1,:])/self.dx
        
        return curl_x, curl_y

    def curl_E(self, E_x, E_y):
        """Returns the in-plane curl of the vectors E"""
        curl = np.zeros((self.nx+2, self.ny+2))
        curl[1:-1,1:-1] += (E_x[1:-1,1:] - E_x[1:-1,:-1])/self.dy
        curl[1:-1,1:-1] -= (E_y[1:,1:-1] - E_y[:-1,1:-1])/self.dx

        return curl
        

    def lap_points(self, points):
        """Calculates laplacian of a quantity saved on the y ribs, using staggered as appropriate"""
        lap_p = np.zeros((self.nx+1, self.ny+1))
        #x direction
        lap_p[1:-1,1:-1] += (points[2:,1:-1] -  2*points[1:-1,1:-1] + points[:-2,1:-1])/self.dx**2
        lap_p[1:-1,1:-1] += (points[1:-1,2:] -  2*points[1:-1,1:-1] + points[1:-1,:-2])/self.dy**2
        return lap_p

    def lap_centres(self, points):
        """Calculates laplacian of a quantity saved on the y ribs, using staggered as appropriate"""
        lap_p = np.zeros((self.nx+2, self.ny+2))
        #x direction
        lap_p[1:-1,1:-1] += (points[2:,1:-1] -  2*points[1:-1,1:-1] + points[:-2,1:-1])/self.dx**2
        lap_p[1:-1,1:-1] += (points[1:-1,2:] -  2*points[1:-1,1:-1] + points[1:-1,:-2])/self.dy**2
        return lap_p

    def div_E(self, E_x, E_y):
        """Returns the in-plane curl of the vectors E"""
        div = np.zeros((self.nx+1, self.ny+1))
        div[:,:] += (E_x[1:,:] - E_x[:-1,:])/self.dx
        div[:,:] += (E_y[:,1:] - E_y[:,:-1])/self.dy

        return div
    
    def grad(self, phi):
        """Returns gradients of phi (at grid points)"""
        grad_x = np.zeros((self.nx+2, self.ny+1))
        grad_y = np.zeros((self.nx+1, self.ny+2))
        grad_x[1:-1,:] = (phi[1:,:] - phi[:-1,:])/self.dx
        grad_y[:,1:-1] = (phi[:,1:] - phi[:,:-1])/self.dy
        
        return grad_x, grad_y
        
        
class FT():
    
    """Creates various quantities related to the Fourier transform (mutliplication matrix etc.)"""
    def __init__(self,grid):
        self.grid = grid
        pass
        
        
    def point_transform(self, rhs_points):
        """2D transform in at the grid points"""
        rhs_transform = ifft2(rhs_points)  
        
        j = np.arange(self.grid.nx+1)
        k = np.arange(self.grid.ny+1)

        d2Xdx2 = (2*np.cos(2*np.pi*j/(self.grid.nx+1)) - 2)/self.grid.dx**2
        d2Ydy2 = (2*np.cos(2*np.pi*k/(self.grid.ny+1)) - 2)/self.grid.dy**2
        
        d2 = d2Xdx2.reshape((self.grid.nx+1,1)) + d2Ydy2
        
        d2[0,0] = 1.0   #Avoid problems with the zero value. All this does is add a constant anyway so set to zero in the transform
        
        a_twiddle = rhs_transform/d2
        
        a_twiddle[0,0] = 0.0
        
        phi = (fft2(a_twiddle)).real
        
        return phi

    def centre_transform(self, rhs_centres):
        """2D transform in at the grid centres"""
        rhs_transform = ifft2(rhs_centres)  
        
        j = np.arange(self.grid.nx+2)
        k = np.arange(self.grid.ny+2)

        d2Xdx2 = (2*np.cos(2*np.pi*j/(self.grid.nx+2)) - 2)/self.grid.dx**2
        d2Ydy2 = (2*np.cos(2*np.pi*k/(self.grid.ny+2)) - 2)/self.grid.dy**2
        
        d2 = d2Xdx2.reshape((self.grid.nx+2,1)) + d2Ydy2
        
        d2[0,0] = 1.0   #Avoid problems with the zero value. All this does is add a constant anyway so set to zero in the transform
        
        a_twiddle = rhs_transform/d2
        
        a_twiddle[0,0] = 0.0
        
        G = (fft2(a_twiddle)).real
        
        return G

class compute_electrics():
    
    def __init__(self, run, init_number, omega = 0., start = 0, end = 500):

        mag_dt = 0.5
        t_ahead = 50.0 #Take the derivate from some time ahead. To be determined kind of randomly

        grid = Grid(run)  #Establish grid (on new scales, not 192)
        
        omega_init = omega
        data_directory = './magnetograms/'
        
        if len(sys.argv) > 1:
            run = int(sys.argv[1])
        else:
            run = 0
        
        import_resolution = 128
        
        paras = np.loadtxt('parameters/variables%03d.txt' % run)
        
        if not os.path.exists('efields'):
            os.mkdir('efields')

        if not os.path.exists('efields/%03d' % init_number):
            os.mkdir('efields/%03d' % init_number)

        hlog_exists = False
        #See if an existing log exists
        try:
            hdiff_log = np.load('hdiffs.npy')
            hlog = np.load('hs.npy')

            hlog = gaussian_filter(hlog, 20.0)

            if len(hdiff_log) > 400:
                hlog_exists = True

            hdiff_log = (hlog[1:] - hlog[:-1])/mag_dt

            hdiff_log = np.concatenate((hdiff_log, np.array([0])))

            #plt.plot(hlog)

            #plt.plot(hdiff_log)
            #plt.show()
        except:
            pass

        if not hlog_exists:
            hdiff_log = []   #Derivateives of the in-plane helicity, so it doesn't need to be done lots
            hlog = []

        for snap in range(start, end):
            t_int = int(t_ahead/mag_dt)

            #print('Importing fields', snap, 'and', snap + 1)
            bfield_fname = '%s%04d.nc' % (data_directory, snap)
            efield_fname = '%s%03d/%04d.nc' % ('./efields/', init_number, snap)
        
            try:
                data = netcdf_file(bfield_fname, 'r', mmap=False)
                #print('File', bfield_fname, 'found')
        
            except:
                #print('File', bfield_fname, 'not found')
                continue
        
            bfield1 = data.variables['bz'][:]
            
            bfield_fname = '%s%04d.nc' % (data_directory, snap + 1)
        
            try:
                data = netcdf_file(bfield_fname, 'r', mmap=False)
                #print('File', bfield_fname, 'found')
        
            except:
                #print('File', bfield_fname, 'not found')
                continue
        
            bfield2 = data.variables['bz'][:]
        
            bfield1 = gaussian_filter(bfield1, sigma = 1.0)
            bfield2 = gaussian_filter(bfield2, sigma = 1.0)
            
            #Filter the magnetic fields a bit to stop instabilities
            
            #plt.pcolormesh(bfield2-bfield1)
            #plt.show()
            
            diff_init = bfield2 - bfield1   #Difference between the magnetic field at import resolution
            
            X, Y = np.meshgrid(grid.xc, grid.yc, indexing = 'ij')
        
            diff_fn = RegularGridInterpolator((grid.xc_import[1:-1], grid.yc_import[1:-1]), diff_init, bounds_error = False, method = 'linear', fill_value = None)
        
            diff = diff_fn((X,Y))   #Difference now interpolated to the new grid
                
            ft = FT(grid)
            G = ft.centre_transform(-diff)
            ex, ey = grid.curl_inplane(G)
            curl_test = grid.curl_E(ex, ey)
        
            #print('Curl test', np.max(np.abs(curl_test[1:-1,1:-1] + diff[1:-1,1:-1])))
        
            #Define distribution of the divergence of the electric field. 
            #Following Cheung and De Rosa, just proportional to the vertical field
            if True:   #'Twist' proportional to magnetic field

                X, Y = np.meshgrid(grid.xs, grid.ys, indexing = 'ij')

                bf_fn = RegularGridInterpolator((grid.xc_import[1:-1], grid.yc_import[1:-1]), 0.5 * (bfield1 + bfield2), bounds_error = False, method = 'linear', fill_value = None)

                bf = bf_fn((X,Y))   #Difference now interpolated to the new grid

                D = omega*(bf)

                #print('OMEGA', omega)
                # plt.pcolormesh(bf)
                # plt.colorbar()
                # plt.show()
            elif True:
                #Pattern proportional to magnetic field but overall magnitude with rate of change in helicity.
                #Could also just try overall helicity as well?

                #Calculate rate of change in Helicity over these timesteps
                if not hlog_exists:
                    bfield_fname = '%s%04d.nc' % (data_directory, snap)

                    try:
                        data = netcdf_file(bfield_fname, 'r', mmap=False)
                    except:
                        continue

                    bx = np.swapaxes(data.variables['bx'][:],0,1)
                    by = np.swapaxes(data.variables['by'][:],0,1)
                    bz = np.swapaxes(data.variables['bz'][:],0,1)

                    mag_nx = bz.shape[0]; mag_ny = bz.shape[1]
                    mag_dx = (grid.xs[-1] - grid.xs[0])/mag_nx
                    mag_dy = (grid.ys[-1] - grid.ys[0])/mag_ny

                    hfield1 = self.compute_inplane_helicity(grid, bx, by, bz)

                    hsum1 = np.sum(np.abs(hfield1)*mag_dx*mag_dy)

                    bfield_fname = '%s%04d.nc' % (data_directory, snap + 1)

                    data.close()
                    try:
                        data = netcdf_file(bfield_fname, 'r', mmap=False)
                    except:
                        continue

                    bx = np.swapaxes(data.variables['bx'][:],0,1)
                    by = np.swapaxes(data.variables['by'][:],0,1)
                    bz = np.swapaxes(data.variables['bz'][:],0,1)

                    data.close()

                    hfield2 = self.compute_inplane_helicity(grid, bx, by, bz)

                    hsum2 = np.sum(np.abs(hfield2)*mag_dx*mag_dy)

                    dhdt = (hsum2 - hsum1)/mag_dt

                    print('Helicity change', dhdt, hsum1, hsum2)

                    hdiff_log.append(dhdt)
                    hlog.append(hsum1)

                else:
                    try:
                        dhdt = hdiff_log[snap+ t_int]
                    except:
                        dhdt = 0.

                X, Y = np.meshgrid(grid.xs, grid.ys, indexing = 'ij')

                bf_fn = RegularGridInterpolator((grid.xc_import[1:-1], grid.yc_import[1:-1]), 0.5 * (bfield1 + bfield2), bounds_error = False, method = 'linear', fill_value = None)

                bf = bf_fn((X,Y))   #Difference now interpolated to the new grid

                bfield_fname = '%s%04d.nc' % (data_directory, snap)

                D = omega_init*bf #omega_init*dhdt*(bf)

            else:   #'Twist' proportional to in-plane helicity
                bfield_fname = '%s%04d.nc' % (data_directory, snap)

                try:
                    data = netcdf_file(bfield_fname, 'r', mmap=False)
                except:
                    continue
                bx = np.swapaxes(data.variables['bx'][:],0,1)
                by = np.swapaxes(data.variables['by'][:],0,1)
                bz = np.swapaxes(data.variables['bz'][:],0,1)
                hfield1 = self.compute_inplane_helicity(grid, bx, by, bz)

                bfield_fname = '%s%04d.nc' % (data_directory, snap + 1)

                try:
                    data = netcdf_file(bfield_fname, 'r', mmap=False)

                except:
                    continue

                bx = np.swapaxes(data.variables['bx'][:],0,1)
                by = np.swapaxes(data.variables['by'][:],0,1)
                bz = np.swapaxes(data.variables['bz'][:],0,1)
                hfield2 = self.compute_inplane_helicity(grid, bx, by, bz)

                X, Y = np.meshgrid(grid.xs, grid.ys, indexing = 'ij')

                hf_fn = RegularGridInterpolator((grid.xc_import[1:-1], grid.yc_import[1:-1]), 0.5*(hfield1 + hfield2), bounds_error = False, method = 'linear', fill_value = None)

                hf = hf_fn((X,Y))   #Difference now interpolated to the new grid

                D = omega*hf

                # plt.pcolormesh(hf)
                # plt.colorbar()
                # plt.show()

            div_test = grid.div_E(ex, ey)
            phi = ft.point_transform(-div_test + D)
            correct_x, correct_y = grid.grad(phi)
            ex += correct_x
            ey += correct_y
        
            div_test = grid.div_E(ex, ey)
            #print('Div Test', np.max(np.abs(div_test[1:-1,1:-1] - D[1:-1,1:-1])))
        
            curl_test = grid.curl_E(ex, ey)
            
            #print('Overall', np.max(np.abs(curl_test[1:-1,1:-1] + diff[1:-1,1:-1])))
        
            #Swap sign of E, as that seems to be the way forward.
            ex = -ex
            ey = -ey

            if False:
                plt.pcolormesh(ey)
                plt.savefig('plots/ey%d.png' % snap)
                plt.close()
                
                plt.pcolormesh(ex)
                plt.savefig('plots/ex%d.png' % snap)
                plt.show()

            fid = netcdf_file(efield_fname, 'w')
            fid.createDimension('xs', grid.nx+1)
            fid.createDimension('ys', grid.ny+1)
            fid.createDimension('xc', grid.nx+2)
            fid.createDimension('yc', grid.ny+2)
        
            vid = fid.createVariable('xs', 'd', ('xs',))
            vid[:] = grid.xs
            vid = fid.createVariable('ys', 'd', ('ys',))
            vid[:] = grid.ys
        
            vid = fid.createVariable('xc', 'd', ('xc',))
            vid[:] = grid.xc
            vid = fid.createVariable('yc', 'd', ('yc',))
            vid[:] = grid.yc
        
            if False:
                plt.pcolormesh(grid.xc, grid.yc, bf.T)

                #plt.pcolormesh(grid.xc, grid.ys, ex.T)
                #plt.pcolormesh(grid.xs, grid.yc, ey.T)

                plt.show()
            #Transposes are necessary as it's easier to flip here than in Fortran
            #Still doesn't seem to like it -- sums are all correct but going in with the wrong direction.

            vid = fid.createVariable('ex', 'd', ('ys','xc'))
            vid[:] = np.swapaxes(ex, 0, 1)
            vid = fid.createVariable('ey', 'd', ('yc','xs'))
            vid[:] = np.swapaxes(ey, 0, 1)
        
            fid.close()

        if len(hdiff_log) > 400:
            np.save('hdiffs.npy', np.array(hdiff_log))
            np.save('hs.npy', np.array(hlog))

    def compute_inplane_helicity(self, grid, bx, by, bz):
        #Need to average these to grid centres to get the FFT to work
        bx0 = 0.5*(bx[:,1:] + bx[:,:-1])
        by0 = 0.5*(by[1:,:] + by[:-1,:])
        bz0 = bz[:,:]

        def norm2d(vec):
            mag = np.linalg.norm(vec)
            if (mag > 0.0):
                v = vec/mag
            else:
                v = np.array([0, 0])
            return np.array([v[0],v[1],0.0])

        def getFrequencyMatrix(ncells,spacing):
            freqlist1da =np.roll(np.linspace(-ncells[0]/2,ncells[0]/2-1,ncells[0]),round(ncells[0]/2))/(ncells[0]*spacing[0])
            freqlist1db =np.roll(np.linspace(-ncells[1]/2,ncells[1]/2-1,ncells[1]),round(ncells[1]/2))/(ncells[1]*spacing[1])
            return np.array([np.array([np.array([2.0*np.pi*freqlist1da[i],2.0*np.pi*freqlist1db[j]]) for j in range(len(freqlist1db))]) for i  in range(len(freqlist1da))]);

        #def curla(ax, ay, az):
            #Outputs the curl of a (which is averaged to grid points so a bit messy)

        #Find in -plane vector potential in the winding gauge

        fm = getFrequencyMatrix([bz.shape[0], bz.shape[1]],[grid.dx, grid.dy]);
        # make the basis

        kparr = np.array([np.array([norm2d(fm[i][j]) for j in range(len(fm[0]))]) for i  in range(len(fm))]);
        kperp = np.array([np.array([np.array([-kparr[i][j][1],kparr[i][j][0],0.0]) for j in range(len(fm[0]))]) for i  in range(len(fm))])
        # note in the k matrix below the k=0 element is set to one so we can divide by it.
        k = np.array([np.array([1.0 if i==j==0 else np.linalg.norm(fm[i][j]) for i in range(len(fm))]) for j  in range(len(fm[0]))]).T

        nx = bz.shape[0]; ny = bz.shape[1]
        aftx = np.zeros([bz.shape[0],bz.shape[1]],dtype=np.complex128)
        afty = np.zeros([bz.shape[0],bz.shape[1]],dtype=np.complex128)
        aftz = np.zeros([bz.shape[0],bz.shape[1]],dtype=np.complex128)

        fbx = fft2(bx0[:,:]); fby = fft2(by0[:,:]); fbz = fft2(bz0[:,:])

        akperp = -1j*fbz/k
        ## fix i =j  element
        akw = 1j*(-(kparr[:,:,1])*fbx + (kparr[:,:,0])*fby)/k
        ## fix i =j  element
        aftx[:,:] = akperp*kperp[:,:,0]
        afty[:,:] = akperp*kperp[:,:,1]
        aftz[:,:] = akperp*kperp[:,:,2]+akw

        ax0 = ifft2(aftx[:,:])
        ay0 = ifft2(afty[:,:])
        az0 = ifft2(aftz[:,:])
        ax0 = np.real(ax0)
        ay0 = np.real(ay0)
        az0 = np.real(az0)

        ax = np.zeros((nx, ny+1))
        ay = np.zeros((nx+1, ny))

        ax[:,1:-1] = 0.5*(ax0[:,1:] + ax0[:,:-1])
        ay[1:-1,:] = 0.5*(ay0[1:,:] + ay0[:-1,:])

        bz_test = (ay[1:,:] - ay[:-1,:])/grid.dx - (ax[:,1:] - ax[:,:-1])/grid.dy
        #This vector potential should be reasonably OK... Need code to test though

        #Should be proportional to the magnetic field strength, so this helicity requires a square root. I'm pretty sure the scaling is OK here...
        hfield = np.sqrt(np.abs(ax0*bx0 + ay0*by0 + az0*bz0))

        return hfield

