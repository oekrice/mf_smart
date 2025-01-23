!*******************************************************************************
PROGRAM main
!*******************************************************************************
! Main fortran file. Use this to initialise parameters, grid, produce/read in initial conditions and then run the code
!*******************************************************************************
    USE shared_data
    USE init
    USE mpi_tools
    USE evolve
    USE output
    !USE output, ONLY: save_snap, print_array, diagnostics
    IMPLICIT NONE

    INTEGER:: nstart, nend, init_mag, mag_interval
    REAL(num):: tmags
    !REAL(num):: mag_interval
    ! Put some of the major variables in here - things that can be changed occasionally but not in a series of runs
    cfl  = 0.1
    mf_delta = 1e-4

    ! Import the parameters and set up the grid
    CALL initialise()

    if (.true.) then
    if (hamilton_flag < 0.5) then
        data_directory_root = '/extra/tmp/trcn27/mf3d/'
    else
        data_directory_root = '/nobackup/trcn27/mf3d0/'
    end if

    if (run_number < 10) then
        write (data_directory, "(A23,A2,I1,A1)") data_directory_root, '00',int(run_number), "/"
    else if (run_number < 100) then
        write (data_directory, "(A23,A1,I2,A1)") data_directory_root, '0', int(run_number), "/"
    else
        write (data_directory, "(A23,I3,A1)") data_directory_root, int(run_number), "/"
    end if

    !CALL update_surface_flows(0)

    if (proc_num == 0) print*, 'Set up, running from', mag_min, 'to', mag_max
    t = tstart
    !Adjust so this aligns with a magnetogram.
    !CALL save_snap(0)
    !CALL diagnostics(0)

    !Adjust so magnetic inputs are precisely at timesteps

    if (nmags > 499) then
        tmags = 0.5_num
    else
        tmags = tmax/nmags !Time in between magnetic field inputs
    end if
    !print*, '1', nt, dt

    !print*, 'tmags', tmags, tmags/dt

    !Ensure a whole number of timesteps goes into this time
    mag_interval = int(tmags/dt) + 1
    nt = (int(tmags/dt) + 1)
    dt = tmags/nt


    !print*, tmags/dt, int(tmags/dt) + 1
    !print*, '2', nt, dt, tmags/dt, tmax/dt

    !Set up for ONE mag step
    nstart = mag_min*nt

    nend = mag_max*nt

    if (mag_min == 0) then
        CALL save_snap(0)
        CALL export_magnetogram(0)
    end if

    CALL import_surface_electric(mag_min, 1.0_num/tmags)

    t = nstart*dt
    do n = nstart, nend-1  ! Actually run the code

        if (MOD(n, nt) == 0) then   !Is a multiple -- import new magnetogram

        CALL import_surface_electric(int(n/nt), 1.0_num/tmags)

        end if

        CALL timestep()  !Does everything except the actual timestep (for diagnostic reasons)

        !if (MOD(n, int(nt/int(ndiags-1))) == 0) then   ! Save a snapshot (prints a message as well)
            !if (proc_num == 0) print*, 'Step', n, 'at time', t

            !CALL diagnostics(int(n/(nt/(ndiags-1))))  !WILL BREAK FOR tstart != 0. Should read in such things maybe?

            !print*, 'Max all currents', maxval(abs(jx(0:nx+1, 0:ny,0:nz))), maxval(abs(jy(0:nx, 0:ny+1,0:nz))), maxval(abs(jz(0:nx, 0:ny,0:nz+1)))
        !end if

        if (MOD(n, nt) == nt - 1) then   !Is a multiple -- do the export
            CALL export_magnetogram(int((n+1)/nt))
            CALL save_snap(int((n+1)/nt))

        end if
        ax = ax - dt*ex
        ay = ay - dt*ey
        az = az - dt*ez

        t = t + dt

    end do

    if (proc_num == -1) then
        print*, 'Open Flux', proc_num, z_rank, sum(abs(bz(1:nx,1:ny,nz)))
        print*, 'Max. currents', proc_num, sum(abs(jx(2:nx-2,2:ny-2,2:nz-1))), &
        sum(abs(jy(2:nx-2,2:ny-2,2:nz-1))), sum(abs(jz(2:nx-2,2:ny-2,2:nz-1)))
    end if
    !CALL diagnostics(int(n/(nt/(ndiags-1))))
    if (proc_num == 0) print*, 'Step', n, 'at time', t

    !CALL diagnostics(ndiags-1)
    end if
    if (proc_num == 0 .and. mag_max == 500) print*, 'Fortran code completed sucessfully. Carry on.'
    CALL mpi_finalize(ierr)
    STOP

END PROGRAM main
