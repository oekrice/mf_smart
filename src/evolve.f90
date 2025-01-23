!*******************************************************************************

MODULE evolve
!*******************************************************************************
! Initialise the simulation. Read in/compute the initial condition (python code for that I think?), read in the parameters and all will be well
!*******************************************************************************
    USE shared_data
    USE mpi_tools
    USE boundary
    USE pressure
    USE netcdf

    !USE output
    IMPLICIT NONE

!*******************************************************************************
CONTAINS

SUBROUTINE timestep()

    !New timestep with a quasi-implicit bit. Using a predictor for the magnfield of the new bit but not the old
    CALL calculate_magnetic()

    !CALL check_solenoidal()
    CALL calculate_jp()

    CALL calculate_current()

    CALL j_to_gridpts()
    CALL b_to_gridpts()

    CALL calculate_velocity()

    CALL calculate_pressure()

    CALL calculate_electric()

    CALL MPI_Barrier(comm,ierr)  !Wait for t to be broadcast everywhere.

END SUBROUTINE timestep


SUBROUTINE calculate_magnetic()

    IMPLICIT NONE

    !Determine the interior points from the vector potential.
    !No boundary conditions here (do that in a bit)

    ! INTERIOR POINTS (DON'T USE INFORMATION FROM A THAT DOESN'T EXIST)

    bx(0:nx, 1:ny,1:nz) = (az(0:nx,1:ny,1:nz) - az(0:nx, 0:ny-1,1:nz))/dy - (ay(0:nx,1:ny,1:nz) - ay(0:nx,1:ny,0:nz-1))/dz

    by(1:nx, 0:ny,1:nz) = (ax(1:nx,0:ny,1:nz) - ax(1:nx, 0:ny,0:nz-1))/dz - (az(1:nx,0:ny,1:nz) - az(0:nx-1,0:ny,1:nz))/dx

    bz(1:nx, 1:ny,0:nz) = (ay(1:nx,1:ny,0:nz) - ay(0:nx-1, 1:ny,0:nz))/dx - (ax(1:nx,1:ny,0:nz) - ax(1:nx,0:ny-1,0:nz))/dy

    CALL bfield_mpi
    CALL magnetic_boundary

    if (n == 0) bz_surf_reference(0:nx+1,0:ny+1) = bz(0:nx+1,0:ny+1,0)  !Save reference lower boundary field to stop annoying instabilities due to lack of upwinding

END SUBROUTINE calculate_magnetic

SUBROUTINE calculate_current()

    IMPLICIT NONE
    !Determine the current from the magnetic field (after boundary conditions etc.)

    jx(0:nx+1, 0:ny,0:nz) = (bz(0:nx+1,1:ny+1,0:nz) - bz(0:nx+1, 0:ny,0:nz))/dy - (by(0:nx+1,0:ny,1:nz+1) - by(0:nx+1,0:ny,0:nz))/dz

    jy(0:nx, 0:ny+1,0:nz) = (bx(0:nx,0:ny+1,1:nz+1) - bx(0:nx, 0:ny+1,0:nz))/dz - (bz(1:nx+1,0:ny+1,0:nz) - bz(0:nx,0:ny+1,0:nz))/dx

    jz(0:nx, 0:ny,0:nz+1) = (by(1:nx+1,0:ny,0:nz+1) - by(0:nx, 0:ny,0:nz+1))/dx - (bx(0:nx,1:ny+1,0:nz+1) - bx(0:nx,0:ny,0:nz+1))/dy

    jx(0:nx+1, 0:ny,0:nz) = jx(0:nx+1, 0:ny,0:nz) - jpx(0:nx+1, 0:ny,0:nz)
    jy(0:nx, 0:ny+1,0:nz) = jy(0:nx, 0:ny+1,0:nz) - jpy(0:nx, 0:ny+1,0:nz)

END SUBROUTINE calculate_current

SUBROUTINE j_to_gridpts
    !Averages the current field to raw grid points
    !Should only need to average in one direction
    IMPLICIT NONE
    jx1(0:nx,0:ny,0:nz) = 0.5_num*(jx(1:nx+1,0:ny,0:nz) + jx(0:nx,0:ny,0:nz))
    jy1(0:nx,0:ny,0:nz) = 0.5_num*(jy(0:nx,1:ny+1,0:nz) + jy(0:nx,0:ny,0:nz))
    jz1(0:nx,0:ny,0:nz) = 0.5_num*(jz(0:nx,0:ny,1:nz+1) + jz(0:nx,0:ny,0:nz))

END SUBROUTINE j_to_gridpts

SUBROUTINE b_to_gridpts
    !Averages the magnetic field to raw grid points
    !Need to average in two dimensions
    IMPLICIT NONE
    bx1(0:nx,0:ny,0:nz) = 0.25_num*(bx(0:nx,0:ny,0:nz) + bx(0:nx,1:ny+1,0:nz) + bx(0:nx,0:ny,1:nz+1) + bx(0:nx,1:ny+1,1:nz+1))
    by1(0:nx,0:ny,0:nz) = 0.25_num*(by(0:nx,0:ny,0:nz) + by(1:nx+1,0:ny,0:nz) + by(0:nx,0:ny,1:nz+1) + by(1:nx+1,0:ny,1:nz+1))
    bz1(0:nx,0:ny,0:nz) = 0.25_num*(bz(0:nx,0:ny,0:nz) + bz(1:nx+1,0:ny,0:nz) + bz(0:nx,1:ny+1,0:nz) + bz(1:nx+1,1:ny+1,0:nz))

END SUBROUTINE b_to_gridpts

SUBROUTINE calculate_velocity
    !Calculates the magnetofrictional velocity
    IMPLICIT NONE
    b2 = bx1**2 + by1**2 + bz1**2 !B squared

    nu(:,:,:) = nu0

    if (abs(mf_delta) < 1e-10) then !No softening
        soft = b2
    else !Softening.
        soft = b2 + mf_delta*exp(-b2/mf_delta)
    end if

    vx = nu*(jy1*bz1 - jz1*by1)/soft
    vy = nu*(jz1*bx1 - jx1*bz1)/soft
    vz = nu*(jx1*by1 - jy1*bx1)/soft

    if (z_down < 0) then
        vx(:,:,0) = 0.0_num; vy(:,:,0) = 0.0_num; vz(:,:,0) = 0.0_num
    end if

END SUBROUTINE calculate_velocity

SUBROUTINE add_boundary_flows()

    !Adds twisting directly onto the electric field. Only affects ez (for the jet model, at least)
    IMPLICIT NONE

    if (z_down < 0 .and. shearfact > 1d-6) then

      vx(0:nx,0:ny,0) = shearfact*surf_vx(0:nx,0:ny)
      vy(0:nx,0:ny,0) = shearfact*surf_vy(0:nx,0:ny)
      vz(0:nx,0:ny,0) = shearfact*surf_vz(0:nx,0:ny)

    end if

END SUBROUTINE add_boundary_flows


SUBROUTINE check_solenoidal()
    !Checks the solenoidal condition by calculating the divergence of all raw grid cells
    IMPLICIT NONE

    real(num), dimension(:,:,:):: div(0:nx+1,0:ny+1,0:nz+1)
    div = 0.0_num
    div(0:nx+1,0:ny+1,0:nz+1) = div(0:nx+1,0:ny+1,0:nz+1) + dx*dy*(bz(0:nx+1,0:ny+1,-1:nz) - bz(0:nx+1,0:ny+1,0:nz+1))
    div(0:nx+1,0:ny+1,0:nz+1) = div(0:nx+1,0:ny+1,0:nz+1) + dy*dz*(bx(-1:nx,0:ny+1,0:nz+1) - bx(0:nx+1,0:ny+1,0:nz+1))
    div(0:nx+1,0:ny+1,0:nz+1) = div(0:nx+1,0:ny+1,0:nz+1) + dx*dz*(by(0:nx+1,-1:ny,0:nz+1) - by(0:nx+1,0:ny+1,0:nz+1))

    print*, 'Max divergence', maxval(abs(div(0:nx+1,0:ny+1,0:nz+1)))

END SUBROUTINE check_solenoidal

SUBROUTINE calculate_electric()

    !Calculates the electric field - resistivity, magnetofriction and boundary effects
    IMPLICIT NONE

    ex = 0.0; ey = 0.0; ez = 0.0

    if (eta > 0) then
        !Determine the current from the magnetic field (after boundary conditions etc.)
        ex(1:nx, 0:ny,0:nz) = ex(1:nx, 0:ny,0:nz) + eta*(jx(1:nx, 0:ny,0:nz))! - jpx(1:nx, 0:ny,0:nz))
        ey(0:nx, 1:ny,0:nz) = ey(0:nx, 1:ny,0:nz) + eta*(jy(0:nx, 1:ny,0:nz))! - jpy(0:nx, 1:ny,0:nz))
        ez(0:nx, 0:ny,1:nz) = ez(0:nx, 0:ny,1:nz) + eta*jz(0:nx, 0:ny,1:nz)
    end if

    !Add shearing (if necessary) directly onto this (averaged) field
    CALL add_boundary_flows()

    ex1 = vz*by1 - vy*bz1
    ey1 = vx*bz1 - vz*bx1
    ez1 = vy*bx1 - vx*by1

    !Average to Ribs (interior only):
    ex(1:nx,0:ny,0:nz) = ex(1:nx,0:ny,0:nz)  + 0.5_num*(ex1(0:nx-1,0:ny,0:nz) + ex1(1:nx,0:ny,0:nz))
    ey(0:nx,1:ny,0:nz) = ey(0:nx,1:ny,0:nz)  + 0.5_num*(ey1(0:nx,0:ny-1,0:nz) + ey1(0:nx,1:ny,0:nz))
    ez(0:nx,0:ny,1:nz) = ez(0:nx,0:ny,1:nz)  + 0.5_num*(ez1(0:nx,0:ny,0:nz-1) + ez1(0:nx,0:ny,1:nz))

    !Add outflow (if necessary) directly onto this field
    if (voutfact > 0) then
    ex(1:nx,0:ny,0:nz) = ex(1:nx,0:ny,0:nz) + voutx(1:nx,0:ny,0:nz)*by(1:nx,0:ny,0:nz)
    ey(0:nx,1:ny,0:nz) = ey(0:nx,1:ny,0:nz) - vouty(0:nx,1:ny,0:nz)*bx(0:nx,1:ny,0:nz)
    end if
    
    !Add electric field loaded in from elsewhere
    if (z_rank == 0) then
        ex(1:nx,0:ny,0) = surf_ex(1:nx,0:ny)
        ey(0:nx,1:ny,0) = surf_ey(0:nx,1:ny)
        
    end if


END SUBROUTINE calculate_electric

SUBROUTINE import_surface_electric(flow_number, dt_fact)

    !Imports the velocity field from the imported DAVE magnetogram
    IMPLICIT NONE

    INTEGER:: flow_number

    CHARACTER(LEN =64):: electric_filename
    CHARACTER(LEN = 4):: flow_id
    CHARACTER(LEN = 4):: run_id

    INTEGER:: ncid, vid
    REAL(num):: dt_fact

    if (flow_number < 499) then
        write (flow_id,'(I4.4)') flow_number
        write (run_id,'(I3.3)') init_number

        electric_filename = trim("./efields/"//trim(run_id)//'/'//trim(flow_id)//'.nc')

        call try(nf90_open(trim(electric_filename), nf90_nowrite, ncid))

        call try(nf90_inq_varid(ncid, 'ex', vid))
        call try(nf90_get_var(ncid, vid, surf_ex(0:nx+1,0:ny), &
        start = (/x_rank*nx+1,y_rank*ny+1/),count = (/nx+2,ny+1/)))

        call try(nf90_inq_varid(ncid, 'ey', vid))
        call try(nf90_get_var(ncid, vid, surf_ey(0:nx,0:ny+1), &
        start = (/x_rank*nx+1,y_rank*ny+1/),count = (/nx+1,ny+2/)))

        call try(nf90_close(ncid))

        surf_ex = surf_ex*dt_fact
        surf_ey = surf_ey*dt_fact

    else
        surf_ex = 0.0_num
        surf_ey = 0.0_num
    end if

END SUBROUTINE import_surface_electric

SUBROUTINE update_surface_flows(flow_number)

    !Imports the velocity field from the imported DAVE magnetogram
    IMPLICIT NONE

    INTEGER:: flow_number

    CHARACTER(LEN =64):: velocity_filename
    CHARACTER(LEN = 4):: flow_id
    !INTEGER:: ncid, vid

    write (flow_id,'(I4.4)') flow_number
    velocity_filename = trim("./magnetograms/velocity"//trim(flow_id)//'.nc')

!     call try(nf90_open(trim(velocity_filename), nf90_nowrite, ncid))
!
!     call try(nf90_inq_varid(ncid, 'vx', vid))
!     call try(nf90_get_var(ncid, vid, surf_vx(0:nx,0:ny), &
!     start = (/x_rank*nx+1,y_rank*ny+1/),count = (/nx+1,ny+1/)))
!
!     call try(nf90_inq_varid(ncid, 'vy', vid))
!     call try(nf90_get_var(ncid, vid, surf_vy(0:nx,0:ny), &
!     start = (/x_rank*nx+1,y_rank*ny+1/),count = (/nx+1,ny+1/)))
!
!     call try(nf90_inq_varid(ncid, 'vz', vid))
!     call try(nf90_get_var(ncid, vid, surf_vz(0:nx,0:ny), &
!     start = (/x_rank*nx+1,y_rank*ny+1/),count = (/nx+1,ny+1/)))
!
!     if (proc_num == 0) print*, 'Surface velocity loaded, fname: ', velocity_filename

    !call try(nf90_inq_varid(ncid, 'bx', vid))
    !call try(nf90_get_var(ncid, vid, bx(0:nx,1:ny,0), &
    !start = (/x_rank*nx+1,y_rank*ny+1/),count = (/nx+1,ny/)))

    !call try(nf90_inq_varid(ncid, 'by', vid))
    !call try(nf90_get_var(ncid, vid, by(1:nx,0:ny,0), &
    !start = (/x_rank*nx+1,y_rank*ny+1/),count = (/nx,ny+1/)))

    !call try(nf90_inq_varid(ncid, 'bz', vid))
    !call try(nf90_get_var(ncid, vid, bz(1:nx,1:ny,0), &
    !start = (/x_rank*nx+1,y_rank*ny+1/),count = (/nx,ny/)))

    !if (proc_num == 0) print*, 'Surface magnetic field loaded'

    !call try(nf90_close(ncid))




END SUBROUTINE update_surface_flows

subroutine try(status)
    ! Catch error in reading netcdf fild.
    INTEGER, INTENT(IN):: status

    if (status /= NF90_noerr) THEN
        PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
        call mpi_abort(comm, ierr)
    end if

end subroutine try





















!*******************************************************************************
END MODULE evolve
!*******************************************************************************
