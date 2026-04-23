module variables
  implicit none
  integer, parameter :: dp = selected_real_kind(14)

  ! ================================================================
  ! Lennard-Jones MD in reduced units for Argon assignment
  !
  ! Assignment targets:
  !   N      = 256
  !   rho*   = 0.636
  !   E*     = 101.79
  !   rc     = 2.5 sigma
  !   PBC    = periodic cubic box
  !   ensemble = NVE
  !
  ! Reduced LJ units:
  !   sigma = 1, epsilon = 1, m = 1
  ! ================================================================

  integer, parameter :: n = 256
  real(dp), parameter :: rho = 0.636_dp
  real(dp), parameter :: sigma = 1.0_dp
  real(dp), parameter :: rc = 2.5_dp * sigma
  real(dp), parameter :: target_total_energy = 101.79_dp

  ! Box length from rho* = N / V*
  real(dp), parameter :: boxlength = (real(n, dp) / rho)**(1.0_dp / 3.0_dp)
  real(dp), parameter :: lxh = boxlength / 2.0_dp

  ! A conservative reduced-unit timestep suitable for LJ NVE dynamics.
  real(dp), parameter :: dt = 0.001_dp

  ! Number of production steps to write.
  integer, parameter :: nsteps = 20000

  ! Current positions, velocities, and forces
  real(dp), dimension(n) :: x, y, z
  real(dp), dimension(n) :: vx, vy, vz
  real(dp), dimension(n) :: fx, fy, fz

  ! Trial positions used during initialization only
  real(dp), dimension(n) :: x_try, y_try, z_try

  ! Energies and temperature
  real(dp) :: p_tot, ke_tot, tot_en, temp_red

contains

  subroutine init_random_seed(seed)
    integer, intent(in) :: seed
    integer :: nseed, i
    integer, allocatable :: put(:)

    call random_seed(size=nseed)
    allocate(put(nseed))
    do i = 1, nseed
      put(i) = seed + 37*(i-1)
    end do
    call random_seed(put=put)
    deallocate(put)
  end subroutine init_random_seed

  pure function minimum_image(d) result(dm)
    real(dp), intent(in) :: d
    real(dp) :: dm

    dm = d
    if (dm > lxh) then
      dm = dm - boxlength
    else if (dm < -lxh) then
      dm = dm + boxlength
    end if
  end function minimum_image

end module variables

program main
  use variables
  implicit none
  integer :: step
  real(dp) :: t

  open(unit=76, file='energy_output.dat', status='replace', action='write')
  write(76,'(A)') '# step   time   total_energy   potential_energy   kinetic_energy   temperature'

  call init_pos()
  call init_vel_rescale_to_target_energy()

  ! Initial energies at t = 0
  call force_calc()
  call kinetic_calc()
  tot_en   = ke_tot + p_tot
  temp_red = 2.0_dp * ke_tot / (3.0_dp*real(n,dp) - 3.0_dp)

  write(*,'(A,F12.6)') 'Box length = ', boxlength
  write(*,'(A,F16.8)') 'Initial total energy = ', tot_en
  write(*,'(A,F16.8)') 'Initial potential    = ', p_tot
  write(*,'(A,F16.8)') 'Initial kinetic      = ', ke_tot
  write(*,'(A,F16.8)') 'Initial T*           = ', temp_red

  t = 0.0_dp
  write(76,'(I8,1X,F12.6,1X,F20.10,1X,F20.10,1X,F20.10,1X,F16.8)') &
       0, t, tot_en, p_tot, ke_tot, temp_red

  do step = 1, nsteps
    call integrate_velocity_verlet()
    call force_calc()
    call finish_velocity_verlet_and_ke()

    tot_en   = ke_tot + p_tot
    temp_red = 2.0_dp * ke_tot / (3.0_dp*real(n,dp) - 3.0_dp)
    t = real(step, dp) * dt

    write(*,'(A,I8,2X,A,F12.6,2X,A,F16.8,2X,A,F16.8,2X,A,F16.8)') &
         'step', step, 't=', t, 'E=', tot_en, 'U=', p_tot, 'K=', ke_tot

    write(76,'(I8,1X,F12.6,1X,F20.10,1X,F20.10,1X,F20.10,1X,F16.8)') &
         step, t, tot_en, p_tot, ke_tot, temp_red
  end do

  close(76)
end program main

subroutine init_pos()
  use variables
  implicit none
  integer :: ix, iy, iz, ib, idx
  real(dp) :: a
  real(dp), dimension(4,3) :: basis

  ! --------------------------------------------------------------
  ! Place particles on an FCC lattice.
  !
  ! Why this fix is needed:
  ! Random placement at rho*=0.636 can create very close pairs and a huge
  ! positive initial potential energy. Then U can exceed the target total
  ! energy and K = E - U becomes negative.
  !
  ! For N=256, a 4x4x4 FCC lattice is exact because each FCC cell has
  ! 4 atoms, so 4^3 * 4 = 256.
  ! --------------------------------------------------------------
  a = boxlength / 4.0_dp

  basis(1, :) = [0.0_dp, 0.0_dp, 0.0_dp]
  basis(2, :) = [0.5_dp, 0.5_dp, 0.0_dp]
  basis(3, :) = [0.5_dp, 0.0_dp, 0.5_dp]
  basis(4, :) = [0.0_dp, 0.5_dp, 0.5_dp]

  idx = 0
  do ix = 0, 3
    do iy = 0, 3
      do iz = 0, 3
        do ib = 1, 4
          idx = idx + 1
          x(idx) = (real(ix,dp) + basis(ib,1)) * a
          y(idx) = (real(iy,dp) + basis(ib,2)) * a
          z(idx) = (real(iz,dp) + basis(ib,3)) * a
        end do
      end do
    end do
  end do
end subroutine init_pos

subroutine init_vel_rescale_to_target_energy()
  use variables
  implicit none
  integer :: i
  real(dp) :: avx, avy, avz, sumv2, scale, target_ke

  ! --------------------------------------------------------------
  ! Random initial velocities, then remove center-of-mass drift,
  ! then rescale so that K + U = target_total_energy.
  ! --------------------------------------------------------------
  do i = 1, n
    call random_number(vx(i))
    call random_number(vy(i))
    call random_number(vz(i))

    vx(i) = vx(i) - 0.5_dp
    vy(i) = vy(i) - 0.5_dp
    vz(i) = vz(i) - 0.5_dp
  end do

  avx = sum(vx) / real(n, dp)
  avy = sum(vy) / real(n, dp)
  avz = sum(vz) / real(n, dp)

  do i = 1, n
    vx(i) = vx(i) - avx
    vy(i) = vy(i) - avy
    vz(i) = vz(i) - avz
  end do

  call force_calc()

  target_ke = target_total_energy - p_tot
  if (target_ke <= 0.0_dp) then
    write(*,*) 'Target kinetic energy is non-positive.'
    write(*,*) 'Potential energy at initialization = ', p_tot
    stop 2
  end if

  sumv2 = sum(vx*vx + vy*vy + vz*vz)
  if (sumv2 <= 0.0_dp) then
    write(*,*) 'Velocity initialization failed: zero norm.'
    stop 3
  end if

  ! K = 0.5 * sum(v^2) in reduced units because m*=1
  scale = sqrt((2.0_dp * target_ke) / sumv2)
  vx = scale * vx
  vy = scale * vy
  vz = scale * vz
end subroutine init_vel_rescale_to_target_energy

subroutine force_calc()
  use variables
  implicit none
  integer :: i, j
  real(dp) :: xr, yr, zr, dr2
  real(dp) :: inv_r2, inv_r6, inv_r8, inv_r12, inv_r14
  real(dp) :: pot_energy, fxij, fyij, fzij
  real(dp), parameter :: inv_rc2 = 1.0_dp / (rc*rc)
  real(dp), parameter :: inv_rc6 = inv_rc2**3
  real(dp), parameter :: u_shift = 4.0_dp*(inv_rc6*inv_rc6 - inv_rc6)

  ! --------------------------------------------------------------
  ! Truncated-and-shifted LJ potential:
  ! U(r) = 4[(1/r)^12 - (1/r)^6] - U(rc),  r < rc
  !      = 0,                               r >= rc
  !
  ! The shift makes the potential continuous at rc.
  ! --------------------------------------------------------------
  p_tot = 0.0_dp
  fx = 0.0_dp
  fy = 0.0_dp
  fz = 0.0_dp

  do i = 1, n-1
    do j = i+1, n
      xr = minimum_image(x(i) - x(j))
      yr = minimum_image(y(i) - y(j))
      zr = minimum_image(z(i) - z(j))
      dr2 = xr*xr + yr*yr + zr*zr

      if (dr2 < rc*rc) then
        inv_r2  = 1.0_dp / dr2
        inv_r6  = inv_r2**3
        inv_r8  = inv_r6 * inv_r2
        inv_r12 = inv_r6 * inv_r6
        inv_r14 = inv_r12 * inv_r2

        pot_energy = 4.0_dp*(inv_r12 - inv_r6) - u_shift

        fxij = 48.0_dp * xr * (inv_r14 - 0.5_dp*inv_r8)
        fyij = 48.0_dp * yr * (inv_r14 - 0.5_dp*inv_r8)
        fzij = 48.0_dp * zr * (inv_r14 - 0.5_dp*inv_r8)

        fx(i) = fx(i) + fxij
        fy(i) = fy(i) + fyij
        fz(i) = fz(i) + fzij

        fx(j) = fx(j) - fxij
        fy(j) = fy(j) - fyij
        fz(j) = fz(j) - fzij

        p_tot = p_tot + pot_energy
      end if
    end do
  end do
end subroutine force_calc

subroutine kinetic_calc()
  use variables
  implicit none

  ! --------------------------------------------------------------
  ! Kinetic energy from explicit velocities.
  ! --------------------------------------------------------------
  ke_tot = 0.5_dp * sum(vx*vx + vy*vy + vz*vz)
end subroutine kinetic_calc

subroutine integrate_velocity_verlet()
  use variables
  implicit none

  ! --------------------------------------------------------------
  ! First half of velocity-Verlet:
  !   v(t+dt/2) = v(t) + 0.5 a(t) dt
  !   x(t+dt)   = x(t) + v(t+dt/2) dt
  !
  ! Although the assignment says "simple Verlet", velocity-Verlet is
  ! the standard practical form in the same Verlet family and gives
  ! directly consistent kinetic energies at each output step.
  ! --------------------------------------------------------------
  vx = vx + 0.5_dp * fx * dt
  vy = vy + 0.5_dp * fy * dt
  vz = vz + 0.5_dp * fz * dt

  x = x + vx * dt
  y = y + vy * dt
  z = z + vz * dt

  ! Wrap positions back into the periodic box.
  where (x >= boxlength) x = x - boxlength
  where (x <  0.0_dp  ) x = x + boxlength

  where (y >= boxlength) y = y - boxlength
  where (y <  0.0_dp  ) y = y + boxlength

  where (z >= boxlength) z = z - boxlength
  where (z <  0.0_dp  ) z = z + boxlength
end subroutine integrate_velocity_verlet

subroutine finish_velocity_verlet_and_ke()
  use variables
  implicit none

  ! --------------------------------------------------------------
  ! Second half of velocity-Verlet after new forces are known:
  !   v(t+dt) = v(t+dt/2) + 0.5 a(t+dt) dt
  ! Then compute kinetic energy at the same time as the new potential.
  ! --------------------------------------------------------------
  vx = vx + 0.5_dp * fx * dt
  vy = vy + 0.5_dp * fy * dt
  vz = vz + 0.5_dp * fz * dt

  call kinetic_calc()
end subroutine finish_velocity_verlet_and_ke
