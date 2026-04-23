module variables
  implicit none
  integer, parameter :: dp = selected_real_kind(14)

  ! ------------------------------------------------------------
  ! Rahman (1964) liquid argon state point in Lennard-Jones units
  ! Paper parameters:
  !   N = 864
  !   epsilon/kB = 120 K
  !   sigma = 3.4 A
  !   cutoff = 2.25 sigma
  !   density = 1.374 g/cm^3  -> rho* ~= 0.80718
  !   T = 94.4 K              -> T*   ~= 0.78667
  !   dt = 1e-14 s            -> dt*  ~= 0.00463
  ! ------------------------------------------------------------

  integer, parameter :: n = 864
  real(dp), parameter :: sigma = 1.0_dp
  real(dp), parameter :: epsilon = 1.0_dp
  real(dp), parameter :: mass = 1.0_dp
  real(dp), parameter :: rho = 0.80718_dp
  real(dp), parameter :: rc = 2.25_dp * sigma
  real(dp), parameter :: rc2 = rc*rc
  real(dp), parameter :: boxlength = (real(n,dp)/rho)**(1.0_dp/3.0_dp)
  real(dp), parameter :: halfbox = 0.5_dp*boxlength
  real(dp), parameter :: dt = 0.00463_dp
  real(dp), parameter :: tstar_target = 94.4_dp/120.0_dp
  integer, parameter :: nsteps = 5000
  integer, parameter :: thermo_every = 1

  ! Truncated-and-shifted LJ potential: U_shifted(rc)=0.
  real(dp), parameter :: ucut = 4.0_dp*((1.0_dp/rc)**12 - (1.0_dp/rc)**6)

  real(dp), dimension(n) :: x, y, z
  real(dp), dimension(n) :: vx, vy, vz
  real(dp), dimension(n) :: fx, fy, fz
  real(dp), dimension(n) :: fx_new, fy_new, fz_new

  real(dp) :: potential_energy, kinetic_energy, total_energy, temperature

contains

  subroutine init_random_seed(seed)
    integer, intent(in) :: seed
    integer :: nseed, i
    integer, allocatable :: put(:)

    call random_seed(size=nseed)
    allocate(put(nseed))
    do i = 1, nseed
      put(i) = seed + 97*(i-1)
    end do
    call random_seed(put=put)
    deallocate(put)
  end subroutine init_random_seed

  subroutine minimum_image(dx, dy, dz)
    real(dp), intent(inout) :: dx, dy, dz

    if (dx >  halfbox) dx = dx - boxlength
    if (dx < -halfbox) dx = dx + boxlength
    if (dy >  halfbox) dy = dy - boxlength
    if (dy < -halfbox) dy = dy + boxlength
    if (dz >  halfbox) dz = dz - boxlength
    if (dz < -halfbox) dz = dz + boxlength
  end subroutine minimum_image

  subroutine wrap_positions()
    integer :: i
    do i = 1, n
      x(i) = modulo(x(i), boxlength)
      y(i) = modulo(y(i), boxlength)
      z(i) = modulo(z(i), boxlength)
    end do
  end subroutine wrap_positions

  subroutine init_positions_random()
    ! Random initial positions with a mild hard-core exclusion.
    ! This follows the paper more closely than an FCC start, while
    ! avoiding catastrophic overlaps that blow up the LJ repulsion.
    integer :: i, j, tries
    real(dp) :: rx, ry, rz, dx, dy, dz, r2
    logical :: overlap
    real(dp), parameter :: rmin = 0.90_dp
    real(dp), parameter :: rmin2 = rmin*rmin

    call init_random_seed(1964)

    do i = 1, n
      tries = 0
      do
        tries = tries + 1
        if (tries > 200000) then
          print *, 'Failed to place particle ', i, ' without overlap.'
          stop 1
        end if

        call random_number(rx)
        call random_number(ry)
        call random_number(rz)

        rx = rx * boxlength
        ry = ry * boxlength
        rz = rz * boxlength

        overlap = .false.
        do j = 1, i-1
          dx = rx - x(j)
          dy = ry - y(j)
          dz = rz - z(j)
          call minimum_image(dx, dy, dz)
          r2 = dx*dx + dy*dy + dz*dz
          if (r2 < rmin2) then
            overlap = .true.
            exit
          end if
        end do

        if (.not. overlap) then
          x(i) = rx
          y(i) = ry
          z(i) = rz
          exit
        end if
      end do
    end do
  end subroutine init_positions_random

  subroutine init_velocities_target_temperature()
    ! Draw random velocities, remove COM drift, and rescale to target T*.
    integer :: i
    real(dp) :: rx, ry, rz, vcmx, vcmy, vcmz, sumv2, scale

    do i = 1, n
      call random_number(rx)
      call random_number(ry)
      call random_number(rz)
      vx(i) = rx - 0.5_dp
      vy(i) = ry - 0.5_dp
      vz(i) = rz - 0.5_dp
    end do

    vcmx = sum(vx)/real(n,dp)
    vcmy = sum(vy)/real(n,dp)
    vcmz = sum(vz)/real(n,dp)

    do i = 1, n
      vx(i) = vx(i) - vcmx
      vy(i) = vy(i) - vcmy
      vz(i) = vz(i) - vcmz
    end do

    sumv2 = sum(vx*vx + vy*vy + vz*vz)
    if (sumv2 <= 0.0_dp) then
      print *, 'Velocity initialization failed.'
      stop 2
    end if

    ! In reduced LJ units with m*=1:
    ! T* = sum(v_i^2) / (3N - 3)
    scale = sqrt(tstar_target*(3.0_dp*real(n,dp) - 3.0_dp)/sumv2)

    do i = 1, n
      vx(i) = scale * vx(i)
      vy(i) = scale * vy(i)
      vz(i) = scale * vz(i)
    end do
  end subroutine init_velocities_target_temperature

  subroutine compute_forces(xin, yin, zin, fxout, fyout, fzout, pe)
    real(dp), intent(in)  :: xin(n), yin(n), zin(n)
    real(dp), intent(out) :: fxout(n), fyout(n), fzout(n)
    real(dp), intent(out) :: pe

    integer :: i, j
    real(dp) :: dx, dy, dz, r2, invr2, invr6, invr12, fij_over_r

    fxout = 0.0_dp
    fyout = 0.0_dp
    fzout = 0.0_dp
    pe = 0.0_dp

    do i = 1, n-1
      do j = i+1, n
        dx = xin(i) - xin(j)
        dy = yin(i) - yin(j)
        dz = zin(i) - zin(j)
        call minimum_image(dx, dy, dz)
        r2 = dx*dx + dy*dy + dz*dz

        if (r2 < rc2) then
          invr2 = 1.0_dp/r2
          invr6 = invr2*invr2*invr2
          invr12 = invr6*invr6

          ! Force on i from j: F = 48 r^-14 - 24 r^-8, times vector r_ij
          fij_over_r = 48.0_dp*invr12*invr2 - 24.0_dp*invr6*invr2

          fxout(i) = fxout(i) + fij_over_r*dx
          fyout(i) = fyout(i) + fij_over_r*dy
          fzout(i) = fzout(i) + fij_over_r*dz

          fxout(j) = fxout(j) - fij_over_r*dx
          fyout(j) = fyout(j) - fij_over_r*dy
          fzout(j) = fzout(j) - fij_over_r*dz

          pe = pe + 4.0_dp*(invr12 - invr6) - ucut
        end if
      end do
    end do
  end subroutine compute_forces

  subroutine compute_kinetic_and_temperature()
    kinetic_energy = 0.5_dp*sum(vx*vx + vy*vy + vz*vz)
    temperature = (2.0_dp*kinetic_energy)/(3.0_dp*real(n,dp) - 3.0_dp)
    total_energy = kinetic_energy + potential_energy
  end subroutine compute_kinetic_and_temperature

  subroutine velocity_verlet_step()
    integer :: i

    do i = 1, n
      x(i) = x(i) + vx(i)*dt + 0.5_dp*fx(i)*dt*dt
      y(i) = y(i) + vy(i)*dt + 0.5_dp*fy(i)*dt*dt
      z(i) = z(i) + vz(i)*dt + 0.5_dp*fz(i)*dt*dt
    end do

    call wrap_positions()

    do i = 1, n
      vx(i) = vx(i) + 0.5_dp*fx(i)*dt
      vy(i) = vy(i) + 0.5_dp*fy(i)*dt
      vz(i) = vz(i) + 0.5_dp*fz(i)*dt
    end do

    call compute_forces(x, y, z, fx_new, fy_new, fz_new, potential_energy)

    do i = 1, n
      vx(i) = vx(i) + 0.5_dp*fx_new(i)*dt
      vy(i) = vy(i) + 0.5_dp*fy_new(i)*dt
      vz(i) = vz(i) + 0.5_dp*fz_new(i)*dt
    end do

    fx = fx_new
    fy = fy_new
    fz = fz_new

    call compute_kinetic_and_temperature()
  end subroutine velocity_verlet_step

end module variables

program rahman_argon_1964
  use variables
  implicit none

  integer :: step
  real(dp) :: time

  open(unit=10, file='rahman_energy_output.dat', status='replace', action='write')
  write(10,'(a)') '# step time total_energy potential_energy kinetic_energy temperature'

  call init_positions_random()
  call init_velocities_target_temperature()
  call compute_forces(x, y, z, fx, fy, fz, potential_energy)
  call compute_kinetic_and_temperature()

  print *, 'Rahman (1964) style LJ Argon run'
  print *, 'N              = ', n
  print *, 'rho*           = ', rho
  print *, 'boxlength      = ', boxlength
  print *, 'rc*            = ', rc
  print *, 'dt*            = ', dt
  print *, 'target T*      = ', tstar_target
  print *, 'initial PE      = ', potential_energy
  print *, 'initial KE      = ', kinetic_energy
  print *, 'initial total E = ', total_energy
  print *, 'initial T*      = ', temperature

  time = 0.0_dp
  write(10,'(i8,1x,f12.6,1x,f18.10,1x,f18.10,1x,f18.10,1x,f18.10)') &
       0, time, total_energy, potential_energy, kinetic_energy, temperature

  do step = 1, nsteps
    call velocity_verlet_step()
    time = real(step,dp)*dt

    if (mod(step, thermo_every) == 0) then
      write(10,'(i8,1x,f12.6,1x,f18.10,1x,f18.10,1x,f18.10,1x,f18.10)') &
           step, time, total_energy, potential_energy, kinetic_energy, temperature
    end if
  end do

  close(10)

  print *, 'Done. Wrote energy_output.dat'
end program rahman_argon_1964
