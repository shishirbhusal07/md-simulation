
module md_params
  implicit none
  integer, parameter :: dp = selected_real_kind(14, 200)

  ! Rahman (1964)-style Lennard-Jones Argon state point in reduced units
  integer, parameter :: n = 864
  real(dp), parameter :: rho = 0.80718_dp
  real(dp), parameter :: sigma = 1.0_dp
  real(dp), parameter :: rc = 2.25_dp
  real(dp), parameter :: boxlength = (real(n,dp)/rho)**(1.0_dp/3.0_dp)
  real(dp), parameter :: halfbox = 0.5_dp*boxlength
  real(dp), parameter :: dt = 0.00463_dp
  real(dp), parameter :: target_t = 94.4_dp/120.0_dp

  ! Run control
  integer, parameter :: n_equil = 5000
  integer, parameter :: n_prod  = 15000
  integer, parameter :: sample_every = 10
  integer, parameter :: trajectory_every = 20

  ! g(r) histogram
  integer, parameter :: nbins_gr = 300
  real(dp), parameter :: dr_gr = halfbox / real(nbins_gr, dp)

  ! Soft initialization
  real(dp), parameter :: min_sep_init = 0.90_dp

  ! Shifted LJ potential at rc
  real(dp), parameter :: u_shift = 4.0_dp*((1.0_dp/rc)**12 - (1.0_dp/rc)**6)

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

  pure subroutine minimum_image(dx, dy, dz)
    real(dp), intent(inout) :: dx, dy, dz
    if (dx >  halfbox) dx = dx - boxlength
    if (dx < -halfbox) dx = dx + boxlength
    if (dy >  halfbox) dy = dy - boxlength
    if (dy < -halfbox) dy = dy + boxlength
    if (dz >  halfbox) dz = dz - boxlength
    if (dz < -halfbox) dz = dz + boxlength
  end subroutine minimum_image

  pure subroutine wrap_position(x, y, z)
    real(dp), intent(inout) :: x, y, z
    x = x - boxlength*floor(x/boxlength)
    y = y - boxlength*floor(y/boxlength)
    z = z - boxlength*floor(z/boxlength)
  end subroutine wrap_position

  subroutine init_positions_random(x, y, z)
    real(dp), intent(out) :: x(n), y(n), z(n)
    integer :: i, j, tries
    real(dp) :: rx, ry, rz, dx, dy, dz, rij

    call init_random_seed(1964)

    do i = 1, n
       tries = 0
place_particle: do
          tries = tries + 1
          if (tries > 100000) then
             write(*,*) 'Failed to place particle ', i
             stop 1
          end if

          call random_number(rx); call random_number(ry); call random_number(rz)
          x(i) = rx*boxlength
          y(i) = ry*boxlength
          z(i) = rz*boxlength

          do j = 1, i-1
             dx = x(i) - x(j)
             dy = y(i) - y(j)
             dz = z(i) - z(j)
             call minimum_image(dx, dy, dz)
             rij = sqrt(dx*dx + dy*dy + dz*dz)
             if (rij < min_sep_init) cycle place_particle
          end do
          exit place_particle
       end do place_particle
    end do
  end subroutine init_positions_random

  subroutine init_velocities_temperature(vx, vy, vz)
    real(dp), intent(out) :: vx(n), vy(n), vz(n)
    integer :: i
    real(dp) :: r1, r2, ga1, ga2
    real(dp) :: avx, avy, avz, kin, scale

    ! Gaussian velocities via Box-Muller
    i = 1
    do while (i <= n)
       call random_number(r1); call random_number(r2)
       r1 = max(r1, 1.0e-12_dp)
       ga1 = sqrt(-2.0_dp*log(r1))*cos(2.0_dp*acos(-1.0_dp)*r2)
       ga2 = sqrt(-2.0_dp*log(r1))*sin(2.0_dp*acos(-1.0_dp)*r2)
       vx(i) = ga1
       if (i < n) vx(i+1) = ga2

       call random_number(r1); call random_number(r2)
       r1 = max(r1, 1.0e-12_dp)
       ga1 = sqrt(-2.0_dp*log(r1))*cos(2.0_dp*acos(-1.0_dp)*r2)
       ga2 = sqrt(-2.0_dp*log(r1))*sin(2.0_dp*acos(-1.0_dp)*r2)
       vy(i) = ga1
       if (i < n) vy(i+1) = ga2

       call random_number(r1); call random_number(r2)
       r1 = max(r1, 1.0e-12_dp)
       ga1 = sqrt(-2.0_dp*log(r1))*cos(2.0_dp*acos(-1.0_dp)*r2)
       ga2 = sqrt(-2.0_dp*log(r1))*sin(2.0_dp*acos(-1.0_dp)*r2)
       vz(i) = ga1
       if (i < n) vz(i+1) = ga2

       i = i + 2
    end do

    ! Remove center-of-mass drift
    avx = sum(vx)/real(n,dp)
    avy = sum(vy)/real(n,dp)
    avz = sum(vz)/real(n,dp)
    vx = vx - avx
    vy = vy - avy
    vz = vz - avz

    ! Rescale to target reduced temperature
    kin = 0.5_dp*sum(vx*vx + vy*vy + vz*vz)
    scale = sqrt((1.5_dp*(real(n,dp)-1.0_dp)*target_t)/kin)
    vx = scale*vx
    vy = scale*vy
    vz = scale*vz
  end subroutine init_velocities_temperature

  subroutine compute_forces(x, y, z, fx, fy, fz, potential)
    real(dp), intent(in)  :: x(n), y(n), z(n)
    real(dp), intent(out) :: fx(n), fy(n), fz(n)
    real(dp), intent(out) :: potential
    integer :: i, j
    real(dp) :: dx, dy, dz, r2, r, invr2, invr6, invr12, fij_over_r

    fx = 0.0_dp; fy = 0.0_dp; fz = 0.0_dp
    potential = 0.0_dp

    do i = 1, n-1
       do j = i+1, n
          dx = x(i)-x(j)
          dy = y(i)-y(j)
          dz = z(i)-z(j)
          call minimum_image(dx,dy,dz)
          r2 = dx*dx + dy*dy + dz*dz
          if (r2 < rc*rc) then
             r = sqrt(r2)
             invr2 = 1.0_dp/r2
             invr6 = invr2*invr2*invr2
             invr12 = invr6*invr6

             potential = potential + 4.0_dp*(invr12 - invr6) - u_shift

             ! Force = 48 r^-14 - 24 r^-8, written as (fij_over_r)*r_vec
             fij_over_r = 48.0_dp*invr2*(invr12 - 0.5_dp*invr6)
             fx(i) = fx(i) + fij_over_r*dx
             fy(i) = fy(i) + fij_over_r*dy
             fz(i) = fz(i) + fij_over_r*dz
             fx(j) = fx(j) - fij_over_r*dx
             fy(j) = fy(j) - fij_over_r*dy
             fz(j) = fz(j) - fij_over_r*dz
          end if
       end do
    end do
  end subroutine compute_forces

  subroutine velocity_verlet_step(x, y, z, vx, vy, vz, fx, fy, fz, potential)
    real(dp), intent(inout) :: x(n), y(n), z(n), vx(n), vy(n), vz(n)
    real(dp), intent(inout) :: fx(n), fy(n), fz(n)
    real(dp), intent(out)   :: potential
    real(dp) :: fx_new(n), fy_new(n), fz_new(n)
    integer :: i

    do i = 1, n
       x(i) = x(i) + vx(i)*dt + 0.5_dp*fx(i)*dt*dt
       y(i) = y(i) + vy(i)*dt + 0.5_dp*fy(i)*dt*dt
       z(i) = z(i) + vz(i)*dt + 0.5_dp*fz(i)*dt*dt
       call wrap_position(x(i), y(i), z(i))
       vx(i) = vx(i) + 0.5_dp*fx(i)*dt
       vy(i) = vy(i) + 0.5_dp*fy(i)*dt
       vz(i) = vz(i) + 0.5_dp*fz(i)*dt
    end do

    call compute_forces(x, y, z, fx_new, fy_new, fz_new, potential)

    do i = 1, n
       vx(i) = vx(i) + 0.5_dp*fx_new(i)*dt
       vy(i) = vy(i) + 0.5_dp*fy_new(i)*dt
       vz(i) = vz(i) + 0.5_dp*fz_new(i)*dt
    end do

    fx = fx_new; fy = fy_new; fz = fz_new
  end subroutine velocity_verlet_step

  pure function kinetic_energy(vx, vy, vz) result(ke)
    real(dp), intent(in) :: vx(n), vy(n), vz(n)
    real(dp) :: ke
    ke = 0.5_dp*sum(vx*vx + vy*vy + vz*vz)
  end function kinetic_energy

  pure function temperature_from_ke(ke) result(temp)
    real(dp), intent(in) :: ke
    real(dp) :: temp
    temp = 2.0_dp*ke/(3.0_dp*real(n-1,dp))
  end function temperature_from_ke

  subroutine rescale_velocities(vx, vy, vz, t_target)
    real(dp), intent(inout) :: vx(n), vy(n), vz(n)
    real(dp), intent(in)    :: t_target
    real(dp) :: ke, tcur, scale
    ke = kinetic_energy(vx, vy, vz)
    tcur = temperature_from_ke(ke)
    if (tcur <= 0.0_dp) return
    scale = sqrt(t_target/tcur)
    vx = scale*vx
    vy = scale*vy
    vz = scale*vz
  end subroutine rescale_velocities

  subroutine accumulate_gr(x, y, z, gr_hist)
    real(dp), intent(in)    :: x(n), y(n), z(n)
    real(dp), intent(inout) :: gr_hist(nbins_gr)
    integer :: i, j, bin
    real(dp) :: dx, dy, dz, rij

    do i = 1, n-1
       do j = i+1, n
          dx = x(i)-x(j)
          dy = y(i)-y(j)
          dz = z(i)-z(j)
          call minimum_image(dx,dy,dz)
          rij = sqrt(dx*dx + dy*dy + dz*dz)
          if (rij < halfbox) then
             bin = int(rij/dr_gr) + 1
             if (bin >= 1 .and. bin <= nbins_gr) then
                ! Count both i->j and j->i so normalization is simple.
                gr_hist(bin) = gr_hist(bin) + 2.0_dp
             end if
          end if
       end do
    end do
  end subroutine accumulate_gr

  subroutine write_gr(gr_hist, n_samples)
    real(dp), intent(in) :: gr_hist(nbins_gr)
    integer, intent(in) :: n_samples
    integer :: b
    real(dp) :: r_lower, r_upper, r_mid, shell_vol, ideal_count, gval
    real(dp) :: volume, number_density
    integer :: unitno

    volume = boxlength**3
    number_density = real(n,dp)/volume

    open(newunit=unitno, file='gr_output.dat', status='replace', action='write')
    write(unitno,'(a)') '# r_mid  g(r)  raw_counts'
    do b = 1, nbins_gr
       r_lower = (b-1)*dr_gr
       r_upper = b*dr_gr
       r_mid   = 0.5_dp*(r_lower + r_upper)
       shell_vol = (4.0_dp/3.0_dp)*acos(-1.0_dp)*(r_upper**3 - r_lower**3)

       ! Expected ideal-gas count per sample in this shell:
       ! N * rho * shell_vol, because histogram counts each pair twice.
       ideal_count = real(n,dp) * number_density * shell_vol
       if (n_samples > 0 .and. ideal_count > 0.0_dp) then
          gval = gr_hist(b) / (real(n_samples,dp) * ideal_count)
       else
          gval = 0.0_dp
       end if
       write(unitno,'(3(1x,es20.10))') r_mid, gval, gr_hist(b)
    end do
    close(unitno)
  end subroutine write_gr


  subroutine write_xyz_frame(unitno, step, time, x, y, z)
    integer, intent(in) :: unitno, step
    real(dp), intent(in) :: time
    real(dp), intent(in) :: x(n), y(n), z(n)
    integer :: i

    write(unitno,'(i0)') n
    write(unitno,'(a,i0,a,es20.10,a,f18.8)') 'step=', step, ' time=', time, ' boxlength=', boxlength
    do i = 1, n
       write(unitno,'(a,3(1x,f18.10))') 'Ar', x(i), y(i), z(i)
    end do
  end subroutine write_xyz_frame

end module md_params


program rahman_lj_argon_with_gr
  use md_params
  implicit none

  real(dp) :: x(n), y(n), z(n), vx(n), vy(n), vz(n), fx(n), fy(n), fz(n)
  real(dp) :: potential, ke, te, temp, time
  real(dp) :: gr_hist(nbins_gr)
  integer :: step, unit_energy, unit_xyz, n_gr_samples

  call init_positions_random(x, y, z)
  call init_velocities_temperature(vx, vy, vz)
  call compute_forces(x, y, z, fx, fy, fz, potential)

  ke = kinetic_energy(vx, vy, vz)
  te = ke + potential
  temp = temperature_from_ke(ke)

  write(*,'(a)') 'Rahman (1964) style LJ Argon run with g(r)'
  write(*,'(a,i0)') ' N              = ', n
  write(*,'(a,f18.8)') ' rho*           = ', rho
  write(*,'(a,f18.8)') ' boxlength      = ', boxlength
  write(*,'(a,f18.8)') ' rc*            = ', rc
  write(*,'(a,es18.8)') ' dt*            = ', dt
  write(*,'(a,f18.8)') ' target T*      = ', target_t
  write(*,'(a,f18.8)') ' initial PE     = ', potential
  write(*,'(a,f18.8)') ' initial KE     = ', ke
  write(*,'(a,f18.8)') ' initial total E= ', te
  write(*,'(a,f18.8)') ' initial T*     = ', temp
  write(*,'(a,i0)') ' equil steps    = ', n_equil
  write(*,'(a,i0)') ' prod steps     = ', n_prod
  write(*,'(a,i0)') ' traj every     = ', trajectory_every

  open(newunit=unit_energy, file='energy_output.dat', status='replace', action='write')
  write(unit_energy,'(a)') '# step  time  total_energy  potential_energy  kinetic_energy  temperature'
  open(newunit=unit_xyz, file='trajectory.xyz', status='replace', action='write')

  ! -------------------------
  ! Equilibration:
  ! keep rescaling to target T* so the dense liquid settles
  ! before the production NVE run.
  ! -------------------------
  time = 0.0_dp
  do step = 1, n_equil
     call velocity_verlet_step(x, y, z, vx, vy, vz, fx, fy, fz, potential)
     call rescale_velocities(vx, vy, vz, target_t)
     ke = kinetic_energy(vx, vy, vz)
     te = ke + potential
     temp = temperature_from_ke(ke)
     time = time + dt
  end do

  ! Recompute forces after the final rescale so the production run
  ! starts from a clean state.
  call compute_forces(x, y, z, fx, fy, fz, potential)

  gr_hist = 0.0_dp
  n_gr_samples = 0

  ! -------------------------
  ! Production:
  ! pure NVE data collection
  ! -------------------------
  do step = 0, n_prod
     ke = kinetic_energy(vx, vy, vz)
     te = ke + potential
     temp = temperature_from_ke(ke)
     write(unit_energy,'(i10,1x,5(es22.12,1x))') step, time, te, potential, ke, temp

     if (mod(step, sample_every) == 0) then
        call accumulate_gr(x, y, z, gr_hist)
        n_gr_samples = n_gr_samples + 1
     end if

     if (mod(step, trajectory_every) == 0) then
        call write_xyz_frame(unit_xyz, step, time, x, y, z)
     end if

     if (step < n_prod) then
        call velocity_verlet_step(x, y, z, vx, vy, vz, fx, fy, fz, potential)
        time = time + dt
     end if
  end do

  close(unit_energy)
  close(unit_xyz)
  call write_gr(gr_hist, n_gr_samples)

  write(*,'(a,i0)') ' g(r) samples accumulated = ', n_gr_samples
  write(*,'(a)') ' Wrote energy_output.dat, gr_output.dat, and trajectory.xyz'

end program rahman_lj_argon_with_gr
