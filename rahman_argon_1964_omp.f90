! ============================================================
! rahman_argon_1964_omp.f90
!
! Molecular dynamics simulation of liquid Argon using a
! Lennard-Jones potential, reproducing the state point from
! Rahman (1964) Phys. Rev. 136, A405.
!
! Features added over the baseline serial code:
!   - Verlet neighbor list: reduces force evaluation from O(N^2)
!     per step to O(N * avg_neighbors), rebuilt only when needed
!   - OpenMP parallelism: force loop, list build, g(r), and
!     integrator half-kicks are all parallelized across threads
!
! All quantities are in Lennard-Jones reduced units:
!   length in sigma, energy in epsilon, time in tau = sigma*sqrt(m/epsilon)
! ============================================================

module md_params

  ! omp_lib provides OpenMP runtime functions such as
  ! omp_get_wtime() (wall-clock timer) and omp_get_max_threads()
  use omp_lib
  implicit none

  ! dp: double-precision kind parameter.
  ! selected_real_kind(14, 200) requests at least 14 decimal digits
  ! of precision and a decimal exponent range of at least 200.
  integer, parameter :: dp = selected_real_kind(14, 200)

  ! ----------------------------------------------------------
  ! Physical system parameters (Rahman 1964 state point)
  ! ----------------------------------------------------------

  ! Number of atoms — 864 = 6^3 * 4 unit cells of an FCC lattice,
  ! a standard size used in the original Rahman paper.
  integer,  parameter :: n         = 864

  ! Reduced number density rho* = N * sigma^3 / V, matching liquid
  ! Argon near its triple point in LJ reduced units.
  real(dp), parameter :: rho       = 0.80718_dp

  ! LJ length scale sigma in reduced units (always 1 by definition).
  real(dp), parameter :: sigma     = 1.0_dp

  ! Cutoff radius rc: pair interactions are truncated and shifted
  ! to zero beyond this distance. rc = 2.25 sigma is the value
  ! used by Rahman (1964).
  real(dp), parameter :: rc        = 2.25_dp

  ! Simulation box side length, derived from N and rho* via V = N/rho.
  ! Cubic root gives the side length of the periodic cubic box.
  real(dp), parameter :: boxlength = (real(n,dp)/rho)**(1.0_dp/3.0_dp)

  ! Half the box length — used in the minimum image convention to
  ! decide which periodic image of j is closest to i.
  real(dp), parameter :: halfbox   = 0.5_dp*boxlength

  ! Integration time step in reduced units. dt = 0.00463 tau
  ! corresponds to ~10 fs for real Argon.
  real(dp), parameter :: dt        = 0.00463_dp

  ! Target reduced temperature T* = kT/epsilon = 94.4 K / 120 K.
  ! 120 K is the LJ epsilon/k_B for Argon; 94.4 K is the Rahman state point.
  real(dp), parameter :: target_t  = 94.4_dp/120.0_dp

  ! ----------------------------------------------------------
  ! Run control parameters
  ! ----------------------------------------------------------

  ! Number of equilibration steps: system is driven toward target_t
  ! by velocity rescaling at every step before production begins.
  integer, parameter :: n_equil          = 5000

  ! Number of production steps: pure NVE (constant N, V, E) dynamics
  ! from which thermodynamic averages are collected.
  integer, parameter :: n_prod           = 15000

  ! g(r) is accumulated every sample_every production steps.
  integer, parameter :: sample_every     = 10

  ! A full XYZ snapshot of all atom positions is written every
  ! trajectory_every steps to trajectory.xyz for visualization.
  integer, parameter :: trajectory_every = 20

  ! ----------------------------------------------------------
  ! Radial distribution function g(r) histogram parameters
  ! ----------------------------------------------------------

  ! Number of equally-spaced radial bins covering [0, halfbox].
  integer,  parameter :: nbins_gr = 300

  ! Width of each g(r) bin in reduced units.
  real(dp), parameter :: dr_gr    = halfbox / real(nbins_gr, dp)

  ! ----------------------------------------------------------
  ! Initialization parameters
  ! ----------------------------------------------------------

  ! Minimum allowed separation between any two atoms during random
  ! placement. Prevents overlapping cores that would cause divergent
  ! forces at the start of the simulation.
  real(dp), parameter :: min_sep_init = 0.90_dp

  ! ----------------------------------------------------------
  ! Shifted Lennard-Jones potential value at the cutoff rc.
  ! Subtracting u_shift from every pair interaction ensures U(rc) = 0,
  ! removing the discontinuity in energy at the cutoff.
  ! u_shift = 4 * [ (1/rc)^12 - (1/rc)^6 ]
  ! ----------------------------------------------------------
  real(dp), parameter :: u_shift = 4.0_dp*((1.0_dp/rc)**12 - (1.0_dp/rc)**6)

  ! ----------------------------------------------------------
  ! Verlet neighbor list parameters
  ! ----------------------------------------------------------

  ! Skin thickness added to rc when building the neighbor list.
  ! Neighbors within rc + r_skin are stored; the list is valid until
  ! any atom drifts more than r_skin/2 from its position at last build.
  ! Larger r_skin means fewer rebuilds but more neighbors per atom.
  real(dp), parameter :: r_skin    = 0.30_dp

  ! List cutoff radius: all atoms within r_list of atom i are stored
  ! as potential neighbors. r_list > rc ensures no true neighbor is missed
  ! between consecutive list rebuilds.
  real(dp), parameter :: r_list    = rc + r_skin

  ! Squared list cutoff — used in distance checks to avoid sqrt().
  real(dp), parameter :: r_list2   = r_list * r_list

  ! Squared force cutoff — same trick for the inner force loop.
  real(dp), parameter :: rc2       = rc * rc

  ! Maximum number of neighbors stored per atom. For liquid Argon at
  ! rho*=0.807 and r_list=2.55, ~100-120 neighbors are typical; 350
  ! provides a safe upper bound with room for density fluctuations.
  integer,  parameter :: max_neigh = 350

  ! ----------------------------------------------------------
  ! Module-level Verlet list storage arrays
  ! ----------------------------------------------------------

  ! nlist(i): number of neighbors currently stored for atom i.
  integer  :: nlist(n)

  ! vlist(k,i): index of the k-th neighbor of atom i.
  ! Stored as (max_neigh, n) — column-major Fortran layout — so that
  ! all neighbors of a single atom i (vlist(1:nlist(i), i)) are
  ! contiguous in memory, giving cache-friendly inner-loop access.
  integer  :: vlist(max_neigh, n)

  ! Reference positions at the time of the last list build.
  ! Used by need_rebuild() to measure how far each atom has drifted.
  real(dp) :: x0_vl(n), y0_vl(n), z0_vl(n)

contains

  ! ============================================================
  ! init_random_seed
  ! Seeds Fortran's intrinsic PRNG deterministically using the
  ! supplied integer seed. Each element of the seed array is
  ! set to seed + 37*(i-1), spreading the seed across all
  ! required seed-array slots to avoid correlated sequences.
  ! ============================================================
  subroutine init_random_seed(seed)
    integer, intent(in) :: seed       ! user-supplied base seed value
    integer :: nseed, i
    integer, allocatable :: put(:)

    ! Query how many integers the compiler's PRNG requires as seed.
    call random_seed(size=nseed)

    ! Allocate an array of that size to hold the seed values.
    allocate(put(nseed))

    ! Fill with a deterministic sequence derived from the base seed.
    do i = 1, nseed
       put(i) = seed + 37*(i-1)
    end do

    ! Install the constructed seed array into the PRNG state.
    call random_seed(put=put)

    deallocate(put)
  end subroutine init_random_seed

  ! ============================================================
  ! minimum_image
  ! Applies the minimum image convention to a displacement vector
  ! (dx, dy, dz) under periodic boundary conditions.
  ! If a component exceeds +halfbox it is folded back by one box
  ! length; if it is less than -halfbox it is extended by one box
  ! length. This gives the shortest-image distance between two atoms.
  ! Declared 'pure' so the compiler can call it freely from within
  ! OpenMP parallel regions without side-effect concerns.
  ! ============================================================
  pure subroutine minimum_image(dx, dy, dz)
    real(dp), intent(inout) :: dx, dy, dz

    ! Fold x-component into [-halfbox, +halfbox]
    if (dx >  halfbox) dx = dx - boxlength
    if (dx < -halfbox) dx = dx + boxlength

    ! Fold y-component
    if (dy >  halfbox) dy = dy - boxlength
    if (dy < -halfbox) dy = dy + boxlength

    ! Fold z-component
    if (dz >  halfbox) dz = dz - boxlength
    if (dz < -halfbox) dz = dz + boxlength
  end subroutine minimum_image

  ! ============================================================
  ! wrap_position
  ! Maps a single atom's coordinates back into the primary periodic
  ! box [0, boxlength) after the integrator has moved it outside.
  ! Uses floor() to handle displacements of more than one box length,
  ! which can occur during initialization or large-step moves.
  ! ============================================================
  pure subroutine wrap_position(x, y, z)
    real(dp), intent(inout) :: x, y, z

    ! Subtract the appropriate integer multiple of boxlength
    ! to bring x back into [0, boxlength).
    x = x - boxlength*floor(x/boxlength)
    y = y - boxlength*floor(y/boxlength)
    z = z - boxlength*floor(z/boxlength)
  end subroutine wrap_position

  ! ============================================================
  ! init_positions_random
  ! Places all N atoms at random positions inside the simulation box
  ! using a rejection algorithm: each candidate position is checked
  ! against all previously placed atoms and rejected if any pair
  ! distance is less than min_sep_init (0.9 sigma). This avoids
  ! hard-core overlaps that would produce huge initial forces.
  ! The PRNG is seeded with 1964 (year of Rahman's paper) for
  ! reproducibility.
  ! ============================================================
  subroutine init_positions_random(x, y, z)
    real(dp), intent(out) :: x(n), y(n), z(n)  ! output: atom positions
    integer :: i, j, tries
    real(dp) :: rx, ry, rz, dx, dy, dz, rij

    ! Seed PRNG for reproducible initial configuration.
    call init_random_seed(1964)

    ! Loop over each atom to be placed.
    do i = 1, n
       tries = 0

       ! Named loop so we can 'cycle' back to it on rejection.
place_particle: do
          tries = tries + 1

          ! Guard against infinite loops in pathologically dense systems.
          if (tries > 100000) then
             write(*,*) 'Failed to place particle ', i
             stop 1
          end if

          ! Draw three uniform random numbers in [0,1) and scale to box.
          call random_number(rx); call random_number(ry); call random_number(rz)
          x(i) = rx*boxlength
          y(i) = ry*boxlength
          z(i) = rz*boxlength

          ! Check against all already-placed atoms j < i.
          do j = 1, i-1
             dx = x(i) - x(j)
             dy = y(i) - y(j)
             dz = z(i) - z(j)

             ! Apply minimum image so we measure the nearest-image distance.
             call minimum_image(dx, dy, dz)
             rij = sqrt(dx*dx + dy*dy + dz*dz)

             ! If too close to atom j, reject and try a new position.
             if (rij < min_sep_init) cycle place_particle
          end do

          ! Passed all overlap checks — accept this position.
          exit place_particle
       end do place_particle
    end do
  end subroutine init_positions_random

  ! ============================================================
  ! init_velocities_temperature
  ! Assigns initial velocities drawn from a Maxwell-Boltzmann
  ! distribution at target_t using the Box-Muller transform,
  ! which converts pairs of uniform random numbers into pairs of
  ! Gaussian random numbers. After generation:
  !   1. Centre-of-mass velocity is zeroed (no net momentum drift).
  !   2. Velocities are rescaled so the instantaneous kinetic
  !      temperature exactly equals target_t.
  ! ============================================================
  subroutine init_velocities_temperature(vx, vy, vz)
    real(dp), intent(out) :: vx(n), vy(n), vz(n)  ! output: atom velocities
    integer :: i
    real(dp) :: r1, r2, ga1, ga2
    real(dp) :: avx, avy, avz, kin, scale

    ! --- Box-Muller transform ---
    ! Each call produces two Gaussian-distributed values (ga1, ga2)
    ! from two uniform values (r1, r2), assigned to consecutive atoms.
    ! Stride of 2 means we fill N velocities in N/2 iterations.
    i = 1
    do while (i <= n)
       ! --- vx component ---
       call random_number(r1); call random_number(r2)
       ! Clamp r1 away from zero to avoid log(0).
       r1 = max(r1, 1.0e-12_dp)
       ! Box-Muller: ga1 and ga2 are independent N(0,1) samples.
       ga1 = sqrt(-2.0_dp*log(r1))*cos(2.0_dp*acos(-1.0_dp)*r2)
       ga2 = sqrt(-2.0_dp*log(r1))*sin(2.0_dp*acos(-1.0_dp)*r2)
       vx(i) = ga1
       if (i < n) vx(i+1) = ga2   ! assign second Gaussian to next atom

       ! --- vy component (fresh pair of uniforms) ---
       call random_number(r1); call random_number(r2)
       r1 = max(r1, 1.0e-12_dp)
       ga1 = sqrt(-2.0_dp*log(r1))*cos(2.0_dp*acos(-1.0_dp)*r2)
       ga2 = sqrt(-2.0_dp*log(r1))*sin(2.0_dp*acos(-1.0_dp)*r2)
       vy(i) = ga1
       if (i < n) vy(i+1) = ga2

       ! --- vz component ---
       call random_number(r1); call random_number(r2)
       r1 = max(r1, 1.0e-12_dp)
       ga1 = sqrt(-2.0_dp*log(r1))*cos(2.0_dp*acos(-1.0_dp)*r2)
       ga2 = sqrt(-2.0_dp*log(r1))*sin(2.0_dp*acos(-1.0_dp)*r2)
       vz(i) = ga1
       if (i < n) vz(i+1) = ga2

       ! Advance by 2 since each Box-Muller call fills two atoms.
       i = i + 2
    end do

    ! --- Remove centre-of-mass drift ---
    ! Compute mean velocity in each direction and subtract it.
    ! This ensures total linear momentum is zero, preventing the
    ! entire system from drifting across the box over time.
    avx = sum(vx)/real(n,dp)
    avy = sum(vy)/real(n,dp)
    avz = sum(vz)/real(n,dp)
    vx = vx - avx;  vy = vy - avy;  vz = vz - avz

    ! --- Rescale to exact target temperature ---
    ! KE = 0.5 * sum(v^2) in reduced units.
    kin = 0.5_dp*sum(vx*vx + vy*vy + vz*vz)

    ! From the equipartition theorem: KE = (3/2)*(N-1)*T* (using N-1
    ! degrees of freedom because COM drift has been removed).
    ! scale = sqrt(T*_target / T*_current).
    scale = sqrt((1.5_dp*(real(n,dp)-1.0_dp)*target_t)/kin)
    vx = scale*vx;  vy = scale*vy;  vz = scale*vz
  end subroutine init_velocities_temperature

  ! ============================================================
  ! build_verlet_list
  ! Constructs a full neighbor list: for each atom i, stores ALL
  ! atoms j (j /= i) within the extended cutoff r_list = rc + r_skin.
  ! 'Full' means both (i,j) and (j,i) are stored, so in compute_forces
  ! each thread can independently compute all forces on atom i without
  ! needing to write to atom j — eliminating OpenMP race conditions.
  !
  ! The build itself is parallelized: each thread independently fills
  ! the neighbor list for its assigned atoms i. Since each thread writes
  ! to vlist(:,i) and nlist(i) for its own i values only, there are no
  ! write conflicts between threads.
  !
  ! After building, reference positions x0_vl/y0_vl/z0_vl are saved
  ! so need_rebuild() can later measure cumulative atomic displacements.
  ! ============================================================
  subroutine build_verlet_list(x, y, z)
    real(dp), intent(in) :: x(n), y(n), z(n)  ! current atom positions
    integer  :: i, j, cnt
    real(dp) :: dx, dy, dz, rij2

    ! Parallelize over atom i; SCHEDULE(dynamic,16) distributes chunks
    ! of 16 iterations dynamically to balance load if neighbor counts
    ! vary across the box (e.g., near density fluctuations).
    ! All local variables (j, dx, dy, dz, rij2, cnt) are thread-private.
    !$OMP PARALLEL DO PRIVATE(j, dx, dy, dz, rij2, cnt) SCHEDULE(dynamic, 16)
    do i = 1, n
      cnt = 0  ! neighbor counter for atom i, reset for each i

      ! Scan every other atom j as a candidate neighbor of i.
      do j = 1, n
        if (j == i) cycle   ! skip self-interaction

        ! Displacement vector from j to i.
        dx = x(i) - x(j)
        dy = y(i) - y(j)
        dz = z(i) - z(j)

        ! Apply minimum image convention to get shortest-image displacement.
        call minimum_image(dx, dy, dz)

        ! Squared distance — avoids an unnecessary sqrt for the list check.
        rij2 = dx*dx + dy*dy + dz*dz

        ! If j is within the list radius, record it as a neighbor of i.
        if (rij2 < r_list2) then
          cnt = cnt + 1
          ! Guard against exceeding the fixed array bound max_neigh.
          if (cnt <= max_neigh) vlist(cnt, i) = j
        end if
      end do

      ! Store the final neighbor count for atom i, capped at max_neigh
      ! in case the array bound was exceeded (should not occur with
      ! max_neigh = 350 for these parameters).
      nlist(i) = min(cnt, max_neigh)
    end do
    !$OMP END PARALLEL DO

    ! Save current positions as the reference for drift detection.
    ! The next call to need_rebuild() will measure displacement from these.
    x0_vl = x;  y0_vl = y;  z0_vl = z
  end subroutine build_verlet_list

  ! ============================================================
  ! need_rebuild
  ! Determines whether the Verlet neighbor list must be rebuilt.
  ! The list is guaranteed valid as long as no atom has moved more
  ! than r_skin/2 since the last build — if any atom has drifted
  ! further, a formerly non-neighboring atom may have entered rc,
  ! meaning the list is stale. Returns .true. to trigger a rebuild.
  !
  ! Comparison is done on squared displacements to avoid sqrt calls.
  ! ============================================================
  logical function need_rebuild(x, y, z)
    real(dp), intent(in) :: x(n), y(n), z(n)  ! current atom positions
    real(dp) :: max_disp2, dx, dy, dz, disp2
    integer  :: i

    ! Threshold: squared half-skin distance. If max displacement^2
    ! exceeds this, the list safety guarantee is violated.
    real(dp), parameter :: half_skin2 = (0.5_dp*r_skin)**2

    max_disp2 = 0.0_dp

    ! Find the maximum squared displacement of any atom from its
    ! position at the time of the last list build.
    do i = 1, n
      ! Displacement since last build (raw Cartesian, no PBC needed here
      ! because we only care about the magnitude of drift, not direction).
      dx = x(i) - x0_vl(i)
      dy = y(i) - y0_vl(i)
      dz = z(i) - z0_vl(i)
      disp2 = dx*dx + dy*dy + dz*dz

      ! Track the maximum across all atoms.
      if (disp2 > max_disp2) max_disp2 = disp2
    end do

    ! Return .true. if the worst-case atom has drifted beyond r_skin/2.
    need_rebuild = (max_disp2 > half_skin2)
  end function need_rebuild

  ! ============================================================
  ! compute_forces
  ! Evaluates the total Lennard-Jones force on every atom and the
  ! total potential energy, using the pre-built Verlet neighbor list.
  !
  ! LJ pair potential (shifted):
  !   U(r) = 4*[ r^-12 - r^-6 ] - u_shift    for r < rc
  !          0                                  for r >= rc
  !
  ! LJ pair force on atom i from atom j:
  !   f_ij = [ 48*r^-14 - 24*r^-8 ] * r_vec
  !        = (48/r^2) * [ r^-12 - 0.5*r^-6 ] * r_vec
  !
  ! Full-list strategy: since both (i,j) and (j,i) are in the list,
  ! each thread computes ALL forces on its atom i without touching atom j.
  ! This is thread-safe by construction.  The potential is counted only
  ! for pairs where j > i (via the guard) to avoid double-counting.
  ! OpenMP REDUCTION(+:potential) gives each thread a private accumulator
  ! that is summed into the shared variable at the end of the parallel region.
  ! ============================================================
  subroutine compute_forces(x, y, z, fx, fy, fz, potential)
    real(dp), intent(in)  :: x(n), y(n), z(n)         ! atom positions
    real(dp), intent(out) :: fx(n), fy(n), fz(n)      ! output: forces
    real(dp), intent(out) :: potential                 ! output: total PE
    integer  :: i, k, j
    real(dp) :: dx, dy, dz, r2, invr2, invr6, invr12, fij_over_r

    ! Initialize total potential to zero before accumulation.
    potential = 0.0_dp

    ! Parallelize the outer loop over atom i.
    ! Each thread handles a static contiguous block of i-values.
    ! PRIVATE: each thread has its own copies of loop indices and
    !          intermediate LJ variables — no shared state in the arithmetic.
    ! REDUCTION(+:potential): each thread accumulates PE privately;
    !          OpenMP adds all thread-local sums at the barrier.
    !$OMP PARALLEL DO PRIVATE(k, j, dx, dy, dz, r2, invr2, invr6, invr12, fij_over_r) &
    !$OMP             REDUCTION(+:potential) SCHEDULE(static)
    do i = 1, n
      ! Initialize force on atom i to zero at the start of each i-iteration.
      ! This is thread-safe because only one thread ever processes a given i.
      fx(i) = 0.0_dp;  fy(i) = 0.0_dp;  fz(i) = 0.0_dp

      ! Loop over all neighbors of atom i stored in the Verlet list.
      do k = 1, nlist(i)
        j = vlist(k, i)   ! retrieve the k-th neighbor index of atom i

        ! Displacement vector from atom j to atom i.
        dx = x(i) - x(j)
        dy = y(i) - y(j)
        dz = z(i) - z(j)

        ! Shortest-image displacement under periodic boundary conditions.
        call minimum_image(dx, dy, dz)

        ! Squared distance — avoids sqrt for the cutoff check.
        r2 = dx*dx + dy*dy + dz*dz

        ! Only compute interaction if pair is within the force cutoff rc.
        ! Neighbors in (rc, r_list] are stored but contribute no force.
        if (r2 < rc2) then

          ! Precompute inverse powers of r for the LJ expressions.
          invr2  = 1.0_dp/r2        ! 1/r^2
          invr6  = invr2*invr2*invr2 ! 1/r^6
          invr12 = invr6*invr6       ! 1/r^12

          ! Accumulate potential energy only once per pair (j > i guard).
          ! The full list stores both directions, so without this guard
          ! every pair would be counted twice.
          if (j > i) potential = potential + 4.0_dp*(invr12 - invr6) - u_shift

          ! Magnitude of force divided by r: F(r)/r = 48/r^2 * (1/r^12 - 0.5/r^6).
          ! Multiplying by the displacement vector dx/dy/dz gives the
          ! Cartesian force components via f = (F/r) * r_vec.
          fij_over_r = 48.0_dp*invr2*(invr12 - 0.5_dp*invr6)

          ! Accumulate force on atom i from atom j.
          ! Note: we do NOT update fx(j) here (full-list approach).
          ! When j is processed in its own loop iteration (with i as neighbor),
          ! the equal and opposite force will be added to fx(j) naturally.
          fx(i) = fx(i) + fij_over_r*dx
          fy(i) = fy(i) + fij_over_r*dy
          fz(i) = fz(i) + fij_over_r*dz
        end if
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine compute_forces

  ! ============================================================
  ! velocity_verlet_step
  ! Advances positions and velocities by one time step dt using the
  ! velocity Verlet (Störmer-Verlet) algorithm, which is time-reversible
  ! and symplectic (conserves phase-space volume), making it well-suited
  ! for NVE MD.
  !
  ! Algorithm (two half-kicks with a full position update between):
  !   Step 1 — position update + first half velocity kick:
  !     r(t+dt) = r(t) + v(t)*dt + 0.5*F(t)*dt^2
  !     v(t+dt/2) = v(t) + 0.5*F(t)*dt
  !   Step 2 — recompute forces at new positions:
  !     F(t+dt) = compute_forces( r(t+dt) )
  !   Step 3 — second half velocity kick:
  !     v(t+dt) = v(t+dt/2) + 0.5*F(t+dt)*dt
  !
  ! Both the position loop and the velocity loops are OpenMP-parallelized
  ! over atoms with static scheduling (equal work per atom).
  ! ============================================================
  subroutine velocity_verlet_step(x, y, z, vx, vy, vz, fx, fy, fz, potential)
    real(dp), intent(inout) :: x(n), y(n), z(n)        ! positions (updated in place)
    real(dp), intent(inout) :: vx(n), vy(n), vz(n)     ! velocities (updated in place)
    real(dp), intent(inout) :: fx(n), fy(n), fz(n)     ! forces at time t (input)
    real(dp), intent(out)   :: potential                ! PE at new positions (output)

    ! Temporary arrays to hold forces at the new positions r(t+dt).
    real(dp) :: fx_new(n), fy_new(n), fz_new(n)
    integer :: i

    ! --- Step 1: position update and first half-kick ---
    ! Each atom's update is independent, so the loop is trivially parallel.
    !$OMP PARALLEL DO SCHEDULE(static)
    do i = 1, n
       ! Full position update: r(t+dt) = r(t) + v(t)*dt + 0.5*F(t)*dt^2
       x(i) = x(i) + vx(i)*dt + 0.5_dp*fx(i)*dt*dt
       y(i) = y(i) + vy(i)*dt + 0.5_dp*fy(i)*dt*dt
       z(i) = z(i) + vz(i)*dt + 0.5_dp*fz(i)*dt*dt

       ! Fold updated position back into the primary box [0, boxlength).
       call wrap_position(x(i), y(i), z(i))

       ! First half-kick: v(t+dt/2) = v(t) + 0.5*F(t)*dt
       vx(i) = vx(i) + 0.5_dp*fx(i)*dt
       vy(i) = vy(i) + 0.5_dp*fy(i)*dt
       vz(i) = vz(i) + 0.5_dp*fz(i)*dt
    end do
    !$OMP END PARALLEL DO

    ! --- Step 2: recompute forces at r(t+dt) ---
    ! compute_forces uses the Verlet list; the list validity check and
    ! rebuild are handled in the main loop before calling this subroutine.
    call compute_forces(x, y, z, fx_new, fy_new, fz_new, potential)

    ! --- Step 3: second half-kick ---
    ! v(t+dt) = v(t+dt/2) + 0.5*F(t+dt)*dt
    !$OMP PARALLEL DO SCHEDULE(static)
    do i = 1, n
       vx(i) = vx(i) + 0.5_dp*fx_new(i)*dt
       vy(i) = vy(i) + 0.5_dp*fy_new(i)*dt
       vz(i) = vz(i) + 0.5_dp*fz_new(i)*dt
    end do
    !$OMP END PARALLEL DO

    ! Replace old forces with the newly computed forces so the next
    ! velocity_verlet_step call has F(t) = F(t+dt) of this call.
    fx = fx_new;  fy = fy_new;  fz = fz_new
  end subroutine velocity_verlet_step

  ! ============================================================
  ! kinetic_energy
  ! Returns total kinetic energy KE = 0.5 * sum_i |v_i|^2
  ! in reduced units. Declared 'pure' (no side effects) so the
  ! compiler can inline and optimize it freely.
  ! ============================================================
  pure function kinetic_energy(vx, vy, vz) result(ke)
    real(dp), intent(in) :: vx(n), vy(n), vz(n)
    real(dp) :: ke
    ! Sum of squared speeds across all atoms, scaled by 0.5 (mass = 1 in LJ units).
    ke = 0.5_dp*sum(vx*vx + vy*vy + vz*vz)
  end function kinetic_energy

  ! ============================================================
  ! temperature_from_ke
  ! Converts kinetic energy to instantaneous reduced temperature T*
  ! using the equipartition theorem:
  !   KE = (3/2) * (N-1) * T*
  ! N-1 degrees of freedom are used (not 3N) because the
  ! centre-of-mass velocity has been zeroed, removing 3 DOF.
  ! ============================================================
  pure function temperature_from_ke(ke) result(temp)
    real(dp), intent(in) :: ke
    real(dp) :: temp
    temp = 2.0_dp*ke/(3.0_dp*real(n-1,dp))
  end function temperature_from_ke

  ! ============================================================
  ! rescale_velocities
  ! Implements simple velocity rescaling (Berendsen-like thermostat
  ! in the limit of instantaneous coupling): multiplies all velocities
  ! by sqrt(T*_target / T*_current) so the next temperature measurement
  ! will read exactly t_target. Used only during equilibration.
  ! ============================================================
  subroutine rescale_velocities(vx, vy, vz, t_target)
    real(dp), intent(inout) :: vx(n), vy(n), vz(n)  ! velocities to rescale
    real(dp), intent(in)    :: t_target              ! desired reduced temperature

    real(dp) :: ke, tcur, scale

    ! Compute current kinetic energy and instantaneous temperature.
    ke   = kinetic_energy(vx, vy, vz)
    tcur = temperature_from_ke(ke)

    ! Guard against divide-by-zero if all velocities are zero.
    if (tcur <= 0.0_dp) return

    ! Scale factor: sqrt(T*_target / T*_current).
    scale = sqrt(t_target/tcur)

    ! Apply uniform rescaling to all velocity components.
    vx = scale*vx;  vy = scale*vy;  vz = scale*vz
  end subroutine rescale_velocities

  ! ============================================================
  ! accumulate_gr
  ! Adds the pair distance counts from the current configuration
  ! to the running g(r) histogram gr_hist.
  !
  ! The radial distribution function g(r) is defined so that
  ! rho * g(r) * 4*pi*r^2 dr is the average number of atoms in
  ! a shell [r, r+dr] around a given atom.
  !
  ! We use the Verlet list (with j > i guard) to visit each pair
  ! exactly once, then add 2.0 to the bin so that both the i->j
  ! and j->i contributions are counted — matching the normalization
  ! in write_gr which treats each atom as a centre.
  !
  ! The outer loop over i is OpenMP-parallelized; gr_hist is an
  ! array reduction so each thread accumulates a private copy of
  ! the histogram and they are summed at the end of the parallel region.
  ! ============================================================
  subroutine accumulate_gr(x, y, z, gr_hist)
    real(dp), intent(in)    :: x(n), y(n), z(n)         ! current positions
    real(dp), intent(inout) :: gr_hist(nbins_gr)         ! running histogram (updated)
    integer  :: i, k, j, bin
    real(dp) :: dx, dy, dz, rij

    ! Array REDUCTION(+:gr_hist) requires OpenMP 4.5+.
    ! Each thread has a private zero-initialized copy of gr_hist;
    ! the copies are added together after the parallel region.
    !$OMP PARALLEL DO PRIVATE(k, j, dx, dy, dz, rij, bin) &
    !$OMP             REDUCTION(+:gr_hist) SCHEDULE(static)
    do i = 1, n
      do k = 1, nlist(i)
        j = vlist(k, i)

        ! Count each pair only once (j > i); the factor 2.0 below
        ! accounts for both the i-centred and j-centred contributions.
        if (j <= i) cycle

        ! Shortest-image displacement between i and j.
        dx = x(i)-x(j);  dy = y(i)-y(j);  dz = z(i)-z(j)
        call minimum_image(dx, dy, dz)
        rij = sqrt(dx*dx + dy*dy + dz*dz)

        ! Only bin distances up to halfbox; beyond that the minimum
        ! image sphere is not fully contained in the box, so g(r)
        ! would be biased by image artifacts.
        if (rij < halfbox) then
          ! Determine which bin this distance falls in (1-indexed).
          bin = int(rij/dr_gr) + 1

          ! Safety check to keep within array bounds.
          if (bin >= 1 .and. bin <= nbins_gr) then
            ! Add 2.0: counts both the (i,j) and (j,i) contribution
            ! so normalization in write_gr uses N as the number of centres.
            gr_hist(bin) = gr_hist(bin) + 2.0_dp
          end if
        end if
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine accumulate_gr

  ! ============================================================
  ! write_gr
  ! Normalizes the accumulated g(r) histogram and writes the result
  ! to 'gr_output.dat' with columns: r_mid, g(r), raw_counts.
  !
  ! Normalization:
  !   g(r) = hist(b) / [ n_samples * N * rho * shell_vol(b) ]
  !
  !   - n_samples: number of configurations accumulated
  !   - N * rho * shell_vol: expected count in the bin for an
  !     ideal gas (uniform density), i.e., the reference value
  !   - Dividing observed by expected gives g(r) -> 1 at large r
  !     for a disordered liquid.
  ! ============================================================
  subroutine write_gr(gr_hist, n_samples)
    real(dp), intent(in) :: gr_hist(nbins_gr)   ! accumulated raw histogram
    integer,  intent(in) :: n_samples           ! number of sampled frames
    integer :: b
    real(dp) :: r_lower, r_upper, r_mid, shell_vol, ideal_count, gval
    real(dp) :: volume, number_density
    integer :: unitno

    ! Total box volume and number density (atoms per unit volume).
    volume         = boxlength**3
    number_density = real(n,dp)/volume

    open(newunit=unitno, file='gr_output.dat', status='replace', action='write')
    write(unitno,'(a)') '# r_mid  g(r)  raw_counts'

    do b = 1, nbins_gr
       ! Radial boundaries and midpoint of this bin.
       r_lower   = (b-1)*dr_gr
       r_upper   = b*dr_gr
       r_mid     = 0.5_dp*(r_lower + r_upper)

       ! Volume of the spherical shell between r_lower and r_upper.
       ! V_shell = (4/3)*pi*(r_upper^3 - r_lower^3)
       shell_vol = (4.0_dp/3.0_dp)*acos(-1.0_dp)*(r_upper**3 - r_lower**3)

       ! Expected number of pair counts in this bin per sample for
       ! an ideal gas: N * rho * shell_vol (both directions counted).
       ideal_count = real(n,dp) * number_density * shell_vol

       ! Compute normalized g(r), guarding against division by zero.
       if (n_samples > 0 .and. ideal_count > 0.0_dp) then
          gval = gr_hist(b) / (real(n_samples,dp) * ideal_count)
       else
          gval = 0.0_dp
       end if

       ! Write: bin midpoint, normalized g(r), and raw accumulated count.
       write(unitno,'(3(1x,es20.10))') r_mid, gval, gr_hist(b)
    end do
    close(unitno)
  end subroutine write_gr

  ! ============================================================
  ! write_xyz_frame
  ! Appends one snapshot of all atom positions to an open XYZ file.
  ! The XYZ format is:
  !   Line 1: number of atoms
  !   Line 2: comment line (step, simulation time, box length)
  !   Lines 3..N+2: element symbol and x y z coordinates
  ! This format is readable by OVITO, VMD, and most MD visualization tools.
  ! ============================================================
  subroutine write_xyz_frame(unitno, step, time, x, y, z)
    integer,  intent(in) :: unitno   ! file unit number (already open)
    integer,  intent(in) :: step     ! current MD step number
    real(dp), intent(in) :: time     ! current simulation time in reduced units
    real(dp), intent(in) :: x(n), y(n), z(n)   ! atom positions
    integer :: i

    ! First line: atom count (required by XYZ format).
    write(unitno,'(i0)') n

    ! Second line: metadata comment — step index, time, and box size.
    write(unitno,'(a,i0,a,es20.10,a,f18.8)') &
         'step=', step, ' time=', time, ' boxlength=', boxlength

    ! Atom lines: element symbol 'Ar' followed by x, y, z coordinates.
    do i = 1, n
       write(unitno,'(a,3(1x,f18.10))') 'Ar', x(i), y(i), z(i)
    end do
  end subroutine write_xyz_frame

end module md_params


! ============================================================
! Main program: rahman_lj_argon_omp
!
! Orchestrates the full MD run:
!   1. Initialize positions (random, overlap-free) and velocities
!      (Maxwell-Boltzmann at target_t).
!   2. Build the Verlet neighbor list and compute initial forces.
!   3. Equilibration loop (NVT): run velocity_verlet_step with
!      velocity rescaling at each step to drive the system to target_t.
!   4. Production loop (NVE): run without thermostat; accumulate
!      g(r) histogram and write energy/trajectory output.
!   5. Normalize and write g(r), report wall-clock time and
!      number of neighbor list rebuilds.
! ============================================================
program rahman_lj_argon_omp
  use md_params    ! imports all parameters, arrays, and subroutines
  implicit none

  ! --- Atom state arrays ---
  ! Positions, velocities, and forces for all N atoms (three Cartesian components each).
  real(dp) :: x(n), y(n), z(n)       ! positions
  real(dp) :: vx(n), vy(n), vz(n)    ! velocities
  real(dp) :: fx(n), fy(n), fz(n)    ! forces

  ! --- Thermodynamic observables ---
  real(dp) :: potential   ! total LJ potential energy
  real(dp) :: ke          ! total kinetic energy
  real(dp) :: te          ! total energy (KE + PE)
  real(dp) :: temp        ! instantaneous reduced temperature T*
  real(dp) :: time        ! simulation clock in reduced time units

  ! --- g(r) accumulation ---
  real(dp) :: gr_hist(nbins_gr)   ! running histogram of pair distances

  ! --- Control and I/O variables ---
  integer  :: step          ! current MD step counter
  integer  :: unit_energy   ! file unit for energy_output.dat
  integer  :: unit_xyz      ! file unit for trajectory.xyz
  integer  :: n_gr_samples  ! number of configurations sampled for g(r)
  integer  :: n_rebuilds    ! total number of Verlet list rebuilds (diagnostic)

  ! --- Timing ---
  real(dp) :: t_start, t_end   ! wall-clock times from omp_get_wtime()

  ! Record wall-clock start time using the OpenMP high-resolution timer.
  t_start = omp_get_wtime()

  ! --- Initialization ---
  ! Place atoms randomly with minimum separation min_sep_init to avoid overlaps.
  call init_positions_random(x, y, z)

  ! Assign velocities from a Maxwell-Boltzmann distribution at target_t,
  ! with COM drift removed and exact temperature scaling applied.
  call init_velocities_temperature(vx, vy, vz)

  ! Build the initial Verlet neighbor list from the starting positions.
  ! This must be done before the first force evaluation.
  call build_verlet_list(x, y, z)

  ! Compute initial forces using the freshly built neighbor list.
  call compute_forces(x, y, z, fx, fy, fz, potential)

  ! Compute initial thermodynamic quantities for the startup printout.
  ke   = kinetic_energy(vx, vy, vz)
  te   = ke + potential
  temp = temperature_from_ke(ke)

  ! --- Print run header ---
  write(*,'(a)')        'Rahman (1964) LJ Argon — Verlet list + OpenMP'
  write(*,'(a,i0)')     ' Threads        = ', omp_get_max_threads()
  write(*,'(a,i0)')     ' N              = ', n
  write(*,'(a,f18.8)')  ' rho*           = ', rho
  write(*,'(a,f18.8)')  ' boxlength      = ', boxlength
  write(*,'(a,f18.8)')  ' rc*            = ', rc
  write(*,'(a,f18.8)')  ' r_skin         = ', r_skin
  write(*,'(a,es18.8)') ' dt*            = ', dt
  write(*,'(a,f18.8)')  ' target T*      = ', target_t
  write(*,'(a,f18.8)')  ' initial PE     = ', potential
  write(*,'(a,f18.8)')  ' initial KE     = ', ke
  write(*,'(a,f18.8)')  ' initial total E= ', te
  write(*,'(a,f18.8)')  ' initial T*     = ', temp
  write(*,'(a,i0)')     ' equil steps    = ', n_equil
  write(*,'(a,i0)')     ' prod steps     = ', n_prod

  ! --- Open output files ---
  ! energy_output.dat: one row per step with time, total E, PE, KE, T*.
  open(newunit=unit_energy, file='energy_output.dat', status='replace', action='write')
  write(unit_energy,'(a)') '# step  time  total_energy  potential_energy  kinetic_energy  temperature'

  ! trajectory.xyz: successive XYZ frames for visualization.
  open(newunit=unit_xyz, file='trajectory.xyz', status='replace', action='write')

  ! ============================================================
  ! Equilibration phase (NVT with velocity rescaling)
  ! Runs n_equil steps during which velocities are rescaled to
  ! target_t after every step. This drives the system from its
  ! (possibly far-from-equilibrium) initial state into a thermally
  ! equilibrated liquid configuration before data collection begins.
  ! ============================================================
  n_rebuilds = 0
  time = 0.0_dp

  do step = 1, n_equil
     ! Check whether any atom has drifted more than r_skin/2 since
     ! the last list build; if so, rebuild before advancing the step.
     if (need_rebuild(x, y, z)) then
        call build_verlet_list(x, y, z)
        n_rebuilds = n_rebuilds + 1   ! track rebuild frequency for diagnostics
     end if

     ! Advance positions and velocities by one dt using velocity Verlet.
     call velocity_verlet_step(x, y, z, vx, vy, vz, fx, fy, fz, potential)

     ! Rescale velocities so instantaneous T* = target_t (NVT thermostat).
     call rescale_velocities(vx, vy, vz, target_t)

     ! Recompute thermodynamic observables after the rescale.
     ke   = kinetic_energy(vx, vy, vz)
     te   = ke + potential
     temp = temperature_from_ke(ke)

     ! Advance simulation clock.
     time = time + dt
  end do

  ! After equilibration, recompute forces from the clean post-equilibration
  ! state to ensure fx/fy/fz are consistent with the current positions
  ! before entering the NVE production run.
  call compute_forces(x, y, z, fx, fy, fz, potential)

  ! Reset g(r) histogram and sample counter — only production data enters g(r).
  gr_hist      = 0.0_dp
  n_gr_samples = 0

  ! ============================================================
  ! Production phase (NVE — no thermostat)
  ! Runs n_prod steps of constant-energy MD.  Energy, g(r), and
  ! trajectory data are written to disk during this phase.
  ! The loop starts at step=0 so that the post-equilibration state
  ! is recorded as the first data point before any integration occurs.
  ! ============================================================
  do step = 0, n_prod

     ! --- Compute and record thermodynamic observables ---
     ke   = kinetic_energy(vx, vy, vz)
     te   = ke + potential
     temp = temperature_from_ke(ke)

     ! Write one row per step: step, time, E_total, E_pot, E_kin, T*.
     write(unit_energy,'(i10,1x,5(es22.12,1x))') step, time, te, potential, ke, temp

     ! --- Accumulate g(r) every sample_every steps ---
     ! Averaging over many configurations reduces statistical noise
     ! in the radial distribution function.
     if (mod(step, sample_every) == 0) then
        call accumulate_gr(x, y, z, gr_hist)
        n_gr_samples = n_gr_samples + 1
     end if

     ! --- Write XYZ trajectory frame every trajectory_every steps ---
     ! Saving every frame would be expensive in storage; subsampling
     ! still gives smooth animations while keeping file sizes manageable.
     if (mod(step, trajectory_every) == 0) then
        call write_xyz_frame(unit_xyz, step, time, x, y, z)
     end if

     ! --- Advance to the next step (skip on the final step) ---
     ! The guard (step < n_prod) ensures we record data at step n_prod
     ! before exiting without performing a spurious extra integration step.
     if (step < n_prod) then
        ! Rebuild Verlet list if any atom has drifted beyond r_skin/2.
        if (need_rebuild(x, y, z)) then
           call build_verlet_list(x, y, z)
           n_rebuilds = n_rebuilds + 1
        end if

        ! Propagate the system by one time step dt.
        call velocity_verlet_step(x, y, z, vx, vy, vz, fx, fy, fz, potential)

        ! Advance simulation clock by one time step.
        time = time + dt
     end if
  end do

  ! --- Finalize output ---
  close(unit_energy)   ! flush and close energy_output.dat
  close(unit_xyz)      ! flush and close trajectory.xyz

  ! Normalize the accumulated g(r) histogram and write gr_output.dat.
  call write_gr(gr_hist, n_gr_samples)

  ! Record wall-clock end time.
  t_end = omp_get_wtime()

  ! --- Final run statistics ---
  write(*,'(a,i0)')    ' g(r) samples   = ', n_gr_samples
  ! A rebuild count near n_total/30 is healthy; far fewer suggests
  ! r_skin is too large; far more suggests it is too small.
  write(*,'(a,i0)')    ' list rebuilds  = ', n_rebuilds
  ! Wall time includes initialization, equilibration, and production.
  write(*,'(a,f12.2)') ' Wall time (s)  = ', t_end - t_start
  write(*,'(a)')       ' Wrote energy_output.dat, gr_output.dat, and trajectory.xyz'

end program rahman_lj_argon_omp
