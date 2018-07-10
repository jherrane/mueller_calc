module setup
   use common
   use io
   use mie
   use translations
   use projection

   implicit none
contains

!****************************************************************************80

   subroutine init_geometry(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      type(data_struct), dimension(:), allocatable :: sphere
      real(dp) :: ka, maxrad, maxrad_sph, vol
      integer :: i, Nspheres, sph, max_level, max_sph
      complex(dp) :: k

      if (matrices%Tmat == 0) print *, 'Order of basis functions  =', mesh%order
      if (mesh%M_ex > 3 .or. mesh%M_ex < 1) then
         print *, 'ERROR: order should be 1, 2 or 3'
         stop
      end if
      if (mesh%order > 1 .or. mesh%M_ex < 0) then
         print *, 'ERROR: Expansion order should be 0 or 1 '
         stop
      end if

      print *, 'Reading mesh...'
      call read_mesh(matrices, mesh) ! io
      call int_points(mesh)

      if (matrices%Tmat == 0) print *, 'Initializing FFT... '
      if (allocated(mesh%nodes)) deallocate (mesh%nodes, mesh%etopol_box, mesh%tetras)
      call build_grid2(mesh) ! geometry
      call build_box(mesh) ! geometry
      call tetras_in_cubes(mesh) ! geometry

      if (matrices%Tmat == 0) then
         print *, '   Grid size            = ', (/mesh%Nx, mesh%Ny, mesh%Nz/)
         print *, '   Delta grid           = ', real(mesh%delta)
         print *, '   Number of cubes      = ', mesh%N_cubes
         print *, '   Delta cubes          = ', real(mesh%box_delta)
         print *, '   Elems. in cube (max) = ', mesh%N_tet_cube
         print *, 'Done'
      end if

      if (matrices%Tmat == 0) call construct_projectors(matrices, mesh) ! projection

! Check whether using maxval is good or not
      do i = 1, matrices%bars
         ka = mesh%ki(i)*(dble(maxval([mesh%Nx, mesh%Ny, mesh%Nz])) &
                          *mesh%delta)/2.0d0
         matrices%Nmaxs(i) = truncation_order(ka)
      end do

   end subroutine init_geometry

!****************************************************************************80

   subroutine construct_projectors(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, M_box
      real(dp) :: k_orig

      k_orig = mesh%k ! Keep the original just in case
      M_box = size(mesh%etopol_box, 1)

      if (.not. allocated(matrices%rhs)) then
         if (mesh%order == 0) then
            allocate (matrices%rhs(3*mesh%N_tet))
            allocate (matrices%x(3*mesh%N_tet))
            allocate (matrices%Ax(3*mesh%N_tet))

            allocate (matrices%listS(M_box, mesh%N_tet, matrices%bars))
            allocate (matrices%listSx(M_box, mesh%N_tet, matrices%bars))
            allocate (matrices%listSy(M_box, mesh%N_tet, matrices%bars))
            allocate (matrices%listSz(M_box, mesh%N_tet, matrices%bars))
            allocate (matrices%listindS(M_box, mesh%N_tet, matrices%bars))
         end if
         if (mesh%order == 1) then
            allocate (matrices%rhs(4*3*mesh%N_tet))
            allocate (matrices%x(4*3*mesh%N_tet))
            allocate (matrices%Ax(4*3*mesh%N_tet))

            allocate (matrices%listS(M_box, 4*mesh%N_tet, matrices%bars))
            allocate (matrices%listSx(M_box, 4*mesh%N_tet, matrices%bars))
            allocate (matrices%listSy(M_box, 4*mesh%N_tet, matrices%bars))
            allocate (matrices%listSz(M_box, 4*mesh%N_tet, matrices%bars))
            allocate (matrices%listindS(M_box, 4*mesh%N_tet, matrices%bars))
         end if
      end if
      print *, 'Constructing ', trim(mesh%projector), '-projectors... '

      do i = 1, size(mesh%ki)
         if (allocated(matrices%S)) then
            deallocate (matrices%S, matrices%Sx, matrices%Sy, matrices%Sz, matrices%indS)
         end if
         mesh%k = mesh%ki(i)
         call print_bar(i, size(mesh%ki))
         if (mesh%order == 1) call pfft_projection_lin(matrices, mesh)
         if (mesh%order == 0) call pfft_projection_const(matrices, mesh)
         matrices%listS(:, :, i) = matrices%S
         matrices%listSx(:, :, i) = matrices%Sx
         matrices%listSy(:, :, i) = matrices%Sy
         matrices%listSz(:, :, i) = matrices%Sz
         matrices%listindS(:, :, i) = matrices%indS
      end do

      mesh%k = k_orig

   end subroutine construct_projectors

!****************************************************************************80

   subroutine update_projections(matrices, mesh, i)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i

      matrices%S = matrices%listS(:, :, i)
      matrices%Sx = matrices%listSx(:, :, i)
      matrices%Sy = matrices%listSy(:, :, i)
      matrices%Sz = matrices%listSz(:, :, i)
      matrices%indS = matrices%listindS(:, :, i)

   end subroutine update_projections

!****************************************************************************80

   subroutine allocate_inc_wave(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      complex(dp), dimension(:), allocatable :: a_in, b_in
      integer :: Nmax, i, las, nm

      Nmax = maxval(matrices%Nmaxs)

      allocate (matrices%as((Nmax + 1)**2 - 1, matrices%bars))
      allocate (matrices%bs((Nmax + 1)**2 - 1, matrices%bars))

      las = (Nmax + 1)*(2*Nmax + 1)*(2*Nmax + 3)/3 - 1

      allocate (matrices%rotDs(las, matrices%bars))
      allocate (matrices%indDs(las, 2, matrices%bars))
      allocate (matrices%rotD90s(las, matrices%bars))
      allocate (matrices%indD90s(las, 2, matrices%bars))
      allocate (matrices%rotXs(las, matrices%bars))
      allocate (matrices%indXs(las, 2, matrices%bars))
      allocate (matrices%rotYs(las, matrices%bars))
      allocate (matrices%indYs(las, 2, matrices%bars))

      do i = 1, matrices%bars
         if (allocated(a_in)) deallocate (a_in, b_in)
         Nmax = matrices%Nmaxs(i)
         nm = (Nmax + 1)**2 - 1
         allocate (a_in(nm))
         allocate (b_in(nm))
         call planewave(Nmax, dcmplx(mesh%ki(i)), a_in, b_in)

         matrices%as(1:nm, i) = a_in
         matrices%bs(1:nm, i) = b_in
      end do

   end subroutine allocate_inc_wave

!****************************************************************************80

   subroutine update_inc_wave(matrices, mesh, i, a_nm, b_nm, a_nm90, b_nm90)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      complex(dp), dimension(:), allocatable :: a_in, b_in, a90, b90, &
                                                a, b, a_nm, b_nm, a_nm90, b_nm90
      complex(dp), dimension(:, :), allocatable :: Taa, Tab, Tba, Tbb
      complex(8), dimension(:), allocatable :: rotD, rotD90
      integer, dimension(:, :), allocatable :: indD, indD90
      integer :: Nmax, i, las, nm

      matrices%E0 = matrices%E_rel(i)*matrices%E*matrices%E0hat
      matrices%E90 = matrices%E_rel(i)*matrices%E*matrices%E90hat
      mesh%k = mesh%ki(i)

      Nmax = matrices%Nmaxs(i)
      las = (Nmax + 1)*(2*Nmax + 1)*(2*Nmax + 3)/3 - 1
      nm = (Nmax + 1)**2 - 1
      rotD = matrices%rotDs(1:las, i)
      indD = matrices%indDs(1:las, :, i)
      rotD90 = matrices%rotD90s(1:las, i)
      indD90 = matrices%indD90s(1:las, :, i)

      Taa = matrices%Taai(1:nm, 1:nm, i)
      Tab = matrices%Tabi(1:nm, 1:nm, i)
      Tba = matrices%Tbai(1:nm, 1:nm, i)
      Tbb = matrices%Tbbi(1:nm, 1:nm, i)

      a_in = matrices%E_rel(i)*matrices%E*matrices%as(1:nm, i)
      b_in = matrices%E_rel(i)*matrices%E*matrices%bs(1:nm, i)

      a = sparse_matmul(rotD, indD, a_in, nm)
      b = sparse_matmul(rotD, indD, b_in, nm)
      a90 = sparse_matmul(rotD90, indD90, a_in, nm)
      b90 = sparse_matmul(rotD90, indD90, b_in, nm)

      a_nm = matmul(Taa, a) + matmul(Tab, b)
      b_nm = matmul(Tbb, b) + matmul(Tba, a)
      a_nm90 = matmul(Taa, a90) + matmul(Tab, b90)
      b_nm90 = matmul(Tbb, b90) + matmul(Tba, a90)

   end subroutine update_inc_wave

!****************************************************************************80

   subroutine rot_setup(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      real(dp), dimension(3, 3) ::  RT, RT90

      RT = transpose(matrices%R)
      RT90 = matmul(RT, matrices%R90_init)

      matrices%Rk = matrices%R
      matrices%khat = matmul(RT, [0d0, 0d0, 1d0])
      matrices%E0hat = dcmplx(matmul(RT, [1d0, 0d0, 0d0]))
      matrices%E90hat = dcmplx(matmul(RT, [0d0, 1d0, 0d0]))

      call gen_rotations(matrices, mesh, RT, RT90)

   end subroutine rot_setup

!****************************************************************************80

   subroutine gen_rotations(matrices, mesh, R, R90)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: las, Nmax, i
      complex(8), dimension(:), allocatable :: rotD, rotD90, rotX, rotY
      integer, dimension(:, :), allocatable :: indD, indD90, indX, indY
      real(dp), dimension(3, 3) :: R, R90, Rx, Ry

      Rx = reshape([0d0, 0d0, 1d0, 0d0, 1d0, 0d0, -1d0, 0d0, 0d0], [3, 3])
      Ry = reshape([1d0, 0d0, 0d0, 0d0, 0d0, 1d0, 0d0, -1d0, 0d0], [3, 3])

      do i = 1, matrices%bars
         Nmax = matrices%Nmaxs(i)
         las = (Nmax + 1)*(2*Nmax + 1)*(2*Nmax + 3)/3 - 1

         if (allocated(rotD)) deallocate (rotD, indD, rotD90, indD90, rotX, indX, rotY, indY)
         allocate (rotD(las))
         allocate (indD(las, 2))
         allocate (rotD90(las))
         allocate (indD90(las, 2))
         allocate (rotX(las))
         allocate (indX(las, 2))
         allocate (rotY(las))
         allocate (indY(las, 2))

         call sph_rotation_sparse_gen(mat2euler(R), Nmax, rotD, indD)
         call sph_rotation_sparse_gen(mat2euler(R90), Nmax, rotD90, indD90)
         call sph_rotation_sparse_gen(mat2euler(Rx), Nmax, rotX, indX)
         call sph_rotation_sparse_gen(mat2euler(Ry), Nmax, rotY, indY)

         matrices%rotDs(1:las, i) = rotD
         matrices%indDs(1:las, :, i) = indD
         matrices%rotD90s(1:las, i) = rotD90
         matrices%indD90s(1:las, :, i) = indD90
         matrices%rotXs(1:las, i) = rotX
         matrices%indXs(1:las, :, i) = indX
         matrices%rotYs(1:las, i) = rotY
         matrices%indYs(1:las, :, i) = indY
      end do

   end subroutine gen_rotations

!****************************************************************************80

   subroutine int_points(mesh)
      type(mesh_struct) :: mesh
      real(dp), dimension(:, :), allocatable :: P
      real(dp), dimension(:), allocatable :: w
      real(dp) :: d, d2
      integer :: i1

      d = 0d0
      do i1 = 1, size(mesh%coord, 2)
         d2 = vlen(mesh%coord(:, i1))
         if (d2 > d) then
            d = d2
         end if
      end do
      d = 1.05d0*d

      call sample_points(P, w, 20, 20)
      if (.not. allocated(mesh%P)) then
         allocate (mesh%P(size(P, 1), size(P, 2)), mesh%w(size(w)))
      end if

! These points should enclose the object
      mesh%P = P*d !* mesh%a
      mesh%w = w*d**2 !* mesh%a**2

   end subroutine int_points

!****************************************************************************80

   subroutine polarization(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      real(dp) :: k0(3), k1(3), R(3, 3)

      k0 = [0d0, 0d0, 1d0]
      k1 = matrices%khat
      R = rotate_a_to_b(k0, k1)

      matrices%E0hat = dcmplx(matmul(R, [1d0, 0d0, 0d0]))
      matrices%E90hat = dcmplx(matmul(R, [0d0, 1d0, 0d0]))

   end subroutine polarization

!****************************************************************************80
! Setup the wavelength band and all matrices inv A in it
   subroutine wavelength_band(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      real(dp) :: lmax
      real(dp), dimension(:), allocatable :: c, diff, absdif
      integer :: i, n

      n = matrices%bars
      allocate (c(n))
      allocate (diff(n))
      allocate (absdif(n))

      call linspace(dlog10(matrices%lambda1), dlog10(matrices%lambda2), n, c)
      do i = 1, n
         c(i) = 10**c(i)
      end do

! Calculate the final wavelengths and the corresponding wave numbers
      do i = 1, n
         mesh%ki(i) = 2d0*pi/c(i)
      end do
      matrices%E_rel = 1d0

   end subroutine wavelength_band


!****************************************************************************80

   subroutine diagonalize_inertia(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, imin, imax, imid, negs

! Test if I is "almost diagonal", which causes diasym to do funny stuff
      if (dabs(mesh%I(1, 1)*mesh%I(2, 2)*mesh%I(3, 3)) > 1d6*dabs(mesh%I(1, 2)*mesh%I(1, 3)*mesh%I(2, 3))) then
         print *, 'FALLBACK MODE IN DIAGONALIZATION OF INERTIA TENSOR'
         matrices%Ip = 0d0
         imin = minloc([mesh%I(1, 1), mesh%I(2, 2), mesh%I(3, 3)], 1)
         imax = maxloc([mesh%I(1, 1), mesh%I(2, 2), mesh%I(3, 3)], 1)
         if (imin == 1) then
            if (imax == 2) then
               imid = 3
               matrices%P = reshape([1d0, 0d0, 0d0, 0d0, 0d0, 1d0, 0d0, 1d0, 0d0], [3, 3])
            end if
            if (imax == 3) then
               imid = 2
               matrices%P = eye(3)
            end if
         else if (imin == 2) then
            if (imax == 1) then
               imid = 3
               matrices%P = reshape([0d0, 1d0, 0d0, 0d0, 0d0, 1d0, 1d0, 0d0, 0d0], [3, 3])
            end if
            if (imax == 3) then
               imid = 1
               matrices%P = reshape([0d0, 1d0, 0d0, 1d0, 0d0, 0d0, 0d0, 0d0, 1d0], [3, 3])
            end if
         else
            if (imax == 1) then
               imid = 2
               matrices%P = reshape([0d0, 0d0, 1d0, 0d0, 1d0, 1d0, 1d0, 0d0, 0d0], [3, 3])
            end if
            if (imax == 2) then
               imid = 1
               matrices%P = reshape([0d0, 0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 1d0, 0d0], [3, 3])
            end if
         end if
         matrices%Ip(1) = minval([mesh%I(1, 1), mesh%I(2, 2), mesh%I(3, 3)])
         matrices%Ip(3) = maxval([mesh%I(1, 1), mesh%I(2, 2), mesh%I(3, 3)])
         matrices%Ip(2) = mesh%I(imid, imid)
      else
! P will be the rotation matrix between principal axes and laboratory
! axes, so that diag(Ip) = P'*I*P
         matrices%P = mesh%I
         call diasym(matrices%P, matrices%Ip)
      end if

! Probably the choice in DDSCAT: force eigenvectors to have mostly positive
! components
      do i = 1, 3
         negs = 0
         if (matrices%P(1, i) < 0d0) negs = negs + 1
         if (matrices%P(2, i) < 0d0) negs = negs + 1
         if (matrices%P(3, i) < 0d0) negs = negs + 1
         if (negs >= 2) matrices%P(:, i) = -matrices%P(:, i)
      end do

      matrices%I = 0d0
      matrices%I_inv = 0d0

      forall (i=1:3) matrices%I(i, i) = matrices%Ip(i)
      forall (i=1:3) matrices%I_inv(i, i) = 1d0/matrices%Ip(i)
      
      mesh%alpha = matrices%Ip/(2d0/5d0*mesh%rho*mesh%V*mesh%a**2)

   end subroutine diagonalize_inertia

!****************************************************************************80
! Mass parameters for VIE-mesh
! output: mesh%V = volume of the mesh
!         mesh%mass = mass of the mesh
!         mesh%CM = coordinates of the center of mass
!         mesh%I = The moment of inertia tensor
!****************************************************************************80
   subroutine vie_params(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      real(dp) :: rho, V, totV, mass, totMass, detJ, &
                  a, b, c, ap, bp, cp
      integer  :: i1
      real(dp), dimension(3) :: p0, p1, p2, p3, COM, CM
      real(dp), dimension(4) :: x, y, z
      real(dp), dimension(3, 3) :: I, E

      rho = mesh%rho ! currently density not in tetrahedral element data

      I = 0.d0
      totV = 0.d0
      totMass = 0.d0
      CM = 0.d0

      do i1 = 1, mesh%N_tet
         ! Get vertices
         p0 = mesh%coord(:, mesh%etopol(1, i1))
         p1 = mesh%coord(:, mesh%etopol(2, i1))
         p2 = mesh%coord(:, mesh%etopol(3, i1))
         p3 = mesh%coord(:, mesh%etopol(4, i1))
         COM = [(p0(1) + p1(1) + p2(1) + p3(1))/4d0, &
                (p0(2) + p1(2) + p2(2) + p3(2))/4d0, &
                (p0(3) + p1(3) + p2(3) + p3(3))/4d0]

         V = abs(dot_product((p0 - p3), &
                             crossRR((p1 - p3), (p2 - p3))))/6d0

         x = [p0(1) - COM(1), p1(1) - COM(1), &
              p2(1) - COM(1), p3(1) - COM(1)]
         y = [p0(2) - COM(2), p1(2) - COM(2), &
              p2(2) - COM(2), p3(2) - COM(2)]
         z = [p0(3) - COM(3), p1(3) - COM(3), &
              p2(3) - COM(3), p3(3) - COM(3)]

         detJ = (x(2) - x(1))*((y(3) - y(1))*(z(4) &
                                              - z(1)) - (y(4) - y(1))*(z(3) - z(1))) - &
                (x(3) - x(1))*((y(2) - y(1))*(z(4) &
                                              - z(1)) - (y(4) - y(1))*(z(2) - z(1))) + &
                (x(4) - x(1))*((y(2) - y(1))*(z(3) &
                                              - z(1)) - (y(3) - y(1))*(z(2) - z(1)))

         a = rho*detJ*diagterm(y, z)/60d0
         b = rho*detJ*diagterm(x, z)/60d0
         c = rho*detJ*diagterm(x, y)/60d0

         ap = rho*detJ*offterm(y, z)/120d0
         bp = rho*detJ*offterm(x, z)/120d0
         cp = rho*detJ*offterm(x, y)/120d0

         E = reshape([a, -bp, -cp, -bp, b, -ap, -cp, -ap, c], &
                     [3, 3])
         mass = rho*V

         totV = totV + V
         totMass = totMass + mass
         CM = CM + mass*COM
         I = I + E + mass*(dot_product(COM, COM)*eye(3) - &
                           real_outer_product(COM, COM))

      end do

      mesh%V = totV
      mesh%mass = totMass
      mesh%CM = CM/totMass
      mesh%I = I

      mesh%coord(1, :) = mesh%coord(1, :) - mesh%CM(1)
      mesh%coord(2, :) = mesh%coord(2, :) - mesh%CM(2)
      mesh%coord(3, :) = mesh%coord(3, :) - mesh%CM(3)
      matrices%CM = mesh%CM ! Save real CM if it happens to be not near origin
      mesh%CM = dble([0d0, 0d0, 0d0])
      matrices%x_CM = mesh%CM

   end subroutine vie_params

!****************************************************************************80

   function diagterm(y, z) result(a)

      real(dp)                 :: a

      real(dp), dimension(4) :: y, &
                                z

      a = (y(1)*y(1) &
           + y(1)*y(2) + y(2)*y(2) &
           + y(1)*y(3) + y(2)*y(3) + y(3)*y(3) &
           + y(1)*y(4) + y(2)*y(4) + y(3)*y(4) + y(4)*y(4) &
           + z(1)*z(1) &
           + z(1)*z(2) + z(2)*z(2) &
           + z(1)*z(3) + z(2)*z(3) + z(3)*z(3) &
           + z(1)*z(4) + z(2)*z(4) + z(3)*z(4) + z(4)*z(4))

   end function diagterm

!****************************************************************************80

   function offterm(y, z) result(ap)
      real(dp)                 :: ap
      real(dp), dimension(4) :: y, z

      ap = (2*y(1)*z(1) &
            + y(2)*z(1) &
            + y(3)*z(1) &
            + y(4)*z(1) &
            + y(1)*z(2) &
            + 2*y(2)*z(2) &
            + y(3)*z(2) &
            + y(4)*z(2) &
            + y(1)*z(3) &
            + y(2)*z(3) &
            + 2*y(3)*z(3) &
            + y(4)*z(3) &
            + y(1)*z(4) &
            + y(2)*z(4) &
            + y(3)*z(4) &
            + 2*y(4)*z(4))

   end function offterm

end module setup
