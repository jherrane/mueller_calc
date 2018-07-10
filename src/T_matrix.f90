module T_matrix
   use transformation_matrices
   use precorrection
   use gmres_module
   use setup
   use common

   implicit none

contains

!****************************************************************************80
! The main calculation routine for the T-matrices. Uses the celebrated
! JVIE-methodology.
   subroutine calc_T(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, ii, nm, sz

      sz = size(mesh%ki)
      if (matrices%singleT == 1) sz = 1

      write (*, '(3(A,I0))') ' Construct matrices for ', sz, ' wavelengths...'
      if (matrices%singleT == 1) then
         write (*, '(2(A, F8.3))') ' Wavelength now is ', &
            2d0*pi/mesh%ki(matrices%whichbar)/1d-6, ' um.'
      end if
      do i = 1, sz
         if (size(mesh%params, 2) > 1) mesh%param = mesh%params(:, matrices%whichbar)
         ii = i
         if (matrices%singleT == 1) then
            ii = matrices%whichbar
            write (*, '(3(A,I0))') ' Step ', ii
         else
            write (*, '(3(A,I0))') ' Step ', ii, '/', sz, ''
         end if
         nm = (matrices%Nmaxs(ii) + 1)**2 - 1
         mesh%k = mesh%ki(ii)

         call allocate_T(matrices, mesh, ii)

         call update_projections(matrices, mesh, ii)

         if (allocated(matrices%Fg)) deallocate (matrices%Fg)
         if (allocated(matrices%sp_mat)) deallocate (matrices%sp_mat, matrices%sp_ind)
         call build_G(matrices, mesh)

         print *, ' Constructing the sparse part of the system matrix'
         if (mesh%order == 0) then
            call compute_near_zone_interactions_const(matrices, mesh)
         end if
         if (mesh%order == 1) then
            call compute_near_zone_interactions_lin(matrices, mesh)
         end if

         print *, ' Compute T-matrix...'
         call compute_T_matrix(matrices, mesh, matrices%Nmaxs(ii), matrices%Taa, &
                               matrices%Tab, matrices%Tba, matrices%Tbb)

         matrices%Taai(1:nm, 1:nm, ii) = matrices%Taa
         matrices%Tabi(1:nm, 1:nm, ii) = matrices%Tab
         matrices%Tbai(1:nm, 1:nm, ii) = matrices%Tba
         matrices%Tbbi(1:nm, 1:nm, ii) = matrices%Tbb
      end do

   end subroutine calc_T

!****************************************************************************80
! Allocate space for the currently used T-matrix, thus deallocation if
! necessary.
   subroutine allocate_T(matrices, mesh, i)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: Nmax, i

      if (allocated(matrices%Taa)) then
         deallocate (matrices%Taa, matrices%Tba, matrices%Tab, matrices%Tbb)
      end if

      Nmax = matrices%Nmaxs(i)
      allocate (matrices%Taa((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1))
      allocate (matrices%Tab((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1))
      allocate (matrices%Tba((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1))
      allocate (matrices%Tbb((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1))

   end subroutine allocate_T

!****************************************************************************80
! Allocate memory for the collection of T-matrices. Wastes space as the
! largest T-matrix usually is much larger than the others. Must be allocated
! before anything else is done with T-matrices.
   subroutine allocate_Ti(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: Nmax, i, ind1, ind2, nm
      
      do i = 1,size(matrices%Nmaxs)
         nm = (matrices%Nmaxs(i)+1)**2-1
         T_size = T_size + nm**2
      end do
      
      Nmax = maxval(matrices%Nmaxs)

      allocate (matrices%Taai((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1, matrices%bars))
      allocate (matrices%Tabi((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1, matrices%bars))
      allocate (matrices%Tbai((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1, matrices%bars))
      allocate (matrices%Tbbi((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1, matrices%bars))

      matrices%Taai(:, :, :) = dcmplx(0.0, 0.0)
      matrices%Tabi(:, :, :) = dcmplx(0.0, 0.0)
      matrices%Tbai(:, :, :) = dcmplx(0.0, 0.0)
      matrices%Tbbi(:, :, :) = dcmplx(0.0, 0.0)

   end subroutine allocate_Ti

!****************************************************************************80
! Śolve the T-matrix via the VIE method and accelerated GMRES.
   subroutine compute_T_matrix(matrices, mesh, Nmax, Taa, Tab, Tba, Tbb)
      type(data) :: matrices
      type(mesh_struct) :: mesh      
      real(dp) :: k
      integer :: Nmax, nm

      complex(dp) :: mat(3*mesh%N_tet, 2*((Nmax + 1)**2 - 1))
      complex(dp) :: T_mat(2*((Nmax + 1)**2 - 1), 2*((Nmax + 1)**2 - 1))
      complex(dp), dimension((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1) :: Taa, Tab, Tba, Tbb
      complex(dp) :: sc

      k = mesh%k

! Compute transformations
      call vswf2constant(mesh, dcmplx(k), Nmax, mat)

      do nm = 1, size(mat, 2)
         matrices%rhs = mat(:, nm)
! Redirect terminal output to trash unless debug mode is on. Maybe
! not very portable...
         if (debug == 0) open (unit=6, file="/dev/null", form="formatted")
         call gmres(matrices, mesh)
         T_mat(:, nm) = matmul(transpose(conjg(mat)), matrices%x)
         if (debug == 0) then
            open (unit=6, file="/dev/stdout", form="formatted")
         else
            print *, nm, '/', size(mat, 2)
         end if
         call print_bar(nm, size(mat, 2))

      end do

      nm = (Nmax + 1)**2 - 1
      sc = dcmplx(0.0d0, k**3.0d0)
      Taa = T_mat(1:nm, 1:nm)*sc
      Tab = T_mat(1:nm, nm + 1:2*nm)*sc
      Tba = T_mat(nm + 1:2*nm, 1:nm)*sc
      Tbb = T_mat(nm + 1:2*nm, nm + 1:2*nm)*sc

   end subroutine compute_T_matrix

!****************************************************************************80
! Return the precalculated fields a and b, given electric field amplitude E and
! the current wavelength index ii
   subroutine incident_fields(matrices, mesh, E, a_in, b_in, ii)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      real(dp) :: E
      complex(dp), dimension(:), allocatable :: a_in, b_in
      integer :: nm, las, ii, Nmax

      Nmax = matrices%Nmaxs(ii)
      las = (Nmax + 1)*(2*Nmax + 1)*(2*Nmax + 3)/3 - 1
      nm = (Nmax + 1)**2 - 1

      if(allocated(a_in)) deallocate(a_in, b_in)
      allocate (a_in(nm), b_in(nm))

      a_in = E*matrices%as(1:nm, ii)
      b_in = E*matrices%bs(1:nm, ii)

   end subroutine incident_fields

!****************************************************************************80
! Solve the scattered fields p and q, given electric field amplitude E and
! the current wavelength index ii
   subroutine scattered_fields(matrices, mesh, E, p, q, p90, q90, ii)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      real(dp) :: E
      complex(dp), dimension(:), allocatable :: a_in, b_in, a90, b90, &
                                                a, b, p, q, p90, q90, ptemp, qtemp, p90temp, q90temp
      complex(dp), dimension(:, :), allocatable :: Taa, Tab, Tba, Tbb
      complex(8), dimension(:), allocatable :: rotD, rotD90, rbak
      integer, dimension(:, :), allocatable :: indD, indD90, ibak
      integer :: nm, las, ii, Nmax

      Nmax = matrices%Nmaxs(ii)
      las = (Nmax + 1)*(2*Nmax + 1)*(2*Nmax + 3)/3 - 1
      nm = (Nmax + 1)**2 - 1

      if (.not. allocated(p)) allocate (p(nm), q(nm), p90(nm), q90(nm))

      allocate (rbak(las), ibak(las, 2))
      allocate (a(nm), b(nm), a90(nm), b90(nm), ptemp(nm), qtemp(nm), p90temp(nm), q90temp(nm))
      allocate (rotD(las), indD(las, 2), rotD90(las), indD90(las, 2))
      allocate (Taa(nm, nm), Tab(nm, nm), Tba(nm, nm), Tbb(nm, nm))

      rotD = matrices%rotDs(1:las, ii)
      indD = matrices%indDs(1:las, :, ii)
      rotD90 = matrices%rotD90s(1:las, ii)
      indD90 = matrices%indD90s(1:las, :, ii)

      Taa = matrices%Taai(1:nm, 1:nm, ii)
      Tab = matrices%Tabi(1:nm, 1:nm, ii)
      Tba = matrices%Tbai(1:nm, 1:nm, ii)
      Tbb = matrices%Tbbi(1:nm, 1:nm, ii)

      a_in = E*matrices%as(1:nm, ii)
      b_in = E*matrices%bs(1:nm, ii)

      a = sparse_matmul(rotD, indD, a_in, nm)
      b = sparse_matmul(rotD, indD, b_in, nm)
      a90 = sparse_matmul(rotD90, indD90, a_in, nm)
      b90 = sparse_matmul(rotD90, indD90, b_in, nm)

      ptemp = matmul(Taa, a) + matmul(Tab, b)
      qtemp = matmul(Tbb, b) + matmul(Tba, a)
      p90temp = matmul(Taa, a90) + matmul(Tab, b90)
      q90temp = matmul(Tbb, b90) + matmul(Tba, a90)

      call sph_rotation_sparse_gen(mat2euler(matrices%Rk), Nmax, rbak, ibak)

      p = sparse_matmul(rbak, ibak, ptemp, nm)
      q = sparse_matmul(rbak, ibak, qtemp, nm)
      p90 = sparse_matmul(rbak, ibak, p90temp, nm)
      q90 = sparse_matmul(rbak, ibak, q90temp, nm)

   end subroutine scattered_fields

end module T_matrix
