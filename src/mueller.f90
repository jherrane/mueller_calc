module mueller
   use common
   use T_matrix
   use gaussquad
   implicit none

contains

!****************************************************************************80

   subroutine scattering_extinction_matrices(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, j, ind_ave, ind2_ave, ind_al, ind2_al, ind3, ii, N_size, N_incangle, halton_init
      real(dp) :: K(4, 4), Cext_al, Csca_al, Cext_ave, Csca_ave
      complex(dp) :: Kevals(4), Kevecs(4, 4)
      real(dp), dimension(:, :), allocatable :: S_ave, S_al, K_ave, K_al, &
      SS_ave, KK_ave, SS_al, KK_al
      real(dp), dimension(N_ia) :: inc_angles
      CHARACTER(LEN=80) :: mueller_al, extinction_al, mueller_av, extinction_av

      mueller_al = 'out/M'//trim(matrices%out)//'-al'
      extinction_al = 'out/E'//trim(matrices%out)//'-al'
      mueller_av = 'out/M'//trim(matrices%out)//'-ave'
      extinction_av = 'out/E'//trim(matrices%out)//'-ave'

      call linspace(ia_range(1), ia_range(2), N_ia, inc_angles)

      allocate (SS_ave(N_points*size(mesh%ki, 1), 20))
      allocate (SS_al(N_points*size(mesh%ki, 1)*(size(inc_angles)), 20))
      allocate (KK_ave(size(mesh%ki, 1), 20))
      allocate (KK_al(size(mesh%ki, 1)*(size(inc_angles)), 20))
      SS_ave = 0d0
      KK_ave = 0d0
      SS_al = 0d0
      KK_ave = 0d0
      allocate (S_al(N_points, 18), K_al(1, 18))

      ind_ave = 0
      ind2_ave = 0
      ind_al = 0
      ind2_al = 0
      ind3 = 0

      do N_size = size(mesh%ki, 1), 1, -1
! Choose wavelength
         mesh%k = mesh%ki(N_size)
         
         do N_incangle = 0, size(inc_angles, 1)      
            if(N_incangle==0) then
               write(*, '(2(A,I0))') '  Computing average Mueller matrix ', 1+size(mesh%ki,1)-N_size, &
               '/', size(mesh%ki,1)
               call mueller_ave(matrices, mesh, N_points, points, N_size, S_ave, &
                  K_ave, Csca_ave, Cext_ave)
               do i = 1,N_points
                  ind_ave = ind_ave + 1
                  SS_ave(ind_ave, 1) = 1+size(mesh%ki,1)-N_size
                  SS_ave(ind_ave, 2) = 0
                  SS_ave(ind_ave, 3) = i
                  SS_ave(ind_ave, 4:19) = S_ave(i, 3:18)
                  SS_ave(ind_ave, 20) = Csca_ave
               end do
               ind2_ave = ind2_ave + 1
               KK_ave(ind2_ave, 1) = 1+size(mesh%ki,1)-N_size
               KK_ave(ind2_ave, 2) = 0
               KK_ave(ind2_ave, 3) = 1
               KK_ave(ind2_ave, 4:19) = K_ave(1, 3:18)
               KK_ave(ind2_ave, 20) = Cext_ave
            else
               ind3 = ind3 + 1
               write(*, '(2(A,I0))') '  Computing aligned Mueller matrix ', &
               ind3, '/', size(mesh%ki, 1)*size(inc_angles, 1)
               if(matrices%xi_in>0d0) then
                  call mueller_align(matrices, mesh, N_points, points, N_size, inc_angles(N_incangle), &
                  S_al, K_al, Csca_al, Cext_al,matrices%xi_in)  
               else
                  call mueller_align(matrices, mesh, N_points, points, N_size, inc_angles(N_incangle), &
                  S_al, K_al, Csca_al, Cext_al)  
               end if
               do i = 1, N_points
                  ind_al = ind_al + 1
                  SS_al(ind_al, 1) = 1+size(mesh%ki,1)-N_size
                  SS_al(ind_al, 2) = N_incangle
                  SS_al(ind_al, 3) = i
                  SS_al(ind_al, 4:19) = S_al(i, 3:18)
                  SS_al(ind_al, 20) = Csca_al
               end do
               ind2_al = ind2_al + 1
               KK_al(ind2_al, 1) = 1+size(mesh%ki,1)-N_size
               KK_al(ind2_al, 2) = N_incangle
               KK_al(ind2_al, 3) = 1
               KK_al(ind2_al, 4:19) = K_al(1, 3:18)
               KK_al(ind2_al, 20) = Cext_al
            end if
         end do
      end do
      call write_RT_matrix(SS_al, mueller_al, 1)
      call write_RT_matrix(KK_al, extinction_al, 2)
      call write_RT_matrix(SS_ave, mueller_av, 1)
      call write_RT_matrix(KK_ave, extinction_av, 2)

   end subroutine scattering_extinction_matrices

!****************************************************************************80

   subroutine mueller_ave(matrices, mesh, N_points, points, ii, SS, KK, Csca_out, Cext_out)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, j, ii, N_avgs, halton_init, nm, Nmax, N_points
      real(dp) :: E, vec(3), k_sph(3)
      real(dp) :: Cext, Cabs, Csca, Csca_out, Cext_out
      real(dp), dimension(:, :), allocatable :: S, SS, K, KK, points
      complex(dp), dimension(:), allocatable :: a_in, b_in
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      N_avgs = 720 ! Number of averaging directions
      if(.not. allocated(SS)) allocate(SS(N_points,18), KK(1,18))
      SS = 0d0
      KK = 0d0
      Csca = 0d0 
      Cext = 0d0 

      matrices%R = eye(3)
      halton_init = 0
      E = matrices%E
      Nmax = matrices%Nmaxs(ii)

      do i = 1, N_avgs
         vec(1) = 1d0
         vec(2) = acos(2d0*halton_seq(halton_init + i, 2) - 1d0)
         vec(3) = halton_seq(halton_init + i, 3)*2d0*pi
         matrices%khat = [dsin(vec(2))*dcos(vec(3)), dsin(vec(2))*dsin(vec(3)), dcos(vec(2))]

         matrices%khat = -matrices%khat/vlen(matrices%khat)
         matrices%R = rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0])

         k_sph = cart2sph(matrices%khat)

         call rot_setup(matrices, mesh)
         call incident_fields(matrices, mesh, E, a_in, b_in, ii)
         call scattered_fields(matrices, mesh, E, p, q, p90, q90, ii)

         if (allocated(S)) deallocate (S, K)

         call scattering_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), N_points, points, S)
         call extinction_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), k_sph(2), k_sph(3), K)
         call cross_sections(p, q, a_in, b_in, dcmplx(mesh%k), Nmax, Cext, Csca, Cabs)

         SS = SS + S/N_avgs
         KK = KK + K/N_avgs
         Csca_out = Csca_out + Csca/N_avgs
         Cext_out = Cext_out + Cext/N_avgs
         call print_bar(i,N_avgs)
      end do

   end subroutine mueller_ave

!****************************************************************************80

   subroutine mueller_align(matrices, mesh, N_points, points, ii, psi, &
      SS, KK, Csca_out, Cext_out, xi_in)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, j, ii, ind, N_points, Nmax, NB, Nxi
      real(dp) :: E, vec(3), phi, omega(3), k_sph(3), &
      aproj(3), theta, xi, B(3), phi_B, psi
      real(dp) :: Cext, Cabs, Csca, Csca_out, Cext_out
      real(dp), dimension(3,3) :: RR, Rxi, Qt, R0
      real(dp), optional :: xi_in
      real(dp), dimension(:, :), allocatable :: S, SS, K, KK, points
      complex(dp), dimension(:), allocatable :: a_in, b_in
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      NB = 36 ! Number of points to calculate the perfect orientations

      SS = 0d0
      KK = 0d0
      Csca = 0d0 
      Cext = 0d0 

      B = [sin(psi), 0d0, cos(psi)]
      Nxi = 20 ! Number of precession averaging directions
      if(.NOT. present(xi_in)) then
         xi = 0d0
         Nxi = 1
      else
         xi = xi_in
      end if 

      Rxi = R_aa([0d0,1d0,0d0], xi)

! Rotation axis for precession averaging
      B = [sin(psi), 0d0, cos(psi)]
      E = matrices%E

      Nmax = matrices%Nmaxs(ii)
      ind = 0
      do j = 1, Nxi
         phi_B = dble(j - 1)*2d0*pi/dble(Nxi)
         R0 = matmul(Rxi,rotate_a_to_b(matrices%P(:, 3), B))
         R0 = matmul(R_aa(B,phi_B),R0)
         omega = matmul(R0, matrices%P(:,3))
         Qt = matmul(R0, matrices%P)
         aproj = [Qt(1, 3), Qt(2, 3), 0d0]
         aproj = aproj/vlen(aproj)
         phi = dacos(aproj(1))
         do i = 1, NB
            ind = ind + 1
            theta = dble(i - 1)*2d0*pi/NB
            RR = matmul(transpose(R_aa(omega, theta)), transpose(R0))

            matrices%khat = -matmul(RR, [0d0, 0d0, 1d0])
            matrices%khat = matrices%khat/vlen(matrices%khat)
            k_sph = cart2sph(matrices%khat)

            matrices%R = rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0])
            matrices%R = matmul(matrices%R,R_aa(matrices%khat, phi))

            call rot_setup(matrices, mesh)
            call incident_fields(matrices, mesh, E, a_in, b_in, ii)
            call scattered_fields(matrices, mesh, E, p, q, p90, q90, ii)

            if (allocated(S)) deallocate (S, K)

            call scattering_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), N_points, points, S)
            call extinction_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), k_sph(2), k_sph(3), K)
            call cross_sections(p, q, a_in, b_in, dcmplx(mesh%k), Nmax, Cext, Csca, Cabs)

            SS = SS + S/NB/Nxi
            KK = KK + K/NB/Nxi
            Csca_out = Csca_out + Csca/NB/Nxi
            Cext_out = Cext_out + Cext/NB/Nxi
            call print_bar(ind,Nxi*NB)
         end do 
      end do

   end subroutine mueller_align

end module mueller
