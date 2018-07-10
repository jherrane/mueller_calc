module mueller
   use T_matrix
   use gaussquad
   implicit none

contains

!****************************************************************************80

   subroutine scattering_extinction_matrices(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, j, ind, ind2, ind3, ii, N_size, N_ia, halton_init
      real(dp) :: inc_angles(2), K(4, 4), Cext_al, Csca_al, Cext_ave, Csca_ave
      complex(dp) :: Kevals(4), Kevecs(4, 4)
      real(dp), dimension(:, :), allocatable :: S_ave, S_al, SS, K_ave, K_al, KK, points
      CHARACTER(LEN=80) :: mueller_out, extinction_out

      mueller_out = 'out/M'//trim(matrices%out)
      extinction_out = 'out/E'//trim(matrices%out)

      inc_angles = [90d0, 180d0]

      allocate (SS(Ntheta*Nphi*size(mesh%ki, 1)*(1+size(inc_angles)), 20))
      allocate (KK(size(mesh%ki, 1)*(1+size(inc_angles)), 20))
      SS = 0d0
      KK = 0d0
      allocate (S_al(Ntheta*Nphi, 18), K_al(1, 18))

      ind = 0
      ind2 = 0
      ind3 = 0
      do N_size = size(mesh%ki, 1), 1, -1
! Choose wavelength
         mesh%k = mesh%ki(N_size)
         
         do N_ia = 0, size(inc_angles, 1)      
            if(N_ia==0) then
               write(*, '(2(A,I0))') '  Computing average Mueller matrix ', 1+size(mesh%ki,1)-N_size, &
               '/', size(mesh%ki,1)
               call mueller_ave(matrices, mesh, Ntheta, 1, N_size, S_ave, &
                  K_ave, Csca_ave, Cext_ave)
               do i = 1,Ntheta
                  ind = ind + 1
                  SS(ind, 1) = 1+size(mesh%ki,1)-N_size
                  SS(ind, 2) = 0
                  SS(ind, 3) = i
                  SS(ind, 4:19) = S_ave(i, 3:18)
                  SS(ind, 20) = Csca_ave
               end do
               ind2 = ind2 + 1
               KK(ind2, 1) = 1+size(mesh%ki,1)-N_size
               KK(ind2, 2) = N_ia
               KK(ind2, 3) = 1
               KK(ind2, 4:19) = K_ave(1, 3:18)
               KK(ind2, 20) = Cext_ave
            else
               ind3 = ind3 + 1
               write(*, '(2(A,I0))') '  Computing aligned Mueller matrix ', &
               ind3, '/', size(mesh%ki, 1)*size(inc_angles, 1)
               call mueller_align(matrices, mesh, Ntheta, Nphi, N_size, inc_angles(N_ia), &
                  S_al, K_al, Csca_al, Cext_al)  
               do i = 1, Ntheta*Nphi
                  ind = ind + 1
                  SS(ind, 1) = 1+size(mesh%ki,1)-N_size
                  SS(ind, 2) = N_ia
                  SS(ind, 3) = i
                  SS(ind, 4:19) = S_al(i, 3:18)
                  SS(ind, 20) = Csca_al
               end do
               ind2 = ind2 + 1
               KK(ind2, 1) = 1+size(mesh%ki,1)-N_size
               KK(ind2, 2) = N_ia
               KK(ind2, 3) = 1
               KK(ind2, 4:19) = K_al(1, 3:18)
               KK(ind2, 20) = Cext_al
            end if
         end do
      end do
      call write_RT_matrix(SS, mueller_out, 1)
      call write_RT_matrix(KK, extinction_out, 2)

   end subroutine scattering_extinction_matrices

!****************************************************************************80

   subroutine mueller_ave(matrices, mesh, N_theta, N_phi, ii, SS, KK, Csca_out, Cext_out)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, ii, N_avgs, halton_init, nm, Nmax, N_theta, N_phi
      real(dp) :: E, vec(3), k_sph(3)
      real(dp) :: Cext, Cabs, Csca, Csca_out, Cext_out
      real(dp), dimension(:, :), allocatable :: S, SS, K, KK, points
      complex(dp), dimension(:), allocatable :: a_in, b_in
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      N_avgs = 720 ! Number of averaging directions
      if(.not. allocated(SS)) allocate(SS(N_theta*N_phi,18), KK(1,18))
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

         call mueller_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), N_theta, N_phi, S)
         call extinction_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), k_sph(2), k_sph(3), K)
         call cross_sections(p, q, a_in, b_in, dcmplx(mesh%k), Nmax, Cext, Csca, Cabs)
         
         SS = SS + S/N_avgs
         KK = KK + K/N_avgs
         Csca_out = Csca_out + Csca/N_avgs
         Cext_out = Cext_out + Cext/N_avgs
      end do

   end subroutine mueller_ave

!****************************************************************************80

   subroutine mueller_align(matrices, mesh, N_theta, N_phi, ii, psi, &
      SS, KK, Csca_out, Cext_out, xi_in)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, j, ii, ind, N_theta, N_phi, Nmax, NB, Nxi
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

            call mueller_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), N_theta, N_phi, S)
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
