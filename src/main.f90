program main
   use common
   use io
   use mueller
   use setup

   implicit none
   type(data) :: matrices
   type(mesh_struct) :: mesh

   call splash('v1.0')

   call check_paramsfile(matrices)
   call read_params(matrices, mesh)
   call read_arguments(matrices, mesh)
   call wavelength_band(matrices, mesh)

   call init_geometry(matrices, mesh)

   call allocate_Ti(matrices, mesh)

   if (matrices%Tmat == 1 .AND. file_exists(matrices%tname)) then
      call read_T(matrices, mesh)
      write (*, '(A, 20F6.3)') ' Wavelengths in um: ', 2d6*pi/mesh%ki
   else
      if (matrices%singleT == 1) then
         call T_empty(matrices, mesh)
      end if
      write (*, '(A, 20F6.3)') ' Wavelengths in um: ', 2d6*pi/mesh%ki
      call calc_T(matrices, mesh)
      call write_T(matrices, mesh)
   end if

   call polarization(matrices, mesh)
   call allocate_inc_wave(matrices, mesh)
   call vie_params(matrices, mesh)

   matrices%x_CM = mesh%CM
   call diagonalize_inertia(matrices, mesh)
   matrices%k_orig = matrices%khat
   matrices%E0_orig = real(matrices%E0hat)
   matrices%E90_orig = real(matrices%E90hat)

   call scattering_extinction_matrices(matrices, mesh)

end program main
