module io
   use common
   use h5io

   implicit none

contains
! MY ROUTINES
! CLEARTEXT READ
! CLEARTEXT WRITE

! MY ROUTINES *****************************************************************
!****************************************************************************80

   subroutine splash(v)
      character(4) :: v ! Version
      print *, '******************************************************************************'
      print *, '**                                                                          **'
      print *, '**                      Mueller matrix calculator ',v,'                     **'
      print *, '**                                                                          **'
      print *, '******************************************************************************'
      print *, ''
      call curr_time()

   end subroutine splash

!****************************************************************************80

   subroutine curr_time()
      character(8)  :: date, fmt, fmt2
      character(10) :: time
      character(5)  :: zone
      character(3)  :: hh, mm, ss, ms
      integer, dimension(8) :: values

      fmt = '(I2.2)'
      fmt2 = '(I3.3)'
      call date_and_time(date, time, zone, values)
      write (hh, fmt) values(5)
      write (mm, fmt) values(6)
      write (ss, fmt) values(7)
      write (ms, fmt2) values(8)
      print'(3(A, I0),8A)', ' Time: ', values(3), '.', values(2), '.', values(1), ' -- ', &
         trim(hh), ':', trim(mm), ':', trim(ss), '.', trim(ms)

   end subroutine curr_time

!****************************************************************************80

   subroutine print_bar(i1, Nmax)
      integer :: i1, Nmax, k
      real(dp) :: r
      character(len=1) :: bar, back

      back = char(8)
      bar = '='
      r = 100d0/Nmax
! print the percentage and the bar without line change, then reset cursor
      if (floor(i1*r) - floor((i1 - 1)*r) > 0) then
         write (6, '(2x,1i3,1a1,2x,1a1,256a1)', advance='no') &
            ceiling(r*i1), '%', '|', (bar, k=1, 50*i1/Nmax)
         write (6, '(256a1)', advance='no') (back, k=1, (50*i1/Nmax) + 9)
         if (i1 == Nmax) then
            write (*, *) ''
         end if
      end if

   end subroutine print_bar

!****************************************************************************80

   subroutine print_mat(mat, matname)
      real(dp), dimension(:, :) :: mat
      character(len=*) :: matname
      integer :: i, m

      write (*, '(2A)') trim(matname), ' = '
      m = size(mat, 1)
      do i = 1, m
         write (*, '(10F8.4)') mat(i, :)
      end do

   end subroutine print_mat

! CLEARTEXT READ **************************************************************
!****************************************************************************80
   subroutine check_paramsfile(matrices)
      type(data) :: matrices
! Subroutine to read input parameters
      integer ::  i
      character(len=80) :: arg

      do i = 1, command_argument_count(), 2
         call get_command_argument(i, arg)

         select case (arg)
         case ('-p', '--paramsfile')
            call get_command_argument(i + 1, matrices%paramsfile)
            write (*, '(2A)') ' Paramsfile: ', trim(matrices%paramsfile)
         end select
      end do
   end subroutine check_paramsfile

!****************************************************************************80
! Subroutine to read input parameters
   subroutine read_arguments(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i
      character(len=80) :: arg_name, arg

      matrices%singleT = 0

      do i = 1, command_argument_count(), 2
         call get_command_argument(i, arg_name)

         select case (arg_name)

         case ('-d', '--debug')
            call get_command_argument(i, arg)
            debug = 1
         case ('-m', '--mesh')
            call get_command_argument(i + 1, arg)
            mesh%meshname = arg
            write (*, '(2A)') ' Mesh: ', trim(mesh%meshname)
         case ('-T', '--Tmat')
            call get_command_argument(i + 1, arg)
            matrices%tname = arg
            write (*, '(2A)') ' T-matrix: ', trim(matrices%tname)
         case ('-o', '--out')
            call get_command_argument(i + 1, arg)
            matrices%out = arg
            write (*, '(2A)') ' Log: ', 'out/log'//trim(matrices%out)
         case ('--mueller_mode')
            call get_command_argument(i + 1, arg)
            matrices%mueller_mode = trim(arg)
            write (*, '(2A)') ' Mueller mode: ', trim(matrices%mueller_mode)
         case ('-p', '--paramsfile') ! Skip, as already read at check_paramsfile
         case ('-w', '--wavelen')
            call get_command_argument(i + 1, arg)
            read (arg, *) matrices%whichbar
         case ('-S', '--singleT')
            call get_command_argument(i + 1, arg)
            read (arg, *) matrices%singleT
         case ('-x', '--xi')
            call get_command_argument(i + 1, arg)
            read (arg, *) matrices%xi_in
         case ('-P', '--points')
            call get_command_argument(i + 1, arg)
            read (arg, *) points_file

         case ('-h', '--help')
            write (*, '(A)') ' Commands        Value       Description'
            write (*, '(A)') '---------------------------------------'
            write (*, '(A)') ' -d --debug                  Print more detailed info'
            write (*, '(A)') ' -m --mesh       mesh.h5     Mesh geometry'
            write (*, '(A)') ' -T --Tmat       T.h5        T-matrix file'
            write (*, '(A)') ' -o --out                    Output file identifier'
            write (*, '(A)') ' -p --paramsfile params.in   Read input parameters from file'
            write (*, '(A)') ' -w --wavelen    0           Choose wavelength from the T-matrix'
            write (*, '(A)') ' -S --singleT    0           Calculate only one T-matrix, of wb'
            write (*, '(A)') ' -x --xi         0.0         Precession angle about B'
            write (*, '(A)') ' -P --points     points      File for Mueller angular grid'
            stop
         case default
            print'(a,a,/)', 'Unrecognized command-line option: ', arg_name
            stop
         end select
      end do

      if (matrices%singleT == 1 .AND. matrices%whichbar == 0) then
         write (*, '(A)') ' Warning: Single T-matrix chosen but all bars used. Now whichbar = 1'
         matrices%whichbar = 1
      end if

   end subroutine read_arguments

!****************************************************************************80
! Subroutine to read a input file for nifty usage
   subroutine read_params(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
! Input related variables
      character(len=150) :: buffer, label
      real(dp) :: temp, tempii
      integer :: pos
      integer, parameter :: fh = 15
      integer :: ios = 0
      integer :: line = 0

! Control file variables
      real(dp) :: khat(3), lambda1, lambda2
      integer :: whichbar
      character(len=4) :: R0 = 'i'
      character(len=8) :: mode

      open (fh, file=matrices%paramsfile)

! ios is negative if an end of record condition is encountered or if
! an endfile condition was detected.  It is positive if an error was
! detected.  ios is zero otherwise.

      do while (ios == 0)
         read (fh, '(A)', iostat=ios) buffer
         if (ios == 0) then
            line = line + 1

            ! Find the first instance of whitespace.  Split label and data.
            pos = scan(buffer, '    ')
            label = buffer(1:pos)
            buffer = buffer(pos + 1:)

            select case (label)
            case ('a')
               read (buffer, *, iostat=ios) mesh%a
            case ('polarization')
               read (buffer, *, iostat=ios) matrices%polarization
            case ('bars')
               read (buffer, *, iostat=ios) matrices%bars
            case ('Tmat')
               read (buffer, *, iostat=ios) matrices%Tmat
            case ('lambda1')
               read (buffer, *, iostat=ios) lambda1
               matrices%lambda1 = lambda1*1.d-9
            case ('lambda2')
               read (buffer, *, iostat=ios) lambda2
               matrices%lambda2 = lambda2*1.d-9
            case ('xi_in')
               read (buffer, *, iostat=ios) matrices%xi_in
            case ('N_ia')
               read (buffer, *, iostat=ios) N_ia
            case ('ia_range')
               read (buffer, *, iostat=ios) ia_range
            case ('whichbar')
               read (buffer, *, iostat=ios) whichbar
               if (matrices%whichbar == 0) matrices%whichbar = whichbar
            case default
               !print *, 'Skipping invalid label at line', line

            end select
         end if
      end do

      close (fh)
      matrices%R90_init = reshape(dble([0d0, 1d0, 0d0, -1d0, 0d0, 0d0, 0d0, 0d0, 1d0]), [3, 3])
! Convert angles to radians
      ia_range = ia_range*pi/180d0
      allocate (mesh%ki(matrices%bars))
      allocate (matrices%E_rel(matrices%bars))
      allocate (matrices%Nmaxs(matrices%bars))

   end subroutine read_params

!****************************************************************************80

   subroutine read_points(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      real(dp), dimension(:,:), allocatable :: input

      N_points = get_last_line_no(points_file)
      allocate(points(2,N_points), input(3,N_points))
      open(1, file = points_file)
      read(1,*) input
      close(1)
      points = input(2:3, :)
      N_theta = calc_uniques(points)

   end subroutine read_points

!****************************************************************************80

   function get_last_line_no(fle) result(lineno)
      character(len=80) :: fle
      integer :: lineno
      integer :: ios = 0
      integer :: nbline = 0

      open (11, file=trim(fle))
      do while (.true.)
         read (11, *, iostat=ios) ! ios should have been declared as an integer
         if (ios > 0) then
            stop'problem somewhere (Called from io:get_last_line_no)'
         else if (ios < 0) then ! end of file is reached
            exit
         else
            nbline = nbline + 1
         end if
      end do
      close (11)
      lineno = nbline

   end function get_last_line_no

!****************************************************************************80

   subroutine read_log(matrices, mesh, no)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer, parameter :: fh = 15
      integer :: line

! Control file variables
      real(dp) :: wlast(3), t_tresh, t1, t2, tf, w_av(3), Rab(3, 3)
      integer :: firstlineno, lastlineno, no, numlines, i, n1, n2
      real(dp) :: x(3), v(3), w(3), J(3), F(3), N(3), t, R(9)

      numlines = no
      lastlineno = get_last_line_no(matrices%out)
      if(numlines>=lastlineno-23) numlines = lastlineno-23
      firstlineno = lastlineno-numlines
      open (fh, file='out/log'//trim(matrices%out))
      do i = 1, lastlineno-1
         read (fh, *)
      end do

      read (fh, *) line, x, v, w, J, N, F, t, R

      wlast = w/dsqrt(w(1)**2 + w(2)**2 + w(3)**2)

      close (fh)

      t_tresh = 4d0*pi/dsqrt(w(1)**2 + w(2)**2 + w(3)**2)
      tf = t
      t1 = t - t_tresh
      t2 = t - 10*t_tresh

      w_av = 0d0
      open (fh, file=trim(matrices%out))
      do i = 1, lastlineno
         if (i >= firstlineno) then
            read (fh, *) line, x, v, w, J, N, F, t, R
            w_av = w_av + w
            if (t > t2) then
               ! print*, line
               t2 = tf
            end if
            if (t > t1) then
               ! print*, line
               t1 = tf
            end if
            if (i == firstlineno) then
               ! print*, line
               n1 = line
            end if
         else
            read (fh, *)
         end if
      end do
      n2 = line
      w_av = w_av/(lastlineno - firstlineno + 1)
      ! print*, 'w_av = ', w_av
      close (fh)
      Rab = rotate_a_to_b(w_av, [1d0, 0d0, w_av(3)])
      ! call print_mat(Rab,'R_ab')
      matrices%R_al = Rab
      allocate (matrices%RRR(3, 3, n2 - n1 + 1))
      allocate (matrices%www(3, n2 - n1 + 1))

      open (fh, file=trim(matrices%out))
      do i = 1, lastlineno
         if (i > firstlineno) then
            read (fh, *) line, x, v, w, J, N, F, t, R
            matrices%RRR(:, :, i - firstlineno) = reshape(R, [3, 3])
            matrices%www(:, i - firstlineno) = w
         else
            read (fh, *)
         end if
      end do
      close (fh)

   end subroutine read_log

!****************************************************************************80

   subroutine read_mueller(A, fname)
      real(dp), allocatable, intent(out) :: A(:, :)
      character(len=80), intent(in) :: fname
      integer, dimension(2) :: dims ! Dataset dimensions

      integer     :: i

      dims = (/181, 18/)
      allocate (A(dims(1) - 1, dims(2)))

      open (unit=1, file=trim(fname))
      read (1, *)
      do i = 2, dims(1)
         read (1, '(2F8.4,16ES14.4)') A(i, :)
      end do

      close (1)

   end subroutine read_mueller

! CLEARTEXT WRITE *************************************************************
!****************************************************************************80

   subroutine write_array(A, fname)
      real(dp), intent(in) :: A(:, :)
      character(len=8), intent(in) :: fname
      integer(HSIZE_T), dimension(2) :: dims ! Dataset dimensions
      integer     :: i

      dims = [int(size(A, 1), 8), int(size(A, 2), 8)]

      open (unit=1, file=trim(fname), ACTION="write", STATUS="replace")

      do i = 1, int(dims(1))
         write (1, '(1024ES14.6)') A(i, :)
      end do

      close (1)

   end subroutine write_array

!****************************************************************************80

   subroutine write_mueller(matrices, A)
      type(data) :: matrices
      real(dp), intent(in) :: A(:, :)
      character(len=120) :: fname
      integer(HSIZE_T), dimension(2) :: dims ! Dataset dimensions
      integer     :: i

      fname = 'out/mueller-'//trim(matrices%mueller_mode)//trim(matrices%out)
      dims = [int(size(A, 1), 8), int(size(A, 2), 8)]

      open (unit=1, file=trim(fname), ACTION="write", STATUS="replace")
      write (1, '(18A)') '   phi ', '    theta ', '      S11  ', '         S12  ', &
         '         S13  ', '         S14  ', '         S21  ', '         S22  ', '         S23  ', &
         '         S24  ', '         S31  ', '         S32  ', '         S33  ', '         S34  ', &
         '         S41  ', '         S42  ', '         S43  ', '         S44  '

      do i = 1, int(dims(1))
         write (1, '(2F8.4,16ES14.6)') A(i, :)
      end do

      close (1)

   end subroutine write_mueller

!****************************************************************************80

   subroutine write_RT_matrix(A, fname, type)
      real(dp), intent(in) :: A(:, :)
      character(len=80), intent(in) :: fname
      integer(HSIZE_T), dimension(2) :: dims ! Dataset dimensions

      integer     :: i, type

      dims = [int(size(A, 1), 8), int(size(A, 2), 8)]
      open (unit=1, file=trim(fname), ACTION="write", STATUS="replace")

      if (type == 1) then
         write (1, '(20A)') ' N_size ', '     N_ia ', '    N_pts ', '      S11  ', '         S12  ', &
            '         S13  ', '         S14  ', '         S21  ', '         S22  ', '         S23  ', &
            '         S24  ', '         S31  ', '         S32  ', '         S33  ', '         S34  ', &
            '         S41  ', '         S42  ', '         S43  ', '         S44  ', '        Csca  '
      else if (type == 2) then
         write (1, '(20A)') ' N_size ', '     N_ia ', '    N_pts ', '      K11  ', '         K12  ', &
            '         K13  ', '         K14  ', '         K21  ', '         K22  ', '         K23  ', &
            '         K24  ', '         K31  ', '         K32  ', '         K33  ', '         K34  ', &
            '         K41  ', '         K42  ', '         K43  ', '         K44  ', '        Cext  '
      end if

      do i = 1, int(dims(1))
         write (1, '(3I8,19ES14.6)') int(A(i, 1:3)), A(i, 4:20)
      end do

      close (1)

   end subroutine write_RT_matrix

!****************************************************************************80

end module
