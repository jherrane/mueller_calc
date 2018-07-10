module sparse_mat
use common
use geometry
implicit none

contains

subroutine allocate_Acorr(sp_mat,sp_ind, mesh)
type (mesh_struct) :: mesh
complex, dimension(:,:,:), allocatable :: sp_mat
integer, dimension(:,:), allocatable :: sp_ind

integer :: t_cube, b_cube, tests(mesh%N_tet_cube), bases(mesh%N_tet_cube)
integer :: i1, t1, b1, tet_T, tet_B, las, M_cube
integer, dimension(:), allocatable :: near_cubes

M_cube = (2*mesh%near_zone + 1)**3
allocate(near_cubes(M_cube))

las = 1
do t_cube = 1, mesh%N_cubes

   tests = mesh%tetras(:,t_cube)   
   near_cubes = find_near_cubes(mesh,t_cube)

   do i1 = 1,M_cube
      if(tests(1)==0) exit
      b_cube = near_cubes(i1)      
      if(b_cube == 0) exit
  
      bases = mesh%tetras(:,b_cube) 

      do t1 = 1,mesh%N_tet_cube
         tet_T = tests(t1)
         if(tet_T == 0) exit 
         
         do b1 = 1,mesh%N_tet_cube    
            tet_B = bases(b1)
            if(tet_B == 0) exit
            if(tet_B <= tet_T) then
            las = las + 1
            end if
         end do

      end do
   end do

end do


if(mesh%order == 0) then
   allocate(sp_mat(3,3,las-1), sp_ind(2,las-1))
end if
if(mesh%order == 1) then
   allocate(sp_mat(12,12,las-1), sp_ind(2,las-1))
end if


end subroutine allocate_Acorr

subroutine allocate_Acorr2(sp_mat,sp_ind, mesh, val, ind)
type (mesh_struct) :: mesh
complex, dimension(:,:,:), allocatable :: sp_mat
integer, dimension(:,:), allocatable :: sp_ind
integer, dimension(:), allocatable :: val, ind

integer :: t_cube, b_cube, tests(mesh%N_tet_cube), bases(mesh%N_tet_cube)
integer :: near_cubes(27), i1, t1, b1, tet_T, tet_B, las, las2, las3
integer :: T_nodes(4)
double precision :: T_coord(3,4)

las = 1
do t_cube = 1, mesh%N_cubes

   tests = mesh%tetras(:,t_cube)   
   near_cubes = find_near_cubes(mesh,t_cube)

   do i1 = 1,27
      if(tests(1)==0) exit
      b_cube = near_cubes(i1)      
      if(b_cube == 0) exit
  
      bases = mesh%tetras(:,b_cube) 

      do t1 = 1,mesh%N_tet_cube
         tet_T = tests(t1)
         if(tet_T == 0) exit 

         do b1 = 1,mesh%N_tet_cube    
            tet_B = bases(b1)
            if(tet_B == 0) exit
            las = las + 1
         end do

      end do
   end do

end do

allocate(sp_mat(12,12,las-1), sp_ind(2,las-1))
allocate(val(las-1))
allocate(ind(mesh%N_tet))
ind(:) = 0
las2 = 1
las3 = 1

do tet_T = 1, mesh%N_tet

   T_nodes = mesh%etopol(:,tet_T)
   T_coord = mesh%coord(:,T_nodes)
   t_cube = point_in_cube(mesh, sum(T_coord,2)/4.0)
   near_cubes = find_near_cubes(mesh,t_cube)

   las = 0
   do i1 = 1,27
      b_cube = near_cubes(i1)      
      if(b_cube == 0) exit
  
      bases = mesh%tetras(:,b_cube) 
       
      do b1 = 1,mesh%N_tet_cube    
         tet_B = bases(b1)
         if(tet_B == 0) exit

         val(las2) = tet_B
         las = las + 1
         las2 = las2 + 1
      end do

        
       
   end do

   ind(tet_T) =  las
     
end do



end subroutine allocate_Acorr2


function matmul_Acorr(val,ind,x) result(Ax)
complex, dimension(:,:,:) :: val
integer, dimension(:,:) :: ind
double complex, dimension(:) :: x
double complex, dimension(:), allocatable :: Ax
integer :: N,M, i1, a,b

N = size(x)
M = size(ind,2)
allocate(Ax(N))
Ax(:) = dcmplx(0.0,0.0)
!$omp parallel num_threads(nthreads) default(private) &
!$omp shared(x, val, ind, Ax,M)
!$omp do reduction(+:Ax)
do i1 = 1,M
   a = 12*(ind(1,i1)-1)
   b = 12*(ind(2,i1)-1)
   Ax(a+1:a+12) =   Ax(a+1:a+12) + matmul(dcmplx(val(:,:,i1)),x(b+1:b+12))
 
   if(a .ne. b) then
      Ax(b+1:b+12) =   Ax(b+1:b+12) + matmul(dcmplx(transpose(val(:,:,i1))),x(a+1:a+12))
   end if

end do

!$omp end do
!$omp end parallel

end function matmul_Acorr

function matmul_Acorr_const(val,ind,x) result(Ax)
complex, dimension(:,:,:) :: val
integer, dimension(:,:) :: ind
double complex, dimension(:) :: x
double complex, dimension(:), allocatable :: Ax
integer :: N,M, i1, a,b

N = size(x)
M = size(ind,2)
allocate(Ax(N))
Ax(:) = dcmplx(0.0,0.0)
!$omp parallel num_threads(nthreads) default(private) &
!$omp shared(x, val, ind, Ax,M)
!$omp do reduction(+:Ax)
do i1 = 1,M
   a = 3*(ind(1,i1)-1)
   b = 3*(ind(2,i1)-1)
   Ax(a+1:a+3) =   Ax(a+1:a+3) + matmul(dcmplx(val(:,:,i1)),x(b+1:b+3)) 

   if(a .ne. b) then
      Ax(b+1:b+3) =   Ax(b+1:b+3) + matmul(dcmplx(transpose(val(:,:,i1))),x(a+1:a+3))
   end if
end do

!$omp end do
!$omp end parallel

end function matmul_Acorr_const

end module sparse_mat
