module sparse
!use omp_lib
use common
implicit none



contains 
!_______________________________________________________________

function complex_sparse_matmul(val,ind, x) result(Ax)
complex, dimension(:,:) :: val
integer, dimension(:,:) :: ind
double complex, dimension(:) :: x
double complex, dimension(:), allocatable :: Ax

integer ::  rows, cols, i, j, N

N=size(ind,2)
allocate(Ax(N))
Ax(:) = dcmplx(0.0,0.0)


rows=size(ind,2)
cols=size(ind,1)


!$omp parallel num_threads(nthreads) default(private) &
!$omp firstprivate(cols, rows)  &
!$omp shared(Ax, x, val, ind)
!$omp do schedule(guided,32) 

do i = 1,rows
   do j = 1, cols
      if(ind(j,i) > 0) then       
         Ax(i) = Ax(i) + dcmplx(val(j,i)) * x(ind(j,i)) 
      end if
   end do
end do

!$omp end do
!$omp end parallel

end function complex_sparse_matmul

!_______________________________________________________________

function sparse_T_matmul(val,ind, x, N) result(Ax)
double complex, dimension(:,:) :: val
integer, dimension(:,:) :: ind
double complex, dimension(:) :: x
double complex, dimension(:), allocatable :: Ax
integer :: N

integer ::  rows, cols, i, j, row_t,t1,t2,rate

allocate(Ax(N))

Ax(:) = dcmplx(0.0,0.0)
rows=size(ind,2)
cols=size(ind,1)

!!$omp parallel num_threads(nthreads) default(private) &
!!$omp firstprivate(cols, rows)  &
!!$omp shared(x, val, ind, Ax)

!!$omp do reduction(+:Ax)
do i = 1,rows  
   do j = 1, cols   
      row_t = ind(j,i)     
      Ax(row_t) = Ax(row_t) + val(j,i) * x(i)
   end do
end do
!!$omp end do

!!$omp end parallel


end function sparse_T_matmul

!_____________________________________________________________

function sparse_T_matmul_mpi(val,ind,S_tet, x, M1v, M2v) result(Ax)
double complex, dimension(:,:) :: val
integer, dimension(:,:) :: ind
integer, dimension(:) :: S_tet
double complex, dimension(:) :: x
double complex, dimension(:), allocatable :: Ax
integer :: M1v, M2v, MM

integer ::  rows, cols, i, j, row_t,t1,t2,rate

allocate(Ax(M2v-M1v+1))
MM=M1v-1

Ax(:) = dcmplx(0.0,0.0)
rows=size(ind,2)
cols=size(ind,1)

!!$omp parallel num_threads(nthreads) default(private) &
!!$omp firstprivate(cols, rows)  &
!!$omp shared(x, val, ind, Ax)

!!$omp do reduction(+:Ax)
do i = 1,rows  
   do j = 1, cols   
      row_t = ind(j,i)     
      if(row_t .ne. 0) then
         Ax(row_t-MM) = Ax(row_t-MM) + val(j,i) * x(S_tet(i))
         !print*,  row_t-MM
      end if
   end do
end do
!!$omp end do

!!$omp end parallel


end function sparse_T_matmul_mpi



!_______________________________________________________________


function sparse_matmul2(val,ind, x) result(Ax)
double complex, dimension(:,:) :: val
integer, dimension(:,:) :: ind
double complex, dimension(:) :: x
double complex, dimension(:), allocatable :: Ax

integer :: N, rows, cols, i, j

N=size(ind,2)
allocate(Ax(N))
Ax(:) = dcmplx(0.0,0.0)
rows=size(ind,2)
cols=size(ind,1)

!$omp parallel num_threads(nthreads) default(none) &
!$omp firstprivate(cols, rows)  &
!$omp private(i,j) &
!$omp shared(Ax,x, val, ind)
!$omp do schedule(guided,32)
do i = 1,rows
   do j = 1, cols
      !if(ind(j,i) > 0) then
         Ax(i) = Ax(i) + val(j,i) * x(ind(j,i))
      !end if
   end do
end do
!$omp end do
!$omp end parallel

end function sparse_matmul2
!____________________________________________________________!

function sparse_eps_matmul(eps,val,ind, x, Nbasis) result(Ax)
double complex, dimension(:,:) :: val
integer, dimension(:,:) :: ind
double complex, dimension(:) :: x, eps
double complex, dimension(:), allocatable :: Ax

double complex :: eri
integer :: N, rows, cols, i, j,tt, Nbasis

N=size(ind,2)
allocate(Ax(N))
Ax = dcmplx(0.0,0.0)
rows=size(ind,2)
cols=size(ind,1)

!$omp parallel num_threads(nthreads) default(none) &
!$omp firstprivate(cols, rows, Nbasis)  &
!$omp private(i,j,eri,tt) &
!$omp shared(Ax,x, val, ind,eps)
!$omp do 
do i = 1,rows
   !eri = eps((i-1)/Nbasis+1) - dcmplx(1.0, 0.0)
   do j = 1, cols     
      !Ax(i) = Ax(i) + val(j,i) * eri * x(ind(j,i))
      Ax(i) = Ax(i) + val(j,i) * x(ind(j,i))
   end do
end do
!$omp end do
!$omp end parallel

end function sparse_eps_matmul

!_____________________________________________________________!


!____________________________________________________________!


function sparse2dense_vec(sp_vec, sp_ind, N) result(vec)
double complex, dimension(:) :: sp_vec
integer, dimension(:) :: sp_ind
integer :: N
double complex, dimension(:), allocatable :: vec
allocate(vec(N))
vec(:) = dcmplx(0.0,0.0) 
vec(sp_ind) = sp_vec

end function sparse2dense_vec


!____________________________________________________________!

subroutine sparse2dense_vec2(sp_vec, sp_ind, vec) 
double complex, dimension(:) :: sp_vec
integer, dimension(:) :: sp_ind

double complex, dimension(:) :: vec

vec(:) = dcmplx(0.0,0.0) 
vec(sp_ind) = sp_vec

end subroutine sparse2dense_vec2

!--------------------------------------------------------------------

function sparse_T_matmul2(val,ind, x, N) result(Ax)
double complex, dimension(:,:) :: val
integer, dimension(:,:) :: ind
double complex, dimension(:) :: x
double complex, dimension(:), allocatable :: Ax, Ax_local
integer :: N

integer ::  rows, cols, i, j, row_t,t1,t2,rate

allocate(Ax(N),Ax_local(N))

Ax(:) = dcmplx(0.0,0.0)
Ax_local(:) = dcmplx(0.0,0.0)

rows=size(ind,2)
cols=size(ind,1)

!$omp parallel num_threads(nthreads) default(private) &
!$omp firstprivate(cols, rows, Ax_local)  &
!$omp shared(x, val, ind)
!$omp do 
do i = 1,rows  
   do j = 1, cols   
      row_t = ind(j,i)     
      Ax_local(row_t) = Ax_local(row_t) + val(j,i) * x(i)
   end do
end do
!$omp end do
!$omp critical
Ax = Ax + Ax_local
!$omp end critical
!$omp end parallel


end function sparse_T_matmul2


end module sparse
