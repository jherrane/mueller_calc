module precorrection 
use common  
use integrals
use integration_points
use possu
use geometry
use sparse !lin
use sparse_mat

implicit none

contains


function compute_pG(mesh, t_cube, b_cube) result(G)
implicit none
type (mesh_struct) :: mesh

integer, intent(in) :: t_cube, b_cube

double complex, dimension(:,:), allocatable :: G 
integer :: N, i1, i2
double precision :: r(3), rp(3)

N = size(mesh%etopol_box,1)

allocate(G(N,N))
G(:,:) = dcmplx(0.0,0.0)

do i1 = 1,N
   r = mesh%nodes(:,mesh%etopol_box(i1,t_cube))
   do i2 = 1,N
      if(mesh%etopol_box(i1,t_cube) .ne. mesh%etopol_box(i2,b_cube)) then
         rp = mesh%nodes(:,mesh%etopol_box(i2,b_cube))
         G(i1,i2) = Gr(r,rp,mesh%k)
      end if
   end do
end do

end function compute_pG


subroutine precorrect_const(matrices, tet_T, tet_B, G, corr1, corr2)
  type (data), intent(in) :: matrices
double complex, dimension(:,:), intent(in) :: G
double complex, dimension(3,3) :: corr1, corr2 
integer :: tet_T, tet_B
double complex :: a, vec_x(size(G,1)), vec_y(size(G,1)), vec_z(size(G,1))


corr1(:,:) = dcmplx(0.0, 0.0)
corr2(:,:) = dcmplx(0.0, 0.0)

a = dot_product(conjg(matmul(G,matrices%S(:,tet_B))), matrices%S(:,tet_T))

corr1(1,1) = a
corr1(2,2) = a
corr1(3,3) = a

vec_x = matmul(G,matrices%Sx(:,tet_B))
vec_y = matmul(G,matrices%Sy(:,tet_B))
vec_z = matmul(G,matrices%Sz(:,tet_B))

corr2(1,1) = dot_product(conjg(vec_x), matrices%Sx(:,tet_T))
corr2(1,2) = dot_product(conjg(vec_y), matrices%Sx(:,tet_T))
corr2(1,3) = dot_product(conjg(vec_z), matrices%Sx(:,tet_T))

corr2(2,1) = dot_product(conjg(vec_x), matrices%Sy(:,tet_T))
corr2(2,2) = dot_product(conjg(vec_y), matrices%Sy(:,tet_T))
corr2(2,3) = dot_product(conjg(vec_z), matrices%Sy(:,tet_T))

corr2(3,1) = dot_product(conjg(vec_x), matrices%Sz(:,tet_T))
corr2(3,2) = dot_product(conjg(vec_y), matrices%Sz(:,tet_T))
corr2(3,3) = dot_product(conjg(vec_z), matrices%Sz(:,tet_T))

end subroutine precorrect_const

!_______________________________________________________________________

subroutine precorrect_lin(matrices, tet_T, tet_B, G, corr1, corr2)
  type (data), intent(in) :: matrices
double complex, dimension(:,:), intent(in) :: G
double complex, dimension(12,12) :: corr1, corr2 
integer :: tet_T, tet_B, i1, i2
double complex :: a, vec_x(size(G,1)), vec_y(size(G,1)), vec_z(size(G,1))


corr1(:,:) = dcmplx(0.0, 0.0)
corr2(:,:) = dcmplx(0.0, 0.0)


do i1 = 1,4
   do i2 = 1,4
      a = dot_product(conjg(matmul(G,matrices%S(:,4*(tet_B-1)+i1))), matrices%S(:,4*(tet_T-1)+i2))
      corr1(i2,i1) = a
      corr1(4+i2,4+i1) = a
      corr1(8+i2,8+i1) = a
   end do
end do


do i1 = 1,4

   vec_x = matmul(G,matrices%Sx(:,4*(tet_B-1)+i1))
   vec_y = matmul(G,matrices%Sy(:,4*(tet_B-1)+i1))
   vec_z = matmul(G,matrices%Sz(:,4*(tet_B-1)+i1))

   do i2 = 1,4
      corr2(i2,i1) = dot_product(conjg(vec_x), matrices%Sx(:,4*(tet_T-1)+i2))
      corr2(i2,4+i1) = dot_product(conjg(vec_y), matrices%Sx(:,4*(tet_T-1)+i2))
      corr2(i2,8+i1) = dot_product(conjg(vec_z), matrices%Sx(:,4*(tet_T-1)+i2))

      corr2(4+i2,i1) = dot_product(conjg(vec_x), matrices%Sy(:,4*(tet_T-1)+i2))
      corr2(4+i2,4+i1) = dot_product(conjg(vec_y), matrices%Sy(:,4*(tet_T-1)+i2))
      corr2(4+i2,8+i1) = dot_product(conjg(vec_z), matrices%Sy(:,4*(tet_T-1)+i2))

      corr2(8+i2,i1) = dot_product(conjg(vec_x), matrices%Sz(:,4*(tet_T-1)+i2))
      corr2(8+i2,4+i1) = dot_product(conjg(vec_y), matrices%Sz(:,4*(tet_T-1)+i2))
      corr2(8+i2,8+i1) = dot_product(conjg(vec_z), matrices%Sz(:,4*(tet_T-1)+i2))

   end do
end do

end subroutine precorrect_lin

!***********************************************************************

subroutine compute_near_zone_interactions_const(matrices,mesh)
type (mesh_struct), intent(in) :: mesh
type (data) :: matrices

integer(kind=8) :: mem
integer :: t_cube, b_cube, i1, t1, b1, tet_T, tet_B
integer :: T_nodes(4), B_nodes(4), bases(mesh%N_tet_cube), tests(mesh%N_tet_cube)
integer, dimension(:), allocatable :: near_cubes
double precision ::  vol_T, vol_B, B_coord(3,4), T_coord(3,4)
double complex :: mat_block(3,3), corr(3,3), corr2(3,3), ele_param
double complex, dimension(:,:), allocatable :: Gcorr
double precision :: T_rot(3,6), T_dN(3,4), B_rot(3,6), B_dN(3,4)
double precision, dimension(:,:), allocatable :: P0_tet, P0_tri
double precision, dimension(:), allocatable :: w0_tet, w0_tri
double complex :: alok(3,3), alok1(3,3), alok2(3,3)

integer :: blok, M_cube

M_cube = (2*mesh%near_zone + 1)**3

allocate(Gcorr(size(mesh%etopol_box,1),size(mesh%etopol_box,1)))

allocate(near_cubes(M_cube))

corr(:,:) = dcmplx(0.0,0.0)
corr2(:,:) = dcmplx(0.0,0.0)

call inttetra(P0_tet,w0_tet,5)
call inttri(P0_tri,w0_tri,5)

call allocate_Acorr(matrices%sp_mat,matrices%sp_ind, mesh)


matrices%sp_mat(:,:,:) = cmplx(0.0,0.0)
matrices%sp_ind(:,:) = 0

mem = sizeof(matrices%sp_mat)/1024/1024 + &
sizeof(matrices%sp_ind)/1024/1024 + &
sizeof(matrices%S)/1024/1024 * 4 + &
sizeof(matrices%indS)/1024/1024 + &
sizeof(matrices%Fg)/1024/1024*2 + &
(5+mesh%maxit)*8*2*3*size(mesh%etopol,2)/1024/1024 + &
sizeof(mesh%coord)/1024/1024 + &
sizeof(mesh%etopol)/1024/1024 + &
sizeof(mesh%etopol_box)/1024/1024 + &
sizeof(mesh%nodes)/1024/1024 + &
sizeof(mesh%tetras)/1024/1024

print*, 'Estimated memory usage:', mem, 'Mb'

blok=1

do t_cube = 1, mesh%N_cubes
  
   tests = mesh%tetras(:,t_cube)
   
   near_cubes = find_near_cubes(mesh,t_cube)

   do i1 = 1,M_cube  
      if(tests(1)==0) exit
      b_cube = near_cubes(i1)      
      if(b_cube == 0) exit
      
      bases = mesh%tetras(:,b_cube)     
      Gcorr = compute_pG(mesh, t_cube, b_cube)


      do t1 = 1,mesh%N_tet_cube
         tet_T = tests(t1)
         if(tet_T == 0) exit 

         T_nodes =  mesh%etopol(:,tet_T)
         T_coord = mesh%coord(:,T_nodes)
         vol_T = tetra_volume(T_coord) 

         call gradshape(T_rot,T_dN,T_coord)

         ele_param = mesh%param(tet_T)      

         do b1 = 1,mesh%N_tet_cube    
            tet_B = bases(b1)

            if(tet_B == 0) exit
            if(tet_B >= tet_T) then
            B_nodes =  mesh%etopol(:,tet_B)
            B_coord = mesh%coord(:,B_nodes)           
            vol_B = tetra_volume(B_coord) 

            call gradshape(B_rot,B_dN,B_coord)

            !_________-Correction term____________________________________
          
            call precorrect_const(matrices, tet_T, tet_B, Gcorr, corr, corr2)
            corr = corr * sqrt(vol_T * vol_B)
            corr2 = corr2 * sqrt(vol_T * vol_B)

            !_____________________________________________________________ 

            alok(:,:) = dcmplx(0.0,0.0)
            if(tet_T == tet_B) then 
               alok(1,1) = 1.0 / (ele_param-1.0) 
               alok(2,2) = 1.0 / (ele_param-1.0) 
               alok(3,3) = 1.0 / (ele_param-1.0) 
              
               !alok(1,1) = ele_param / (ele_param-1.0) 
               !alok(2,2) = ele_param / (ele_param-1.0) 
               !alok(3,3) = ele_param / (ele_param-1.0) 
            end if
          
          

            call integrate_V_V_G(alok1,mesh,P0_tet, w0_tet,P0_tet, w0_tet,T_coord, B_coord) 
            call integrate_dV_dV_G(alok2,mesh,P0_tri,w0_tri,P0_tri,w0_tri ,T_coord,B_coord)
           
            !call integrate_nx_dV_dV_G(alok2,mesh,P0_tri,w0_tri,P0_tri,w0_tri ,T_coord,B_coord)
           
            ! call integrate_dV_V_gradG(alok2,mesh,P0_tet,w0_tet,P0_tri,w0_tri ,T_coord,B_coord)
  
          

            mat_block = alok + 1.0/sqrt(vol_B*vol_T) * &
(-mesh%k**2 * (alok1-corr) + (alok2 - corr2) )

            !mat_block = alok + 1.0/sqrt(vol_B*vol_T) * &
!(-mesh%k**2 * (-corr) + (-alok2 - corr2) )

           
            matrices%sp_mat(:,:,blok) = cmplx(mat_block)
            matrices%sp_ind(1,blok) = tet_T
            matrices%sp_ind(2,blok) = tet_B

            blok=blok+1
        
         end if
         end do  
      end do       
   end do  
end do


end subroutine compute_near_zone_interactions_const

!*************************************************
!*
!***********************************************

subroutine compute_near_zone_interactions_lin(matrices,mesh)

type (mesh_struct), intent(in) :: mesh
type (data) :: matrices

integer(kind = 8) :: mem
integer :: t_cube, b_cube, i1, t1, b1, tet_T, tet_B
integer :: T_nodes(4), B_nodes(4), bases(mesh%N_tet_cube), tests(mesh%N_tet_cube)
integer, dimension(:), allocatable :: near_cubes
double precision ::  vol_T, vol_B,  B_coord(3,4), T_coord(3,4)
double complex :: mat_block(12,12), corr(12,12), corr2(12,12), ele_param
double complex, dimension(:,:), allocatable :: Gcorr
double precision :: T_rot(3,6), T_dN(3,4), B_rot(3,6), B_dN(3,4)
double precision, dimension(:,:), allocatable :: P0_tet, P0_tri
double precision, dimension(:), allocatable :: w0_tet, w0_tri
double complex :: alok(12,12), alok1(12,12), alok2(12,12), alok3(12,12), alok4(12,12), alok5(12,12)
double complex :: GNN(4,4)
integer :: ind(12,2), blok, M_cube

M_cube = (2*mesh%near_zone + 1)**3

allocate(Gcorr(size(mesh%etopol_box,1),size(mesh%etopol_box,1)))

allocate(near_cubes(M_cube))

corr(:,:) = dcmplx(0.0,0.0)
corr2(:,:) = dcmplx(0.0,0.0)

call inttetra(P0_tet,w0_tet,5)
!call inttetra2(P0_tet,w0_tet,4)
call inttri(P0_tri,w0_tri,5)

call allocate_Acorr(matrices%sp_mat,matrices%sp_ind, mesh)

matrices%sp_mat(:,:,:) = cmplx(0.0,0.0)
matrices%sp_ind(:,:) = 0

mem = sizeof(matrices%sp_mat)/1024/1024 + &
sizeof(matrices%sp_ind)/1024/1024 + &
sizeof(matrices%S)/1024/1024 * 4 + &
sizeof(matrices%indS)/1024/1024 + &
sizeof(matrices%Fg)/1024/1024*2 + &
(5+mesh%maxit)*8*2*12*size(mesh%etopol,2)/1024/1024 + &
sizeof(mesh%coord)/1024/1024 + &
sizeof(mesh%etopol)/1024/1024 + &
sizeof(mesh%etopol_box)/1024/1024 + &
sizeof(mesh%nodes)/1024/1024 + &
sizeof(mesh%tetras)/1024/1024

print*, 'Estimated memory usage:', mem, 'Mb'

blok=1


do t_cube = 1, mesh%N_cubes
  
   tests = mesh%tetras(:,t_cube)
   
   near_cubes = find_near_cubes(mesh,t_cube)

   do i1 = 1,M_cube  
      if(tests(1)==0) exit
      b_cube = near_cubes(i1)      
      if(b_cube == 0) exit
      
      bases = mesh%tetras(:,b_cube)     
      Gcorr = compute_pG(mesh, t_cube, b_cube)


      do t1 = 1,mesh%N_tet_cube
         tet_T = tests(t1)
         if(tet_T == 0) exit 

         T_nodes =  mesh%etopol(:,tet_T)
         T_coord = mesh%coord(:,T_nodes)
         vol_T = tetra_volume(T_coord) 

         call gradshape(T_rot,T_dN,T_coord)

         ele_param = mesh%param(tet_T)      

         do b1 = 1,mesh%N_tet_cube    
            tet_B = bases(b1)

            if(tet_B == 0) exit
            if(tet_B >= tet_T) then
            B_nodes =  mesh%etopol(:,tet_B)
            B_coord = mesh%coord(:,B_nodes)           
            vol_B = tetra_volume(B_coord) 

            call gradshape(B_rot,B_dN,B_coord)

            !_________-Correction term____________________________________
          
            call precorrect_lin(matrices, tet_T, tet_B, Gcorr, corr, corr2)
            corr = corr * sqrt(vol_T * vol_B)
            corr2 = corr2 * sqrt(vol_T * vol_B)

            !_____________________________________________________________ 

            alok(:,:) = dcmplx(0.0,0.0)
            if(tet_T == tet_B) then    
              call integrate_NN(alok, P0_tet, w0_tet, T_coord)
            end if
          
            call integrate_V_V_GNN(GNN, mesh, P0_tet, w0_tet, T_coord, B_coord)
            call integrate_dV_dV_GNN(alok2, mesh, P0_tri, w0_tri, T_coord, B_coord)
            call integrate_dV_V_GNN(alok4, mesh, P0_tet, w0_tet, P0_tri, w0_tri, T_coord, B_coord, B_dN)
            call integrate_V_dV_GNN(alok5, mesh, P0_tet, w0_tet, P0_tri, w0_tri, T_coord, B_coord, T_dN)
            call aloks(alok1,alok3,ind, GNN, tet_T, tet_B, T_dN, B_dN)

                  
            !mat_block = 1.0/sqrt(vol_B*vol_T) * (alok + (ele_param-1.0) * &
!(-mesh%k**2 * (alok1-corr) + alok2 + alok3 - alok4 - alok5 - corr2))

            mat_block = 1.0/sqrt(vol_B*vol_T) * (alok/(ele_param-1.0) +  &
(-mesh%k**2 * (alok1-corr) + (alok2 + alok3 - alok4 - alok5) - corr2))

           
            matrices%sp_mat(:,:,blok) = cmplx(mat_block)
            matrices%sp_ind(1,blok) = tet_T
            matrices%sp_ind(2,blok) = tet_B

            blok=blok+1
          
         end if          
         end do  
      end do       
   end do  
end do


end subroutine compute_near_zone_interactions_lin


end module precorrection
