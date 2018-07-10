module integrals
use singularity_subtraction_N
use singularity_subtraction
implicit none

contains


!*************************************************************************

subroutine integrate_NN(alok, P0_tet, w0_tet, T_coord)
double precision, dimension(:,:), intent(in) :: P0_tet
double precision, dimension(:), intent(in) :: w0_tet
double precision, intent(in) ::  T_coord(3,4)
double complex, dimension(12,12) :: alok

double complex, dimension(4,4) :: intN
double precision :: Pt(3,size(w0_tet)), wt(size(w0_tet))
double precision :: shape_tet(4,size(w0_tet))
integer :: i1

shape_tet(1,:) = 1-P0_tet(1,:)-P0_tet(2,:) -P0_tet(3,:) 
shape_tet(2,:) = P0_tet(1,:)
shape_tet(3,:) = P0_tet(2,:)
shape_tet(4,:) = P0_tet(3,:)

call linmap_tet(Pt, wt, T_coord, P0_tet, w0_tet)
alok(:,:) = dcmplx(0.0,0.0)
intN(:,:) = dcmplx(0.0,0.0)
do i1 = 1,size(wt)

   intN(:,1) = intN(:,1) + shape_tet(1,i1) * shape_tet(:,i1) * wt(i1)
   intN(:,2) = intN(:,2) + shape_tet(2,i1) * shape_tet(:,i1) * wt(i1)
   intN(:,3) = intN(:,3) + shape_tet(3,i1) * shape_tet(:,i1) * wt(i1)
   intN(:,4) = intN(:,4) + shape_tet(4,i1) * shape_tet(:,i1) * wt(i1)

end do

alok(1:4,1:4) = intN
alok(5:8,5:8) = intN
alok(9:12,9:12) = intN

end subroutine integrate_NN

subroutine integrate_V_V_GNN(intN, mesh, P0_tet, w0_tet, T_coord, B_coord)
type (mesh_struct), intent(in) :: mesh
double precision, dimension(:,:), intent(in) :: P0_tet
double precision, dimension(:), intent(in) :: w0_tet
double precision, intent(in) ::  T_coord(3,4), B_coord(3,4)

integer :: i1

double precision :: Pb(3,size(w0_tet)), wb(size(w0_tet))
double precision :: Pt(3,size(w0_tet)), wt(size(w0_tet))
double precision :: shape_tet(4,size(w0_tet)), rf(3)
 
double complex :: intN(4,4), I(4) 

shape_tet(1,:) = 1-P0_tet(1,:)-P0_tet(2,:) -P0_tet(3,:) 
shape_tet(2,:) = P0_tet(1,:)
shape_tet(3,:) = P0_tet(2,:)
shape_tet(4,:) = P0_tet(3,:)

call linmap_tet(Pt, wt, T_coord, P0_tet, w0_tet)
call linmap_tet(Pb, wb, B_coord, P0_tet, w0_tet)

intN(:,:) = dcmplx(0,0)

do i1 = 1,size(wt)
   rf = Pt(:,i1)
   I = singularity_subtraction_int_V_GN(rf,mesh%k, Pb, wb, B_coord,shape_tet)

   intN(1,:) = intN(1,:) + I * shape_tet(1,i1) * wt(i1)
   intN(2,:) = intN(2,:) + I * shape_tet(2,i1) * wt(i1)
   intN(3,:) = intN(3,:) + I * shape_tet(3,i1) * wt(i1)
   intN(4,:) = intN(4,:) + I * shape_tet(4,i1) * wt(i1)

end do

end subroutine integrate_V_V_GNN


!************************************************************************

subroutine integrate_dV_V_GNN(alok, mesh, P0_tet, w0_tet, P0_tri, w0_tri, T_coord, B_coord, B_dN)
type (mesh_struct), intent(in) :: mesh
double precision, dimension(:,:), intent(in) :: P0_tet, P0_tri
double precision, dimension(:), intent(in) :: w0_tet, w0_tri
double precision, intent(in) ::  T_coord(3,4), B_coord(3,4), B_dN(3,4)
double complex :: alok(12,12)

integer :: i1, i2, face_t, T_nodes(3)
integer, parameter :: ind(12)=[2,3,4,1,4,3,1,2,4,1,3,2]

double precision :: Pb(3,size(w0_tet)), wb(size(w0_tet))
double precision :: Pt(3,size(w0_tri)), wt(size(w0_tri))
double precision :: shape_tet(4,size(w0_tet)),shape_tri(3,size(w0_tri)), rf(3)
double precision :: T_tri_coord(3,3), T_nvec(3), dd(12)
 
double complex :: intN(4), I, II(12) 

shape_tet(1,:) = 1-P0_tet(1,:)-P0_tet(2,:) -P0_tet(3,:) 
shape_tet(2,:) = P0_tet(1,:)
shape_tet(3,:) = P0_tet(2,:)
shape_tet(4,:) = P0_tet(3,:)

shape_tri(1,:) = 1-P0_tri(1,:)-P0_tri(2,:) 
shape_tri(2,:) = P0_tri(1,:)
shape_tri(3,:) = P0_tri(2,:)

dd = [B_dN(1,:),B_dN(2,:),B_dN(3,:)]

call linmap_tet(Pb, wb, B_coord, P0_tet, w0_tet)

II(:) = dcmplx(0.0,0.0)

do face_t = 1,4
   T_nodes(1) = ind(3*(face_t-1)+1)
   T_nodes(2) = ind(3*(face_t-1)+2)
   T_nodes(3) = ind(3*(face_t-1)+3)

   T_tri_coord = T_coord(:,T_nodes)
   T_nvec = tri_n_vectors(T_tri_coord)
   call linmap_tri(Pt, wt, T_tri_coord, P0_tri, W0_tri)

   intN(:) = dcmplx(0.0,0.0)
   do i1 = 1,size(wt)
      rf = Pt(:,i1)
      I = sum(singularity_subtraction_int_V_GN(rf,mesh%k, Pb, wb, B_coord,shape_tet))      
      intN(T_nodes) = intN(T_nodes) + I * shape_tri(:,i1) * wt(i1)
   end do

   II(1:4) = II(1:4) + intN*T_nvec(1)
   II(5:8) = II(5:8) + intN*T_nvec(2)
   II(9:12) = II(9:12) + intN*T_nvec(3)

end do

do i2 = 1,12
   alok(i2,:) = II(i2)*dd
end do

end subroutine integrate_dV_V_GNN




!*****************************************************************************
!
!
!
!****************************************************************************

subroutine integrate_dV_dV_GNN(alok, mesh, P0_tri, w0_tri, T_coord, B_coord)
type (mesh_struct), intent(in) :: mesh
double precision, dimension(:,:), intent(in) :: P0_tri
double precision, dimension(:), intent(in) :: w0_tri
double precision, intent(in) :: T_coord(3,4), B_coord(3,4)

integer :: T_nodes(3), B_nodes(3), face_b, face_t, i1

double precision :: Pt(3,size(w0_tri)), wt(size(w0_tri))
double precision :: Pb(3,size(w0_tri)), wb(size(w0_tri)), shape_tri(3,size(w0_tri))
integer, parameter :: ind(12)=[2,3,4,1,4,3,1,2,4,1,3,2]
double precision :: T_tri_coord(3,3), B_tri_coord(3,3)
double precision :: T_nvec(3), B_nvec(3), rf(3) 
double complex :: Ix(4), Iy(4), Iz(4), IIx(4,4), IIy(4,4), IIz(4,4), alok(12,12)
double complex :: intN(3)

shape_tri(1,:) = 1.0 - P0_tri(1,:) - P0_tri(2,:)
shape_tri(2,:) = P0_tri(1,:)
shape_tri(3,:) = P0_tri(2,:)

alok(:,:) = dcmplx(0.0,0.0)

do face_t = 1,4
   T_nodes(1) = ind(3*(face_t-1)+1)
   T_nodes(2) = ind(3*(face_t-1)+2)
   T_nodes(3) = ind(3*(face_t-1)+3)

   !T_tri_coord = tetra_face(T_coord,face_t)
   T_tri_coord = T_coord(:,T_nodes)
   T_nvec = tri_n_vectors(T_tri_coord)
   call linmap_tri(Pt, wt, T_tri_coord, P0_tri, W0_tri)

   IIx(:,:) = dcmplx(0.0,0.0)
   IIy(:,:) = dcmplx(0.0,0.0)
   IIz(:,:) = dcmplx(0.0,0.0)

   do i1 = 1,size(wt)
      rf = Pt(:,i1)

      !________________________________________________________
      Ix(:) = dcmplx(0.0,0.0)
      Iy(:) = dcmplx(0.0,0.0)
      Iz(:) = dcmplx(0.0,0.0)

      do face_b = 1,4
      
         B_nodes(1) = ind(3*(face_b-1)+1)
         B_nodes(2) = ind(3*(face_b-1)+2)
         B_nodes(3) = ind(3*(face_b-1)+3)

         !B_tri_coord = tetra_face(B_coord,face_b)
         B_tri_coord = B_coord(:,B_nodes)
         B_nvec = tri_n_vectors(B_tri_coord)
         call linmap_tri(Pb, wb, B_tri_coord, P0_tri, W0_tri)
         
         intN = singularity_subtraction_int_S_GN(rf, B_tri_coord, B_nvec, mesh%k, Pb, wb,shape_tri)
     
         Ix(B_nodes) = Ix(B_nodes) + intN*B_nvec(1)
         Iy(B_nodes) = Iy(B_nodes) + intN*B_nvec(2)
         Iz(B_nodes) = Iz(B_nodes) + intN*B_nvec(3)
 
      end do
      !_________________________________________________________

      IIx(T_nodes(1),:) = IIx(T_nodes(1),:) + Ix * shape_tri(1,i1) * wt(i1)
      IIx(T_nodes(2),:) = IIx(T_nodes(2),:) + Ix * shape_tri(2,i1) * wt(i1)
      IIx(T_nodes(3),:) = IIx(T_nodes(3),:) + Ix * shape_tri(3,i1) * wt(i1)

      IIy(T_nodes(1),:) = IIy(T_nodes(1),:) + Iy * shape_tri(1,i1) * wt(i1)
      IIy(T_nodes(2),:) = IIy(T_nodes(2),:) + Iy * shape_tri(2,i1) * wt(i1)
      IIy(T_nodes(3),:) = IIy(T_nodes(3),:) + Iy * shape_tri(3,i1) * wt(i1)

      IIz(T_nodes(1),:) = IIz(T_nodes(1),:) + Iz * shape_tri(1,i1) * wt(i1)
      IIz(T_nodes(2),:) = IIz(T_nodes(2),:) + Iz * shape_tri(2,i1) * wt(i1)
      IIz(T_nodes(3),:) = IIz(T_nodes(3),:) + Iz * shape_tri(3,i1) * wt(i1)

   end do

   alok(1:4,1:4) =  alok(1:4,1:4) + IIx * T_nvec(1)
   alok(1:4,5:8) = alok(1:4,5:8) + IIy * T_nvec(1)
   alok(1:4,9:12) = alok(1:4,9:12) + IIz * T_nvec(1)

   alok(5:8,1:4) = alok(5:8,1:4) + IIx * T_nvec(2)
   alok(5:8,5:8) = alok(5:8,5:8) + IIy * T_nvec(2)
   alok(5:8,9:12) = alok(5:8,9:12)+ IIz * T_nvec(2)

   alok(9:12,1:4) = alok(9:12,1:4) + IIx * T_nvec(3)
   alok(9:12,5:8) = alok(9:12,5:8) + IIy * T_nvec(3)
   alok(9:12,9:12) = alok(9:12,9:12) + IIz * T_nvec(3)

end do

end subroutine integrate_dV_dV_GNN

!**********************************************************************


subroutine integrate_V_dV_GNN(alok, mesh, P0_tet, w0_tet, P0_tri, w0_tri, T_coord, B_coord, T_dN)
type (mesh_struct), intent(in) :: mesh
double precision, dimension(:,:), intent(in) :: P0_tri, P0_tet
double precision, dimension(:), intent(in) :: w0_tri, w0_tet
double precision, intent(in) :: T_coord(3,4), B_coord(3,4), T_dN(3,4)

integer :: B_nodes(3), face_b, i1, i2

double precision :: Pt(3,size(w0_tet)), wt(size(w0_tet)), shape_tet(4,size(w0_tet))
double precision :: Pb(3,size(w0_tri)), wb(size(w0_tri)), shape_tri(3,size(w0_tri))
integer, parameter :: ind(12)=[2,3,4,1,4,3,1,2,4,1,3,2]
double precision :: B_tri_coord(3,3)
double precision :: B_nvec(3), rf(3), dd(12) 
double complex :: I(12), II(12), alok(12,12)
double complex :: intN(3)

shape_tri(1,:) = 1.0 - P0_tri(1,:) - P0_tri(2,:)
shape_tri(2,:) = P0_tri(1,:)
shape_tri(3,:) = P0_tri(2,:)

shape_tet(1,:) = 1-P0_tet(1,:)-P0_tet(2,:) -P0_tet(3,:) 
shape_tet(2,:) = P0_tet(1,:)
shape_tet(3,:) = P0_tet(2,:)
shape_tet(4,:) = P0_tet(3,:)

dd = [T_dN(1,:),T_dN(2,:),T_dN(3,:)]

alok(:,:) = dcmplx(0.0,0.0)
II(:) = dcmplx(0.0,0.0)


call linmap_tet(Pt, wt, T_coord, P0_tet, w0_tet)

do i1 = 1,size(wt)

   rf = Pt(:,i1)

   I(:) = dcmplx(0.0,0.0)
  

   do face_b = 1,4
      
      B_nodes(1) = ind(3*(face_b-1)+1)
      B_nodes(2) = ind(3*(face_b-1)+2)
      B_nodes(3) = ind(3*(face_b-1)+3)

      !B_tri_coord = tetra_face(B_coord,face_b)
      B_tri_coord = B_coord(:,B_nodes)
      B_nvec = tri_n_vectors(B_tri_coord)
      call linmap_tri(Pb, wb, B_tri_coord, P0_tri, W0_tri)
         
      intN = singularity_subtraction_int_S_GN(rf, B_tri_coord, B_nvec, mesh%k, Pb, wb,shape_tri)
     
      I(B_nodes) = I(B_nodes) + intN*B_nvec(1)
      I(B_nodes+4) = I(B_nodes+4) + intN*B_nvec(2)
      I(B_nodes+8) = I(B_nodes+8) + intN*B_nvec(3)
 
   end do
     

   II = II + I * wt(i1)

end do

do i2 = 1,12
   alok(i2,:) = dd(i2) * II 
end do 

end subroutine integrate_V_dV_GNN








!**************************************************************************
! local matrices
!
!****************************************************************************

subroutine aloks(alok,alok2,ind, GNN, test_T, basis_T, T_dN, B_dN)
double complex :: alok(12,12), alok2(12,12)  
integer :: ind(12,2), i2
double complex :: GNN(4,4)
double precision :: T_dN(3,4), B_dN(3,4), ddT(12),ddB(12)
integer :: test_T, basis_T, i1, nt, nb, ct, cb

alok = dcmplx(0.0,0.0)
alok2 = dcmplx(0.0,0.0)

alok(1:4,1:4) = GNN
alok(5:8,5:8) = GNN
alok(9:12,9:12) = GNN

do i1 = 1,12
   ind(i1,1) = 12*(test_T-1)+i1
   ind(i1,2) = 12*(basis_T-1)+i1 
end do

!do cb = 1,3 
!   do nb = 1,4
!      do ct = 1,3 
!         do nt = 1,4
!            alok2(4*(ct-1)+nt,4*(cb-1)+nb ) = T_dN(ct,nt) * B_dN(cb,nb)
!         end do
!      end do
!   end do
!end do


ddT = [T_dN(1,:),T_dN(2,:),T_dN(3,:)]
ddB = [B_dN(1,:),B_dN(2,:),B_dN(3,:)]

do i2 = 1,12
   alok2(i2,:) = ddT(i2)*ddB 
end do

alok2 = alok2 * sum(GNN)

end subroutine



!*************************************************************************

subroutine integrate_V_V_G(intN2, mesh, P0_tet, w0_tet,P02_tet, w02_tet , T_coord, B_coord)
type (mesh_struct), intent(in) :: mesh
double precision, dimension(:,:), intent(in) :: P0_tet, P02_tet
double precision, dimension(:), intent(in) :: w0_tet, w02_tet
double precision, intent(in) ::  T_coord(3,4), B_coord(3,4)

integer :: i1
double precision :: Pt(3,size(w02_tet)), wt(size(w02_tet))
double precision :: Pb(3,size(w0_tet)), wb(size(w0_tet))

double precision :: rf(3), bt(3,3) 
double complex :: intN, intN2(3,3) 



bt(:,:) = 0
bt(1,1) = 1
bt(2,2) = 1
bt(3,3) = 1

call linmap_tet(Pt, wt, T_coord, P02_tet, w02_tet)
call linmap_tet(Pb, wb, B_coord, P0_tet, w0_tet)

intN2(:,:) = dcmplx(0,0)

do i1 = 1,size(wt)

   rf = Pt(:,i1)
   intN = singularity_subtraction_int_V_G(rf,mesh%k, Pb, wb, B_coord)
   intN2 = intN2 + bt * intN * wt(i1)
  
end do

end subroutine integrate_V_V_G

!*****************************************************************************
!
!
!
!****************************************************************************

subroutine integrate_dV_dV_G(intN2, mesh, P0_tri, w0_tri, P02_tri, w02_tri, T_coord, B_coord)
type (mesh_struct), intent(in) :: mesh
double precision, dimension(:,:), intent(in) :: P0_tri, P02_tri
double precision, dimension(:), intent(in) :: w0_tri, w02_tri
double precision, intent(in) :: T_coord(3,4), B_coord(3,4)

integer :: face_t, face_b, i1
double precision :: Pt(3,size(w02_tri)), wt(size(w02_tri))
double precision :: Pb(3,size(w0_tri)), wb(size(w0_tri))

double precision :: T_tri_coord(3,3), B_tri_coord(3,3), bt(3,3)
double precision :: rf(3), T_nvec(3), B_nvec(3) 
double complex :: intN, intN2(3,3)

intN2(:,:) = dcmplx(0,0)

do face_t = 1,4

  T_tri_coord = tetra_face(T_coord,face_t)
  T_nvec = tri_n_vectors(T_tri_coord)
  call linmap_tri(Pt, wt, T_tri_coord, P02_tri, W02_tri)
  
  do face_b = 1,4
     B_tri_coord = tetra_face(B_coord,face_b)
     B_nvec = tri_n_vectors(B_tri_coord)
     call linmap_tri(Pb, wb, B_tri_coord, P0_tri, W0_tri)
   

     bt(:,1) = T_nvec * B_nvec(1)
     bt(:,2) = T_nvec * B_nvec(2)
     bt(:,3) = T_nvec * B_nvec(3)

     do i1 = 1, size(wt)
           
        rf = Pt(:,i1)
        intN = singularity_subtraction_int_S_G(rf, B_tri_coord, B_nvec, mesh%k, Pb, wb)
        intN2 = intN2 + bt * intN * wt(i1)
     end do
   
  end do

end do

end subroutine integrate_dV_dV_G

!**************************************************************************
!
!
!****************************************************************************

subroutine integrate_dV_V_gradG(intN2, mesh, P0_tet, w0_tet, P02_tri, w02_tri, T_coord, B_coord)
type (mesh_struct), intent(in) :: mesh
double precision, dimension(:,:), intent(in) :: P0_tet, P02_tri
double precision, dimension(:), intent(in) :: w0_tet, w02_tri
double precision, intent(in) :: T_coord(3,4), B_coord(3,4)

integer :: face_t, i1
double precision :: Pt(3,size(w02_tri)), wt(size(w02_tri))
double precision :: Pb(3,size(w0_tet)), wb(size(w0_tet))

double precision :: T_tri_coord(3,3), B_tri_coord(3,3)
double precision :: rf(3), T_nvec(3)
double complex :: intN(3), intN2(3,3), bt(3,3)

intN2(:,:) = dcmplx(0,0)

do face_t = 1,4

  T_tri_coord = tetra_face(T_coord,face_t)
  T_nvec = tri_n_vectors(T_tri_coord)
  call linmap_tri(Pt, wt, T_tri_coord, P02_tri, W02_tri)
   
  call linmap_tet(Pb, wb, B_coord, P0_tet, W0_tet)
  
  do i1 = 1, size(wt)
     
     rf = Pt(:,i1)
     intN = singularity_subtraction_int_V_gradG(rf,mesh%k, Pb, wb, B_coord)

     bt(:,1) = T_nvec * intN(1)
     bt(:,2) = T_nvec * intN(2)
     bt(:,3) = T_nvec * intN(3)

     intN2 = intN2 + bt * wt(i1)

  end do
   
end do

end subroutine integrate_dV_V_gradG

!*********************************************************************

subroutine integrate_nx_dV_dV_G(intN2, mesh, P0_tri, w0_tri, P02_tri, w02_tri, T_coord, B_coord)
type (mesh_struct), intent(in) :: mesh
double precision, dimension(:,:), intent(in) :: P0_tri, P02_tri
double precision, dimension(:), intent(in) :: w0_tri, w02_tri
double precision, intent(in) :: T_coord(3,4), B_coord(3,4)

integer :: face_t, face_b, i1
double precision :: Pt(3,size(w02_tri)), wt(size(w02_tri))
double precision :: Pb(3,size(w0_tri)), wb(size(w0_tri))

double precision :: T_tri_coord(3,3), B_tri_coord(3,3), bt(3,3)
double precision :: rf(3), T_nvec(3), B_nvec(3) 
double complex :: intN, intN2(3,3)

intN2(:,:) = dcmplx(0,0)

do face_t = 1,4

  T_tri_coord = tetra_face(T_coord,face_t)
  T_nvec = tri_n_vectors(T_tri_coord)
  call linmap_tri(Pt, wt, T_tri_coord, P02_tri, W02_tri)
  
  do face_b = 1,4
     B_tri_coord = tetra_face(B_coord,face_b)
     B_nvec = tri_n_vectors(B_tri_coord)
     call linmap_tri(Pb, wb, B_tri_coord, P0_tri, W0_tri)
   

     bt(1,1) = T_nvec(3)*B_nvec(3) + (-T_nvec(2)) * (-B_nvec(2)) 
     bt(2,1) = T_nvec(1)* (-B_nvec(2)) 
     bt(3,1) = -T_nvec(1) * (B_nvec(3)) 

     bt(1,2) = -T_nvec(2)*B_nvec(1)
     bt(2,2) = T_nvec(3) * B_nvec(3) + T_nvec(1) * B_nvec(1)
     bt(3,2) = -T_nvec(1)*B_nvec(3)

     bt(1,3) = T_nvec(3) * (-B_nvec(1))
     bt(2,3) = -T_nvec(3) * (B_nvec(2))
     bt(3,3) = T_nvec(2) * B_nvec(2) + T_nvec(1) * B_nvec(1)
    
     do i1 = 1, size(wt)
           
        rf = Pt(:,i1)
        intN = singularity_subtraction_int_S_G(rf, B_tri_coord, B_nvec, mesh%k, Pb, wb)
        intN2 = intN2 + bt * intN * wt(i1)
     end do
   
  end do

end do

end subroutine integrate_nx_dV_dV_G




end module integrals

