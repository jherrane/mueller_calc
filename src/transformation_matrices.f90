module transformation_matrices
use common
use mie

implicit none

contains

subroutine vswf2constant(mesh, k, Nmax, mat)
type (mesh_struct) :: mesh
integer :: Nmax, inout
complex(dp), dimension(3,(Nmax+1)**2-1) :: MM_nm, NN_nm
complex(dp), dimension(3*mesh%N_tet, 2*((Nmax+1)**2 -1)) :: mat
complex(dp) :: F(3), G(3), k
integer :: tet, ind1(3), nm
real(dp) :: T_coord(3,4), r(3), vol


nm = (Nmax+1)**2 -1

do tet = 1, mesh%N_tet

   T_coord = mesh%coord(:,mesh%etopol(:,tet))
   r = (T_coord(:,1) +T_coord(:,2) + T_coord(:,3) + T_coord(:,4))/4.0

   call calc_MN(MM_nm, NN_nm, Nmax, k, r, 0)

   vol= tetra_volume(T_coord)

   ind1 = [3*(tet-1)+1,3*(tet-1)+2,3*(tet-1)+3]

   mat(ind1,1:nm) = MM_nm * sqrt(vol)
   mat(ind1,nm+1:2*nm) = NN_nm * sqrt(vol)

end do


end subroutine vswf2constant

end module transformation_matrices
