! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@ This source file is part of the elsA software owned by Airbus, Safran and ONERA    @@
! @@ Any copy or distribution, total or partial, of this source file should strictly    @@
! @@ observe the commitments described in the Cooperation Agreement of 29 May 2015      @@
! @@ between Airbus, Safran and ONERA for elsA development.                             @@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine f_computemorse_fromelmts(MorseVtxIdx, MorseVtx, MorseNbf, ngon, nEC, nFac, nVtxFAF)
  !$ use OMP_LIB
  implicit none
  ! -----------------------------------------------------------------
  ! -> In
  integer,intent(in) :: nEC                                          ! Size of NGon
  integer,intent(in) :: nFac                                         ! Number of faces
  integer,intent(in) :: nVtxFAF  ! (ForAllFaces)                     ! Total number of vertex in NGon
  !
  integer,dimension(0:nEC-1),intent(in) :: ngon                      ! NGon node
  ! -> In/Out
  integer,dimension(0:nFac     ),intent(inout) :: MorseVtxIdx        ! Morse Index [0:nFac]
  integer,dimension(0:nVtxFAF-1),intent(inout) :: MorseVtx           ! Morse Vtx [0:nVtxFAF-1]
  integer,dimension(0:nFac-1   ),intent(inout) :: MorseNbf           ! Temporary array
  ! -> Local
  integer :: nf
  integer :: NbVtxOnFace                                             ! Nombre de vertex de la face courante
  integer :: indvtx
  integer :: cptmor
  integer :: cptngo
  ! -----------------------------------------------------------------
  ! -> Begin
  !
  ! > Init
  cptmor = 0
  cptngo = 0
  !
  ! I/ Loop on face and build MorseVtx, MorseVtxIdx is temporal in this Loop
  do nf=0,nFac-1

    ! >>> Read NGON Node : NbVtx [ Vtx1, Vtx2, Vtx3, ...]
    NbVtxOnFace = ngon(cptngo)

    ! >>> Store it !
    MorseNbf(nf) = NbVtxOnFace

    ! >>> Determine the number of vertex in NGon : It's the maximum of NGon numbering
    do indvtx=0, NbVtxOnFace-1
        MorseVtx(cptmor + indvtx) = ngon(cptngo + indvtx + 1)
    end do

    ! >>> Go to next elements in NGon
    cptngo  = cptngo  + NbVtxOnFace + 1

    ! >>> Go to next elements in Morse
    cptmor  = cptmor  + NbVtxOnFace

  end do
  ! > Manage Boundary
  ! MorseNbfNew[newf] = NbVtxOnFace;
  ! MorseNbfOld[nf  ] = NbVtxOnFace;


  ! II/ Make the full index array for morse and manage the last elements
  MorseVtxIdx(0) = 0
  do nf=1,nFac
    !
    MorseVtxIdx(nf) = MorseVtxIdx(nf-1) + MorseNbf(nf-1);
    !
  end do

  ! > Debug [Verbose]
  ! write(*,*) "### nFaces [calculated/FromNGON] : [", MorseVtxIdx(nFac), " / ", nVtxFAF, " ]"

  ! -> End
end subroutine f_computemorse_fromelmts



subroutine f_preparemorse_fromelmts(nVtxFAF, nVtx, ngon, nEC, nFac)
  !$ use OMP_LIB
  implicit none
  ! -----------------------------------------------------------------
  ! -> In
  integer,intent(in) :: nEC                                          ! Size of NGon
  integer,intent(in) :: nFac                                         ! Number of faces
  !
  integer,dimension(0:nEC-1),intent(in) :: ngon                      ! NGon node
  ! -> In/Out
  integer,intent(inout) :: nVtxFAF                                   ! Total number of vertex in NGon
  integer,intent(inout) :: nVtx                                      ! Total number of vertex
  ! -> Local
  integer :: nf
  integer :: NbVtxOnFace                                             ! Nombre de vertex de la face courante
  integer :: indvtx
  integer :: cptloc
  ! -----------------------------------------------------------------
  ! -> Begin
  ! > Init
  cptloc  = 0
  nVtx    = 0
  nVtxFAF = 0
  !
  do nf=0,nFac-1

    ! >>> Read NGON Node : NbVtx [ Vtx1, Vtx2, Vtx3, ...]
    NbVtxOnFace = ngon(cptloc)

    ! >>> Determine the number of vertex in NGon : It's the maximum of NGon numbering
    do indvtx=0, NbVtxOnFace-1
        nVtx = max(nVtx, ngon(cptloc+indvtx+1))
    end do

    ! >>> Stack in
    nVtxFAF = nVtxFAF + NbVtxOnFace

    ! >>> Go to next elements in NGon
    cptloc = cptloc + NbVtxOnFace + 1

  end do
  ! -> End
end subroutine f_preparemorse_fromelmts
