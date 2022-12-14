program msh2scrip

  use netcdf

  implicit none

  integer, parameter :: grid_rank = 1    ! mesh
  integer, parameter :: nvert = 3        ! triangles

  integer :: i,ie,n,id,rc,ncid,dim2(2),dim1(1)
  integer :: idimid,jdimid,kdimid

  integer, dimension(grid_rank) :: gdims

  character(len=20) :: vname

  character(len=200) :: dirsrc, dirout
  character(len=200) :: mshfname, scpfname
  character(len=200) :: fname
  character(len=120) :: chead

  integer :: unum,nn,ne,ecnt
  integer :: nelements, nnodes, nvalid
  integer :: i1,i2,i3,i4,i5,i6,node1,node2,node3

  real(kind=8) :: nlon, nlat, ndpt
  real(kind=8)   , allocatable :: nodeCoords(:,:)
  integer(kind=4), allocatable :: ownedNodes(:,:)
  real(kind=8)   , allocatable :: elemCornerCoordX(:,:)
  real(kind=8)   , allocatable :: elemCornerCoordY(:,:)
  real(kind=8)   , allocatable :: elemCoordX(:)
  real(kind=8)   , allocatable :: elemCoordY(:)
  integer(kind=4), allocatable :: mask(:)

  real(kind=8) :: xcnrs(3), ycnrs(3), xcen, ycen
  !-------------------------------------------------------------------

  dirsrc = '/scratch1/NCEPDEV/nems/Denise.Worthen/WORK/WaveIn_unstr/'
  !mshfname='globa_1deg.msh'
  mshfname='globa_1deg_no_land.msh'

  !dirsrc = '/scratch1/NCEPDEV/stmp2/Ali.Abdolali/Source/Aug7Dev/regtests/ww3_tp2.17/input/'  
  !mshfname = 'inlet.msh'

  dirout = './'
  fname = trim(dirsrc)//trim(mshfname)
  open(newunit=unum,file=trim(fname), form='formatted', status='old')

  read(unum,*)chead
  read(unum,*)chead
  read(unum,*)chead
  read(unum,*)chead
  read(unum,*)nnodes
  print *,'number nodes = ',nnodes

  ! read the node coords
  allocate(nodeCoords(1:2,1:nnodes))
  do nn = 1,nnodes
     read(unum,*)i1,nlon,nlat,ndpt
     !print *,i1,nlon,nlat,ndpt
     nodeCoords(1,nn) = nlon
     nodeCoords(2,nn) = nlat
  end do
  print '(a,2f8.3)','longitude range : ',minval(nodeCoords(1,:)),maxval(nodeCoords(1,:))

  read(unum,*)chead
  read(unum,*)chead
  read(unum,*)nelements
  print *,'number elements = ',nelements

  allocate(mask(1:nelements))
  allocate(ownedNodes(1:3,1:nelements))
  mask = 1
  ownedNodes = -1
  ecnt = 0
  do ne = 1,nelements
     read(unum,*)i1,i2,chead
     if (i2 /= 15) then
        backspace(unum)
        read(unum,*)i1,i2,i3,i4,i5,i6,node1,node2,node3
        ecnt = ecnt+1
        ownedNodes(1,ecnt) = node1
        ownedNodes(2,ecnt) = node2
        ownedNodes(3,ecnt) = node3
        !print *,i1,i2,ecnt,node1,node2,node3
     end if
  end do
  close(unum)
  print '(a,2i8)','node range ',minval(ownedNodes),maxval(ownedNodes)

  nvalid = ecnt
  print *,'number elements not on open boundary = ',nvalid
  allocate(elemCornerCoordX(1:3,1:nvalid))
  allocate(elemCornerCoordY(1:3,1:nvalid))
  allocate(elemCoordX(1:nvalid))
  allocate(elemCoordY(1:nvalid))
  elemCornerCoordX = -999.0
  elemCornerCoordY = -999.0
  elemCoordX = -999.0
  elemCoordY = -999.0

  do ne = 1,nvalid
     node1 = ownedNodes(1,ne)
     node2 = ownedNodes(2,ne)
     node3 = ownedNodes(3,ne)

     elemCornerCoordX(1,ne) = nodeCoords(1,node1)
     elemCornerCoordX(2,ne) = nodeCoords(1,node2)
     elemCornerCoordX(3,ne) = nodeCoords(1,node3)

     elemCornerCoordY(1,ne) = nodeCoords(2,node1)
     elemCornerCoordY(2,ne) = nodeCoords(2,node2)
     elemCornerCoordY(3,ne) = nodeCoords(2,node3)
  end do
  
  ! calculate center
  do ne = 1,nvalid
     do i = 1,3
        xcnrs(i) = elemCornerCoordX(i,ne)
        ycnrs(i) = elemCornerCoordY(i,ne)
     end do
     call calc_center(xcnrs,ycnrs,xcen,ycen)

     elemCoordX(ne) = xcen
     elemCoordY(ne) = ycen
  end do
  print '(a,2f8.2)','lon range ',minval(elemCoordX),maxval(elemCoordX)

  ! arbitrary points
  ne = nvalid/2
  print *,elemCornerCoordX(1,ne),elemCornerCoordY(1,ne),elemCornerCoordX(2,ne),elemCornerCoordY(2,ne),&
       elemCornerCoordX(3,ne),elemCornerCoordY(3,ne)
  print *,elemCoordX(ne),elemCoordY(ne)

  ne = nvalid-15
  print *,elemCornerCoordX(1,ne),elemCornerCoordY(1,ne),elemCornerCoordX(2,ne),elemCornerCoordY(2,ne),&
       elemCornerCoordX(3,ne),elemCornerCoordY(3,ne)
  print *,elemCoordX(ne),elemCoordY(ne)

  !-------------------------------------------------------------------

  fname = trim(dirout)//trim(mshfname)//'.SCRIP.nc'
  print *,'creating '//trim(fname)

  rc = nf90_create(trim(fname), nf90_64bit_offset, ncid)
  rc = nf90_def_dim(ncid, 'grid_size',     nvalid, idimid)
  rc = nf90_def_dim(ncid, 'grid_corners',   nvert, jdimid)
  rc = nf90_def_dim(ncid, 'grid_rank',  grid_rank, kdimid)

  !grid_dims
  dim1(:) = (/kdimid/)
  rc = nf90_def_var(ncid, 'grid_dims', nf90_int, dim1, id)
  ! mask
  dim1(:) = (/idimid/)
  rc = nf90_def_var(ncid, 'grid_imask', nf90_int, dim1, id)
  rc = nf90_put_att(ncid, id,     'units',      'unitless')

  ! centers
  dim1(:) =  (/idimid/)
  vname = 'grid_center_lon'
  rc = nf90_def_var(ncid, vname, nf90_double, dim1, id)
  rc = nf90_put_att(ncid, id,          'units',   'degrees')
  rc = nf90_put_att(ncid, id,  'standard_name', 'longitude')
  vname = 'grid_center_lat'
  rc = nf90_def_var(ncid, vname, nf90_double, dim1, id)
  rc = nf90_put_att(ncid, id,          'units',  'degrees')
  rc = nf90_put_att(ncid, id,  'standard_name', 'latitude')

  ! corners
  dim2(:) =  (/jdimid,idimid/)
  vname = 'grid_corner_lon'
  rc = nf90_def_var(ncid, vname, nf90_double, dim2, id)
  rc = nf90_put_att(ncid, id,          'units',   'degrees')
  rc = nf90_put_att(ncid, id,  'standard_name', 'longitude')
  vname = 'grid_corner_lat'
  rc = nf90_def_var(ncid, vname, nf90_double, dim2, id)
  rc = nf90_put_att(ncid, id,          'units',  'degrees')
  rc = nf90_put_att(ncid, id,  'standard_name', 'latitude')
  rc = nf90_enddef(ncid)

  rc = nf90_inq_varid(ncid,  'grid_dims',        id)
  rc = nf90_put_var(ncid,             id,     gdims)
  rc = nf90_inq_varid(ncid, 'grid_imask',        id)
  rc = nf90_put_var(ncid,             id,      mask(1:nvalid))

  rc = nf90_inq_varid(ncid,  'grid_center_lon',         id)
  rc = nf90_put_var(ncid,                   id, elemCoordX(1:nvalid))
  rc = nf90_inq_varid(ncid,  'grid_center_lat',         id)
  rc = nf90_put_var(ncid,                   id, elemCoordY(1:nvalid))

  rc = nf90_inq_varid(ncid,  'grid_corner_lon',               id)
  rc = nf90_put_var(ncid,                   id, elemCornerCoordX(:,1:nvalid))
  rc = nf90_inq_varid(ncid,  'grid_corner_lat',               id)
  rc = nf90_put_var(ncid,                   id, elemCornerCoordY(:,1:nvalid))
  rc = nf90_close(ncid)

end program msh2scrip

subroutine calc_center(xs,ys,xc,yc)

  real(kind=8), dimension(3), intent(in)  :: xs,ys
  real(kind=8)              , intent(out) :: xc,yc

  !local variables
  real(kind=8) :: xtmp(1:3)
  real(kind=8) :: diff12, diff13, diff23,maxdiff
  logical      :: xlon

  ! determine if nodes are discontinuous across 0E
  ! or 180E
  xlon = .false.
  diff12 = abs(xs(1) - xs(2))
  diff13 = abs(xs(1) - xs(3))
  diff23 = abs(xs(2) - xs(3))
  maxdiff = max(max(diff12,diff13),diff23)
  if (maxdiff .ge. 180.0) xlon = .true.

  ! find the triangle center
  yc = sum(ys)/3.0
  if (.not. xlon) then
     xc = sum(xs)/3.0
  else
     !print *,diff12,diff13,diff23
     xtmp = xs
     where(xtmp .ge. 180.0) xtmp = xtmp - 360.0
     !print *,xtmp
     xc = sum(xtmp)/3.0
     !print *,xc
     if (xc .lt. 0.0)xc = xc + 360.0
  end if

end subroutine calc_center
