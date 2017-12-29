subroutine makeHeader(nvertices, ntiles, hfilenameOrig, hfilenameDest)
	!Write nvertices and ntiles to taco_parameters.h file for TACO run
	character(len=*), intent(in) :: nvertices, ntiles
	character(len=*), intent(in) :: hfilenameOrig, hfilenameDest
	
	character(len=1000) :: line
	
	open(unit=10, file=hfilenameOrig, action="read", status="old")
	open(unit=11, file=hfilenameDest, action="write", status="replace")
	ReadLines: do
		read(10, '(A)', end=20) line
		if (index(line, 'nvertices') /= 0) then
			write(11,*) "integer, parameter :: nvertices = ", nvertices
		else if (index(line, 'npairs') == 0 .and. index(line, 'ntiles') /= 0) then
			write(11,*) "integer, parameter :: ntiles = ", ntiles
		else
			write(11,*) line
		endif
	end do ReadLines
	20 continue
	close(10)
	close(11)
end subroutine makeHeader

!subroutine makeAsteroidInputFile(fileName, nvertices, ntiles, vertices, normals, centroids, horizonmap, areas, tiles, totalmass, principalmoments, nvispairs)
!	real*8, dimension(0:2,0:nvertices-1), intent(in) :: vertices
!   	real*8, dimension(0:2,0:ntiles-1), intent(in) :: normals,centroids
!   	real*8, dimension(0:nhorizon-1,0:ntiles-1), intent(in) :: horizonmap
!   	real*8, dimension(0:ntiles-1), intent(in) :: areas
!   	real*8, dimension(0:2), intent(in) :: principalmoments
!   	real*8, intent(in) :: totalmass
!   	integer, dimension(0:2,0:ntiles-1), intent(in) :: tiles
!   	integer, dimension(:), allocatable, intent(in) :: vispairs
!   	integer, intent(in) :: nvertices,ntiles,nhorizon,nvispairs
   		
!end subroutine makeAsteroidInputFile

subroutine str2int(str, int, status1)
	implicit none
	character(len=*), intent(in) :: str
	integer, intent(out) :: int
	integer, intent(out) :: status1
	
	read(str, *, iostat=status1) int
end subroutine str2int


subroutine getVerticesTiles(unitNum, filename, vertices, tiles)
	use asteroid_struct
	character(len=20), intent(in) :: filename
	integer, intent(in) :: unitNum
	real*8, dimension(1:asteroid%nvertices,1:3), intent(out) :: vertices
	integer, dimension(1:asteroid%ntiles,1:3), intent(out) :: tiles
	character(len=600) :: line
	character(len=1) :: dataType
	
	integer :: vcount=1, tcount=1
	
	open(unitNum, file=filename, status='old', action='read')
	ReadLines: do
		read (unitNum, '(A)', end=21) line
		dataType = line(1:1)
		if (dataType /= "#" .and. dataType /= '') then
			if (dataType == "v") then
				read(line,*, end=21) dataType, vertices(vcount, 1), vertices(vcount,2), vertices(vcount,3)
				!print*, vcount
				!print*, vertices(vcount, 1)
				!print*, vertices(vcount, 2)
				!print*, vertices(vcount, 3)
				vcount = vcount + 1
			else if (dataType == "f") then
				read(line, *, end=21) dataType, tiles(tcount, 1), tiles(tcount, 2), tiles(tcount, 3)
				tcount = tcount + 1
			end if
		end if
	end do ReadLines
	21 continue
	close(unitNum)	
end subroutine getVerticesTiles

subroutine printMultiDimArrayInt(array, rows, cols)
	integer, intent(in) :: rows, cols
	integer, dimension(0:cols, 0:rows - 1), intent(in) :: array
	integer :: i = 1, j = 1
	
	do i = 1, rows
		do j = 1, cols
			print*, array(i, j)
		end do
	end do	
end subroutine printMultiDimArrayInt

subroutine printMultiDimArrayReal(array, rows, cols)
	integer, intent(in) :: rows, cols
	real*8, dimension(1:rows, 1:cols), intent(in) :: array
	integer :: i, j
	
	do i = 1, rows
		do j = 1, cols
			print*, i
			print*, j
			print*, array(i, j)
		end do
	end do	
end subroutine printMultiDimArrayReal

subroutine getTileVertices(vertices, tileSet, Point1, Point2, Point3, nvertices)
	integer, intent(in) :: nvertices
	real*8, dimension(1:nvertices,1:3), intent(in) :: vertices
	integer, dimension(1:3), intent(in) :: tileSet
	
	real*8, dimension(1:3), intent(out) :: Point1, Point2, Point3
	
	integer :: xRowNum, yRowNum, zRowNum
	
	xRowNum = tileSet(1)
	yRowNum = tileSet(2)
	zRowNum = tileSet(3)
	
	Point1(1) = vertices(xRowNum, 1)
	Point1(2) = vertices(yRowNum, 1)
	Point1(3) = vertices(zRowNum, 1)
	
	Point2(1) = vertices(xRowNum, 2)
	Point2(2) = vertices(yRowNum, 2)
	Point2(3) = vertices(zRowNum, 2)
	
	Point3(1) = vertices(xRowNum, 3)
	Point3(2) = vertices(yRowNum, 3)
	Point3(3) = vertices(zRowNum, 3) 
end subroutine getTileVertices

subroutine getCentroids(centroids)
	use calculationsMod
	use asteroid_struct
	
	real*8, dimension(1:asteroid%ntiles, 1:3), intent(out) :: centroids
	real*8, dimension(1:3) :: centroid
	
	integer :: i
	integer, dimension(1:3) :: tileSet
	real*8, dimension(1:3) :: Point1, Point2, Point3
	
	do i = 1, asteroid%ntiles
		tileSet(1) = asteroid%tiles(i,1)
		tileSet(2) = asteroid%tiles(i,2)
		tileSet(3) = asteroid%tiles(i,3)
		
		call getTileVertices(asteroid%vertices, tileSet, Point1, Point2, Point3, asteroid%nvertices)
		call calcCentroids(Point1, Point2, Point3, centroid)
		centroids(i, 1:3) = centroid
		!centroids(i, 2) = centroid(2)
		!centroids(i, 3) = centroid(3)
	end do
end subroutine getCentroids

subroutine getConnectedVertexPairs(pairs)
	use asteroid_struct
	implicit none
	
	integer, dimension (1) :: maxtile,mintile
	integer, dimension(1:npairs+1, 1:2), intent(out) :: pairs
	integer, dimension(1:3) :: tile
	integer :: ii, paircount, l, k, j
	integer, dimension(1:2) :: pair
	logical :: yes
	! First use the list of tiles to get the unique list of connected vertex pairs
	paircount=1
	!print*, nvertices
	!print*, ntiles

	do ii=1,asteroid%ntiles ! do1
	! Sort each tile's vertices in ascending order
		tile=asteroid%tiles(ii,1:3)
		maxtile=maxloc(tile)
		mintile=minloc(tile)
		tile=tile((/mintile(1),6-mintile(1)-maxtile(1),maxtile(1)/))
	! Look at each vertex pair in the tile and add it to the list of
	! pairs if it isn't already there
		do j=1,2 !do3
			do k=j+1,3 !do4
				pair=(/tile(j),tile(k)/)
				yes=.false.
				do l=1,paircount !do5
					yes=yes.or.(pair(1).eq.pairs(l,1) .and.pair(2).eq.pairs(l,2))
				enddo !do5
				if (.not.yes) then
					pairs(paircount,1:2)=pair
					paircount=paircount+1
				endif
			enddo !do4
		enddo !do3
	enddo !do1
end subroutine getConnectedVertexPairs

subroutine getTileHorizon(pairs, horizonmap)
	!use calculationsMod
	use asteroid_struct
	implicit none
	integer, dimension(1:npairs+1, 1:2), intent(in) :: pairs
	real*8, dimension(1:nhorizon,1:ntiles), intent(out) :: horizonmap
	integer :: i
	do i=1,ntiles
		call tilehorizon(i,pairs,horizonmap(1,i), ntiles, npairs, nvertices, nhorizon, &
		pi, asteroid%normals, asteroid%vertices, asteroid%centroids)
		print*, horizonmap(1,i)
	enddo
end subroutine getTileHorizon

subroutine test(a1, b1, crossit)
	!real, dimension(1:3), intent(in) :: a, b
	!real*8, dimension(1:3), intent(in) :: a, b
	real, dimension(1:3), intent(in):: a1, b1
	real, dimension(1:3), intent(out) :: crossit
	print*, a1
	print*, b1
	crossit = (/0.2136159,0.774723606,0.00000/)
	print*, crossit
	!crossit = cshift(a1,shift=1)*cshift(b1,shift=-1)-cshift(b1,shift=1)*cshift(a1,shift=-1)	
end subroutine test

subroutine cross(a,b, crossed)
	real, dimension(1:3), intent(in) :: a, b
	real, dimension(1:3), intent(out) :: crossed
	
	crossed(1) = a(2)*b(3) - a(3)*b(2)
	crossed(2) = a(3) * b(1) - a(1) * b(3)
	crossed(3) = a(1) * b(2) - a(2) * b(1)
end subroutine cross

subroutine tilehorizon(tilenumber,pairs,tilehorizonmap, ntiles, npairs, nvertices, nhorizon, pi, normals, vertices, centroids)
	use calculationsMod
	implicit none

	integer, intent(in) :: tilenumber, ntiles, npairs, nvertices, nhorizon
	real*8, intent(in) :: pi
	real*8, dimension(1:ntiles, 1:3), intent(in):: normals, centroids
	real, dimension(1:ntiles, 1:3) :: normals4
	real*8, dimension(1:nvertices,1:3), intent(in) :: vertices
	integer, dimension(1:npairs+1, 1:2), intent(in) :: pairs
	real*8, dimension(1:ntiles), intent(inout) :: tilehorizonmap
	real*8, dimension(1:nvertices) :: elevation,azimuth
	real*8, dimension(1:3) :: dr,dr2,refvector,direction,direct1, out1, refvec
	real*8, dimension (1) :: maxseg
	real*8 :: sense,az
	real, dimension(1:3) :: refvec1
	real*8, dimension(1:npairs+1) :: el0,el1,az0,az1,segmentelev,temp
	integer, dimension(1:npairs+1,1:2) :: usepairs
	integer, dimension(1:nvertices) :: tileindexnumbers,blank1
	integer, dimension(1:nvertices) :: abovehorizon
	integer, dimension(1:3) :: tile
	integer, dimension(1:2) :: pair
	integer, dimension (1) :: maxtile,mintile
	integer :: nwh,nabovehorizon,nusepairs,ii,l,j,k,i,paircount
	logical, dimension(1:nvertices) :: mask
	logical :: yes

	do i=1,nvertices
		tileindexnumbers(i)=i
		blank1(i)=0
	end do

	i=tilenumber
	print*, i

	! Get azimuth relative to north and elevation of all tile vertices as seen
	! from the centroid of the ith tile.
	print*, normals(i, 1:3)
	normals4(1:ntiles,1:3) = real(normals(1:ntiles,1:3))
	print*, normals4(i, 1:3)
	call test(normals4(i, 1:3),(/0.00,0.00,1.00/),refvec1)
	!call test((/-0.77472360658401374,0.21361593920594180,0.59512315021230833/),(/0.00,0.00,1.00/),refvec1)
	refvec = (/0.2136159,0.774723606,0.00000/)
	refvector=refvec/sqrt(sum(refvec**2))
	do j=1,nvertices
		dr2=vertices(j,1:3)-centroids(i,1:3)
		dr=dr2/sqrt(sum(dr2**2))
		call cross_product(normals(i,1:3),dr, direct1)
		direction=direct1/sqrt(sum(direct1**2))
		azimuth(j)=acos(dot_product(refvector,direction))
		call cross_product(refvector,direction, out1)
		sense=dot_product(normals(i,1:3), out1)
		azimuth(j)=azimuth(j)*sense/abs(sense)
		elevation(j)=pi/2.d0 &
			-acos(dot_product(dr,normals(i,1:3)))
	enddo

	! Now compute the horizon map.
	! We care only about those vertices that are above the horizon.

	usepairs=1
	mask=elevation > 1.d-6
	nabovehorizon=count(mask)
	tilehorizonmap(1:nhorizon)=0.d0
	if (nabovehorizon .gt. 0) then  ! Now handle topography above horizon

		abovehorizon=pack(tileindexnumbers,mask,blank1)
		paircount=1
		do j=1,npairs+1
			yes=.false.
			do l=1,nabovehorizon
				yes=yes.or.(pairs(j,1).eq.abovehorizon(l) .or.pairs(j,2).eq.abovehorizon(l))
			enddo
			if (yes) then
				usepairs(paircount,1:2)=(/pairs(j,1),pairs(j,2)/)
				paircount=paircount+1
			endif
		enddo
		paircount=paircount-1
		nusepairs=paircount

		! At each azimuth, find the highest tile edge.

		el0=-99.d0
		el1=-99.d0
		az0=0.d0
		az1=0.d0
		temp=0.d0
		el0(1:nusepairs)=elevation(usepairs(1:nusepairs,1))
		el1(1:nusepairs)=elevation(usepairs(1:nusepairs,2))
		az0(1:nusepairs)=azimuth(usepairs(1:nusepairs,1))
		az1(1:nusepairs)=azimuth(usepairs(1:nusepairs,2))
		where (az1.lt.az0)
			temp=az1
			az1=az0
			az0=temp
			temp=el1
			el1=el0
			el0=temp
		endwhere
		do j=1,nhorizon
			segmentelev=0.d0
			az=(dble(j)+0.5d0)*2.d0*pi/dble(nhorizon)-pi
			where (el0.gt.-99d0.and.az1-az0 .le. pi .and. (az0.le.az.and.az.le.az1))
				segmentelev=el0+(az-az0)*(el1-el0)/(az1-az0)
			else where (el0.gt.-99.d0.and.az1-az0.gt.pi.and.az.le.az0)
				az1=az1-2.d0*pi
				segmentelev=el1+(az-az1)*(el0-el1)/(az0-az1)
			else where (el0.gt.-99.d0.and.az1-az0.gt.pi.and.az.ge.az1)
				az0=az0+2.*pi
				segmentelev=el1+(az-az1)*(el0-el1)/(az0-az1)
			end where
			maxseg=maxval(segmentelev)
			tilehorizonmap(j)=maxseg(1)
			print*, tilehorizonmap(j)
		enddo		

	endif  ! Done handling topography above horizon

	return
end subroutine tilehorizon

subroutine getNormals(normals)
	use calculationsMod
	use asteroid_struct
	
	real*8, dimension(1:asteroid%ntiles, 1:3), intent(out) :: normals
	real*8, dimension(1:3) :: normal
	
	integer :: i
	integer, dimension(1:3) :: tileSet
	real*8, dimension(1:3) :: Point1, Point2, Point3
	
	do i = 1, asteroid%ntiles
		tileSet(1) = asteroid%tiles(i,1)
		tileSet(2) = asteroid%tiles(i,2)
		tileSet(3) = asteroid%tiles(i,3)
		
		call getTileVertices(asteroid%vertices, tileSet, Point1, Point2, Point3, asteroid%nvertices)
		call calcUnitNormals(Point1, Point2, Point3, normal)
		
		normals(i,1) = normal(1)
		normals(i,2) = normal(2)
		normals(i,3) = normal(3)
		
	end do
end subroutine getNormals

subroutine getAreas(areas)
	use calculationsMod
	use asteroid_struct
	
	real*8, dimension(1:asteroid%ntiles), intent(out) :: areas
	real*8 :: area
	
	integer :: i
	integer, dimension(1:3) :: tileSet
	real*8, dimension(1:3) :: Point1, Point2, Point3
	
	do i = 1, asteroid%ntiles
		tileSet(1) = asteroid%tiles(i,1)
		tileSet(2) = asteroid%tiles(i,2)
		tileSet(3) = asteroid%tiles(i,3)
		
		call getTileVertices(asteroid%vertices, tileSet, Point1, Point2, Point3, asteroid%nvertices)
		call calcTileArea(Point1, Point2, Point3, area)
		
		areas(i) = area
	end do
end subroutine getAreas

subroutine getTotalMass(totalMass)
	use calculationsMod
	use asteroid_struct
	
	real*8, intent(out) :: totalMass
	
	call calcTotalMass(asteroid%centroids, asteroid%normals, asteroid%areas, asteroid%ntiles, totalMass)
	print*, "TotalMass"
	print*, totalMass
end subroutine getTotalMass

!subroutine makeAsteroidFile(outFile, vertices, normals, centroids, areas, tiles, nvertices, ntiles)
subroutine makeAsteroidFile(outFile)
	use asteroid_struct
	!real*8, dimension(1:nvertices), intent(in) :: vertices
	!real*8, dimension(1:3,1:ntiles), intent(in) :: normals, centroids
	!real*8, dimension(1:ntiles), intent(in) :: areas
	!integer, dimension(1:3,1:ntiles), intent(in) :: tiles
	!integer, intent(in) :: nvertices, ntiles
	character*80, intent(in) :: outFile
	
	!variables that are not currently populated by storeObjData, but are necessary for asteroid_structure in TACO.
	!integer :: nhorizon = 50 !dummy variable, value gotten from AAAtaco_paramters.h file
	real*8, dimension(1:nhorizon,1:ntiles) :: horizonmap
	real*8, dimension(0:2) :: principalmoments
	real*8 :: totalmass
	
	

	open(3,file=outFile,form='unformatted',status='new')

	write(3)asteroid%nvertices 	!got it
	write(3)asteroid%ntiles		!got it	
	write(3)asteroid%nhorizon	!got it
	write(3)asteroid%vertices	!got it
	write(3)asteroid%normals	!got it
	write(3)asteroid%centroids	!got it
	write(3)asteroid%horizonmap !got it
	write(3)asteroid%areas		!got it
	write(3)asteroid%tiles		!got it
	write(3)asteroid%totalmass	!got it
	write(3)asteroid%principalmoments
	write(3)asteroid%nvispairs
	if (asteroid%nvispairs.gt.0) write(3)asteroid%vispairs

	close(3)

	return
end subroutine makeAsteroidFile

program storeObjData
	use asteroid_struct
	implicit none
	
	integer :: narg, cptArg, index
	character(len=20) :: name, objFile, vCount='0', fCount='0'
	logical :: fileExist
	logical :: lookForObjFile = .FALSE.
	logical :: vertex = .FALSE.
	logical :: facet = .FALSE.
	character(len=10) :: vrefstr
	character(len=10) :: trefstr
	character(len=100) :: headerIn, headerOut
	real*8, dimension(:,:), allocatable :: vertices !rank 1
	integer, dimension(:,:), allocatable :: tiles !rank 2
	real*8, dimension(:,:), allocatable :: centroids !rank 3
	real*8, dimension(:,:), allocatable :: normals !rank 4
	real*8, dimension(:), allocatable :: areas !rank 5
	integer, dimension(1:npairs+1, 1:2) :: pairs
	real*8, dimension(1:nhorizon,1:ntiles) :: horizonmap
	real*8 :: totalMass
	!integer :: nvertices, ntiles, status1, i
	integer :: status1, i
			
	headerIn = '/mnt/c/Users/kdatin01/Documents/TACOpack/taco_src/taco_parameters.h'
	headerOut = '/mnt/c/Users/kdatin01/Documents/UMD/KellieTACOpack/taco_src/taco_parameters.h'
	
	!Check if arguments are found
	narg = command_argument_count()
	objFile = "NA"
	
	if (narg>0) then
		!loop across options
		do cptArg=1, narg
			call get_command_argument(cptArg, name)
			!print*,name
			select case(adjustl(name))
				!First known args
				case("--obj")
					lookForObjFile=.TRUE. !change logical value
				case ("--v")
					vertex = .TRUE.
				case ("--f")
					facet = .TRUE.
				case default
				
				if (vertex) then
					vCount = adjustl(name)
					!print*, vCount
					vertex = .FALSE.
				else if (facet) then
					fCount = adjustl(name)
					!print*, fCount
					facet = .FALSE.
				else if(lookForObjFile) then
					objFile = adjustl(name) !assign value to objFile
					!print*, objFile
					inquire(file=objFile, exist=fileExist) !check if objFile exists
					if(.not.fileExist) then
						write(*,*) 'file ',objFile, ' not found'
						stop
					lookForObjFile = .FALSE.
					end if
				else
					write(*,*)"Option ",adjustl(name)," unknown"
				end if
			end select
		end do
	end if
	
	if (vCount /= '0' .and. fCount /= '0' .and. objFile /= "NA") then
		!Write nvertices and ntiles to taco_parameters.h file for TACO run
		call makeHeader(vCount, fCount, headerIn, headerOut)
		
		!call str2int(vCount, nvertices, status1)
		!do i=1,3
		!	if (stat(i) /= 0) then
		!		print*, "Failed to convert nverticesStr to int"
		!	end if
		!end do
	
		!call str2int(fCount, ntiles, status1)
		!do i=1,3
		!	if (stat(i) /= 0) then
		!		print*, "Failed to convert nverticesStr to int"
		!	end if
		!end do
		print*, nvertices
		print*, ntiles
		print*, npairs+1
		allocate(vertices(1:nvertices,1:3))
		allocate(tiles(1:ntiles,1:3))
		allocate(centroids(1:ntiles, 1:3))
		allocate(normals(1:ntiles, 1:3))
		allocate(areas(1:ntiles))
		
		asteroid%nvertices = nvertices
		asteroid%ntiles = ntiles
		!print*, asteroid%nvertices
		!print*, asteroid%ntiles
		
		call getVerticesTiles(20, objFile, vertices, tiles)
		asteroid%vertices = vertices
		asteroid%tiles = tiles
		call getCentroids(centroids)
		asteroid%centroids = centroids
		call getNormals(normals)
		asteroid%normals = normals
		call getAreas(areas)
		asteroid%areas = areas
		call getTotalMass(totalMass)
		asteroid%totalmass = totalMass
		call getConnectedVertexPairs(pairs)
		print*, pairs
		call getTileHorizon(pairs, horizonmap)
		print*, horizonmap
		asteroid%horizonmap=horizonmap
		print*, vertices
		print*, tiles
		print*, objFile
	else
		print*, "failure somewhere"
	end if
		
end program storeObjData