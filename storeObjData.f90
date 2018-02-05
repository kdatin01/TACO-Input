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
	integer, dimension(1:npairs, 1:2), intent(out) :: pairs
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
					if (l <= npairs) then
						yes=yes.or.(pair(1).eq.pairs(l,1) .and.pair(2).eq.pairs(l,2))
					end if
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
	use calculationsMod
	use asteroid_struct
	implicit none
	integer, dimension(1:npairs, 1:2), intent(in) :: pairs
	real*8, dimension(1:nhorizon,1:ntiles), intent(out) :: horizonmap
	integer :: i
	do i=1,ntiles
		call tilehorizon(i,pairs,horizonmap(1,i))
	enddo
end subroutine getTileHorizon

subroutine getAxisRatios
	use calculationsMod
	use asteroid_struct
	implicit none
	real*8 :: longAxis, middleAxis,shortAxis
	real*8 :: midToLong, shortToLong
	call calcAxes(longAxis, middleAxis, shortAxis, nvertices, asteroid%vertices)
end subroutine getAxisRatios

!subroutine getVisPairs()
!	use calculationsMod
!	use asteroid_struct
!	implicit none
!	integer:: nvispairs
!	integer, dimension(:), allocatable:: vispairs
	
!	call calcVisTiles(ntiles, nhorizon, asteroid%horizonmap, asteroid%normals, asteroid%centroids, nvispairs, vispairs, &
!	asteroid%vertices, asteroid%tiles)
	
!end subroutine getVisPairs



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
	integer, dimension(1:npairs, 1:2) :: pairs
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
		!print*, nvertices
		!print*, ntiles
		!print*, npairs
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
		call getTileHorizon(pairs, horizonmap)
		asteroid%horizonmap=horizonmap
		call getAxisRatios
		!print*, vertices
		!print*, tiles
		!print*, objFile
	else
		print*, "failure somewhere"
	end if
		
end program storeObjData