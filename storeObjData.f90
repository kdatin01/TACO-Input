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

subroutine str2int(str, int, stat)
	implicit none
	character(len=*), intent(in) :: str
	integer, intent(out) :: int
	integer, intent(out) :: stat
	
	read(str, *, iostat=stat) int
end subroutine str2int


subroutine getVerticesTiles(unitNum, filename, nvertices, ntiles, vertices, tiles)
	character(len=20), intent(in) :: filename
	integer, intent(in) :: nvertices, ntiles, unitNum
	real*8, dimension(1:nvertices,1:3), intent(out) :: vertices
	integer, dimension(1:ntiles,1:3), intent(out) :: tiles
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
	Point1(2) = vertices(xRowNum, 2)
	Point1(3) = vertices(xRowNum, 3)
	
	Point2(1) = vertices(yRowNum, 1)
	Point2(2) = vertices(yRowNum, 2)
	Point2(3) = vertices(yRowNum, 3)
	
	Point3(1) = vertices(zRowNum, 1)
	Point3(2) = vertices(zRowNum, 2)
	Point3(3) = vertices(zRowNum, 3) 
end subroutine getTileVertices

subroutine getCentroids(vertices, tiles, centroids, ntiles, nvertices)
	use calculationsMod
	integer, intent(in) :: ntiles, nvertices
	real*8, dimension(1:nvertices,1:3), intent(in) :: vertices
	integer, dimension(1:ntiles,1:3), intent(in) :: tiles
	
	real*8, dimension(1:ntiles, 1:3), intent(out) :: centroids
	real*8, dimension(1:3) :: centroid
	
	integer :: i
	integer, dimension(1:3) :: tileSet
	real*8, dimension(1:3) :: Point1, Point2, Point3
	
	do i = 1, ntiles
		tileSet(1) = tiles(i,1)
		tileSet(2) = tiles(i,2)
		tileSet(3) = tiles(i,3)
		
		call getTileVertices(vertices, tileSet, Point1, Point2, Point3, nvertices)
		call calcCentroids(Point1, Point2, Point3, centroid)
		centroids(i, 1) = centroid(1)
		centroids(i, 2) = centroid(2)
		centroids(i, 3) = centroid(3)
	end do
end subroutine getCentroids

subroutine getNormals(vertices, tiles, normals, ntiles, nvertices)
	use calculationsMod
	
	integer, intent(in) :: ntiles, nvertices
	real*8, dimension(1:nvertices,1:3), intent(in) :: vertices
	integer, dimension(1:ntiles,1:3), intent(in) :: tiles
	
	real*8, dimension(1:ntiles, 1:3), intent(out) :: normals
	real*8, dimension(1:3) :: normal
	
	integer :: i
	integer, dimension(1:3) :: tileSet
	real*8, dimension(1:3) :: Point1, Point2, Point3
	
	do i = 1, ntiles
		tileSet(1) = tiles(i,1)
		tileSet(2) = tiles(i,2)
		tileSet(3) = tiles(i,3)
		
		call getTileVertices(vertices, tileSet, Point1, Point2, Point3, nvertices)
		print*, Point1
		print*, Point2
		print*, Point3
		call calcUnitNormals(Point1, Point2, Point3, normal)
		
		normals(i,1) = normal(1)
		normals(i,2) = normal(2)
		normals(i,3) = normal(3)
		
	end do
end subroutine getNormals

subroutine getAreas(vertices, tiles, areas, ntiles, nvertices)
	use calculationsMod
	
	integer, intent(in) :: ntiles, nvertices
	real*8, dimension(1:nvertices,1:3), intent(in) :: vertices
	integer, dimension(1:ntiles,1:3), intent(in) :: tiles
	
	real*8, dimension(1:ntiles), intent(out) :: areas
	real*8 :: area
	
	integer :: i
	integer, dimension(1:3) :: tileSet
	real*8, dimension(1:3) :: Point1, Point2, Point3
	
	do i = 1, ntiles
		tileSet(1) = tiles(i,1)
		tileSet(2) = tiles(i,2)
		tileSet(3) = tiles(i,3)
		
		call getTileVertices(vertices, tileSet, Point1, Point2, Point3, nvertices)
		print*, Point1
		print*, Point2
		print*, Point3
		call calcTileArea(Point1, Point2, Point3, area)
		print*, area
		
		areas(i) = area
	end do
end subroutine getAreas

program storeObjData

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
	integer :: nvertices, ntiles, stat, i
			
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
		
		call str2int(vCount, nvertices, stat)
		!do i=1,3
		!	if (stat(i) /= 0) then
		!		print*, "Failed to convert nverticesStr to int"
		!	end if
		!end do
	
		call str2int(fCount, ntiles, stat)
		!do i=1,3
		!	if (stat(i) /= 0) then
		!		print*, "Failed to convert nverticesStr to int"
		!	end if
		!end do
		
		allocate(vertices(1:nvertices,1:3))
		allocate(tiles(1:ntiles,1:3))
		allocate(centroids(1:ntiles, 1:3))
		allocate(normals(1:ntiles, 1:3))
		allocate(areas(1:ntiles))
		
		call getVerticesTiles(20, objFile, nvertices, ntiles, vertices, tiles)
		call getCentroids(vertices, tiles, centroids, ntiles, nvertices)
		call getNormals(vertices, tiles, normals, ntiles, nvertices)
		call getAreas(vertices, tiles, areas, ntiles, nvertices)
		print*, vertices
		print*, tiles
		print*, objFile
	else
		print*, "failure somewhere"
	end if
		
end program storeObjData