
module calculationsMod
!include "KellieTACOpack/taco_src/sethorizon.f90"

contains
subroutine calcCentroids(point1, point2, point3, centroid)
	implicit none
	real*8, dimension(1:3), intent(in) :: point1, point2, point3
	real*8, dimension(1:3), intent(out) :: centroid
	
	centroid=(point1+point2+point3)/3.

end subroutine calcCentroids

subroutine calcTileArea(Point1, Point2, Point3, Area)
	real*8, dimension(1:3), intent(in) :: Point1, Point2, Point3
	real*8, intent(out) :: Area
	
	real*8, dimension(1:3) :: crossed
	real*8 :: magnitude
	
	call cross_product(Point2-Point1,Point3-Point1, crossed)
	
	call calcMagnitude(crossed, magnitude)
	
	Area = 0.5 * magnitude

end subroutine calcTileArea

subroutine cross_product(vector1, vector2, crossit)
	real*8, dimension(3), intent(in) :: vector1, vector2
	real*8, dimension(3) :: crossit

	crossit=cshift(vector1,shift=1)*cshift(vector2,shift=-1)-cshift(vector2,shift=1)*cshift(vector1,shift=-1)
end subroutine cross_product

subroutine cross1(a,b, crossed)
	real*8, dimension(1:3), intent(in) :: a, b
	real*8, dimension(1:3), intent(out) :: crossed
	
	crossed(1) = a(2)*b(3) - a(3)*b(2)
	crossed(2) = a(3) * b(1) - a(1) * b(3)
	crossed(3) = a(1) * b(2) - a(2) * b(1)
	
return
end subroutine cross1

subroutine calcMagnitude(crossed, mag)
	real*8, dimension(1:3), intent(in) :: crossed
	real*8, intent(out) :: mag
	
	mag=sqrt(sum(crossed**2))
end subroutine calcMagnitude

subroutine calcNormals(Point1, Point2, Point3, normals)
	implicit none
	real*8, dimension(1:3), intent(in) :: Point1, Point2, Point3
	
	real*8, dimension(1:3), intent(out) :: normals
	
	call cross_product(Point2-Point1,Point3-Point1, normals)
end subroutine calcNormals

subroutine calcUnitNormals(Point1, Point2, Point3, unitNormals)
	implicit none
	real*8, dimension(1:3), intent(in) :: Point1, Point2, Point3
	
	real*8, dimension(1:3), intent(out) :: unitNormals
	
	real*8, dimension(1:3) :: normals
	real*8 :: mag
	
	call calcNormals(Point1, Point2, Point3, normals)
	
	call calcMagnitude(normals, mag)
	
	unitNormals = normals/mag
end subroutine calcUnitNormals

subroutine calcTotalMass(centroids, normals, areas, ntiles, totalMass)
	! Recalculate total mass
	integer, intent(in) :: ntiles
	real*8, dimension(1:ntiles, 1:3), intent(in) :: normals
	real*8, dimension(1:ntiles, 1:3), intent(in) :: centroids
	real*8, dimension(1:ntiles), intent(in) :: areas
	
	real*8, intent(out) :: totalMass
	real*8, dimension (1:ntiles) :: mass
	integer :: i

	mass=0.d0
	totalMass=0.d0

	do i=1,ntiles
		mass(i)=dot_product(centroids(i,1:3),normals(i,1:3)) * areas(i)/3.d0
		totalMass=totalMass+mass(i)
	enddo
end subroutine calcTotalMass

subroutine tilehorizon(tilenumber,pairs,tilehorizonmap)
	use asteroid_struct
	implicit none

	integer, intent(in) :: tilenumber
	!real*8 :: pi = 3.14159265358979d0
	integer, dimension(1:npairs, 1:2), intent(in) :: pairs
	real*8, dimension(1:ntiles), intent(inout) :: tilehorizonmap
	real*8, dimension(1:nvertices) :: elevation,azimuth
	real*8, dimension(1:3) :: dr,dr2,refvector,direction,direct1, out1, refvec
	real*8, dimension (1) :: maxseg
	real*8 :: sense,az
	real*8, dimension(1:npairs) :: el0,el1,az0,az1,segmentelev,temp
	integer, dimension(1:npairs,1:2) :: usepairs
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

	! Get azimuth relative to north and elevation of all tile vertices as seen
	! from the centroid of the ith tile.
	call cross_product(asteroid%normals(i,1:3),(/0.d0,0.d0,1.d0/), refvec)
	refvector=refvec/sqrt(sum(refvec**2))
	do j=1,nvertices
		dr2=asteroid%vertices(j,1:3)-asteroid%centroids(i,1:3)
		dr=dr2/sqrt(sum(dr2**2))
		call cross_product(asteroid%normals(i,1:3),dr, direct1)
		direction=direct1/sqrt(sum(direct1**2))
		azimuth(j)=acos(dot_product(refvector,direction))
		call cross_product(refvector,direction, out1)
		sense=dot_product(asteroid%normals(i,1:3), out1)
		azimuth(j)=azimuth(j)*sense/abs(sense)
		elevation(j)=pi/2.d0 &
			-acos(dot_product(dr,asteroid%normals(i,1:3)))
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
		do j=1,npairs
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
			elsewhere (el0.gt.-99.d0.and.az1-az0.gt.pi.and.az.le.az0)
				az1=az1-2.d0*pi
				segmentelev=el1+(az-az1)*(el0-el1)/(az0-az1)
			elsewhere (el0.gt.-99.d0.and.az1-az0.gt.pi.and.az.ge.az1)
				az0=az0+2.*pi
				segmentelev=el1+(az-az1)*(el0-el1)/(az0-az1)
			endwhere
			maxseg=maxval(segmentelev)
			tilehorizonmap(j)=maxseg(1)
		enddo		

	endif  ! Done handling topography above horizon

	return
end subroutine tilehorizon

subroutine calcAxes(longAxis, middleAxis, shortAxis, nvertices, vertices)
	real*8 :: longAxis, middleAxis, shortAxis
	integer, intent(in) :: nvertices
	real*8, dimension(1:nvertices,1:3), intent(in) :: vertices
	real*8, dimension(1:3) :: point1, point2, maxPoint1, maxPoint2
	real*8 :: distance, distemp
	integer :: i, j
	distance = 0.d0
	distemp = 0.d0
	do i=1,nvertices
		point1 = vertices(i, 1:3)
		do j=1, nvertices
			point2 = vertices(j, 1:3)
			distemp = sqrt((point2(1)-point1(1))**2+(point2(2)-point1(2))**2+(point2(3)-point1(3))**2)
			if (distemp.gt.distance) then
				maxPoint1 = point1
				maxPoint2 = point2
				distance = distemp
				print*, distance
			endif
		enddo
	enddo
	
	longAxis = distance
	print*, "Long Axis = "
	print*, distance
	
end subroutine calcAxes

!subroutine calcVisTiles(ntiles, nhorizon, horizonmap, normals, centroids, nvispairs, vispairs, vertices, tiles) 
!	implicit none
!	
!	integer, intent(in) :: ntiles, nhorizon
!	real*8, dimension(1:nhorizon,1:ntiles), intent(in) :: horizonmap
!	real*8, dimension(1:ntiles, 1:3), intent(in):: normals, centroids
!	integer, intent(out) :: nvispairs
!	integer, dimension(:), allocatable, intent(out) :: vispairs
!	real*8, dimension(1:3) :: norm0,norm1,vec,rij,rki,rkj,x,v0,v1,verti,vertj,vertk
!	real*8, dimension(1:3) :: normal
!	real*8 :: small,di,dj
!	integer, dimension(:,:), allocatable :: facingpairs,unblockedpairs
!	integer, dimension(:), allocatable :: col0,col1,so0,so1,wh
!	integer, dimension(1:ntiles) :: list0,list1
!	logical, dimension(:), allocatable :: notblocked
!	integer :: i,j,k,ii,nupairs,nfpairs,tile0,tile1,ipair,iblock,ivis
!	integer :: ind00,ind01,ind10,ind11,nwh,nwhi,currentcount
!	logical :: currentinlist
!	! -- sort definitions --
!	real*8, parameter :: aln2i = 1.442695042d0, tiny=1.d-9
!	integer, dimension(:), allocatable :: copy
!	integer :: t
!	integer :: lognb2,m,nn,tt,l ! ,k,j,i
!
!	small=1.d-2
!
!	! Find all pairs of tiles that face each other

!	print*,"Making pair list"
!	nfpairs=0
!	open(4,file='visible.temp')
!	do i=1,ntiles-1
	!	if (mod(i,100) .eq. 0) print*,i
!		if (any(horizonmap(1:nhorizon,i) .gt. 0.d0)) then
!			norm0=normals(i, 1:3)
!			do j=i+1,ntiles
!				vec=centroids(j, 1:3)-centroids(i,1:3)
!				vec=vec/sqrt(sum(vec**2))
!				if (dot_product(vec,norm0) .gt. small) then
!					norm1=normals(j,1:3)
!					if (dot_product(vec,norm1) .lt. -small) then
!						write(4,'(2(1x,i7))')i,j
!						nfpairs=nfpairs+1
!					endif	
!				endif
!			end do
!		endif
!	end do
!	close(4)

	! If there are no visible pairs then exit now.

!			print*,nfpairs," facing pairs"

!	if (nfpairs.eq.0) then
!		nvispairs=0
!		return
!	endif

!	if(allocated(facingpairs)) deallocate(facingpairs)
!	allocate(facingpairs(1:nfpairs,1:2))
!	open(4,file='visible.temp')
!	do i=1,nfpairs
!		read(4,*)facingpairs(i,1),facingpairs(i,2)
!	end do
!	close(4)

	! Exclude pairs whose view of each other is blocked by other tiles

!	print*,"Excluding blocked pairs"

!	if(allocated(notblocked)) deallocate(notblocked)
!	allocate(notblocked(1:nfpairs))
!	notblocked=.true.

	! As a preliminary, pre-sort and pre-index each column of the facing pairs list

!	if(allocated(col0)) deallocate(col0)
!	allocate(col0(1:nfpairs))
!	if(allocated(col1)) deallocate(col1)
!	allocate(col1(1:nfpairs))
!	if(allocated(so0)) deallocate(so0)
!	allocate(so0(1:nfpairs))
!	if(allocated(so1)) deallocate(so1)
!	allocate(so1(1:nfpairs))
!	if(allocated(wh)) deallocate(wh)
!	allocate(wh(1:nfpairs))
!	if(allocated(copy)) deallocate(copy)
!	allocate(copy(1:nfpairs))

!	col0=facingpairs(1:nfpairs, 1)

!	do i=1,nfpairs
!		so0(i)=i
!	end do
!	copy=col0
!	lognb2=int(log(dble(nfpairs))*aln2i+tiny)
!	m=nfpairs
!	do  nn=1,lognb2
!		m=m/2
!		k=nfpairs-m
!		do j=1,k
!			i=j
!3           continue
!			l=i+m
!			if (copy(l-1) .lt. copy(i-1)) then
!				t=copy(i-1)
!				copy(i-1)=copy(l-1)
!				copy(l-1)=t
!				tt=so0(i-1)
!				so0(i-1)=so0(l-1)
!				so0(l-1)=tt
!				i=i-m
!				if (i .ge. 1) goto 3
!			endif
!		end do
!	end do

!	col0=col0(so0)

!	col1=facingpairs(1:nfpairs, 1)

!	do i=1,nfpairs
!		so1(i)=i
!	end do
!	copy=col1
!	lognb2=int(log(dble(nfpairs))*aln2i+tiny)
!	m=nfpairs
!	do  nn=1,lognb2
!        m=m/2
!        k=nfpairs-m
!        do j=1,k
!			i=j
!4           continue
!            l=i+m
!            if (copy(l-1) .lt. copy(i-1)) then
!                t=copy(i-1)
!                copy(i-1)=copy(l-1)
!                copy(l-1)=t
!                tt=so1(i-1)
!                so1(i-1)=so1(l-1)
!                so1(l-1)=tt
!                i=i-m
!                if (i .ge. 1) goto 4
!            endif
!        end do
!	end do

!	col1=facingpairs(so1, 1)

!	list0=0
!	list1=0
!	tile0=1
!	tile1=1
!	do ii=1,nfpairs
!        if (col0(ii).ne.tile0) then
!			do while (tile0.ne.col0(ii))
!                list0(tile0)=ii
!                tile0=tile0+1
!            end do
!        endif
!        if (col1(ii).ne.tile1) then
!            do while (tile1.ne.col1(ii))
!                list1(tile1)=ii
!                tile1=tile1+1
!           end do
!        endif
!	end do

	! Loop through pairs looking for blockages

!	do ipair=1,nfpairs
!		if (mod(ipair,100000) .eq. 0) print*,'Done',ipair
!		i=facingpairs(ipair,1)
!		j=facingpairs(ipair,2)
!		rij=centroids(j,1:3)-centroids(i,1:3)
		! Find other facing pairs in which the first tile is either tile of this pair
!		nwhi=1
!		if (i.eq.0) then
!			ind00=1
!		else
!			ind00=list0(i)
!		endif
!		ind01=list0(i)
!		if (ind01.gt.ind00) then
!			nwhi=nwhi+(ind01-ind00)
!			wh(1:nwhi)=so0(ind00:ind01)
!		endif
!		nwh=nwhi
!		if (j.eq.0) then
!			ind00=1
!		else
!			ind00=list0(j)
!		endif
!		ind01=list0(j)
!		if (ind01.gt.ind00) then
!			nwh=nwh+(ind01-ind00)
!			wh(nwhi:nwh)=so0(ind00:ind01)
!		endif

		! Loop through the second tile of these pairs
!		if (nwh.gt.1) then
!			do iblock=1,nwh
!				k=facingpairs(wh(iblock),2)
!				if (.not.(i.eq.facingpairs(wh(iblock),1).and.j.eq.k)) then
!					normal=normals(k,1:3)
!					rki=centroids(i,1:3)-centroids(k,1:3)
!					rkj=centroids(j,1:3)-centroids(k,1:3)
					! Does segment rij pass through this tile?
!					! (a) Are centroids of tiles i and j on opposite sides of tile k?
!					di=dot_product(normal,rki)
!					dj=dot_product(normal,rkj)
!					if (di*dj.lt.0.d0) then
!						! (b) Find point where ray intersects the tile plane
!						x=centroids(k,1:3)+(di*rkj-dj*rki)/(di-dj)
!						! (c) Is the intersection inside the tile?
!						verti=vertices(tiles(2,k),1:3)
!						vertj=vertices(tiles(2,k),1:3)
!						v0=vertj-verti
!						v1=x-verti
!						if (dot_product(normal,cross_product(v0,v1)).gt.0.d0) then
!							vertk=verti
!1							verti=vertj
!							vertj=vertices(tiles(3,k),1:3)
!							v0=vertj-verti
!							v1=x-verti
!							if (dot_product(normal,cross_product(v0,v1)).gt.0.d0) then
!								verti=vertj
!								vertj=vertk
!								v0=vertj-verti
!								v1=x-verti
!								if (dot_product(normal,cross_product(v0,v1)).gt.0.d0) then
!									notblocked(ipair)=.false.
!									goto 99
!								endif
!							endif
!						endif
!					endif
!				endif
!			end do
!		endif
		! Find other facing pairs in which the second tile is either tile of this pair
!		nwhi=1
!		if (i.eq.1) then
!			ind10=1
!		else
!			ind10=list1(i)
!		endif
!		ind11=list1(i)
!		if (ind11.gt.ind10) then
!			nwhi=nwhi+(ind11-ind10)
!			wh(1:nwhi)=so1(ind10:ind11)
!		endif
!		nwh=nwhi
!		if (j.eq.1) then
!			ind10=1
!		else
!1			ind10=list1(j)
!		endif
!		ind11=list1(j)
!		if (ind11.gt.ind10) then
!			nwh=nwh+(ind11-ind10)
!			wh(nwhi:nwh)=so1(ind10:ind11)
!		endif
		! Loop through the first tile of these pairs
!		if (nwh.gt.1) then
!			do iblock=1,nwh
!				k=facingpairs(wh(iblock),2)
!				if (.not.(i.eq.k.and.j.eq.facingpairs(wh(iblock)),2)) then
!					normal=normals(k,1:3)
!					rki=centroids(i,1:3)-centroids(k,1:3)
!					rkj=centroids(j,1:3)-centroids(k,1:3)
!					! Does segment rij pass through this tile?
!					! (a) Are centroids of tiles i and j on opposite sides of tile k?
!					di=dot_product(normal,rki)
!					dj=dot_product(normal,rkj)
!					if (di*dj.lt.0.d0) then
!						! (b) Find point where ray intersects the tile plane
!						x=centroids(k,1:3)+(di*rkj-dj*rki)/(di-dj)
!						! (c) Is the intersection inside the tile?
!						verti=vertices(tiles(k,1),1:3)
!						vertj=vertices(tiles(k,2),1:3)
!						v0=vertj-verti
!						v1=x-verti
!						if (dot_product(normal,cross_product(v0,v1)).gt.0.d0) then
!							vertk=verti
!							verti=vertj
!							vertj=vertices(tiles(k,3),1:3)
!							v0=vertj-verti
!							v1=x-verti
!							if (dot_product(normal,cross_product(v0,v1)).gt.0.d0) then
!								verti=vertj
!								vertj=vertk
!								v0=vertj-verti
!								v1=x-verti
!								if (dot_product(normal,cross_product(v0,v1)).gt.0.d0) then
!									notblocked(ipair)=.false.
!									goto 99
!								endif
!							endif
!						endif
!					endif
!				endif
!			end do
!		endif
!	99 end do
!
!	nupairs=count(notblocked)
!	if(allocated(unblockedpairs)) deallocate(unblockedpairs)
!	allocate(unblockedpairs(1:nupairs,1:2))
!	unblockedpairs(1:nupairs,1)=pack(facingpairs(1:nfpairs,1),notblocked)
!	unblockedpairs(1:nupairs,2)=pack(facingpairs(1:nfpairs,2),notblocked)
!
!	print*,nupairs," unblocked pairs"
!
!	! Finally compress the list of pairs into the paircode array.
!
!	print*,"Compressing list"
!
!	! First time through is just to count the needed array elements
!
!	ipair=1
!	currentinlist=.false.
!	currentcount=0
!	ivis=1
!	do i=1,ntiles-1
!		do j=i+1,ntiles
!			if (ipair.lt.nupairs) then
!				if (currentinlist) then
!					if (i.eq.unblockedpairs(ipair,1).and.j.eq.unblockedpairs(ipair,2)) then
!						ipair=ipair+1
!						currentcount=currentcount+1
!					else 
!						ivis=ivis+1
!						currentcount=1
!						currentinlist=.false.
!					endif
!				else
!					if (i.eq.unblockedpairs(ipair,1).and.j.eq.unblockedpairs(ipair,2)) then
!						ipair=ipair+1
!						ivis=ivis+1
!						currentcount=1
!						currentinlist=.true.
!					else
!						currentcount=currentcount+1
!					endif
!				endif
!			else
!				if (currentinlist) then
!					ivis=ivis+1
!					currentcount=1
!					currentinlist=.false.
!				else
!					currentcount=currentcount+1
!				endif
!			endif
!		end do
!	end do
!	ivis=ivis+1
!
!	! Allocate the array in the asteroid structure
!
!	nvispairs=ivis
!	allocate(vispairs(1:nvispairs))


!	! Second time through, actually build the encoded vispairs array.
!
!	ipair=1
!1	currentinlist=.false.
!	currentcount=0
!	ivis=1
!	do i=1,ntiles-1
!		do j=i+1,ntiles
!			if (ipair.lt.nupairs) then
!				if (currentinlist) then
!					if (i.eq.unblockedpairs(ipair,1).and.j.eq.unblockedpairs(ipair,2)) then
!						ipair=ipair+1
!						currentcount=currentcount+1
!					else 
!						vispairs(ivis)=currentcount
!						ivis=ivis+1
!						currentcount=1
!						currentinlist=.false.
!					endif
!				else
!					if (i.eq.unblockedpairs(ipair,1).and.j.eq.unblockedpairs(ipair,2)) then
!						ipair=ipair+1
!						vispairs(ivis)=currentcount
!						ivis=ivis+1
!						currentcount=1
!						currentinlist=.true.
!					else
!						currentcount=currentcount+1
!					endif
!				endif
!			else
!				if (currentinlist) then
!					vispairs(ivis)=currentcount
!					ivis=ivis+1
!					currentcount=1
!					currentinlist=.false.
!				else
!					currentcount=currentcount+1
!				endif
!			endif
!		end do
!	end do
!	vispairs(ivis)=currentcount

!	return
!end subroutine calcVisTiles
end module