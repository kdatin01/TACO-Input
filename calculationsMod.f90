
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

subroutine tilehorizon1(tilenumber,pairs,tilehorizonmap, ntiles, npairs, nvertices, nhorizon, pi, normals, vertices, centroids)
	implicit none

	integer, intent(in) :: tilenumber, ntiles, npairs, nvertices, nhorizon
	real*8, intent(in) :: pi
	real*8, dimension(1:ntiles, 1:3), intent(in):: normals, centroids
	real*8, dimension(1:nvertices,1:3), intent(in) :: vertices
	integer, dimension(1:npairs, 1:2), intent(in) :: pairs
	real*8, dimension(1:ntiles), intent(inout) :: tilehorizonmap
	real*8, dimension(1:nvertices) :: elevation,azimuth
	real*8, dimension(1:3) :: dr,dr2,refvector,refvec,direction,direct, out1
	real*8, dimension (1) :: maxseg
	real*8 :: sense,az
	real*8, dimension(1:npairs) :: el0,el1,az0,az1,segmentelev,temp
	integer, dimension(1:npairs,1:2) :: usepairs
	integer, dimension(1:nvertices) :: tileindexnumbers,blank
	integer, dimension(1:nvertices) :: abovehorizon
	integer, dimension(1:3) :: tile
	integer, dimension(1:2) :: pair
	integer, dimension (3) :: index
	integer, dimension (1) :: maxtile,mintile
	integer :: nwh,nabovehorizon,nusepairs,ii,l,j,k,i,paircount
	logical, dimension(1:nvertices) :: mask
	logical :: yes

	do i=1,nvertices
		tileindexnumbers(i)=i
		blank(i)=0
	end do

	i=tilenumber

	! Get azimuth relative to north and elevation of all tile vertices as seen
	! from the centroid of the ith tile.

	call cross_product(normals(i, 1:3),(/0.d0,0.d0,1.d0/),refvec)
	refvector=refvec/sqrt(sum(refvec**2))
	do j=1,nvertices
		dr2=vertices(j,1:3)-centroids(i,1:3)
		dr=dr2/sqrt(sum(dr2**2))
		call cross_product(normals(i,1:3),dr, direct)
		direction=direct/sqrt(sum(direct**2))
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

		abovehorizon=pack(tileindexnumbers,mask,blank)
		paircount=1
		do j=1,npairs
			yes=.false.
			do l=1,nabovehorizon
				yes=yes.or.(pairs(1,j).eq.abovehorizon(l) .or.pairs(2,j).eq.abovehorizon(l))
			enddo
			if (yes) then
				usepairs(paircount,1:2)=(/pairs(1,j),pairs(2,j)/)
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
end subroutine tilehorizon1
end module