
module calculationsMod

contains
subroutine calcCentroids(XPoints, YPoints, ZPoints, centroid)
	implicit none
	real*8, dimension(1:3), intent(in) :: XPoints, YPoints, ZPoints
	real*8, dimension(1:3), intent(out) :: centroid
	
	
	centroid(1) = (XPoints(1) + XPoints(2) + XPoints(3))/3
	centroid(2) = (YPoints(1) + YPoints(2) + YPoints(3))/3
	centroid(3) = (ZPoints(1) + ZPoints(2) + ZPoints(3))/3

end subroutine calcCentroids

subroutine calcTileArea(Point1, Point2, Point3, Area)
	real*8, dimension(1:3), intent(in) :: Point1, Point2, Point3
	real*8, intent(out) :: Area
	
	real*8, dimension(1:3) :: crossed
	
	real*8, dimension(1:3) :: vector1, vector2
	real*8 :: magnitude
	call getVector(Point1, Point2, vector1)
	call getVector(Point1, Point3, vector2)
	call cross(vector1, vector2, crossed)
	
	call calcMagnitude(crossed, magnitude)
	
	Area = 0.5 * magnitude
	

end subroutine calcTileArea

subroutine calcMagnitude(vector, mag)
 implicit none
 real*8, dimension(1:3), intent(in) :: vector
 real*8, intent(out) :: mag
 
 real*8 :: sqrSum
 
 sqrSum = vector(1)**2 + vector(2)**2 + vector(3)**2
 mag = sqrt(sqrSum)

end subroutine calcMagnitude

subroutine cross(a,b, crossed)
	real*8, dimension(1:3), intent(in) :: a, b
	real*8, dimension(1:3), intent(out) :: crossed
	
	crossed(1) = a(2)*b(3) - a(3)*b(2)
	crossed(2) = a(3) * b(1) - a(1) * b(3)
	crossed(3) = a(1) * b(2) - a(2) * b(1)
end subroutine cross

subroutine getVector(Point1, Point2, vector)
	real*8, dimension(1:3), intent(in) :: Point1, Point2
	real*8, dimension(1:3), intent(out) :: vector
	
	vector(1) = Point1(1) - Point2(1)
	vector(2) = Point1(2) - Point2(2)
	vector(3) = Point1(3) - Point2(3)
end subroutine getVector

subroutine calcNormals(Point1, Point2, Point3, normals)
	implicit none
	real*8, dimension(1:3), intent(in) :: Point1, Point2, Point3
	
	real*8, dimension(1:3), intent(out) :: normals
	
	real*8, dimension(1:3) :: vector1, vector2
	
	call getVector(Point1, Point2, vector1)
	call getVector(Point1, Point3, vector2)
	call cross(vector1, vector2, normals)
end subroutine calcNormals

subroutine calcUnitNormals(Point1, Point2, Point3, unitNormals)
	implicit none
	real*8, dimension(1:3), intent(in) :: Point1, Point2, Point3
	
	real*8, dimension(1:3), intent(out) :: unitNormals
	
	real*8, dimension(1:3) :: normals
	
	call calcNormals(Point1, Point2, Point3, normals)
	
	unitNormals(1) = normals(1)/sqrt((normals(1)**2) + (normals(2)**2) + (normals(3)**2))
	unitNormals(2) = normals(2)/sqrt((normals(1)**2) + (normals(2)**2) + (normals(3)**2))
	unitNormals(3) = normals(3)/sqrt((normals(1)**2) + (normals(2)**2) + (normals(3)**2))
end subroutine calcUnitNormals

subroutine calcTotalMass(centroids, normals, areas, ntiles, totalMass)
	! Recalculate total mass
	real*8, dimension (1:ntiles) :: mass
	real*8 :: totmass
	integer :: i

	mass=0.d0
	totmass=0.d0

	do i=1,ntiles
		mass(i)=dot_product(centroids(1:3,i),normals(1:3,i)) & *areas(i)/3.d0
		totmass=totmass+mass(i)
	enddo
end subroutine calcTotalMass
end module