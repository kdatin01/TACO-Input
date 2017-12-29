module asteroid_struct

implicit None
include "KellieTACOpack/taco_src/taco_parameters.h"
type asteroid_structure
   real*8, dimension(1:nvertices,1:3) :: vertices
   real*8, dimension(1:ntiles,1:3) :: normals,centroids
   real*8, dimension(0:nhorizon-1,0:ntiles-1) :: horizonmap
   real*8, dimension(0:ntiles-1) :: areas
   real*8, dimension(0:2) :: principalmoments
   real*8 :: totalmass
   integer, dimension(1:ntiles,1:3) :: tiles
   integer, dimension(:), allocatable :: vispairs
   integer :: nvertices,ntiles,nhorizon,nvispairs
end type

type (asteroid_structure) asteroid

contains

integer function write_asteroid(filename)
implicit none
character*80 :: filename

open(3,file=filename,form='unformatted',status='new', err=91)

write(3,err=92)asteroid%nvertices
write(3,err=92)asteroid%ntiles
write(3,err=92)asteroid%nhorizon
write(3,err=92)asteroid%vertices
write(3,err=92)asteroid%normals
write(3,err=92)asteroid%centroids
write(3,err=92)asteroid%horizonmap
write(3,err=92)asteroid%areas
write(3,err=92)asteroid%tiles
write(3,err=92)asteroid%totalmass
write(3,err=92)asteroid%principalmoments
write(3,err=92)asteroid%nvispairs
if (asteroid%nvispairs.gt.0) write(3,err=92)asteroid%vispairs

close(3,err=93)

write_asteroid=0
return

91 write_asteroid=1
return
92 write_asteroid=2
return
93 write_asteroid=3
return

end function




end module asteroid_struct