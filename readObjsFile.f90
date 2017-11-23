module readObjsFile
	implicit none
contains
!*****************************************************************************************
subroutine read_file(unitNum, fileName, numRows, numCols, Array)
	integer, intent(in) :: unitNum
	character(len=*), intent(in) :: fileName
	integer, intent(in) :: numRows, numCols
	real, dimension (1:numRows, 1:numCols), intent(out) :: Array
	
	character(len=300) :: line
	character(len=1) :: dataType
	integer :: i, j
	print*, "got in"
	line = "NA"
	open (unit=unitNum, file=fileName, status='old', action='read')
	ReadComments: do
		read (unitNum, '(A)') line
		print*, line
		if (line (1:1) /= "#") exit ReadComments
	end do ReadComments
	
	backspace (unitNum)
	print*, "got in1"
	do i=1, numRows
		print*, "got in2"
		read (unitNum, *) dataType, (Array (i, j), j=1,numCols)
		print*, "got in3"
	end do
	
	close(unitNum)
	return
end subroutine read_file

subroutine getVertexCount(unitNum, fileName, vCount)
	integer, intent(in) :: unitNum
	character(len=*), intent(in) :: fileName
	integer, intent(out) :: vCount
	
	character(len=300) :: line
	
	vCount = 0
	open (unit=unitNum, file=fileName, status='old', action='read')
	ReadLines: do
		read (unitNum, '(A)') line
		if (line (1:1) /= "#") then
			if (line(1:1) == "v") then
				vCount = vCount + 1
			end if
			if (line(1:1) == '') exit ReadLines
			if (line(1:1) == 'f') exit ReadLines
			if (line(1:1) /= "v") exit ReadLines
		end if
	end do ReadLines
	
	close(unitNum)
	return
end subroutine getVertexCount

subroutine getFacetCount(unitNum, fileName, fCount)
	integer, intent(in) :: unitNum
	character(len=*), intent(in) :: fileName
	integer, intent(out) :: fCount
	
	character(len=300) :: line
	character(len=1) :: dataType
	
	fCount = 0
	line = "NA"
	open (unit=unitNum, file=fileName, status='old', action='read')
	ReadLines: do
		read (unitNum, '(A)') line
		dataType = line(1:1)
		if (dataType /= "#" .and. dataType /= 'v' .and. dataType /= '') then
			if (dataType == "f") then
				fCount = fCount + 1
			end if
		else if (dataType /= "f" .and. fCount /= 0) then
			exit ReadLines
		end if
	end do ReadLines
	
	close(unitNum)
	return
end subroutine getFacetCount

subroutine getArrayCounts(unitNum, fileName, numDataType, dataTypeList, dataTypeCount)
	integer, intent(in) :: unitNum, numDataType
	character(len=*), intent(in) :: fileName
	character(len=*), dimension(1:numDataType), intent(in) :: dataTypeList
	integer, dimension(1:numDataType), intent(out) :: dataTypeCount
	
	character(len=300) :: line
	character(len=1) :: dataType
	integer :: i, loc
	dataTypeCount = 0
	loc = 1000
	print*, dataTypeCount
	
	open(unit=unitNum, file=fileName, status='old', action='read')
	ReadLines: do
		read(unitNum, '(A)', end=20) line
		dataType = line(1:1)
		if(any(dataTypeList==dataType)) then
			do i=1, numDataType
				if (dataTypeList(i) == dataType) then
					loc = i
					exit
				end if
			end do
			if (loc /= 1000) then
				dataTypeCount(loc) = dataTypeCount(loc) + 1
			end if
		end if
	end do ReadLines
	20 continue
	close(unitNum)
	return
end subroutine getArrayCounts
END MODULE readObjsFile