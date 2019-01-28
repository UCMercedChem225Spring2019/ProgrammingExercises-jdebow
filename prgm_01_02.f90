      Program program_02
!
!     The only change needed is to add the second matrix to the code.
!
      implicit none
      integer,parameter::inFileUnitA=10,inFileUnitB=11
      integer::errorFlag,i
      real,dimension(3,3)::matrixInA,matrixInB
      character(len=128)::fileNameA,fileNameB
!
!     We now add the second write and read command.
!
      write(*,*)' What is the name of the input data file?'
      read(*,*) fileNameA
      write(*,*)' What is the name of the second input data file?'
      read(*,*) fileNameB
!
      open(unit=inFileUnitA,file=TRIM(fileNameA),status='old',  &
        iostat=errorFlag)
      if(errorFlag.ne.0) then
        write(*,*)' There was a problem opening the input file.'
        goto 999
      endIf
      do i = 1,3
        read(inFileUnitA,*) matrixInA(1,i),matrixInA(2,i),matrixInA(3,i)
      endDo
      close(inFileUnitA)
!
!     This code will read the second matrix.
!
      open(unit=inFileUnitB,file=TRIM(fileNameB),status='old',  &
        iostat=errorFlag)
      if(errorFlag.ne.0) then
        write(*,*)' There was a problem opening the input file.'
        goto 999
      endIf
      do i = 1,3
        read(inFileUnitB,*) matrixInB(1,i),matrixInB(2,i),matrixInB(3,i)
      endDo
      close(inFileUnitB)
!
!     A second call is necessary to print the second matrix.
!
      call PrintMatrix3x3(matrixInA)
      call PrintMatrix3x3(matrixInB)
!
  999 continue

      End Program program_02


      Subroutine PrintMatrix3x3(matrix)
!
      implicit none
      real,dimension(3,3),intent(in)::matrix
      integer::i
!
 1000 format(3(2x,f5.1))
!
      write(*,*)' Printing Matrix'
!
      write(*,*)matrix
      return
      End Subroutine PrintMatrix3x3
