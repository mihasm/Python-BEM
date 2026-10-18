INCLUDE 'splprocs.f90'
INCLUDE 'epspsi.f90'
INCLUDE 'nacax.f90'
!+
PROGRAM NACA456                                                 ! naca456.f90
! ------------------------------------------------------------------------------
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!  2May01  0.1   RLC   Built first combined program from existing components
!  3May01  0.11  RLC   xxx
!  6Sep01  0.2   RLC   Fixed MakeFinalAirfoil; started ProcessGraphics
!  9Sep01  0.3   RLC   Added FitSymmetric and FitAsymmetric
! 19Sep01  0.4   RLC   Revised spline logic
! 25Sep01  0.5   RLC   Removed the combo CL mean line stuff
! 10Oct01  0.6   RLC   Started using SplineZero as the zero solver
! 15Nov01  0.7   RLC   Changed input variable from tableSize to denCode
! 24Nov01  0.8   RLC   Fixed coding of optional variables
! 27Nov01  0.85  RLC   Reversed order of arguments on MeanLine2
! 05Dec01  0.9   RLC   Replace MakeFinalAirfoil with calls to NACAauxilary
! 10Dec01  0.91  RLC   Worked on making output look nice
! 22Dec01  0.92  RLC   More rearrangement of output
!  4Jan02  1.0   RLC   Final cleanup for release of PDAS 7
! 07Nov10  1.1   RLC   Fixed typo in GetRk1

IMPLICIT NONE

  CHARACTER(LEN=15):: dateTimeStr
  CHARACTER(LEN=*),PARAMETER:: FAREWELL = &
   "End of naca456. Files naca.out, naca.gnu and naca.dbg created."
  INTEGER,PARAMETER:: IN=1, OUT=2, DBG=3, GNU=4   ! unit assignments
  CHARACTER(LEN=*),PARAMETER:: VERSION = "1.1 (7 November 2010)"
!-------------------------------------------------------------------------------
  CALL Welcome()
  CALL ReadDataAndMakeAirfoils()
  WRITE(*,*) FAREWELL
  STOP

CONTAINS

!+
PURE FUNCTION CountInputArray(preset, a) RESULT(n)
! ------------------------------------------------------------------------------
! PURPOSE - Count the number of elements of an array that have been input
!   thru a NAMELIST call. This is to save the user from the error-prone
!   process of counting the number of entries in a list. The array is
!   preset to a very negative value that no one would ever use for real data.
!   This routine simply searches for the first entry in the array that still
!   still has this value.
  REAL,INTENT(IN):: preset   ! the value that was originally placed in each
  REAL,INTENT(IN),DIMENSION(:):: a   ! array to be counted
  INTEGER:: n

  INTEGER:: k
!-------------------------------------------------------------------------------
  DO k=1,SIZE(a)
    IF (a(k) <= preset) EXIT
  END DO
  n=k-1   ! failed on k; k-1 was the last good one

  RETURN
END Function CountInputArray   ! -----------------------------------------------


!+
SUBROUTINE DescribeMeanLine(efu,camber,cmax,xmaxc,a,cl)
! ------------------------------------------------------------------------------
! PURPOSE - Print a description of the mean line in the printed output

  INTEGER,INTENT(IN):: efu   ! external file unit for output
  CHARACTER(LEN=*),INTENT(IN):: camber
  REAL,INTENT(IN):: cmax,xmaxc,a,cl

  CHARACTER(LEN=*),PARAMETER:: FMT2 = &
   "('Two digit mean line,  max camber=', F6.3, '  x of max. camber=',F6.3)"
  CHARACTER(LEN=*),PARAMETER:: FMT3 = &
   "('Three digit mean line, cl=', F7.4, '  x of max. camber=',F6.3)"
  CHARACTER(LEN=*),PARAMETER:: FMT3R = &
   "('Three digit mean line, cl=', F7.4, '  x of max. camber=',F6.3)"
  CHARACTER(LEN=*),PARAMETER:: FMT6 = &
   "('Six-series mean line, a=', F7.4, '   cl=',F7.4)"
  CHARACTER(LEN=*),PARAMETER:: FMT6M = &
   "('Six-series mean line (modified), cl=',F7.4, '  (a=0.8 always for 6M)')"

!-------------------------------------------------------------------------------
  SELECT CASE(camber)
    CASE('0')
      WRITE(efu,*) "No camber, symmetrical section"
    CASE('2')
      WRITE(efu,FMT2) cmax, xmaxc
    CASE('3')
      WRITE(efu,FMT3)  cl,xmaxc
    CASE('3R','3r')
      WRITE(efu,FMT3R) cl,xmaxc
    CASE('6')
      WRITE(efu,FMT6) a,cl
    CASE('6M','6m')
      WRITE(efu,FMT6M) cl
    CASE DEFAULT

  END SELECT

  RETURN
END Subroutine DescribeMeanLine   ! --------------------------------------------

!+
SUBROUTINE DescribeProfile(efu,profile,toc,leIndex,xmaxt)
! ------------------------------------------------------------------------------
! PURPOSE - Print a description of the profile in the printed output

  INTEGER,INTENT(IN):: efu   ! external file unit for output
  CHARACTER(LEN=*),INTENT(IN):: profile
  REAL,INTENT(IN):: toc,leIndex,xmaxt

  CHARACTER(LEN=*),PARAMETER:: FMT4 = &
   "('Four-digit profile, t/c=',F7.4)"
  CHARACTER(LEN=*),PARAMETER:: FMT4M = &
   "('Four-digit-modified profile, t/c=',F7.4/"// &
   "'   l.e. index=',F7.4, '   x of max. thickness=',F7.4)"
  CHARACTER(LEN=*),PARAMETER:: FMT6 = &
   "('Six-series profile, family= ',A, '   t/c=',F7.4)"

!-------------------------------------------------------------------------------
  SELECT CASE(profile)
    CASE('4')
      WRITE(efu,FMT4) toc
    CASE('4M','4m')
      WRITE(efu,FMT4M) toc,leIndex,xmaxt
    CASE DEFAULT
      WRITE(efu,FMT6) profile,toc
  END SELECT
  RETURN
END Subroutine DescribeProfile   ! ------------------------------------------

!+
FUNCTION GetDateTimeStr() RESULT(s)
! ------------------------------------------------------------------------------
! PURPOSE - Return a string with the current date and time
!   This is now standard Fortran. It should work on ALL compilers!
!   You can change the first I2.2 below to I2 if you don't like a leading
!   zero on early morning times, i.e., 6:45 instead of 06:45
  IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER:: MONTH='JanFebMarAprMayJunJulAugSepOctNovDec'
  CHARACTER(LEN=*),PARAMETER:: FMT = '(I2.2,A1,I2.2,I3,A3,I4)'
  CHARACTER(LEN=15):: s
  INTEGER,DIMENSION(8):: v
!-------------------------------------------------------------------------------
  CALL DATE_AND_TIME(VALUES=v)

  WRITE(s,FMT) v(5), ':', v(6), v(3), MONTH(3*v(2)-2:3*v(2)), v(1)
  RETURN
END FUNCTION GetDateTimeStr   ! ------------------------------------------------

!+
SUBROUTINE ReadDataAndMakeAirfoils()
! ------------------------------------------------------------------------------
! PURPOSE - Read data from the input file (UNIT=IN) and define an airfoil.

USE NACAauxilary

  REAL:: a   ! chordwise extent of constant loading
  CHARACTER(LEN=2):: camber ='0'  ! Camber line (2-digit or 3-digit;
                          ! = '0' no camber
                          ! = '2' 2-digit mean line
                          ! = '3' 3-digit mean line
                          ! = '3R' 3-digit reflex mean line
                          ! = '6' 6-series mean line
                          ! = '6M' 6-series (modified) mean line

  REAL:: chord = 1.0   ! model chord length used for listing ordinates
  REAL:: cl            ! design lift coefficient of the camber line
                       ! applies to 3-digit, 3-digit-reflex and
                       !   6-series camber lines
  REAL:: cmax = 0.04   ! maximum camber ordinate-to-chord ratio
  INTEGER:: denCode = 1  ! =0 if specified by user
                         ! =1 for coarse output
                         ! =2 for fine output
                         ! =3 for very fine output
  INTEGER:: errCode
  REAL:: leIndex = 6.0   ! leading-edge index (4-digit modified)
  CHARACTER(LEN=80):: name = ' '   ! title desired on output
  INTEGER:: nTable = 1   ! size of xTable (only used if denCode=0 )
  INTEGER:: nupper, nlower
  REAL,PARAMETER:: PRESET = -1E20
  CHARACTER(LEN=3):: profile = '4'   ! NACA airfoil family 
                           ! = '4' four digit profile
                           ! = '4M' modified-four-digit profile
                           ! = '63' 63 series profile
                           ! = '64' 64 series profile
                           ! = '65' 65 series profile
                           ! = '66' 66 series profile
                           ! = '67' 67 series profile
                           ! = '63A' 63A series profile
                           ! = '64A' 64A series profile
                           ! = '65A' 65A series profile
  REAL:: toc = 0.1   ! thickness-chord ratio of airfoil.
  REAL:: xorigin=0.0, yorigin=0.0  ! leading edge point
  REAL:: xmaxc = 0.4  ! location of maximum camber position
                      ! used for 2-,3- and 3-digit-reflex mean lines
  REAL:: xmaxt = 0.3  ! Nondimensional chordwise location of maximum
                      ! thickness.  Used for 4-digit modified airfoils only.
  REAL,DIMENSION(1000):: xTable  ! user supplied table of x-coordinates
                                 ! only used if denCode=0

  REAL,ALLOCATABLE,DIMENSION(:):: x,ymean,ymeanp,yt,ytp
  REAL,ALLOCATABLE,DIMENSION(:):: xupper,xlower, yupper,ylower
  REAL,ALLOCATABLE,DIMENSION(:):: yu,yl,yup,ylp

  CHARACTER(LEN=*),PARAMETER:: TITLE1 = "COORDINATES OF NACA AIRFOILS"
  CHARACTER(LEN=*),PARAMETER:: TITLE2 = "SYMMETRICAL AIRFOIL DEFINITION"
  CHARACTER(LEN=*),PARAMETER:: TITLE3 = "x         y        dy/dx"
  
  CHARACTER(LEN=*),PARAMETER:: TITLE4 = "THICKNESS          MEAN LINE"
  CHARACTER(LEN=*),PARAMETER:: TITLE5 = &
    "x         y      dy/dx         y      dy/dx"
  CHARACTER(LEN=*),PARAMETER:: TITLE6 = "COMBINED THICKNESS AND CAMBER"
  CHARACTER(LEN=*),PARAMETER:: TITLE7 = &
     "xupper    yupper    xlower    ylower"
  CHARACTER(LEN=*),PARAMETER:: TITLE8 = "INTERPOLATED COORDINATES"
  CHARACTER(LEN=*),PARAMETER:: TITLE9 = "x      yupper   ylower"
  CHARACTER(LEN=*),PARAMETER:: TITLE10 = "SCALED AND TRANSLATED AIRFOIL"

  CHARACTER(LEN=*),PARAMETER:: FMT1 = "('Scaled chord length=',F12.5)"
  CHARACTER(LEN=*),PARAMETER:: FMt2 = &
    "('   leading edge at x=', F12.5, '   y=',F12.5)"

  NAMELIST /NACA/ a, camber, chord, cl, cmax, dencode,leIndex,name, &
    profile, toc, xmaxc,xmaxt,xorigin,yorigin, xTable

!-------------------------------------------------------------------------------
  xTable(:)=PRESET
  READ(UNIT=IN,NML=NACA,IOSTAT=errCode)
  
  IF (errCode < 0) THEN
    WRITE(*,*)   "End of file reading input"
    WRITE(DBG,*) "End of file reading input"
    WRITE(OUT,*) "End of file reading input"
    STOP
  END IF
  
!  IF (errCode > 0) THEN
!    WRITE(*,*)   "Abnormal termination of naca456"
!    WRITE(DBG,*) "Abnormal termination of naca456"
!    WRITE(OUT,*) "Abnormal termination of naca456"
!    STOP
!  END IF

  profile=AdjustL(profile)   ! no leading blanks
  camber=AdjustL(camber)

  WRITE(DBG,*) "Namelist read successfully"
  WRITE(DBG,'(T2,A)') name
  
  IF (denCode == 0) THEN
    nTable=CountInputArray(PRESET, xTable)
  ELSE
    CALL LoadX(denCode,nTable,xTable)   ! sets nTable and xTable
  END IF

  WRITE(DBG,*) "Allocating x,ymean,... with denCode=", denCode
  ALLOCATE (x(nTable),ymean(nTable),ymeanp(nTable), yt(nTable),ytp(nTable))
  ALLOCATE (xUpper(nTable),yUpper(nTable), xLower(nTable),yLower(nTable))

  x=xTable(1:nTable)
  WRITE(DBG,*) "nTable is now ", nTable

  WRITE(DBG,*) "THICKNESS DEFINITION"
  SELECT CASE(AdjustL(profile))
    CASE('4')
      WRITE(DBG,*) "Four-digit profile"
      CALL Thickness4(toc, x,yt,ytp)
    CASE('4M','4m')
      WRITE(DBG,*) "Four-digit-modified profile", leIndex, xmaxt
      CALL Thickness4M(toc,leIndex,xmaxt, x,yt,ytp)
    CASE('63')
      WRITE(DBG,*) "63-series profile"
      CALL Thickness6(1,toc, x,yt,ytp)
    CASE('64')
      WRITE(DBG,*) "64-series profile"
      CALL Thickness6(2,toc, x,yt,ytp)
    CASE('65')
      WRITE(DBG,*) "65-series profile"
      CALL Thickness6(3,toc, x,yt,ytp)
    CASE('66')
      WRITE(DBG,*) "66-series profile"
      CALL Thickness6(4,toc, x,yt,ytp)
    CASE('67')
      WRITE(DBG,*) "67-series profile"
      CALL Thickness6(5,toc, x,yt,ytp)
    CASE('63A','63a')
      WRITE(DBG,*) "63A-series profile"
      CALL Thickness6(6,toc, x,yt,ytp)
    CASE('64A','64a')
      WRITE(DBG,*) "64A-series profile"
      CALL Thickness6(7,toc, x,yt,ytp)
    CASE('65A','65a')
      WRITE(DBG,*) "65A-series profile"
      CALL Thickness6(8,toc, x,yt,ytp)
    CASE DEFAULT
      WRITE(*,*) "Not a valid profile: " // profile
      WRITE(DBG,*) "Not a valid profile: " // profile
      STOP
  END SELECT
  WRITE(DBG,*) "          x         y        dy/dx"
  CALL PrintArraysNumbered(DBG, x,yt,ytp)

  
  SELECT CASE(AdjustL(camber))
    CASE('0')
      WRITE(DBG,*) "No camber"
      ymean=0.0
      ymeanp=0.0
    CASE('2')
      WRITE(DBG,*) "Two digit mean line ", xmaxc, cmax
      CALL MeanLine2(cmax,xmaxc, x,ymean,ymeanp)
     CASE('3')
      WRITE(DBG,*) "Three digit camber line"
      CALL MeanLine3(cl,xmaxc, x,ymean,ymeanp)
    CASE('3R','3r')
      WRITE(DBG,*) "Three digit (reflex) camber line"
      CALL MeanLine3Reflex(cl,xmaxc, x,ymean,ymeanp)
    CASE('6')
      WRITE(DBG,*) "6-series camber line"
      CALL MeanLine6(a,cl, x,ymean,ymeanp)
    CASE('6A','6a','6M','6m')
      WRITE(DBG,*) "6A-series camber line"
      CALL MeanLine6M(cl, x,ymean,ymeanp)
    CASE DEFAULT
      WRITE(*,*) "Not a valid mean line: " // camber
      WRITE(DBG,*) "Not a valid mean line: " // camber
      STOP
  END SELECT
  
!... If the airfoil is symmetrical, the output is simplified...  
  IF (MAXVAL(ABS(ymean))==0.0) THEN
!    WRITE(OUT,'(A,T40,A)') TITLE1,dateTimeStr
    WRITE(OUT,*) name
    CALL DescribeProfile(OUT,profile,toc,leIndex,xmaxt)
    CALL DescribeMeanLine(OUT,camber,cmax,xmaxc,a,cl)
    WRITE(OUT,*) TITLE2
    WRITE(OUT,'(T9,A)') TITLE3
    CALL PrintArraysNumbered(OUT,x,yt,ytp)
    CALL ScaleAirfoil(chord, xOrigin,yOrigin, x,yt, x,-yt)
    WRITE(OUT,*) "End of output for symmetrical airfoil"
    RETURN
  END IF

  WRITE(DBG,*) "MEAN LINE DEFINITION"
  WRITE(DBG,*) "         x          y        dy/dx"
  CALL PrintArraysNumbered(DBG, x,ymean,ymeanp)

!!!  WRITE(OUT,'(A,T40,A)') TITLE1,dateTimeStr
  WRITE(OUT,*) " "
  WRITE(OUT,*) name
  CALL DescribeProfile(OUT,profile,toc,leIndex,xmaxt)
  CALL DescribeMeanLine(OUT,camber,cmax,xmaxc,a,cl)
  WRITE(OUT,*) " "

  WRITE(OUT,'(T24,A)') TITLE4
  WRITE(OUT,'(T13,A)') TITLE5
  CALL PrintArraysNumbered(OUT,x,yt,ytp,ymean,ymeanp)

  CALL CombineThicknessAndCamber(x,yt,ymean,ymeanp, &
        xupper,yupper, xlower,ylower)
  nupper=nTable
  nlower=nupper      
  WRITE(DBG,*) TITLE6
  WRITE(DBG,'(T9,A)') TITLE7
  CALL PrintArraysNumbered(DBG, xupper(1:nupper),yupper(1:nupper), &
                                xlower(1:nlower),ylower(1:nlower))
  WRITE(OUT,'(A)') ACHAR(12)   ! form-feed (new page)
  WRITE(OUT,'(A,T40,A)') TITLE1,dateTimeStr
  WRITE(OUT,*) name
  WRITE(OUT,*) TITLE6
  WRITE(OUT,'(T10,A)') TITLE7
  CALL PrintArraysNumbered(OUT, xupper(1:nupper),yupper(1:nupper), &
                                xlower(1:nlower),ylower(1:nlower))

  ALLOCATE(yu(nTable),yl(nTable), yup(nTable),ylp(nTable))
  CALL InterpolateUpperAndLower(xupper,yupper, xlower,ylower, &
    x,yu,yl,yup,ylp)

  WRITE(OUT,'(A)') ACHAR(12)   ! form-feed (new page)
  WRITE(OUT,'(A,T40,A)') TITLE1,dateTimeStr
  WRITE(OUT,*) name
  WRITE(OUT,*) TITLE8
  WRITE(OUT,'(T13,A)') TITLE9

  CALL PrintArraysNumbered(DBG, x,yu,yup,yl,ylp)
  CALL PrintArraysNumbered(OUT, x,yu,yl)

  IF (xOrigin /= 0.0 .OR. yOrigin /= 0.0 .OR. ABS(chord-1.0) > 1E-6) THEN
    WRITE(OUT,'(A)') ACHAR(12)   ! form-feed (new page)
    WRITE(OUT,'(A,T40,A)') TITLE1,dateTimeStr
    WRITE(OUT,*) name
    WRITE(OUT,FMT1) chord
    WRITE(OUT,FMT2) xOrigin,yOrigin
    WRITE(OUT,*) TITLE10
    WRITE(OUT,'(T10,A)') TITLE7
  END IF
  CALL ScaleAirfoil(chord, xOrigin,yOrigin, x,yu, x,yl)

  DEALLOCATE(ylp,yup,yl,yu,yLower,yUpper,yMean,x)
  WRITE(OUT,*) "End of output for cambered airfoil"
  RETURN
END Subroutine ReadDataAndMakeAirfoils   ! -------------------------------------

!+
SUBROUTINE PrintArrays(efu, x,y, a,b,c,d)
! ------------------------------------------------------------------------------
! PURPOSE - Print arrays, neatly formatted. X and Y are always printed.
!   The next 4 optional variables are printed when present.

  INTEGER,INTENT(IN):: efu   ! unit number of the external file
  REAL,INTENT(IN),DIMENSION(:):: x,y
  REAL,INTENT(IN),DIMENSION(:),OPTIONAL:: a,b,c,d

  CHARACTER(LEN=80):: buffer
  CHARACTER(LEN=*),PARAMETER:: FMT = '(2F10.6)'
  INTEGER:: k
!-------------------------------------------------------------------------------
  buffer=""
  DO k=1,SIZE(x)
    WRITE(buffer(1:20),FMT) x(k),y(k)
    IF (Present(a)) WRITE(buffer(21:30),FMT) a(k)
    IF (Present(b)) WRITE(buffer(31:40),FMT) b(k)
    IF (Present(c)) WRITE(buffer(41:50),FMT) c(k)
    IF (Present(d)) WRITE(buffer(51:60),FMT) d(k)
    WRITE(efu,*) Trim(buffer)
  END DO
  RETURN
END Subroutine PrintArrays   ! -------------------------------------------------

!+
SUBROUTINE PrintArraysNumbered(efu, x,y, a,b,c,d)
! ------------------------------------------------------------------------------
! PURPOSE - Print arrays, neatly formatted. X and Y are always printed.
!   The next 4 optional variables are printed when present.

  INTEGER,INTENT(IN):: efu   ! unit number of the external file
  REAL,INTENT(IN),DIMENSION(:):: x,y
  REAL,INTENT(IN),DIMENSION(:),OPTIONAL:: a,b,c,d

  CHARACTER(LEN=80):: buffer
  CHARACTER(LEN=*),PARAMETER:: FMT = '(2F10.6)'
  INTEGER:: k
!-------------------------------------------------------------------------------
  buffer=""
  DO k=1,SIZE(x)
    WRITE(buffer(1:5),'(I5)') k
    WRITE(buffer(6:25),FMT) x(k),y(k)
    IF (Present(a)) WRITE(buffer(26:35),FMT) a(k)
    IF (Present(b)) WRITE(buffer(36:45),FMT) b(k)
    IF (Present(c)) WRITE(buffer(46:55),FMT) c(k)
    IF (Present(d)) WRITE(buffer(56:65),FMT) d(k)
    WRITE(efu,*) Trim(buffer)
  END DO
  RETURN
END Subroutine PrintArraysNumbered   ! -----------------------------------------

!+
SUBROUTINE ScaleAirfoil(scale,xOrigin,yOrigin, &
  xupper,yupper, xlower,ylower)
! ------------------------------------------------------------------------------
! PURPOSE - Translate and scale the airfoil
!   If no scaling or translation is done, just return
  REAL,INTENT(IN):: scale, xOrigin,yOrigin
  REAL,INTENT(IN),DIMENSION(:):: xupper,yupper,xlower,ylower

  REAL,DIMENSION(SIZE(xupper)):: xu
  REAL,DIMENSION(SIZE(yupper)):: yu
  REAL,DIMENSION(SIZE(xlower)):: xl
  REAL,DIMENSION(SIZE(ylower)):: yl

!-------------------------------------------------------------------------------
  xu=xOrigin + scale*xupper
  yu=yOrigin + scale*yupper   
  xl=xOrigin + scale*xlower
  yl=yOrigin + scale*ylower

  IF (xOrigin /= 0.0 .OR. yOrigin /= 0.0 .OR. ABS(scale-1.0) > 1E-6) &
    CALL PrintArrays(OUT, xu,yu,xl,yl)

  CALL PrintArrays(GNU,xu,yu)      ! plot upper surface
  WRITE(GNU,*) " "                 ! invisible move to start of lower surface
  CALL PrintArrays(GNU,xl,yl)      ! plot lower surface
  CLOSE(UNIT=GNU)                  ! plot is complete
  
  RETURN
END Subroutine ScaleAirfoil   ! ------------------------------------------------

!+
SUBROUTINE Welcome()
! ------------------------------------------------------------------------------
! PURPOSE - Set up all the files
! NOTE - dateTimeStr and all file numbers are global variables

  CHARACTER(LEN=*),PARAMETER:: AUTHOR = &
    "Ralph L. Carmichael, Public Domain Aeronautical Software"
  INTEGER:: errCode
  CHARACTER(LEN=132):: fileName
  CHARACTER(LEN=*),PARAMETER:: GREETING = "COORDINATES OF NACA AIRFOILS"
!-------------------------------------------------------------------------------
  dateTimeStr=GetDateTimeStr()

  WRITE(*,*) GREETING
  WRITE(*,*) "Program naca456, version "//VERSION
  WRITE(*,*) AUTHOR
  DO
    WRITE(*,*) "Enter the name of the input file:"
    READ(*,'(A)',IOSTAT=errCode) fileName
    IF (Len_Trim(fileName)==0) STOP   ! user gives up
    OPEN(UNIT=IN, FILE=fileName, STATUS='OLD', &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT

    OPEN(UNIT=IN, FILE=Trim(fileName)//'.nml', STATUS='OLD', &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT

    OPEN(UNIT=IN, FILE=Trim(fileName)//'.inp', STATUS='OLD', &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT

    WRITE(*,*) "Unable to open " // Trim(fileName) // ". Try again."
  END DO

  INQUIRE(UNIT=IN, NAME=fileName)
  WRITE(*,*) "Reading from "//Trim(fileName)

  OPEN(UNIT=DBG, FILE='naca.dbg', STATUS='REPLACE', ACTION='WRITE')
  WRITE(DBG,*) GREETING
  WRITE(DBG,*) "Created by naca456, version "//VERSION//" on "//dateTimeStr
  WRITE(DBG,*) AUTHOR

  OPEN(UNIT=OUT, FILE='naca.out', STATUS='REPLACE', ACTION='WRITE')
  WRITE(OUT,*) GREETING
  WRITE(OUT,*) "Created by naca456, version "//VERSION//" on "//dateTimeStr
  WRITE(OUT,*) AUTHOR

  OPEN(UNIT=GNU, FILE='naca.gnu', STATUS='REPLACE', ACTION='WRITE')

  RETURN
END Subroutine Welcome   ! -----------------------------------------------------

END Program NACA456   ! ========================================================
