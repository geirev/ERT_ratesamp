program sample
! Samples the measurement perturbation matrix E for the ert-Eclipse applications.
! The program reads the ert/observations/obshist.txt file and the schedulefile defined below.
! Only configured for opr, wpr, gpr observations.
! Predefined well names.

   use m_pseudo1D
   use m_fixsample1D
   use m_set_random_seed2
   use m_random
   implicit none

   integer, parameter :: nrens=100        ! ensemble size
! read observations info assumed stored as follows:
   integer, parameter :: nx=36            ! Number of dates keywords
   integer :: nrwells                     ! number of wells
   integer, parameter :: nrdata=3         ! number of variables per well (OPR, WPR, GPR)
! We will then sample nrwells*nrdata time series of length nx and store these as one realization
! A total of nrens such realizations are computed.

   real :: xlength=35.0                   ! length of time series is 36 months
   real :: cor1=15.0                      ! decorrelation length in months
   character(len=100) :: obshistfile='/home/geve/Dropbox/Statoil/erterr/observations/obshistnew.txt'
   character(len=100) :: schedulefile='/home/geve/Dropbox/Statoil/eclipse/include/history.sch'
   character(len=100) :: controlpath='/home/geve/Dropbox/Statoil/eclipse/Priors/control/'

   integer n1                             ! x-dimension of grid
   integer i,j,k,l,m,ic,j1
   integer idat
   integer iwell
   integer ipos(5)
   real dx


   character(len=8) wellname

   real, allocatable :: B(:,:)                        !  1D samples
   real, allocatable :: C(:,:,:,:)
   real, allocatable :: obs(:,:,:) 
   real, allocatable :: stddev(:,:,:)

   type wellrates
      character(len=4)  wellname
      real opr(nx)
      real gpr(nx)
      real wpr(nx)
   end type
   type(wellrates), allocatable :: rates(:)


   character(len=100) cline
   character(len=4) var
   character(len=4) well
   character(len=4) cens
   real rel_err
   real min_err
   logical ex

   character(len=1) tmp
   integer istat
   integer noise
   character(len=200) string1
   character(len=200) string2
   character(len=200) string3

! Extract rates and well names from schedule file amd dump them in rates.txt
   string1="grep -A1 -e 'DATES' -e 'OPEN RESV' "//trim(schedulefile)// " > rates.txt"
   string2="sed -i -e '/^\/$/d' -e '/--/d' -e 's?/??' -e 's/ OPEN RESV //' -e 's/  */#/g' rates.txt"
   print '(a,a)','string1=',trim(string1)
   istat = system( string1 )
   print '(a,a)','string2=',trim(string2)
   istat = system( string2 )
   istat = system( 'head rates.txt' )




! Optionally generate the list of well names
   inquire(file='wells.txt',exist=ex)
   tmp='y'
   if (ex) then
      write(*,'(a)',advance='no')'Overwrite existing wells.txt file (Y/n)'
      read(*,*)tmp
   endif
   if (tmp=='y' .or. tmp=='Y') then
      string3="cat rates.txt | sed -e ?/'/d? -e '/DATES/d' | cut -f2 -d# | sort -u > wells.txt"
      istat = system( string3 )
   endif

   open(10,file='wells.txt')
      do i=1,1000
         read(10,'(a)',end=100)tmp
      enddo
      100 rewind(10)
      nrwells=i-1
      allocate (rates(nrwells))
      do i=1,nrwells
         read(10,'(a)')rates(i)%wellname
         write(*,'(a)')rates(i)%wellname
      enddo
   close(10)

! Allocating matrices for all wells
   allocate( B(nx, nrwells*nrdata*nrens) )
   allocate( C(nx, nrdata, nrwells, nrens) )
   allocate( stddev(nx,nrdata,nrwells) )
   allocate( obs(nx,nrdata,nrwells) )

   write (*,'(a)',advance='no')'Select type of noise (0-null, 1-pseudo, 2-white, 3-constant): '
   read(*,*)noise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set a variable random seed
   call set_random_seed2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Correction for Gaussian variogram
   cor1=cor1/sqrt(3.0)

   dx=xlength/float(nx-1)
   n1=nint(real(nx)*1.5)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sample all random time series
   select case (noise)
   case (0)
      B=0.0
      print *,'Simulating noise=0.0'
   case (1)
      call pseudo1D(B,nx,nrwells*nrdata*nrens,cor1,dx,n1)
      call fixsample1D(B,nx*nrwells*nrdata,nrens)
      print *,'Simulating pseudo random red noise'

   case (2)
      print *,'Simulating white random noise'
      call random(B,nx*nrwells*nrdata*nrens)

   case (3)
      print *,'Simulating random bias noise'
      B=-1.0
      call random(B(1,1:nrwells*nrdata*nrens),nrwells*nrdata*nrens)
      do i=1,nrwells*nrdata*nrens
         B(2:nx,i)=B(1,i)
      enddo
   end select

! test
   do i=1,nx
      write(*,'(i4,12f10.2)')i,B(i,1:12)
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initializing rates to 0
   do j=1,nrwells
      do i=1,nx
         rates(j)%opr(i)=0.0
         rates(j)%wpr(i)=0.0
         rates(j)%gpr(i)=0.0
      enddo
   enddo

   print *,'B'
! reading rates and times from rates.txt
   open (10,file='rates.txt')
   idat=1
   do i=1,10000
      read(10,'(a)',end=200)cline
      print '(a)',trim(cline)
      if (cline(1:5) == 'DATES') then
         read(10,'(a)',end=200)cline
         idat=idat+1
      else
         ic=0
         do j=1,len_trim(cline)
            if(cline(j:j)=='#') then
               ic=ic+1
               ipos(ic)=j
            endif
         enddo
         wellname(:)=' '
         read(cline(ipos(1)+1:ipos(2)-1),'(a)') wellname(1:ipos(2)-ipos(1)-1)
         iwell=0
         do k=1,nrwells
            if (trim(wellname) == trim(rates(k)%wellname)) then
               iwell=k
            endif
         enddo

         read(cline(ipos(2)+1:ipos(3)-1),*)rates(iwell)%opr(idat)
         read(cline(ipos(3)+1:ipos(4)-1),*)rates(iwell)%wpr(idat)
         read(cline(ipos(4)+1:ipos(5)-1),*)rates(iwell)%gpr(idat)
      endif
   enddo
   200 close(10)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scaling perturbations according to measurement magnitude and error defs in obshistfile
   do i=1,nrwells
      do k=1,nx
         obs(k,1,i) = rates(i)%opr(k)
         obs(k,2,i) = rates(i)%wpr(k)
         obs(k,3,i) = rates(i)%gpr(k)
      enddo
   enddo

   open(10,file=trim(obshistfile),err=300)
   do m=1,100000
      read(10,'(tr21,a3,tr1,a4,tr11,f4.2,tr36,f6.1)',end=300) var,well,rel_err,min_err
      write(*,'(a,a3,a,a4,a,f5.3,a,f8.3)') ' var=',var,' well=',well,' rel_err=',rel_err,' min_err=',min_err
      select case (var)
      case('OPR')
         l=1
      case('WPR')
         l=2
      case('GPR')
         l=3
      end select
      read(well(4:4),'(i1)')i

      do k=1,nx
         if (obs(k,l,i) /= 0.0) then
            stddev(k,l,i)=max(abs(rel_err*obs(k,l,i)),min_err)
         else
            stddev(k,l,i)=0.0
         endif
      enddo
   enddo
   300 close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Saving perturbations for schedule parcing
   C=reshape(B,(/ nx, nrdata, nrwells, nrens /) )
   do j=1,nrens
      cens(:)='    '
      j1=j-1
      if (j1<10)                 write(cens(1:1),'(i1.1)')j1
      if (j1>=10 .and. j1<100)   write(cens(1:2),'(i2.2)')j1
      if (j1>=100 .and. j1<1000) write(cens(1:3),'(i3.3)')j1
      open(11,file=trim(controlpath)//'CONTROL_'//trim(cens)) 
      do i=1,nrwells
         do l=1,nrdata
            do k=1,nx
               C(k,l,i,j)=stddev(k,l,i)*C(k,l,i,j)
               write(11,'(e13.5)')C(k,l,i,j)
            enddo
         enddo
      enddo
      close(11)
   enddo

   open(12,file=trim(controlpath)//'control.txt')
   do i=1,nrwells
      write(12,'(a,a,i4)')rates(i)%wellname,':OPR:',nx
      write(12,'(a,a,i4)')rates(i)%wellname,':WPR:',nx
      write(12,'(a,a,i4)')rates(i)%wellname,':GPR:',nx
      write(12,*)
   enddo
   close(12)

   open(12,file=trim(controlpath)//'wells.txt')
   do i=1,nrwells
      write(12,'(a)')rates(i)%wellname
   enddo
   close(12)

end program
