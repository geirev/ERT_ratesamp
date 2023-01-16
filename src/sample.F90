program sample
! Samples the measurement perturbation matrix E for the ert-Eclipse applications.
! The program reads the ert/observations/obshist.txt file and the schedulefile defined below.
! Only configured for opr, wpr, gpr observations.
! Predefined well names?

   use m_pseudo1D
   use m_fixsample1D
   use m_set_random_seed2
   use m_random
   implicit none

   integer, parameter :: nrens=1000        ! ensemble size
! read observations info assumed stored as follows:
   integer, parameter :: nx=36            ! Number of dates keywords
   integer :: nrwells                     ! number of wells
   integer, parameter :: nrvars=3         ! number of variables per well (OPR, WPR, GPR)
! We will then sample nrwells*nrdata time series of length nx and store these as one realization
! A total of nrens such realizations are computed.

   real :: xlength=35.0                   ! length of time series is 36 months
   real :: cor1=15.0                      ! decorrelation length in months
   character(len=100) :: obshistfile= '/home/AD.NORCERESEARCH.NO/geev/ERT/ert1000/observations/obshistnew.txt'
   character(len=100) :: schedulefile='/home/AD.NORCERESEARCH.NO/geev/Dropbox/Statoil/eclipse/include/history.sch'
   character(len=100) :: priorpath= '/home/AD.NORCERESEARCH.NO/geev/ERT/Prior1000/'
   character(len=200) controlpath


   integer n1                             ! x-dimension of grid
   integer i,j,k,l,m,ic,j1
   integer idat
   integer iwell
   integer ipos(5)
   real dx

   character, parameter :: dq = '"'

   character(len=8) wellname

   type wellrates
      character(len=4)  wellname
      real opr(nx)
      real gpr(nx)
      real wpr(nx)
   end type
   type(wellrates), allocatable :: rates(:),stddev(:)
   type(wellrates), allocatable :: A(:,:)
   real, allocatable :: B(:,:,:)


   character(len=100) cline
   character(len=4) var
   character(len=4) well
   character(len=4) cens
   character(len=4) cnx
   character(len=4) ck
   real rel_err
   real min_err
   logical ex

   character(len=1) tmp
   integer istat
   integer noise
   character(len=1) cnoise
   character(len=200) string1
   character(len=200) string2
   character(len=200) string3
   character(len=200) string4


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
      write(*,'(a)',advance='no')'Overwrite existing wells.txt file (y/n)'
      read(*,*)tmp
   endif
   if (tmp=='y' .or. tmp=='Y') then
      string3="cat rates.txt | sed -e "//dq//"/'/d"//dq//" -e '/DATES/d' | cut -f2 -d# | sort -u > wells.txt"
      print '(a,a)','string2=',trim(string3)
      istat = system( string3 )
   endif

   open(10,file='wells.txt')
      do i=1,1000
         read(10,'(a)',end=100)tmp
      enddo
      100 rewind(10)
      nrwells=i-1
      allocate (rates(nrwells))
      allocate (stddev(nrwells))
      do i=1,nrwells
         read(10,'(a)')rates(i)%wellname
         write(*,'(a)')rates(i)%wellname
      enddo
   close(10)

! Allocating matrices for all wells
   allocate( A(nrwells,nrens) )

   write (*,'(a)',advance='no')'Select type of noise (0-null, 1-red, 2-white, 3-constant): '
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
      do j=1,nrens
      do i=1,nrwells
         A(i,j)%opr(:)=0.0
         A(i,j)%gpr(:)=0.0
         A(i,j)%wpr(:)=0.0
      enddo
      enddo
      cnoise='0'
   case (1)
      print '(a,f8.2)','Simulating pseudo random red noise with decorrelation=',cor1
      allocate(B(nx,nrens,nrvars))
      do i=1,nrwells
         call pseudo1D(B,nx,nrens*nrvars,cor1,dx,n1)
         do j=1,nrens
         do k=1,nx
            A(i,j)%opr(k)=B(k,j,1)
            A(i,j)%wpr(k)=B(k,j,2)
            A(i,j)%gpr(k)=B(k,j,3)
         enddo
         enddo
!            call fixsample1D(A(i,j)%opr,nx,nrens)
      enddo
      cnoise='R'

   case (2)
      print '(a)','Simulating white random noise'
      do j=1,nrens
      do i=1,nrwells
         call random(A(i,j)%opr,nx)
         call random(A(i,j)%wpr,nx)
         call random(A(i,j)%gpr,nx)
      enddo
      enddo
      cnoise='W'

   case (3)
      print '(a)','Simulating random bias noise'
      do j=1,nrens
      do i=1,nrwells
         call random(A(i,j)%opr(1),1)
         call random(A(i,j)%wpr(1),1)
         call random(A(i,j)%gpr(1),1)
         A(i,j)%opr(2:nx)=A(i,j)%opr(1)
         A(i,j)%wpr(2:nx)=A(i,j)%wpr(1)
         A(i,j)%gpr(2:nx)=A(i,j)%gpr(1)
      enddo
      enddo
      cnoise='B'
   end select
   print *,'Sampling completed'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initializing rates to 0
   do i=1,nrwells
      rates(i)%opr(:)=0.0
      rates(i)%wpr(:)=0.0
      rates(i)%gpr(:)=0.0
   enddo

! reading rates and times from rates.txt
   print '(a)','Reading info from rates.txt'
   open (10,file='rates.txt')
   idat=1
   do i=1,10000
      read(10,'(a)',end=200)cline
!      print '(a)',trim(cline)
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
         print '(2i3,tr1,a4,3f10.2)',idat,iwell,rates(iwell)%wellname,rates(iwell)%opr(idat),&
                                          rates(iwell)%wpr(idat),rates(iwell)%gpr(idat)
      endif
   enddo
   200 close(10)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scaling perturbations according to measurement magnitude and error defs in obshistfile
!   do i=1,nrwells
!      do k=1,nx
!         obs(k,1,i) = rates(i)%opr(k)
!         obs(k,2,i) = rates(i)%wpr(k)
!         obs(k,3,i) = rates(i)%gpr(k)
!      enddo
!   enddo

   open(10,file=trim(obshistfile),err=300)
   print '(2a)','Potentially unstable read of: ',trim(obshistfile)
   print '(a)','Check output!'
   do m=1,100000
      read(10,'(tr21,a3,tr1,a4,tr11,f4.2,tr36,f6.1)',end=300) var,well,rel_err,min_err
      write(*,'(a,a3,a,a4,a,f5.3,a,f8.3)') ' var=',var,' well=',well,' rel_err=',rel_err,' min_err=',min_err

      do i=1,nrwells
         if (well == rates(i)%wellname) then
            iwell=i
            exit
         endif
      enddo
      select case (var)
      case('OPR')
         stddev(i)%opr(:)=max(abs(rel_err*rates(iwell)%opr(:)), min_err)
      case('WPR')
         stddev(i)%wpr(:)=max(abs(rel_err*rates(iwell)%wpr(:)), min_err)
      case('GPR')
         stddev(i)%gpr(:)=max(abs(rel_err*rates(iwell)%gpr(:)), min_err)
      end select

   enddo
   300 close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Saving perturbations for schedule parcing
   controlpath(:)=' '
   controlpath=trim(priorpath)//'control'//cnoise//'/'
   string4="mkdir -p "//trim(controlpath)
   print '(a,a)','string4=',trim(string4)
   istat = system( string4 )
   do j=1,nrens
      cens(:)='    '
      j1=j-1
      if (j1<10)                   write(cens(1:1),'(i1.1)')j1
      if (j1>=10   .and. j1<100)   write(cens(1:2),'(i2.2)')j1
      if (j1>=100  .and. j1<1000)  write(cens(1:3),'(i3.3)')j1
      if (j1>=1000 .and. j1<10000) write(cens(1:4),'(i4.4)')j1

      open(11,file=trim(controlpath)//'CONTROL_'//trim(cens))
      do i=1,nrwells
         do l=1,3
         do k=1,nx
            ck(:)='    '
            if (k<10)                  write(ck(1:1),'(i1.1)')k
            if (k>=10   .and. k<100)   write(ck(1:2),'(i2.2)')k
            if (k>=100  .and. k<1000)  write(ck(1:3),'(i3.3)')k
            if (k>=1000 .and. k<10000) write(ck(1:4),'(i4.4)')k
            if (l==1)  write(11,'(4a,f13.3)')trim(rates(i)%wellname),'_WOPR_',ck,' ',stddev(i)%opr(k)*A(i,j)%opr(k)
            if (l==2)  write(11,'(4a,f13.3)')trim(rates(i)%wellname),'_WGPR_',ck,' ',stddev(i)%gpr(k)*A(i,j)%gpr(k)
            if (l==3)  write(11,'(4a,f13.3)')trim(rates(i)%wellname),'_WWPR_',ck,' ',stddev(i)%wpr(k)*A(i,j)%wpr(k)
!            if (l==1)  write(11,'(E13.5)')stddev(i)%opr(k)*A(i,j)%opr(k)
!            if (l==2)  write(11,'(E13.5)')stddev(i)%gpr(k)*A(i,j)%gpr(k)
!            if (l==3)  write(11,'(E13.5)')stddev(i)%wpr(k)*A(i,j)%wpr(k)
         enddo
         enddo
      enddo
      close(11)
   enddo

   open(12,file=trim(controlpath)//'control.txt')
      cnx(:)=' '
      if (nx<10)                   write(cnx(1:1),'(i1.1)')nx
      if (nx>=10   .and. nx<100)   write(cnx(1:2),'(i2.2)')nx
      if (nx>=100  .and. nx<1000)  write(cnx(1:3),'(i3.3)')nx
      if (nx>=1000 .and. nx<10000) write(cnx(1:4),'(i4.4)')nx
      print *,'nx',nx,cnx
      do i=1,nrwells
         write(12,'(a,a,a)')trim(rates(i)%wellname),' OPR ',trim(cnx)
         write(12,'(a,a,a)')trim(rates(i)%wellname),' WPR ',trim(cnx)
         write(12,'(a,a,a)')trim(rates(i)%wellname),' GPR ',trim(cnx)
      enddo
   close(12)

   open(12,file=trim(controlpath)//'wells.txt')
   do i=1,nrwells
      write(12,'(a)')rates(i)%wellname
   enddo
   close(12)

end program
