Program average

integer, parameter:: narg = 13
character*10 listarg(narg)
character*100 fin,dir,fout,fssn(6),HS*1,dataset*4,lev*3
!character(LEN=:),allocatable:: dir1
integer fy,ly,y,ny,len,m,lm,nreg
integer, allocatable:: nc(:,:)
real,allocatable:: arg(:,:,:)
real av(narg), sig(narg)
integer ndm(12)

read(*,'(a)') dir;len=len_trim(dir);print*, dir, len

read(*,'(a)') HS;print *,'HS=',HS
read(*,'(a)') dataset;print *,'d=',dataset
read(*,'(a)') lev;print *,'lev=',lev
read(*,'(a)') fout;print *,'fout=',fout
read(*,'(a)') fssn;print *,'ssn=',ssn

read(*,*) fy
read(*,*) ly;print *,"years=",fy,"-", ly
read(*,*) nreg; print *, "nreg=",nreg
ny=ly-fy+1

write(dir,'(a,"av.",a,".",a,".")') dir(1:len),lev,dataset;len=len_trim(dir)

allocate (arg(narg,ny,12))

do 1 y=fy,ly

if(nreg==0)write(fin,'(a,i4,a)') dir(1:len),y,HS
if(nreg==1)write(fin,'(a,i4,a,".90_225")') dir(1:len),y,HS
if(nreg==2)write(fin,'(a,i4,a,".160_340")') dir(1:len),y,HS
if(nreg==3)write(fin,'(a,i4,a,".180_280")') dir(1:len),y,HS
if(nreg==4)write(fin,'(a,i4,a,".180_250")') dir(1:len),y,HS

print *, 'reading from ',fin
open(10,file=fin,action='read')
read(10,"(a6,a5,a8,7(a,10x),4a)")(listarg(i),i=1,narg)


ndm=[31,28,31,30,31,30,31,31,30,31,30,31]
if(any([1980,1984,1988,1992,1996,2000,2004,2008,2012,2016]==y))ndm(2)=29

do 2 m=1,12
lm=m
2 read(10,1000,end=1) (arg(i,y-fy+1,m),i=1,narg)

1000 format(6x,f5.0,f8.2,7(f10.2,10x),4f10.2)

1 close(10)

open(20,file=fout)
print *, 'writing to ',fout

write(20,'("year",13(6x,a4,7x,"sig"))')(listarg(i),i=1,narg)

do 3 m=1,12
 if(m<lm.or.lm==12)then
  do i=1,narg
   call moment(arg(i,:,m),ny,av(i),sig(i))
  end do
  write(20,1003)m,(av(i),sig(i),i=1,narg)
 endif
3 continue
close(20)
1003 format(i2,26f10.2)

!do i=1,6
!open(21,file=fssn(i))
!if(i==1)then
!write(21,'("year",12(6x,a4))')"Nc","Na","P","dsqP","D"," R","RR","IDsq","ID","IR","IRR"
!else
!write(21,'("year",7(6x,a4))')"Nc","Na","P","dsqP","D"," R","RR"
!endif

!do 4 y=1,ly-fy+1 ;print*, y
!if(i==1)then
!call moment(float(nc(y,6:8)),3,avnc,signc)
!call moment(ncan(y,6:8),3,avncan,signcan)
!call moment(p(y,6:8),3,avp,sigp)
!call moment(c(y,6:8),3,avc,sigc)
!!call moment(dp(y,6:8),3,avdp,sigdp)
!call moment(d(y,6:8),3,avd,sigd)
!call moment(rd(y,6:8),3,avrd,sigrd)
!call moment(rr(y,6:8),3,avrr,sigrr)
!call moment(IDsq(y,6:8),3,avIDsq,sigIDsq)
!call moment(ID(y,6:8),3,avID,sigID)
!call moment(IR(y,6:8),3,avIR,sigIR)
!call moment(IRR(y,6:8),3,avIRR,sigIRR)

! write(21,1004)y+fy-1,avnc,avncan,avp,avc,avd,avrd,avrr,&
!!&avIDsq,avID,avIR,avIRR
!&sum(IDsq(y,6:8)),sum(ID(y,6:8)),sum(IR(y,6:8)),sum(IRR(y,6:8))

!else
!ii=i+3
!avnc=float(nc(y,ii))
!avncan=ncan(y,ii)
!avp=p(y,ii)
!avc=c(y,ii)
!avd=d(y,ii)
!avrd=rd(y,ii)
!avrr=rr(y,ii)

! write(21,1005)y+fy-1,avnc,avncan,avp,avc,avd,avrd,avrr
!endif
!4 continue
!close(21)
!enddo
1004 format(i4,11f10.2)
1005 format(i4,7f10.2)





end

      SUBROUTINE moment(data,n,ave,sdev)
!      SUBROUTINE moment(data,n,ave,adev,sdev,var,skew,curt)
      INTEGER n
      REAL adev,ave,curt,sdev,skew,var,data(n)
      INTEGER j
      REAL p,s,ep
      if(n.le.1)then
       Write(*,'( "n must be at least 2 in moment")')
       read(*,*)
      endif
      s=0.
      do 11 j=1,n
        s=s+data(j)
11    continue
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        adev=adev+abs(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
12    continue
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=sqrt(var)
      if(var.ne.0.)then
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3.
      else
        write(*,'("no skew or kurtosis when zero variance in moment")')
        read(*,*)
      endif
      return
      END
