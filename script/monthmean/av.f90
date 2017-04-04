Program average

character*100 fin,din,fout,fssn(6),HS*1
!character(LEN=:),allocatable:: din1
integer fy,ly,y,ny,len,m,lm,nreg
integer, allocatable:: nc(:,:)
real,allocatable,dimension(:,:):: ncan,p,c,rd,d,rr,n70,n70cont
!!! mn70 - average number of time steps that a cyclone spent >70N
!!! mn70 - average number of time steps that a cyclone continuously stayed >70N

real,allocatable:: IC(:,:),IDp(:,:),IR(:,:),IRR(:,:)
real avnc,avncan,avp,avc,avd,avdr,avdp,avsD,avCarea,avarea
real signc,signcan,sigp,sigc,sigd,sigdp,sigdr,sigsD,sigCarea,sigarea
integer ndm(12)

read(*,'(a)') din;len=len_trim(din)
read(*,'(a)') HS
read(*,'(a)') fout
read(*,'(a)') fssn

read(*,*) fy
read(*,*) ly
read(*,*) nreg
ny=ly-fy+1
allocate (p(ny,12),c(ny,12),d(ny,12),rd(ny,12),rr(ny,12),nc(ny,12),ncan(ny,12),IDsq(ny,12),&
&ID(ny,12),IR(ny,12),IRR(ny,12),n70(ny,12),n70cont(ny,12))

do 1 y=fy,ly

if(nreg==0)write(fin,'(a,"av.",i4,a)') din(1:len),y,HS
if(nreg==1)write(fin,'(a,"av.",i4,a,".90_225")') din(1:len),y,HS
if(nreg==2)write(fin,'(a,"av.",i4,a,".160_340")') din(1:len),y,HS
if(nreg==3)write(fin,'(a,"av.",i4,a,".180_280")') din(1:len),y,HS
if(nreg==4)write(fin,'(a,"av.",i4,a,".180_250")') din(1:len),y,HS

open(10,file=fin,action='read')
read(10,*)

ndm=[31,28,31,30,31,30,31,31,30,31,30,31]
if(any([1980,1984,1988,1992,1996,2000,2004,2008,2012,2016]==y))ndm(2)=29

do 2 m=1,12;print*, m
lm=m
2 read(10,1000,end=1) nc(y-fy+1,m),ncan(y-fy+1,m),p(y-fy+1,m),c(y-fy+1,m),d(y-fy+1,m),rd(y-fy+1,m),rr(y-fy+1,m),&
& n70(y-fy+1,m),n70cont(y-fy+1,m),IC(y-fy+1,m),IDp(y-fy+1,m),IR(y-fy+1,m),IRR(y-fy+1,m)


1000 format(6x,i5,f8.2,7(f10.2,10x),4f10.2)

1 close(10)

open(20,file=fout)
write(20,'("year",13(3x,a7,7x,"sig"))')"Nc","Na","P","dsqP","D","R","RR","N70","N70cont","IDsq","ID","IR","IRR"

do 3 m=1,12 ;print*, m
if(m<lm.or.lm==12)then
call moment(float(nc(:,m)),ny,avnc,signc)
call moment(ncan(:,m),ny,avncan,signcan)
call moment(p(:,m),ny,avp,sigp)
call moment(c(:,m),ny,avc,sigc)
!call moment(dp(:,m),ny,avdp,sigdp)
call moment(d(:,m),ny,avd,sigd)
call moment(rd(:,m),ny,avrd,sigrd)
call moment(rr(:,m),ny,avrr,sigrr)
call moment(n70(:,m),ny,avn70,sign70)
call moment(n79(:,m),ny,avn70cont,sign70cont)
call moment(IDsq(:,m),ny,avIDsq,sigIDsq)
call moment(ID(:,m),ny,avID,sigID)
call moment(IR(:,m),ny,avIR,sigIR)
call moment(IRR(:,m),ny,avIRR,sigIRR)

write(20,1003)m,avnc,signc,avncan,signcan,avp,sigp,avc,sigc,avd,sigd,avrd,sigrd,avrr,sigrr,&
avn70,sign70,avn70cont,sign70cont,&
&avIDsq,sigIDsq,avID,sigID,avIR,sigIR,avIRR,sigIRR
 endif
3 continue
close(20)
1003 format(i2,26f10.2)

open(21,file=fssn)
write(21,'("year",13(6x,a4))')"Nc","Na","P","dsqP","D"," R","RR","N70","N70cont","IDsq","ID","IR","IRR"

do 4 y=1,ly-fy+1 ;print*, y
call moment(float(nc(y,6:8)),3,avnc,signc)
call moment(ncan(y,6:8),3,avncan,signcan)
call moment(p(y,6:8),3,avp,sigp)
call moment(c(y,6:8),3,avc,sigc)
!call moment(dp(y,6:8),3,avdp,sigdp)
call moment(d(y,6:8),3,avd,sigd)
call moment(rd(y,6:8),3,avrd,sigrd)
call moment(rr(y,6:8),3,avrr,sigrr)
call moment(n70(y,6:8),3,avn70,sign70)
call moment(n70cont(y,6:8),3,avn70cont,sign70cont)
call moment(IDsq(y,6:8),3,avIDsq,sigIDsq)
call moment(ID(y,6:8),3,avID,sigID)
call moment(IR(y,6:8),3,avIR,sigIR)
call moment(IRR(y,6:8),3,avIRR,sigIRR)

4 write(21,1004)y+fy-1,avnc,avncan,avp,avc,avd,avrd,avrr,&
!&avIDsq,avID,avIR,avIRR
avn70,avn70cont,&
&sum(IDsq(y,6:8)),sum(ID(y,6:8)),sum(IR(y,6:8)),sum(IRR(y,6:8))
close(21)
1004 format(i4,11f10.2)





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
