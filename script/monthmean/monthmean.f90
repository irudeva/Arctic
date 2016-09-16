Program pdf

character*100 ftrk,fout
integer yr,nreg

integer,allocatable:: cyr(:),m(:),num(:),loc(:)
real, allocatable:: x(:),y(:),p(:),c(:),dp(:),d(:),rd(:),rr(:)
real mp,md,mdp,mrd,sigp,sigrd,sigd,sigdp,mc,sigc,mrr,sigrr,mn70,mn70cont,sign70,sign70cont

integer i,nc,n,nm,ndm(12)
integer ntrk, len,len70,maxlen
integer, allocatable,dimension(:)::n70,n70cont

read (*,'(a)') ftrk
read (*,*) yr
read (*,'(a)') fout
read (*,*) nreg

ndm=[31,28,31,30,31,30,31,31,30,31,30,31]
if(any([1980,1984,1988,1992,1996,2000,2004,2008,2012,2016]==yr))ndm(2)=29

if(yr<1999)then
nyr=yr-1900
elseif(yr>2000) then
nyr=yr-2000
elseif (yr==1999)then
nyr=19
elseif (yr==2000)then
nyr=20
endif

open(30,file=fout)
write(30,'("mmyyyy ncyc ncyc/an"7(6x,a,7x,"sig"),4(6x,a))')"   P"," Dsq","  dp","   R","  RR"," n70","n70c","IDsq","  ID","  IR"," IRR"

do 4 nm=1,12
open (10,file=ftrk,action='read')
open (20,file="dummy")
open (21,file="dummy1")

do 1 i=1,66
1 read(10,*)


nc=0
ntrk=0
do
 ptrk=0
 read(10,'(55x,i4)',end=22) n
 allocate (cyr(n),m(n),x(n),y(n),p(n),c(n),dp(n),rd(n),d(n),loc(n))
 read(10,*)
 read(10,*)
 do 21 i=1,n
 read (10,1000)cyr(i),m(i),x(i),y(i),p(i),c(i),dp(i),rd(i)
1000 format(12x,2i2,28x,6f9.3)
! d(i)=rd(i)*rd(i)*c(i)/4
 if (cyr(i)==nyr.and.n>=5.and.y(i)>=70..and.m(i)==nm)then
 if (nreg==0.or.(nreg==1.and.(x(i)<-135..or.x(i)>90.)).or.&
&(nreg==2.and.(x(i)<-20..or.x(i)>160.)).or.(nreg==3.and.x(i)<-80.).or.(nreg==4.and.(x(i)<-110.)))then
  nc=nc+1
  if (ptrk==0)then
   ntrk=ntrk+1
   ptrk=1
  end if
  write(20,1002) nc,ntrk,i,n, p(i),c(i),dp(i),rd(i)
1002 format(2i7,2i4,4f10.2)
 endif
 endif
 21 continue

 ! write length
 loc=0
 if(ptrk==1)then
  do 23 i=1,n
   if (cyr(i)==nyr.and.y(i)>=70.)then
   if (nreg==0.or.(nreg==1.and.(x(i)<-135..or.x(i)>90.)).or.&
&(nreg==2.and.(x(i)<-20..or.x(i)>160.)).or.(nreg==3.and.x(i)<-80.).or.(nreg==4.and.(x(i)<-110.)))then
    loc(i)=1
   end if
   end if
  23 continue
  len70=count(loc==1)
  len=0
  maxlen=0
  do 24 i=1,n
   if(loc(i)==1)then
    len=len+1
   else
    if (len>=maxlen)maxlen=len
    len=0
   end if
   if(i==n) maxlen=maxval((/maxlen,len/) )
  24 end do
  write(21,1005) ntrk,n,len70,maxlen

  1005 format(4i7)
 end if

 read(10,*)
 deallocate (cyr,x,y,p,c,dp,rd,d,m,loc)

end do

22 close(10)
rewind(20);rewind(21);print*, nm,yr,nc

allocate (num(nc),p(nc),dp(nc),d(nc),rd(nc),c(nc),rr(nc))

if (nc>1)then
do 3 i=1,nc
read (20,1004) num(i),p(i),c(i),dp(i),rd(i)
3 rr(i)=rd(i)**2
1004 format (7x,i7,8x,5f10.2)
close(20,status='delete')
!close(20)

allocate(n70(ntrk),n70cont(ntrk))
do 31 i=1,ntrk
31 read (21,'(14x,2i7)') n70(i),n70cont(i)
close(21,status='delete')
!close(21)


call moment(p,nc,mp,sigp)
!call moment(d,nc,md,sigd)
call moment(dp,nc,mdp,sigdp)
call moment(rd,nc,mrd,sigrd)
call moment(c,nc,mc,sigc)
call moment(rr,nc,mrr,sigrr)
call moment(float(n70),ntrk,mn70,sign70)
call moment(float(n70cont),ntrk,mn70cont,sign70cont)


write(30,1003)nm, yr,nc,float(nc)/ndm(nm),mp,sigp,mc,sigc,mdp,sigdp,mrd,sigrd,mrr,sigrr,mn70,sign70,mn70cont,sign70cont,sum(c),&
& sum(dp),sum(rd),sum(rr)
endif  ! nc>1
4 deallocate(num,p,d,dp,rd,c,rr,n70,n70cont)
close(30)
1003 format(i2,i4,i5,f8.2,20f10.2)

stop

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
