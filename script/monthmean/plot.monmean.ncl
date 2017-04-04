; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

d="ncep"
lev="slp"
npts = 10
col = (/"black","red","orange","blue","paleturquoise","green","olivedrab"/)
mon = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
var = (/"ncyc_anal","slp","dp","radius","n70c"/)
var_title = (/"Ncyc/analysis","av slp","av dp","av radius","av N time steps > 70N"/)
varn = (/2,3,7,9,15/)

y1=1979
y2=2016

ssn=(/"djf","mam","jja","son"/)
fm=(/-1,2,5,8/)
lm=(/1,3,7,10/)


data=new((/y2-y1+1,12/),float)
data@_FillValue = -9999.0
ssndata=new((/y2-y1+1,dimsizes(ssn)/),float)
data@_FillValue = data@_FillValue

pltarry=new((/2,dimsizes(data(:,0))/),typeof(data),data@_FillValue)

 do iv = 0,4
  do loc  = 0,0
   fout = "../../output/trkmonmean/"+var(iv)+"."+lev+"."+d+"."+y1+"-"+y2

   data=data@_FillValue

; ***  Plotting  *******************************************************************
  wks  = gsn_open_wks("ps",fout)

  res                    = True

  ;res@tiMainString       = "Ray Path"
  ;res@tiMainFontHeightF  = 0.013
  ;res@gsnCenterString    = "("+lat(loc)+"N; "+lon(loc)+"E)      period = "+p+"days     root = "+root
  ;res@gsnCenterStringFontHeightF = 0.008

  yr = ispan(y1,y2,1)
  yr1 = ispan(1,y2-y1+1,1)

; *** read data ******************************************************************

     do y = y1,y2

      fin = "../../output/trkmonmean/av."+lev+"."+d+"."+y+"N"
      if (fileexists(fin)) then
      print("Data from "+fin)
       getvar = readAsciiTable (fin,21,"float",1)
       data(y-y1,:dimsizes(getvar(:,varn(iv)))-1) = getvar(:,varn(iv))
; fill the last values with FillValue due to incomplete data
       if(y.eq.2016)then
        data(y-y1,8:)= data@_FillValue
       end if
;seasons
       do sn=0,3
        if(ssn(sn).eq."djf")then
         if(y.ne.y1)then
          ssndata(y-y1,sn)=avg((/data(y-y1-1,11),data(y-y1,0),data(y-y1,1)/))
         else
          ssndata(y-y1,sn)=ssndata@_FillValue
         end if
        else
         ssndata(y-y1,sn)=avg(data(y-y1,fm(sn):lm(sn)))
         ;print (data(y-y1,fm(sn):lm(sn))+" ")
         ;print ("ssn = "+ssndata(y-y1,sn))
        end if
       end do

       delete( getvar)

      end if  ;file exist
      end do  ; yr


;*** XY plot *************************************************************************
print("Start plotting")
;*** monthly
     do m=0,11
; regression
      rc = regline(yr-y1+1, data(:,m))
      rc@units = "6h time step/yr"
      pltarry(0,:)=data(:,m)
      pltarry(1,:)=rc*((yr-y1+1)-rc@xave)+rc@yave

      res@tiMainString = var_title(iv)+"  "+mon(m)+" "+y1+" - "+y2

      plot = gsn_csm_xy(wks,yr,pltarry,res)
;---Drawing the plot will draw the attached polylines and markers.
;      draw(plot)
;      frame(wks)
      end do

;*** seasonal

     do ns=0,3
; regression
      rc = regline(yr-y1+1, ssndata(:,ns))
      rc@units = "6h time step/yr"
      pltarry(0,:)=ssndata(:,ns)
      pltarry(1,:)=rc*((yr-y1+1)-rc@xave)+rc@yave

      res@tiMainString = var_title(iv)+"  "+ssn(ns)+" "+y1+" - "+y2
      if(ssn(ns).eq."djf")then
       res@tiMainString = var_title(iv)+"  "+ssn(ns)+" "+(y1+1)+" - "+y2
      end if

      plot = gsn_csm_xy(wks,yr,pltarry,res)
;---Drawing the plot will draw the attached polylines and markers.
;      draw(plot)
;      frame(wks)
      print ("ok")
      end do


; ***  End plotting  ***************************************************************
  end do ; loc
 end do  ;var

end
