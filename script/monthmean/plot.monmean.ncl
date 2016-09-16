; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

d="ncep"
lev="slp"
npts = 10
col = (/"black","red","orange","blue","paleturquoise","green","olivedrab"/)
dum1 = new(6, graphic)
dum2 = dum1
var = (/"n70c"/)
varn = (/15/)

y1=1979
y2=2015

data=new((/y2-y1+1,12/),float)

 do iv = 0,0
  do loc  = 0,0
   fout = "../../output/trkmonmean/"+var(iv)

   data=0.

; ***  Plotting  *******************************************************************
  wks  = gsn_open_wks("png",fout)

  res                    = True

  ;res@tiMainString       = "Ray Path"
  ;res@tiMainFontHeightF  = 0.013
  ;res@gsnCenterString    = "("+lat(loc)+"N; "+lon(loc)+"E)      period = "+p+"days     root = "+root
  ;res@gsnCenterStringFontHeightF = 0.008

  yr = ispan(y1,y2,1)

; *** read data ******************************************************************

     do y = y1,y2

      fin = "../../output/trkmonmean/av."+lev+"."+d+"."+y+"N"
      if (fileexists(fin)) then
      print("Data from "+fin)
       getvar = readAsciiTable (fin,21,"float",1)
       data(y-y1,:) = getvar(:,varn(iv))

       delete(getvar)

       end if  ; if filein exist
      end do  ; k

;*** XY plot *************************************************************************

     do m=8,8
      plot = gsn_csm_xy(wks,yr,data(:,m),res)
;---Drawing the plot will draw the attached polylines and markers.
;      draw(plot)
;      frame(wks)
      print ("ok")
  end do


; ***  End plotting  ***************************************************************
  end do ; loc
 end do  ;var

end
