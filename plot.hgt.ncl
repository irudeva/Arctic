; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin

;Dates     MM DD
tstart = (/ 8, 1/)
tend   = (/ 8,31/)


;************************************************
; read in atmospheric data
;************************************************
  a = addfile("../DATA/NCEP/hgt.2016.nc","r")
  lev = a->level

  print (lev)

  do il=0,dimsizes(lev)-1
    if (lev(il).eq.300) then
      break
    end if
  end do

  print(il)

  hgt = a->hgt(:,il,:,:)

  time = a->time
  date = cd_calendar(time, 0)
  date2 = cd_calendar(time, -2)

  do t=0,dimsizes(time)-1
   if (date(t,1).eq.tstart(0))then
    if (date(t,2).eq.tstart(1))then
     t1=t
    end if
   end if
   if (date(t,1).eq.tend(0))then
    if (date(t,2).eq.tend(1))then
     t2=t
     break
    end if
   end if
  end do
;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("ps","../DATA/NCEP/hgt.2016")             ; send graphics to PNG file
  gsn_define_colormap(wks, "wgne15")

  res            = True                         ; plot mods desired
  res@gsnPolar   = "NH"                         ; specify the hemisphere
  res@mpMinLatF  = 60                           ; minimum lat to plot
  res@mpFillOn   = False

  res@cnFillOn          = True                  ; color fill
  ;res@cnLevelSpacingF   = 3.                    ; interval spacing
  res@cnLinesOn         = False                 ; no contour lines

  res@gsnSpreadColorEnd     = -3

  res@tiMainString        = "Hgt - " + date2(t1)
  res@gsnLeftString       = "hgt "+lev(il)+" hPa"

  plot = gsn_csm_contour_map_polar(wks,hgt(t1,:,:),res)

  do t=t1+1,t2

; The more efficient way to do animations is to use "setvalues"
; to change the data only.
;
; Use "(/.../)" to avoid re-reading the coordinate arrays, which are
; the same for every timestep.
;
    setvalues plot@data
      "sfDataArray" : hgt(t,:,:)
    end setvalues

; Use setvalues to change the title too.
    setvalues plot
      "tiMainString"  : "Hgt - + date2(t)
    end setvalues

; Draw the plot and advance the frame.
    draw(plot)
    frame(wks)
;    maximize_output(wks,False)

  end do

end
