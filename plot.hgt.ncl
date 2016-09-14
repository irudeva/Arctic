; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin

;Dates     MM DD
tstart = (/ 8, 5/)
tend   = (/ 8, 31/)


;************************************************
; read in atmospheric data
;************************************************
  a = addfile("../DATA/hgt.ncep.2016.nc","r")
  lev = a->level

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
  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                   "Oct","Nov","Dec"/)

  do t=0,dimsizes(time)-1
   if (date(t,1).eq.tstart(0))then
    if (date(t,2).eq.tstart(1))then
     if (date(t,3).eq.0)then
      t1=t
     end if
    end if
   end if
   if (date(t,1).eq.tend(0))then
    if (date(t,2).eq.tend(1))then
     t2=t
     break
    end if
   end if
  end do

  a1 = addfile("../DATA/slp.ncep.2016.nc","r")
  slp = a1->slp(:,:,:)
  slp=slp/100.
;************************************************
; create plot
;************************************************
  wks_type                 = "png"                  ;-- plot output type
  ;wks_type@wkBackgroundColor = "black"              ;-- set workstation background to black (or grey18)
  wks_type@wkWidth         =  2500                  ;-- for presentations
  wks_type@wkHeight        =  2000                  ;-- for presentations

  wks = gsn_open_wks(wks_type,"../output/hgt_slp.2016")             ; send graphics to PNG file
  gsn_define_colormap(wks, "wgne15")

  res            = True                         ; plot mods desired
  res@gsnDraw    =  False                 ;-- don't draw plot yet
  res@gsnFrame   =  False                 ;-- don't advance frame

  res@gsnPaperOrientation = "landscape"
  ;res@gsnMaximize         = True        ; must include w/ Paper Orientation


  res@gsnPolar   = "NH"                         ; specify the hemisphere
  res@mpMinLatF  = 60                           ; minimum lat to plot
  res@mpFillOn   = False

  res@cnFillOn          = True                  ; color fill
  ;res@cnLevelSpacingF   = 3.                    ; interval spacing
  res@cnLinesOn         = False                 ; no contour lines

  res@gsnSpreadColorEnd     = -3

  res@tiMainFontHeightF    =  0.018                  ;-- main title font size
  res@tiMainOffsetYF       =  0.06                  ;-- move title upward
  res@gsnLeftStringFontHeightF       =  0.016                  ;-- main title font size


  res@vpYF                 =  0.95                  ;-- set viewport Y position
  res@vpHeightF            =  0.8                  ;-- set viewport width
  res@vpWidthF             =  0.35                   ;-- set viewport height

  res1                     =  res                   ;-- resource setting for left plot
  res1@vpXF                =  0.049                 ;-- set viewport X position 1. panel
  ;res1@tiMainString        = "Hgt - " + date2(t1)
  res1@gsnLeftString       = "Hgt "+lev(il)+" hPa"


  res2                     =  res                   ;-- resource setting for right plot
  res2@vpXF                =  0.503                 ;-- set viewport X position 2. panel
  ;res2@tiMainString        = "Slp - " + date2(t1)
  res2@gsnLeftString       = "Slp "
  res2@gsnRightString      = "hPa "

  tires1                   =  True                  ;-- text resources title string
  tires1@txFontColor       = "red"                ;-- change to white
  tires1@txFontHeightF     =  0.015                 ;-- text font size
  tires1@txJust            = "CenterCenter"         ;-- text justification



; assign plot array for the data sets displayed on two globes
  plot = new(2,graphic)


  plot(0) = gsn_csm_contour_map_polar(wks,hgt(t1,:,:),res1)
  plot(1) = gsn_csm_contour_map_polar(wks,slp(t1,:,:),res2)

  draw(plot)

  date_str =  sprinti("%0.2iZ ",  tointeger(date(t1,3))) + \
           sprinti("%0.2i ", tointeger(date(t1,2))) + \
           month_abbr(tointeger(date(t1,1))) + " "  + \
           sprinti("%0.4i", tointeger(date(t1,0)))
  gsn_text_ndc(wks,date_str,0.45,0.85,tires1) ;-- center title string
  print(date_str)

  frame(wks)


  do t=t1+1,t2

; The more efficient way to do animations is to use "setvalues"
; to change the data only.
;
; Use "(/.../)" to avoid re-reading the coordinate arrays, which are
; the same for every timestep.
;

;    setvalues plot@data
;      "sfDataArray" : (/a->hgt(t,il,:,:)/)
;    end setvalues
;    setvalues plot@data
;      "sfDataArray" : (/a1->slp(t,:,:)/)
;    end setvalues

    plot(0) = gsn_csm_contour_map_polar(wks,hgt(t,:,:),res1)
    plot(1) = gsn_csm_contour_map_polar(wks,slp(t,:,:),res2)


; Use setvalues to change the title too.
;    setvalues plot
;      "tiMainString"  : "Hgt - + date2(t)
;    end setvalues

    draw(plot)

    date_str =  sprinti("%0.2iZ ",  tointeger(date(t,3))) + \
             sprinti("%0.2i ", tointeger(date(t,2))) + \
             month_abbr(tointeger(date(t,1))) + " "  + \
             sprinti("%0.4i", tointeger(date(t,0)))
    gsn_text_ndc(wks,date_str,0.45,0.85,tires1) ;-- center title string
    print(date_str)


    frame(wks)
;    maximize_output(wks,False)

  end do

end
