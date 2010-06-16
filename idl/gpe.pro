function contourdata, data, x0, y0, z0, phase=phase, vx=vx, vy=vy, vz=vz, $
  dir=dir
  ; Get the contour data if a contour plot is requested.
  if not keyword_set(dir) then begin
    ; Assume contour plot in xy-plane (dir = 'z').
    if keyword_set(phase) then begin
      contourdata = data.phase[*,*,z0]
    endif else if keyword_set(vx) then begin
      contourdata = data.vx[*,*,z0]
    endif else if keyword_set(vy) then begin
      contourdata = data.vy[*,*,z0]
    endif else if keyword_set(vz) then begin
      contourdata = data.vz[*,*,z0]
    endif else begin
      contourdata = data.dens[*,*,z0]
    endelse
  endif else begin
    if dir eq 'x' then begin
      ; Contour plot in yz-plane.
      if keyword_set(phase) then begin
        contourdata = data.phase[x0,*,*]
      endif else if keyword_set(vx) then begin
        contourdata = data.vx[x0,*,*]
      endif else if keyword_set(vy) then begin
        contourdata = data.vy[x0,*,*]
      endif else if keyword_set(vz) then begin
        contourdata = data.vz[x0,*,*]
      endif else begin
        contourdata = data.dens[x0,*,*]
      endelse
    endif else if dir eq 'y' then begin
      ; Contour plot in xz-plane.
      if keyword_set(phase) then begin
        contourdata = data.phase[*,y0,*]
      endif else if keyword_set(vx) then begin
        contourdata = data.vx[*,y0,*]
      endif else if keyword_set(vy) then begin
        contourdata = data.vy[*,y0,*]
      endif else if keyword_set(vz) then begin
        contourdata = data.vz[*,y0,*]
      endif else begin
        contourdata = data.dens[*,y0,*]
      endelse
    endif else if dir eq 'z' then begin
      ; Contour plot in xy-plane.
      if keyword_set(phase) then begin
        contourdata = data.phase[*,*,z0]
      endif else if keyword_set(vx) then begin
        contourdata = data.vx[*,*,z0]
      endif else if keyword_set(vy) then begin
        contourdata = data.vy[*,*,z0]
      endif else if keyword_set(vz) then begin
        contourdata = data.vz[*,*,z0]
      endif else begin
        contourdata = data.dens[*,*,z0]
      endelse
    endif
  endelse

  return, reform(contourdata)
end

pro get_pos, data, x0, y0, z0, xpos=xpos, ypos=ypos, zpos=zpos
  ; Get the plot position for a slice or contour plot.
  x0=(data.nx)/2
  y0=(data.ny)/2
  z0=(data.nz)/2
  if keyword_set(xpos) then begin
    x0 = fix(x0 + xpos*(data.nx)/(data.x(data.nx-1)-data.x(0)))
  endif
  if keyword_set(ypos) then begin
    y0 = fix(y0 + ypos*(data.ny)/(data.y(data.ny-1)-data.y(0)))
  endif
  if keyword_set(zpos) then begin
    z0 = fix(z0 + zpos*(data.nz)/(data.z(data.nz-1)-data.z(0)))
  endif
end

pro get_filenames, digits, i, ifile, idir, iprefix, isuffix, ofile, odir, $
  oprefix
  ; Get in and out filenames.

  ;let 'n' be a string corresponding to the integer 'digits'
  n = string(digits)
  ;format code - 'n' integers with any blanks filled with zeros
  format_str = '(I' + n + '.' + n +')'
  ;the digits in the i'th filename converted to a string
  index = string(format = format_str, i)
  ;input file with all whitespace removed
  ifile = idir + iprefix + strcompress(index, /remove_all) + '.' + isuffix
  ;output file with all whitespace removed
  ofile = odir + oprefix + strcompress(index, /remove_all)
end

pro gpe, nstart, nend, cntr=cntr, c_anim=c_anim, phase=phase, $
  skip=skip, slice=slice, dir=dir, xpos=xpos, ypos=ypos, zpos=zpos, vx=vx, $
  vy=vy, vz=vz, eps=eps, dbl=dbl, xsize=xsize, ysize=ysize, _extra=_extra

;no. digits (including leading 0's) in filename
digits = 7
;no. of start filename
nstart = nstart ;6850L
;no. of end filename
nend = nend ;6850L
;diff. between filenames
default, skip, 1L
;input dir of files
idir = './'
;output dir of files
odir = idir+'images/'
;prefix of input filenames
iprefix = 'dens'
;prefix of output filenames
if keyword_set(cntr) then begin
  if keyword_set(phase) then begin
    oprefix = 'con_phase'
  endif else if keyword_set(vx) then begin
    oprefix = 'con_vx'
  endif else if keyword_set(vy) then begin
    oprefix = 'con_vy'
  endif else if keyword_set(vz) then begin
    oprefix = 'con_vz'
  endif else begin
    oprefix = 'con_dens'
  endelse
endif else begin
  oprefix = 'iso_dens'
endelse
;suffix of input filenames
isuffix = 'dat'
default, xsize, 640
default, ysize, 480

; Rotation matrix.
tmat=[[      0.73090459,     -0.79570162,    -0.034048109,       0.0000000], $
      [      0.36181851,      0.29057828,      0.97630603,       0.0000000], $
      [     -0.70949820,     -0.67152441,      0.46280538,       0.0000000], $
      [       0.0000000,       0.0000000,       0.0000000,       1.0000000]]

;tmat=[[1.0,  0.0, 0.0, 0.0], $
;      [0.0,  0.0, 1.0, 0.0], $
;      [0.0, -1.0, 0.0, 0.0], $
;      [0.0,  0.0, 0.0, 1.0]]

; Contour animation.
if keyword_set(c_anim) then begin
  ;The current graphics device, eg in Linux this is 'X' for X11
  thisDevice = !D.Name
  ;Set the output device to be the Z-buffer and copy colour table from IDL's
  ;internal colour table
  Set_Plot, 'Z', /COPY
  ;Set the Z-buffer resolution, and set the Z-buffer to a standard 2D device
  xres=640
  zres=480
  Device, Set_Resolution=[xres,zres], Z_Buffer=0
  ;Clear any previous contents of the Z-buffer
  Erase
endif

loadct, 39, /silent
;loadct, 3, /silent
TVLCT, r, g, b, /Get

; If doing an animation then check to make sure an output directory exists,
; find the range of the data over all snapshots, and define contour levels.
if nstart ne nend then begin
  if not file_test(odir, /directory) then begin
    print, 'ERROR: Directory '+odir+' does not exist.'
    stop
  endif
  if keyword_set(c_anim) then begin
    print, '*****Finding range of data for contour levels.  Please wait...'
    if keyword_set(dbl) then begin
      varminmax = dblarr(2)
      varminmax = [0.0D, 0.0D]
    endif else begin
      varminmax = fltarr(2)
      varminmax = [0.0, 0.0]
    endelse
    for j = nstart, nend, skip do begin
      print, j
      get_filenames, digits, j, ifile, idir, iprefix, isuffix, ofile, odir, $
        oprefix
      data=read_data(file=ifile, phase=phase, vx=vx, vy=vy, vz=vz, dbl=dbl)
      get_pos, data, x0, y0, z0, xpos=xpos, ypos=ypos, zpos=zpos
      plotdata = contourdata(data, x0, y0, z0, phase=phase, $
        vx=vx, vy=vy, vz=vz, dir=dir)
      tmpminmax = minmax(plotdata)
      if tmpminmax[0] lt varminmax[0] then begin
        varminmax[0] = tmpminmax[0]
      endif
      if tmpminmax[1] gt varminmax[1] then begin
        varminmax[1] = tmpminmax[1]
      endif
    endfor
    print, '*****Done.'
    print, '*****Range of selected data:', varminmax[0], varminmax[1]
    lev = 255
    q = findgen(lev-1)*((abs(varminmax[1]) + abs(varminmax[0])) / $
      (lev-1)) - abs(varminmax[0])
  endif
endif

; Loop over all data files from 'nstart' to 'nend' in steps of 'skip', and do
; the plots.
for i = nstart, nend, skip do begin
  get_filenames, digits, i, ifile, idir, iprefix, isuffix, ofile, odir, $
    oprefix
  print,'******Datafile ',ifile,'******'

  data = read_data(file=ifile, phase=phase, vx=vx, vy=vy, vz=vz, dbl=dbl)
  
  if not keyword_set(cntr) and not keyword_set(slice) then begin
    ; 3D isosurface.

    print, data.t

    vizit, data.dens, $
           tmat = tmat, $
           xsize = xsize, ysize = ysize, $
           bgcolor = [255,255,255], $
           level = 0.75, $
           index = 250, $
           dx = (data.x(data.nx-1)-data.x(0))/data.nx, $
           dy = (data.y(data.ny-1)-data.y(0))/data.ny, $
           dz = (data.z(data.nz-1)-data.z(0))/data.nz, $
           ;xx = data.x, $
           ;yy = data.y, $
           ;zz = data.z, $
           filename = ofile, $
           drawbbox = 1, $
           bboxthick = 1, $ ;6, $
           drawcontent = 1, $
           drawdatabox = 0, $
           eps = eps, $
           _extra = _extra
  endif else if keyword_set(slice) then begin
    ; 1D slice through the data.
    get_pos, data, x0, y0, z0, xpos=xpos, ypos=ypos, zpos=zpos
    ylabel = "!M!!!7w!X!M!!!E2!X"
    if keyword_set(eps) then begin
      !p.font = 0
      set_plot, 'ps'
      ylabel = "|!My!7!X|!E2"
      device, /encapsulated, /isolatin1, filename=ofile+'.eps', $
        xsize=xsize/72., ysize=ysize/72., /inches
    endif
    if not keyword_set(dir) then begin
      ; Assume a slice through the x-direction.
      plotdata=data.dens[*,y0,z0]
      xtitle='!8x!X'
      xrange = [data.x[0],data.x[data.nx-1]]
      xdata = data.x
    endif else begin
      if dir eq 'x' then begin
        ; x-direction.
        plotdata=data.dens[*,y0,z0]
        xtitle='!8x!X'
        xrange = [data.x[0],data.x[data.nx-1]]
        xdata = data.x
      endif else if dir eq 'y' then begin
        ; y-direction.
        plotdata=data.dens[x0,*,z0]
        xtitle='!8y!X'
        xrange = [data.y[0],data.y[data.ny-1]]
        xdata = data.y
      endif else if dir eq 'z' then begin
        ; z-direction.
        plotdata=data.dens[x0,y0,*]
        xtitle='!8z!X'
        xrange = [data.z[0],data.z[data.nz-1]]
        xdata = data.z
      endif
    endelse
    plot, xdata, plotdata, xtitle=xtitle, xrange=xrange, xstyle=1, $
      charsize=1.5, xmargin=7, ymargin=3, _extra=_extra
    xyouts, 0.0, 0.5, ylabel, charsize=1.5, /normal
    if keyword_set(eps) then begin
      device, /close
      !p.font = -1
      set_plot, 'x'
    endif
  endif else begin
    ; 2D contour plot.
    get_pos, data, x0, y0, z0, xpos=xpos, ypos=ypos, zpos=zpos
    plotdata = contourdata(data, x0, y0, z0, phase=phase, $
      vx=vx, vy=vy, vz=vz, dir=dir)
    datarange = minmax(plotdata)
    ; Set the color for the axes.
    color = 0
    ; Use 255 colors for contour levels.
    c_colors = indgen(255)
    print,"Data range:", datarange[0], datarange[1]
    if nstart eq nend then begin
      ; Just do a single plot.
      if keyword_set(eps) then begin
        !p.font = 0
        set_plot, 'ps'
        device, /color, bits=8, /encapsulated, /isolatin1, $
          filename=ofile+'.eps', xsize=xsize/72., ysize=ysize/72., /inches
      endif
      contour, plotdata, data.x, data.y, /iso, /fill, $
        xrange=[data.x[0],data.x[data.nx-1]], xstyle=1, $
        yrange=[data.y[0],data.y[data.ny-1]], ystyle=1, $
        color=color, background=255, nlevels=255, c_colors=c_colors, $
        _extra=_extra 
      colorbar, /vertical, min=datarange[0], max=datarange[1], ncolors=255, $
        color=color
      if keyword_set(eps) then begin
        device, /close
        !p.font = -1
        set_plot, 'x'
      endif
    endif else begin
      ; Use the pre-computed contour levels, so the levels are the same for all
      ; plots.
      contour, plotdata, data.x, data.y, levels=q, /iso, /fill, $
        xrange=[data.x[0],data.x[data.nx-1]], xstyle=1, $
        yrange=[data.y[0],data.y[data.ny-1]], ystyle=1, $
        color=color, background=255, c_colors=c_colors, _extra=_extra
      colorbar, /vertical, min=varminmax[0], max=varminmax[1], ncolors=255, $
        color=color
    endelse
    if keyword_set(c_anim) then begin
      ;return the entire display device area as a byte array
      snapshot = TVRD()
      ;return the RGB values from the internal colour table into TVLCT
      ;setup a 24-bit image as a byte array
      image24 = BytArr(3, xres, zres)
      ;read the colour values from the image into correct part of 'image24'
      image24[0,*,*] = r[snapshot]
      image24[1,*,*] = g[snapshot]
      image24[2,*,*] = b[snapshot]
      ;write the image to a png file
      Write_PNG, ofile+'.png', image24
    endif
  endelse
endfor
end
