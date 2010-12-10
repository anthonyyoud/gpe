function read_data, file=ifile, phase=phase, vx=vx, vy=vy, vz=vz, dbl=dbl, $
  f77=f77

  isdouble = keyword_set(dbl)

  t = isdouble ? 0D : 0.0
  nx = 0L
  ny = 0L
  nz = 0L
  nytot = 0L
  nztot = 0L
  nyprocs = 0L
  nzprocs = 0L
  nprocs = 0L
  jsta = 0L
  jend = 0L
  ksta = 0L
  kend = 0L
  eye = dcomplex(0.0D, 1.0D)

  spawn, 'ls -1d proc* | wc -l', nprocs
  
  for k=0,nprocs(0)-1 do begin
    if k lt 10 then begin
      data_dir = 'proc0'+str(k)+'/'
    endif else begin
      data_dir = 'proc'+str(k)+'/'
    endelse
    file = data_dir+ifile

    if keyword_set(f77) then begin
      openr, 5, file, /F77
    endif else begin
      openr, 5, file
    endelse
    readu, 5, t
    readu, 5, nx, nytot, nztot
    readu, 5, nyprocs, nzprocs
    readu, 5, jsta, jend, ksta, kend
    ny = jend-jsta+1
    nz = kend-ksta+1

    tmp = isdouble ? dcomplexarr(nx,ny,nz) : complexarr(nx,ny,nz)
    if k eq 0 then begin
      psi = isdouble ? dcomplexarr(nx,nytot,nztot) : complexarr(nx,nytot,nztot)
      dens = isdouble ? dblarr(nx,nytot,nztot) : fltarr(nx,nytot,nztot)
      if keyword_set(phase) then begin
        phase = isdouble ? dblarr(nx,nytot,nztot) : fltarr(nx,nytot,nztot)
      endif
      if keyword_set(vx) then begin
        vx = isdouble ? dblarr(nx,nytot,nztot) : fltarr(nx,nytot,nztot)
        dpsidx = isdouble ? dcomplexarr(nx,nytot,nztot) : $
          complexarr(nx,nytot,nztot)
        dpsistardx = isdouble ? dcomplexarr(nx,nytot,nztot) : $
          complexarr(nx,nytot,nztot)
      endif
      if keyword_set(vy) then begin
        vy = isdouble ? dblarr(nx,nytot,nztot) : fltarr(nx,nytot,nztot)
        dpsidy = isdouble ? dcomplexarr(nx,nytot,nztot) : $
          complexarr(nx,nytot,nztot)
        dpsistardy = isdouble ? dcomplexarr(nx,nytot,nztot) : $
          complexarr(nx,nytot,nztot)
      endif
      if keyword_set(vz) then begin
        vz = isdouble ? dblarr(nx,nytot,nztot) : fltarr(nx,nytot,nztot)
        dpsidz = isdouble ? dcomplexarr(nx,nytot,nztot) : $
          complexarr(nx,nytot,nztot)
        dpsistardz = isdouble ? dcomplexarr(nx,nytot,nztot) : $
          complexarr(nx,nytot,nztot)
      endif
    endif
    
    x = isdouble ? dblarr(nx) : fltarr(nx)
    y = isdouble ? dblarr(nytot) : fltarr(nytot)
    z = isdouble ? dblarr(nztot) : fltarr(nztot)
    readu, 5, tmp
    readu, 5, x
    readu, 5, y
    readu, 5, z
    close, 5

    psi[*,jsta:jend,ksta:kend] = tmp[*,*,*]
  endfor

  dens = abs(psi)^2
  data_struct = CREATE_STRUCT( $
    ['t','nx','ny','nz','x','y','z','dens'], $
    t, nx, nytot, nztot, x, y, z, dens)
  if keyword_set(phase) then begin
    phase = atan(psi, /phase)
    phase = phase*(dens ge 1.0D-4)
    data_struct = CREATE_STRUCT(['phase'], phase, data_struct)
  endif
  if keyword_set(vx) then begin
    ;xx = rebin(reform(x,nx,1,1),nx,nytot,nztot)
    dpsidx = derx(x,psi)
    ;dpsidx = deriv(xx,psi)
    dpsistardx = derx(x,conj(psi))
    ;dpsistardx = deriv(xx,conj(psi))
    vx = (dens gt 0.01D) * (-0.5D*eye * (conj(psi)*dpsidx - psi*dpsistardx) / dens)
    data_struct = CREATE_STRUCT(['vx'], double(vx), data_struct)
  endif
  if keyword_set(vy) then begin
    ;yy = rebin(reform(y,nytot,1,1),nx,nytot,nztot)
    dpsidy = dery(y,psi)
    ;dpsidy = deriv(yy,psi)
    dpsistardy = dery(y,conj(psi))
    ;dpsistardy = deriv(yy,conj(psi))
    vy = (dens gt 0.01D) * (-0.5D*eye * (conj(psi)*dpsidy - psi*dpsistardy) / dens)
    data_struct = CREATE_STRUCT(['vy'], double(vy), data_struct)
  endif
  if keyword_set(vz) then begin
    ;zz = rebin(reform(z,nztot,1,1),nx,nytot,nztot)
    dpsidz = derz(z,psi)
    ;dpsidz = deriv(zz,psi)
    dpsistardz = derz(z,conj(psi))
    ;dpsistardz = deriv(zz,conj(psi))
    vz = (dens gt 0.01D) * (-0.5D*eye * (conj(psi)*dpsidz - psi*dpsistardz) / dens)
    data_struct = CREATE_STRUCT(['vz'], double(vz), data_struct)
  endif

  return, data_struct
end
