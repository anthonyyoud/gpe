;
;  $Id: setup_derivs.pro,v 1.1 2006/03/21 14:35:13 mee Exp $
;
;  First derivative d/dx
;  - 6th-order
;  - with ghost cells
;  - on potentially non-equidistant grid
;
;***********************************************************************
pro setup_derivs,object=object,x=x_,y=y_,z=z_,nx=nx_,ny=ny_,nz=nz_,lequidist=lequidist_
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0 
  common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
;
;  Check for degenerate case (no x-extension)
;
  if keyword_set(object) then begin
    x=safe_get_tag(object,"x",default=findgen(7)/6.-0.5)
    y=safe_get_tag(object,"y",default=findgen(7)/6.-0.5)
    z=safe_get_tag(object,"z",default=findgen(7)/6.-0.5)
    nx=safe_get_tag(object,"nx",default=1L)
    ny=safe_get_tag(object,"ny",default=1L)
    nz=safe_get_tag(object,"nz",default=1L)
    lequidist=safe_get_tag(object,"lequidist",default=[1,1,1])
  endif else begin
    default,x_,findgen(7)/6.-0.5
    default,y_,findgen(7)/6.-0.5
    default,z_,findgen(7)/6.-0.5
    default,nx_,1L
    default,ny_,1L
    default,nz_,1L
    default,lequidist_,[1,1,1]
    x=x_
    y=y_
    z=z_
    nx=nx_
    ny=ny_
    nz=nz_
    lequidist=lequidist_
  endelse

  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
;
end
