;
; CLASS tony_polyline  DEFINITION
; $Revision $
;
;
FUNCTION tony_polyline::INIT,npoints=npoints, $
                             points=points,   $
                             radii=radii,     $
                             tangents=tangents, $
                             linkage=linkage,   $
                             normals=normals,   $
                             binormals=binormals, $
                             xyz0=xyz0, $
                             xyz1=xyz1, $
                             periodic=periodic, $
                             _EXTRA=e
  COMPILE_OPT HIDDEN

  if arg_present(npoints) then begin
    self.npoints=npoints
  endif else begin
    if arg_present(points) then self.npoints=(size(points))[1]
  endelse
  if arg_present(points) then self.points=PTR_NEW(points)
  if arg_present(radii) then self.radii=PTR_NEW(radii)
  if arg_present(tangents) then self.tangents=PTR_NEW(tangents)
  if arg_present(linkage) then begin
   self.linkage=PTR_NEW(linkage)
  endif else begin
   self.linkage=PTR_NEW(lonarr(self.npoints,2))
   for i=0,self.npoints-1 do begin
     (*self.linkage)[i,0]=i+1
     (*self.linkage)[i,1]=i-1
   endfor
   (*self.linkage)[0,1]=-1
   (*self.linkage)[self.npoints-1,0]=-1
  endelse

  if arg_present(normals) then self.normals=PTR_NEW(normals)
  if arg_present(binormals) then self.binormals=PTR_NEW(binormals)

  if keyword_set(periodic) then self.periodic=1

  if arg_present(xyz0) then self.xyz0=xyz0
  if arg_present(xyz1) then self.xyz1=xyz1

  self.segments=PTR_NEW(lonarr(npoints))
  (*self.segments)[*]=0

  self.Lxyz=self.xyz1-self.xyz0

  RETURN, 1
END
;
;
;
;
;
PRO tony_polyline::CLEANUP
  COMPILE_OPT HIDDEN
  ret=self->free_points()
  ret=self->free_linkage()
  ret=self->free_normals()
  ret=self->free_binormals() 
  ret=self->free_tangents()
  ret=self->free_radii()
  ret=self->free_segments()
END
;
FUNCTION tony_polyline::calc_tangents,force=force
  COMPILE_OPT HIDDEN

  if (keyword_set(force) or (not self->has_tangents())) then begin
    if not self->free_tangents() then return, 0
    self.tangents=PTR_NEW(fltarr(self.npoints,3))
    for pnt=0,self.npoints-1 do begin
      if ((*self.linkage)[pnt,0] lt 0) then begin 
        ds=(*self.points)[pnt,0:2]-(*self.points)[(*self.linkage)[pnt,1],0:2] 
        (*self.tangents)[pnt,0:2]=ds
      endif else if ((*self.segments)[pnt] ne (*self.segments)[(*self.linkage)[pnt,0]]) then begin
        ds=(*self.points)[pnt,0:2]-(*self.points)[(*self.linkage)[pnt,1],0:2] 
        (*self.tangents)[pnt,0:2]=ds
      endif else if ((*self.linkage)[pnt,1] lt 0) then begin
        ds=(*self.points)[(*self.linkage)[pnt,0],0:2]-(*self.points)[pnt,0:2]
        (*self.tangents)[pnt,0:2]=ds
      endif else if ((*self.segments)[pnt] ne (*self.segments)[(*self.linkage)[pnt,1]]) then begin
        ds=(*self.points)[(*self.linkage)[pnt,0],0:2]-(*self.points)[pnt,0:2]
        (*self.tangents)[pnt,0:2]=ds
      endif else begin
        ds=(*self.points)[(*self.linkage)[pnt,0],0:2]-(*self.points)[(*self.linkage)[pnt,1],0:2] 
        (*self.tangents)[pnt,0:2]=ds/2.
      endelse
    endfor
    (*self.tangents)[*,*]=(*self.tangents)/spread(sqrt(dot2(*self.tangents)),1,3)
  endif

  return, 1
END
; 
FUNCTION tony_polyline::calc_normals,force=force
  COMPILE_OPT HIDDEN

  if (keyword_set(force) or (not self->has_normals())) then begin
    if (not self->has_tangents()) then ret=self->calc_tangents()
    if not self->free_normals() then return, 0
    self.normals=PTR_NEW(fltarr(self.npoints,3))
    for pnt=0,self.npoints-1 do begin
      ds=[0.,0.,0.]
      fix_right=0
      fix_left=0
      leftpnt=pnt
      rightpnt=pnt

      ittr=0L
      while sqrt(dot2(ds)) eq 0. do begin
        if fix_right eq 0 then begin
          if ((*self.linkage)[rightpnt,0] lt 0) then begin
            fix_right=1
          endif else if ((*self.segments)[rightpnt] ne (*self.segments)[(*self.linkage)[rightpnt,0]]) then begin
            fix_right=1
          endif else begin
            rightpnt=(*self.linkage)[rightpnt,0]
          endelse
        endif
        if fix_left eq 0 then begin
          if ((*self.linkage)[leftpnt,1] lt 0) then begin
            fix_left=1
          endif else if ((*self.segments)[leftpnt] ne (*self.segments)[(*self.linkage)[leftpnt,1]]) then begin
            fix_left=1
          endif else begin
            leftpnt=(*self.linkage)[leftpnt,1]
          endelse
        endif
        ds=(*self.tangents)[rightpnt,*]-(*self.tangents)[leftpnt,*] 
        if fix_left ne 0 and fix_right ne 0 then break
        ittr=ittr+1
      endwhile

      if ittr gt 1L then print,"Took ", ittr, "itterations"

      if sqrt(dot2(ds)) eq 0. then begin
        print,sqrt(dot2(ds)),leftpnt,rightpnt,(*self.segments)[leftpnt-2:rightpnt+2]
      endif

      (*self.normals)[pnt,*]=ds
      len=(*self.normals)[pnt,0]^2+(*self.normals)[pnt,1]^2+(*self.normals)[pnt,2]^2
      (*self.normals)[pnt,0:2]=(*self.normals)[pnt,0:2]/len
    endfor
    (*self.normals)[*,0:2]=(*self.normals)[*,0:2]/spread(sqrt(dot2(*self.normals)),1,3)
  endif

  return, 1
END
; 
FUNCTION tony_polyline::calc_binormals,force=force
  COMPILE_OPT HIDDEN

  if (keyword_set(force) or (not self->has_binormals())) then begin
    if (not self->has_tangents()) then ret=self->calc_tangents()
    if (not self->has_normals()) then ret=self->calc_normals()
    if not self->free_binormals() then return, 0
    self.binormals=PTR_NEW(fltarr(self.npoints,3))
    (*self.binormals)[*,0]=(*self.normals)[*,2]*(*self.tangents)[*,1]-(*self.normals)[*,1]*(*self.tangents)[*,2]
    (*self.binormals)[*,1]=(*self.normals)[*,0]*(*self.tangents)[*,2]-(*self.normals)[*,2]*(*self.tangents)[*,0]
    (*self.binormals)[*,2]=(*self.normals)[*,1]*(*self.tangents)[*,0]-(*self.normals)[*,0]*(*self.tangents)[*,1]
    (*self.binormals)[*,0:2]=(*self.binormals)/spread(dot2(*self.binormals),1,3)
    (*self.binormals)[*,0:2]=(*self.binormals)[*,0:2]/spread(sqrt(dot2(*self.binormals)),1,3)
  endif

  return, 1
END
; 
FUNCTION tony_polyline::construct_tube,radius=radius,nphi=nphi,maxds=maxds,untwist=untwist
  default,radius,1.
  default,nphi,8
  default,maxds,1E37

  dphi=2.*!pi/nphi

  if (not self->has_radii()) then ret=self->set_constant_radius(radius)
  if (not self->has_normals()) then ret=self->calc_normals()
  if (not self->has_binormals()) then ret=self->calc_binormals()

  count=0L
  point=0L
  finished=0
  vertices=fltarr(3,nphi*self.npoints)
  triangles=lonarr(nphi*self.npoints*10)
  discs=lonarr(nphi*self.npoints)
  discs[*]=-1L
  free_vertex=0L
  free_triangle=0L
  untwist_offset=0.
  lastpoint=-1L
  while (finished eq 0) do begin
    if count gt 0 then begin
      if (*self.segments)[point] ne (*self.segments)[lastpoint] then begin
        count=0
        untwist_offset=0.
      endif
    endif

    if discs[point] lt 0 then begin
      if lastpoint ge 0 then begin
        thisnormal=(*self.normals)[point,0:2]
        lastnormal=(*self.normals)[lastpoint,0:2]
        crossprod=[0.,0.,0.]
        crossprod[0]=lastnormal[1]*thisnormal[2]-thisnormal[1]*lastnormal[2]
        crossprod[1]=lastnormal[2]*thisnormal[0]-thisnormal[2]*lastnormal[0]
        crossprod[2]=lastnormal[0]*thisnormal[1]-thisnormal[0]*lastnormal[1]

        dotprod=thisnormal[0]*lastnormal[0]+thisnormal[1]*lastnormal[1]+thisnormal[2]*lastnormal[2]

        twist=min([max([dotprod,-1.]),1.])
        twistang=acos(twist)
; if not finite(twist) then begin
;  print,"NAN Found",lastpoint
; endif        
        tangent=(*self.tangents)[lastpoint,0:2]

        direction=crossprod[0]*tangent[0]+crossprod[1]*tangent[1]+crossprod[2]*tangent[2]
        if direction ge 0 then direction=1. else direction=-1.

        untwist_offset=untwist_offset+direction*twistang
        untwist_offset=untwist_offset-(2.*!dpi*round(untwist_offset/!dpi))
        if not keyword_set(untwist) then untwist_offset=0.
        ;untwist_offset=0.
;        print,[twist,acos(twist),untwist_offset, direction]
      endif 

      discs[point]=free_vertex
      for phipnt=0,nphi-1 do begin
        centre=(*self.points)[point,0:2]
        normal=(*self.normals)[point,0:2]
        binormal=(*self.binormals)[point,0:2]
        rad=(*self.radii)[point]
        newpnt=centre+(sin(dphi*phipnt+untwist_offset)*normal+cos(dphi*phipnt+untwist_offset)*binormal)*rad
        vertices[0:2,free_vertex]=newpnt
        free_vertex=free_vertex+1
      endfor 
    endif


    if (count gt 0) then begin
      for phipnt=0,nphi-1 do begin
        leftseg=phipnt mod nphi
        rightseg=(phipnt+1) mod nphi
        triangles[free_triangle]=3
        triangles[free_triangle+1]=discs[lastpoint]+leftseg
        triangles[free_triangle+2]=discs[point]+leftseg
        triangles[free_triangle+3]=discs[lastpoint]+rightseg
        free_triangle=free_triangle+4
        triangles[free_triangle]=3
        triangles[free_triangle+1]=discs[lastpoint]+rightseg
        triangles[free_triangle+2]=discs[point]+leftseg
        triangles[free_triangle+3]=discs[point]+rightseg
        free_triangle=free_triangle+4
      endfor
      if point le 0 then finished=1
    endif

    lastpoint=point
 
    point=(*self.linkage)[point,0] 
    count=count+1
    if (point lt 0) then finished=1
  endwhile

  triangles=triangles[0:free_triangle-1] 
  vertices=vertices[0:2,0:free_vertex-1] 
 
  return,CREATE_STRUCT(['vertices','triangles'],PTR_NEW(vertices),PTR_NEW(triangles))
END
;
FUNCTION tony_polyline::set_constant_radius,radius
  COMPILE_OPT HIDDEN

  ret=self->free_radii()
  self.radii=PTR_NEW(fltarr(self.npoints))
  (*self.radii)[0:self.npoints-1]=radius
  return,1
END
;
FUNCTION tony_polyline::segment_by_maxds,maxds,Lxyz=Lxyz,unfragment=unfragment
  COMPILE_OPT HIDDEN

  default,minsegpnt,5
  default,Lxyz,[1.,1.,1.]
  (*self.segments)[*]=-1L
  if not keyword_set(unfragment) then minsegpnt=0

  segment=0L
  seglength=0L
  segstart=0L

  finished=0
  point=0L
  lastpoint=-1L
  flip=[0,0,0]
  ds=[0.,0.,0.]
  while (finished eq 0) do begin
    if lastpoint ge 0 then begin
      if (*self.segments)[point] ge 0 then begin
        oldsegment=(*self.segments)[point]
        found=where((*self.segments) eq segment)
        if (min(found) ge 0) then begin
          (*self.segments)[found]=oldsegment
        endif 
        segment=segment-1
        break
      endif 

      ds = (*self.points)[point,0:2]-(*self.points)[lastpoint,0:2]
      modds = sqrt(dot2(ds))
      if modds gt maxds then begin
        if seglength gt minsegpnt then begin
          segment=segment+1
          seglength=0L
          segstart=point
        endif else begin 
         presegstart=(*self.linkage)[segstart,1]
         ds = (*self.points)[presegstart,0:2]-(*self.points)[segstart,0:2]
         flip[0:2]=0
         if abs(ds[0]) gt maxds then begin
           if ds[0] gt 0. then flip[0]=-1 else flip[0]=1
         endif 
         if abs(ds[1]) gt maxds then begin
           if ds[1] gt 0. then flip[1]=-1 else flip[1]=1
         endif 
         if abs(ds[2]) gt maxds then begin
           if ds[2] gt 0. then flip[2]=-1 else flip[2]=1
         endif 
         found=segstart
         while found ne point do begin
              delta=flip*1.*Lxyz
              (*self.points)[found,0:2]=((*self.points)[found,0:2])-delta
          found=(*self.linkage)[found,0] 
          endwhile
          ;segment=segment-1

          ds = (*self.points)[point,0:2]-(*self.points)[lastpoint,0:2]
          modds = sqrt(dot2(ds))
          if modds gt maxds then begin
            segment=segment+1
            seglength=0L
            segstart=point
          endif
        endelse
      endif
      if (point eq 0) then break
    endif

    seglength=seglength+1
    (*self.segments)[point]=segment

    lastpoint=point
    point=(*self.linkage)[point,0] 
    if (point lt 0) then break
  endwhile
  return,segment
END
;
FUNCTION tony_polyline::reattach_orphan_segments,xyz0=xyz0,xyz1=xyz1,minsegment=minsegment
  COMPILE_OPT HIDDEN

  default,minsegment,5
  default,xyz0,[0.,0.,0.]
  default,xyz1,[1.,1.,1.]
  Lxyz=xyz1-xyz0

  nsegs=max((*self.segments))+1

  if nsegs eq 1 then return,1
 
  for seg=0,nsegs-1 do begin
    points=where((*self.segments) eq seg)
    if min(points) ge 0 then begin
      segstart=min(points)
      segend=max(points)
      while (*self.segments)[segend  ] eq (*self.segments)[(*self.linkage)[segend  ,0]] do segend  =(*self.linkage)[segend  ,0]
      while (*self.segments)[segstart] eq (*self.segments)[(*self.linkage)[segstart,1]] do segstart=(*self.linkage)[segstart,1]

      segpre =(*self.linkage)[segstart,1]
      segpost=(*self.linkage)[segend,0]

    endif else begin
      points=where((*self.segments) gt seg)
      if min(points) ge 0 then begin
        (*self.segments)[points]=(*self.segments)[points]-1
        nsegs=nsegs-1
      endif
    endelse  
  endfor

  return,1
END
; 
FUNCTION tony_polyline::get_binormals
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.binormals) then return,(*self.binormals)
  return,-1
END
;
FUNCTION tony_polyline::get_segments
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.segments) then return,(*self.segments)
  return,-1
END
;
FUNCTION tony_polyline::get_normals
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.normals) then return,(*self.normals)
  return,-1
END
;
FUNCTION tony_polyline::get_tangents
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.tangents) then return,(*self.tangents)
  return,-1
END
;
FUNCTION tony_polyline::get_points
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.points) then return,(*self.points)
  return,-1
END
;
FUNCTION tony_polyline::get_radii
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.radii) then return,(*self.radii)
  return,-1
END
;
FUNCTION tony_polyline::set_radii,radii
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.radii) then ret=self->free_radii()
  self.radii=PTR_NEW(radii)
  return,-1
END
;
FUNCTION tony_polyline::free_linkage
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.linkage) then PTR_FREE,self.linkage
  return,1
END
; 
FUNCTION tony_polyline::free_points
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.points) then PTR_FREE,self.points
  return,1
END
; 
FUNCTION tony_polyline::free_segments
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.segments) then PTR_FREE,self.segments
  return,1
END
; 
FUNCTION tony_polyline::free_tangents
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.tangents) then PTR_FREE,self.tangents
  return,1
END
; 
FUNCTION tony_polyline::free_segments
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.segments) then PTR_FREE,self.segments
  return,1
END
; 
FUNCTION tony_polyline::free_normals
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.normals) then PTR_FREE,self.normals
  return,1
END
; 
FUNCTION tony_polyline::free_radii
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.radii) then PTR_FREE,self.radii
  return,1
END
; 
FUNCTION tony_polyline::free_binormals
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.binormals) then PTR_FREE,self.binormals
  return,1
END
;
FUNCTION tony_polyline::has_linkage
  COMPILE_OPT HIDDEN

  return, ((self.npoints gt 0) and PTR_VALID(self.linkage))
END
; 
FUNCTION tony_polyline::has_points
  COMPILE_OPT HIDDEN

  return, ((self.npoints gt 0) and PTR_VALID(self.points))
END
; 
FUNCTION tony_polyline::has_tangents
  COMPILE_OPT HIDDEN

  return, ((self.npoints gt 0) and PTR_VALID(self.tangents))
END
; 
FUNCTION tony_polyline::has_normals
  COMPILE_OPT HIDDEN

  return, ((self.npoints gt 0) and PTR_VALID(self.normals))
END
; 
FUNCTION tony_polyline::has_radii
  COMPILE_OPT HIDDEN

  return, ((self.npoints gt 0) and PTR_VALID(self.radii))
END
; 
FUNCTION tony_polyline::has_binormals
  COMPILE_OPT HIDDEN

  return, ((self.npoints gt 0) and PTR_VALID(self.binormals))
END
; 
PRO tony_polyline__DEFINE
  COMPILE_OPT HIDDEN
  struct={tony_polyline, $
            npoints:0L,      $  
            xyz0: fltarr(3), $
            xyz1: fltarr(3), $
            Lxyz: fltarr(3), $
            periodic: 0,     $
            points:PTR_NEW(),   $  
            radii:PTR_NEW(),    $  
            linkage:PTR_NEW(),  $  
            tangents:PTR_NEW(), $  
            normals:PTR_NEW(),  $  
            segments:PTR_NEW(),  $  
            binormals: PTR_NEW()  }
END
;

