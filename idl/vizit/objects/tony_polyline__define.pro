;
; CLASS tony_polyline  DEFINITION
; $Revision $
;
;
FUNCTION tony_polyline_cross,a,b
  v1=a
  v2=b

  v1=reform(v1)
  v2=reform(v2)
;  if (n_elements(v1) < 3) then v1=[v1, 0., 0., 0.]
;  if (n_elements(v2) < 3) then v2=[v2, 0., 0., 0.]
;
;  if (n_elements(v1) > 3) then v1=v1[0:2]
;  if (n_elements(v2) > 3) then v2=v2[0:2]

  switch ((size(v1))[0]) OF
    0: return,0
    1: return, [ v1[1]*v2[2]-v2[1]*v1[2], $
                 v1[2]*v2[0]-v2[2]*v1[0], $
                 v1[0]*v2[1]-v2[0]*v1[1] ]
    2: return, [ v1[*,1]*v2[*,2]-v2[*,1]*v1[*,2], $
                 v1[*,2]*v2[*,0]-v2[*,2]*v1[*,0], $
                 v1[*,0]*v2[*,1]-v2[*,0]*v1[*,1] ]
    3: return, [ v1[*,*,1]*v2[*,*,2]-v2[*,*,1]*v1[*,*,2], $
                 v1[*,*,2]*v2[*,*,0]-v2[*,*,2]*v1[*,*,0], $
                 v1[*,*,0]*v2[*,*,1]-v2[*,*,0]*v1[*,*,1] ]
  endswitch

  return, 0.
END
;
FUNCTION tony_polyline_dot,a,b
  v1=a
  v2=b

  v1=reform(v1)
  v2=reform(v2)
;  help,v1
;  if (n_elements(v1) < 3) then v1=[v1, 0., 0., 0.]
;  if (n_elements(v2) < 3) then v2=[v2, 0., 0., 0.]
;
;  if (n_elements(v1) > 3) then v1=v1[0:2]
;  if (n_elements(v2) > 3) then v2=v2[0:2]

  switch ((size(v1))[0]) OF
    0: return,0
    1: return,v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
    2: return,v1[*,0]*v2[*,0]+v1[*,1]*v2[*,1]+v1[*,2]*v2[*,2]
    3: return,v1[*,*,0]*v2[*,*,0]+v1[*,*,1]*v2[*,*,1]+v1[*,*,2]*v2[*,*,2]
  endswitch
END
;
FUNCTION tony_polyline_dot2,a
  return,tony_polyline_dot(a,a)
END
;
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

  if keyword_set(npoints) then begin
    self.npoints=npoints
  endif else begin
    if keyword_set(points) then self.npoints=(size(points))[1]
  endelse
  if keyword_set(points) then self.points=PTR_NEW(points)
  if keyword_set(radii) then self.radii=PTR_NEW(radii)
  if keyword_set(tangents) then self.tangents=PTR_NEW(tangents)
  if keyword_set(linkage) then begin
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

  if keyword_set(normals) then self.normals=PTR_NEW(normals)
  if keyword_set(binormals) then self.binormals=PTR_NEW(binormals)

  if keyword_set(periodic) then self.periodic=1

  if keyword_set(xyz0) then self.xyz0=xyz0
  if keyword_set(xyz1) then self.xyz1=xyz1

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
    self.tangents=PTR_NEW(dblarr(self.npoints,3))
    for pnt=0L,self.npoints-1 do begin
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
    (*self.tangents)[*,*]=(*self.tangents)/spread(sqrt(tony_polyline_dot2(*self.tangents)),1,3)
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
    for pnt=0L,self.npoints-1 do begin
      ds=[0d,0d,0d]
      fix_right=0
      fix_left=0
      leftpnt=pnt
      rightpnt=pnt

      ittr=0L
      while ((sqrt(tony_polyline_dot2(ds)) eq 0.) and (ittr lt 1000)) do begin
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
        if fix_left ne 0 and fix_right ne 0 then begin
          ds=[0.,0.,0.]
          break
         endif
        ittr=ittr+1
      endwhile
      tangent=(*self.tangents)[pnt,*]
      dstangent=ds
      dstangent[0]=ds[1]*tangent[2]-tangent[1]*ds[2]
      dstangent[1]=tangent[0]*ds[2]-ds[0]*tangent[2]
      dstangent[2]=ds[0]*tangent[1]-tangent[1]*ds[0]
;      if (tony_polyline_dot2(ds) eq 0.) then begin
;        print,"TANGENT == DS"
;      endif
      if ((tony_polyline_dot2(ds) eq 0) or (tony_polyline_dot2(dstangent) eq 0)) then begin
;       print, "Found straight section"
       ds=tony_polyline_cross([1.,0.,0.],tangent)
       if (tony_polyline_dot2(ds) eq 0) then ds=tony_polyline_cross([0.,1.,0.],tangent)
      endif 
;      	if ittr gt 1L then print,"Took ", ittr, " iterations"

      (*self.normals)[pnt,*]=ds
      len=(*self.normals)[pnt,0]^2+(*self.normals)[pnt,1]^2+(*self.normals)[pnt,2]^2
      (*self.normals)[pnt,0:2]=(*self.normals)[pnt,0:2]/len
    endfor
     
    (*self.normals)[*,0:2]=(*self.normals)[*,0:2]/spread(sqrt(tony_polyline_dot2(*self.normals)),1,3)
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
    (*self.binormals)[*,0]=(*self.normals)[*,1]*(*self.tangents)[*,2]-(*self.tangents)[*,1]*(*self.normals)[*,2]
    (*self.binormals)[*,1]=(*self.tangents)[*,0]*(*self.normals)[*,2]-(*self.normals)[*,0]*(*self.tangents)[*,2]
    (*self.binormals)[*,2]=(*self.normals)[*,0]*(*self.tangents)[*,1]-(*self.tangents)[*,1]*(*self.normals)[*,0]
;    (*self.binormals)[*,0:2]=(*self.binormals)/spread(tony_polyline_dot2(*self.binormals),1,3)
    (*self.binormals)[*,0:2]=(*self.binormals)[*,0:2]/spread(sqrt(tony_polyline_dot2(*self.binormals)),1,3)
  endif

  return, 1
END
; 
FUNCTION tony_polyline::construct_tube,radius=radius,nphi=nphi,maxds=maxds,untwist=untwist
  default,radius,1.
  default,nphi,8

  dphi=2.*!pi/nphi

  if (not self->has_radii()) then ret=self->set_constant_radius(radius)
  if (not self->has_normals()) then ret=self->calc_normals()
  if (not self->has_binormals()) then ret=self->calc_binormals()

  count=0L
  point=0L
  looped=0
  finished=0
  vertices=fltarr(3,nphi*self.npoints)
  normals=fltarr(3,nphi*self.npoints)
  triangles=lonarr(nphi*self.npoints*8)
  discs=lonarr(nphi*self.npoints)
  discs[*]=-1L
  free_vertex=0L
  free_triangle=0L
  untwist=0.
  untwist_offset=0L
  lastpoint=-1L
  flip_binormal=1.

  while (finished eq 0) do begin
    if count gt 0 then begin
      if (*self.segments)[point] ne (*self.segments)[lastpoint] then begin
        count=0
        untwist=0.
        untwist_offset=0L
        flip_binormal=1.
      endif
    endif

    if discs[point] lt 0 then begin
      thisnormal=reform((*self.normals)[point,0:2])
      thistangent=reform((*self.tangents)[point,0:2])
      thiscentre=reform((*self.points)[point,0:2])
      thisbinormal=reform((*self.binormals)[point,0:2])
      thisradius=(*self.radii)[point]

      if lastpoint ge 0 then begin
        lasttangent=reform((*self.tangents)[lastpoint,0:2])
        lastnormal=reform((*self.normals)[lastpoint,0:2])
        lastbinormal=reform((*self.binormals)[lastpoint,0:2])
        lastcentre=reform((*self.points)[lastpoint,0:2])

        BN1N2=tony_polyline_dot(lastbinormal,thisnormal)
        N1N2=tony_polyline_dot(lastnormal,thisnormal)
        T1T2=tony_polyline_dot(lasttangent,thistangent)
        modN1N2=sqrt(dot2(thisnormal)*dot2(lastnormal))

 ;       if T1T2 lt 0. then begin
 ;         print,"TUBE DOUBLES BACK!!!"
 ;         theta=0.
 ;       endif 
        theta=acos(min([max([N1N2,-1.]),1.]))
        if ((BN1N2 ge 0.) and (N1N2 ge 0.)) then begin
           untwist=untwist-theta 
        endif else if ((BN1N2 lt 0.) and (N1N2 ge 0.)) then begin
           untwist=untwist+theta 
        endif else if ((BN1N2 lt 0.) and (N1N2 lt 0.)) then begin
           untwist=untwist+theta 
        endif else if ((BN1N2 ge 0.) and (N1N2 lt 0.)) then begin
           untwist=untwist-theta 
        endif
        if not keyword_set(untwist) then untwist=0.

        nseg=round(untwist/dphi)
        untwist=untwist-(nseg*dphi)
        untwist_offset=(nseg+(nphi*2)) mod nphi
      endif 

      discs[point]=free_vertex
      for phipnt=0,nphi-1 do begin
        effectivenormal=cos((dphi*phipnt)+untwist)*thisnormal+sin((dphi*phipnt)+untwist)*thisbinormal
        effectivenormal = effectivenormal/sqrt(dot2(effectivenormal))
        newpnt=thiscentre+(effectivenormal*spread(thisradius,0,3))
        vertices[0:2,free_vertex]=newpnt
        normals[0:2,free_vertex]=effectivenormal
        free_vertex=free_vertex+1
      endfor 
    endif


    tesselate=0
    if (count gt 0) then begin
      tesselate=1
      if (n_elements(maxds) eq 1) then begin
        if (tony_polyline_dot2(lastcentre-thiscentre) gt (maxds^2)) then tesselate=0
      endif
    endif

    if (tesselate eq 1) then begin
      for phipnt=0,nphi-1 do begin
        leftseg=phipnt mod nphi
        rightseg=(phipnt+1) mod nphi
        lastleftseg=(leftseg-untwist_offset+2*nphi) mod nphi
        lastrightseg=(rightseg-untwist_offset+2*nphi) mod nphi
        triangles[free_triangle]=3
        triangles[free_triangle+1]=discs[lastpoint]+lastleftseg
        triangles[free_triangle+2]=discs[point]+leftseg
        triangles[free_triangle+3]=discs[lastpoint]+lastrightseg
        free_triangle=free_triangle+4
        triangles[free_triangle]=3
        triangles[free_triangle+1]=discs[lastpoint]+lastrightseg
        triangles[free_triangle+2]=discs[point]+leftseg
        triangles[free_triangle+3]=discs[point]+rightseg
        free_triangle=free_triangle+4
      endfor
    endif

    if looped eq 1 then begin
      if point le 0 then break
    endif

    lastpoint=point
 
    point=(*self.linkage)[point,0] 
    count=count+1
    if (point lt 0) then finished=1
    looped=1
  endwhile

  triangles=triangles[0:free_triangle-1] 
  vertices=vertices[0:2,0:free_vertex-1] 
  normals=normals[0:2,0:free_vertex-1] 
 
  return,CREATE_STRUCT(['vertices','triangles','normals'],PTR_NEW(vertices),PTR_NEW(triangles),PTR_NEW(normals))
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
      ds = (*self.points)[point,0:2]-(*self.points)[lastpoint,0:2]
      modds = sqrt(tony_polyline_dot2(ds))

      if modds gt maxds then begin
;        if seglength gt minsegpnt then begin
          segment=segment+1
          seglength=0L
          segstart=point
;        endif else begin 
;         presegstart=(*self.linkage)[segstart,1]
;         ds = (*self.points)[presegstart,0:2]-(*self.points)[segstart,0:2]
;         flip[0:2]=0
;         if abs(ds[0]) gt maxds then begin
;           if ds[0] gt 0. then flip[0]=-1 else flip[0]=1
;         endif 
;         if abs(ds[1]) gt maxds then begin
;           if ds[1] gt 0. then flip[1]=-1 else flip[1]=1
;         endif 
;         if abs(ds[2]) gt maxds then begin
;           if ds[2] gt 0. then flip[2]=-1 else flip[2]=1
;         endif 
;         found=segstart
;         while found ne point do begin
;              delta=flip*1.*Lxyz
;              (*self.points)[found,0:2]=((*self.points)[found,0:2])-delta
;          found=(*self.linkage)[found,0] 
;          endwhile
;          ;segment=segment-1
;
;          ds = (*self.points)[point,0:2]-(*self.points)[lastpoint,0:2]
;          modds = sqrt(tony_polyline_dot2(ds))
;          if modds gt maxds then begin
;            segment=segment+1
;            seglength=0L
;            segstart=point
;          endif
;        endelse
      endif else begin
        if (*self.segments)[point] ge 0 then begin
          oldsegment=(*self.segments)[point]
          found=where((*self.segments) eq segment,found_count)
          if (found_count gt 0) then begin
            minfound=min(found)
            (*self.segments)[found]=oldsegment
            if (segment eq 1 ) then (*self.segments)[segstart]=segment
            segment=segment-1
          endif 
          break
        endif 
      endelse
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
FUNCTION tony_polyline::get_linkage
  COMPILE_OPT HIDDEN

  if PTR_VALID(self.linkage) then return,(*self.linkage)
  return,-1
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

