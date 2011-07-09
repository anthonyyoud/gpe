;Copyright 2011 Anthony Youd/Newcastle University
;
;   Licensed under the Apache License, Version 2.0 (the "License");
;   you may not use this file except in compliance with the License.
;   You may obtain a copy of the License at
;
;       http://www.apache.org/licenses/LICENSE-2.0
;
;   Unless required by applicable law or agreed to in writing, software
;   distributed under the License is distributed on an "AS IS" BASIS,
;   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
;   See the License for the specific language governing permissions and
;   limitations under the License.

pro save_vapor_data, ind, nstart, nend, skip, data, odir, phase=phase, $
  vx=vx, vy=vy, vz=vz, append=append, num_levels=num_levels

  ; Set the number of VAPOR data refinement levels.
  default, num_levels, 2

  ; Set the names of the variables to be saved.
  varnames = ['density']
  if keyword_set(phase) then begin
    varnames = [[varnames], 'phase']
  endif
  if keyword_set(vx) then begin
    varnames = [[varnames], 'vx']
  endif
  if keyword_set(vy) then begin
    varnames = [[varnames], 'vy']
  endif
  if keyword_set(vz) then begin
    varnames = [[varnames], 'vz']
  endif
  vdffile = odir + 'gpe.vdf'

  if (ind eq 0) then begin
    ; If this is the first data file to be saved...
    if not keyword_set(append) then begin
      ; If we're not explicitly appending (/append) to an existing VDF file...
      mfd = vdf_create([data.nx,data.ny,data.nz], num_levels)
      timesteps = (nend-nstart)/skip + 1
      vdf_setnumtimesteps, mfd, timesteps
      vdf_setvarnames, mfd, varnames
      extents = [data.x[0], data.y[0], data.z[0], data.x[data.nx-1], $
        data.y[data.ny-1], data.z[data.nz-1]]
      vdf_setextents, mfd, extents
    endif else begin
      ; If we are explicitly appending to an existing VDF file...
      mfd = vdf_create(vdffile)
      ind = ind + vdf_getnumtimesteps(mfd)
      timesteps = vdf_getnumtimesteps(mfd) + (nend-nstart)/skip + 1
      vdf_setnumtimesteps, mfd, timesteps
    endelse
  endif else begin
    ; If this is not the first data file to be saved...
    mfd = vdf_create(vdffile)
  endelse

  ; Write the data.
  vdf_settusertime, mfd, ind, [data.t]
  vdf_write, mfd, vdffile
  dfd = vdc_bufwritecreate(vdffile)
  vdc_openvarwrite, dfd, ind, 'density'
  for k=0,data.nz-1 do begin
    vdc_bufwriteslice, dfd, float(data.dens[*,*,k])
  endfor
  vdc_closevar, dfd
  if keyword_set(phase) then begin
    vdc_openvarwrite, dfd, ind, 'phase'
    for k=0,data.nz-1 do begin
      vdc_bufwriteslice, dfd, float(data.phase[*,*,k])
    endfor
    vdc_closevar, dfd
  endif
  if keyword_set(vx) then begin
    vdc_openvarwrite, dfd, ind, 'vx'
    for k=0,data.nz-1 do begin
      vdc_bufwriteslice, dfd, float(data.vx[*,*,k])
    endfor
    vdc_closevar, dfd
  endif
  if keyword_set(vy) then begin
    vdc_openvarwrite, dfd, ind, 'vy'
    for k=0,data.nz-1 do begin
      vdc_bufwriteslice, dfd, float(data.vy[*,*,k])
    endfor
    vdc_closevar, dfd
  endif
  if keyword_set(vz) then begin
    vdc_openvarwrite, dfd, ind, 'vz'
    for k=0,data.nz-1 do begin
      vdc_bufwriteslice, dfd, float(data.vz[*,*,k])
    endfor
    vdc_closevar, dfd
  endif
  vdc_bufwritedestroy, dfd
  vdf_destroy, mfd
end
