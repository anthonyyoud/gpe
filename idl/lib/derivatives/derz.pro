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

function derz, z, var
  ; First z-derivative.

  nz = n_elements(z)
  nz1 = nz-1

  dz = (z[nz1] - z[0])/double(nz1)

  fz = var

  fz[*,*,0] = -var[*,*,2] + 4D*var[*,*,1] - 3D*var[*,*,0]
  fz[*,*,1:nz-2] = shift(var[*,*,1:nz-2], -1) - shift(var[*,*,1:nz-2], 1)
  fz[*,*,nz1] = var[*,*,nz-3] - 4D*var[*,*,nz-2] + 3D*var[*,*,nz1]

  return, fz / (2D * dz)
end
