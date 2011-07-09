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

function dery, y, var
  ; First y-derivative.

  ny = n_elements(y)
  ny1 = ny-1

  dy = (y[ny1] - y[0])/double(ny1)

  fy = var

  fy[*,0,*] = -var[*,2,*] + 4D*var[*,1,*] - 3D*var[*,0,*]
  fy[*,1:ny-2,*] = shift(var[*,1:ny-2,*], -1) - shift(var[*,1:ny-2,*], 1)
  fy[*,ny1,*] = var[*,ny-3,*] - 4D*var[*,ny-2,*] + 3D*var[*,ny1,*]

  return, fy / (2D * dy)
end
