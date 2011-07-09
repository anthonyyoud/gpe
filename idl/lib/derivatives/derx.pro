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

function derx, x, var
  ; First x-derivative.

  nx = n_elements(x)
  nx1 = nx-1

  dx = (x[nx1] - x[0])/double(nx1)

  fx = var

  fx[0,*,*] = -var[2,*,*] + 4D*var[1,*,*] - 3D*var[0,*,*]
  fx[1:nx-2,*,*] = shift(var[1:nx-2,*,*], -1) - shift(var[1:nx-2,*,*], 1)
  fx[nx1,*,*] = var[nx-3,*,*] - 4D*var[nx-2,*,*] + 3D*var[nx1,*,*]

  return, fx / (2D * dx)
end
