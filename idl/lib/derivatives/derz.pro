; $Id$
;----------------------------------------------------------------------------

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
