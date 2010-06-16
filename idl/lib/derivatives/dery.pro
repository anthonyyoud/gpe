; $Id$
;----------------------------------------------------------------------------

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
