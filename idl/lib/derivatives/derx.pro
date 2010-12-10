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
