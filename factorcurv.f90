double precision function factorcurv (bas, pos)

use multicapa
implicit none
integer bas, pos
factorcurv = (float(bas)-0.5)/(float(pos)-0.5)

return
end

