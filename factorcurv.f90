double precision function factorcurv (base, pos)

use multicapa
implicit none
integer base, pos
factorcurv = (float(base)-0.5)/(float(pos)-0.5)

return
end

