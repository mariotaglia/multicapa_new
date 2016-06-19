double precision function factorcurv (bas, pos)

use multicapa
implicit none
integer bas, pos


select case (abs(curvature))
case (0)
factorcurv = 1.0
case (1)
factorcurv = (float(bas)-0.5)/(float(pos)-0.5)
case (2)
factorcurv = (float(bas)-0.5)/(float(pos)-0.5)**2
end select

return
end

