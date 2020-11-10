!  
!  Written by Leandro Martínez, 2009-2011.
!  Copyright (c) 2009-2011, Leandro Martínez, Jose Mario Martinez,
!  Ernesto G. Birgin.
!  
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
! Function that computes the atom-to-atom component of the objective
! function
!

double precision function fparc(icart,firstjcart)

  use sizes
  use compute_data
  implicit none

  ! SCALAR ARGUMENTS
  integer :: icart,firstjcart

  ! LOCAL SCALARS
  integer :: jcart
  double precision :: a1,a2,a3,datom,tol

  fparc = 0.0d0
  jcart = firstjcart
  do while ( jcart .ne. 0 )
    if(comptype(ibtype(jcart))) then
      if(ibmol(icart).ne.ibmol(jcart).or. &
         ibtype(icart).ne.ibtype(jcart)) then
        a1 = xcart(icart, 1)-xcart(jcart, 1) 
        a2 = xcart(icart, 2)-xcart(jcart, 2) 
        a3 = xcart(icart, 3)-xcart(jcart, 3) 
        datom = a1 * a1 + a2 * a2 + a3 * a3
        tol = (radius(icart)+radius(jcart))**2
        a1 = dmin1(datom - tol, 0.d0)
        fparc = fparc + a1 * a1
        tol = (radius_ini(icart)+radius_ini(jcart))**2
        fdist = dmax1(tol-datom,fdist)
      end if
    end if
    jcart = latomnext(jcart)
  end do

end function fparc

