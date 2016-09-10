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
!
! Compute gradient relative to atom-to-atom distances
!

subroutine gparc(icart,firstjcart)

  use sizes
  use compute_data
  implicit none

  ! SCALAR ARGUMENTS
  integer :: icart,firstjcart

  ! LOCAL SCALARS
  integer :: jcart
  double precision :: a1,a2,a3,datom,dtemp,xdiff,tol

  jcart = firstjcart
  do while ( jcart .ne. 0 )
    if(comptype(ibtype(jcart))) then
      if(ibmol(icart).ne.ibmol(jcart).or. &
        ibtype(icart).ne.ibtype(jcart)) then
        tol = (radius(icart)+radius(jcart))**2
        a1 = xcart(icart, 1)-xcart(jcart, 1) 
        a1 = a1 * a1
        if(a1.lt.tol) then
          a2 = xcart(icart, 2)-xcart(jcart, 2) 
          a2 = a1 + a2 * a2
          if(a2.lt.tol) then
            a3 = xcart(icart, 3)-xcart(jcart, 3)
            datom = a2 + a3 * a3 
            if(datom.lt.tol) then 
              dtemp = 4.d0 * (datom - tol)
              xdiff = dtemp*(xcart(icart,1) - xcart(jcart,1)) 
              gxcar(icart,1)= gxcar(icart,1) + xdiff
              gxcar(jcart,1)= gxcar(jcart,1) - xdiff 
              xdiff = dtemp*(xcart(icart,2) - xcart(jcart,2)) 
              gxcar(icart,2)= gxcar(icart,2) + xdiff
              gxcar(jcart,2)= gxcar(jcart,2) - xdiff 
              xdiff = dtemp*(xcart(icart,3) - xcart(jcart,3)) 
              gxcar(icart,3)= gxcar(icart,3) + xdiff
              gxcar(jcart,3)= gxcar(jcart,3) - xdiff 
            end if
          end if
        end if 
      end if
    end if
    jcart = latomnext(jcart)
  end do

  return
end subroutine gparc

