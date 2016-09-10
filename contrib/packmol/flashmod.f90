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
! Arrays required by the flashsort package. Used only in heuristics, but
! defined here to be allocated dynamically
!

module flashsort

  use sizes
  implicit none
  integer, allocatable :: indflash(:) ! (ntotat)
  integer, allocatable :: lflash(:) ! (ntotat)
  integer :: mflash

end module flashsort

