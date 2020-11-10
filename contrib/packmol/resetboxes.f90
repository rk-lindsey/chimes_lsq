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
! Subroutine resetboxes
!

subroutine resetboxes()
      
  use sizes
  use compute_data, only : nboxes, latomfirst, latomfix
  implicit none

  integer :: i,j,k

  ! Reset boxes

  do i = 1, nboxes(1)
    do j = 1, nboxes(2)
      do k = 1, nboxes(3)
        latomfirst(i,j,k) = latomfix(i,j,k)
      end do
    end do
  end do

  ! Reset margins
      
  do j = 0, nboxes(2)+1
    do k = 0, nboxes(3)+1
      latomfirst(0,j,k) = 0
      latomfirst(nboxes(1)+1,j,k) = 0
    end do
  end do

  do i = 0, nboxes(1)+1
    do k = 0, nboxes(3)+1
      latomfirst(i,0,k) = 0
      latomfirst(i,nboxes(2)+1,k) = 0
    end do
  end do

  do i = 0, nboxes(1)+1
    do j = 0, nboxes(2)+1
      latomfirst(i,j,0) = 0
      latomfirst(i,j,nboxes(3)+1) = 0
    end do
  end do      

  return
end subroutine resetboxes

