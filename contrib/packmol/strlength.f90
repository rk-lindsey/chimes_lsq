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
! Function that determines the length of a string (better than 
! intrinsic "len" because considers tabs as empty characters)
!
function strlength(string)

  implicit none
  integer :: strlength
  character(len=200) :: string
  logical empty_char
  
  strlength = 200
  do while(empty_char(string(strlength:strlength)))
    strlength = strlength - 1
    if ( strlength == 0 ) exit
  end do

end function strlength      

!
! Function that determines if a character is empty (empty, space, or tab)
! (nice suggestion from Ian Harvey -IanH0073- at github)
!

function empty_char(ch)
  character :: ch
  logical empty_char
  empty_char = .false.
  if ( ch == '' .or. &
       ch == achar(9) .or. &
       ch == achar(32) ) then
    empty_char = .true.
  end if
end function empty_char
 
