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

! Routine to print the title 

subroutine title()

  use ahestetic
  write(*,hash3_line)
  write(*,"(' PACKMOL - Packing optimization for the automated generation of', /&
           &' starting configurations for molecular dynamics simulations.', /&
           &' ',/&
           &t62,' Version 16.228 ')")
  write(*,hash3_line)

end subroutine title
