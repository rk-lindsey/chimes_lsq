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
!  Module that contains some ahestetic output definitions
!
module ahestetic

  character(len=13), parameter :: dash1_line = "(  80('-')  )",&
                                  dash2_line = "(/,80('-')  )",&
                                  dash3_line = "(/,80('-'),/)"

  character(len=13), parameter :: hash1_line = "(  80('#')  )",&
                                  hash2_line = "(/,80('#')  )",&
                                  hash3_line = "(/,80('#'),/)"

  character(len=31), parameter :: prog1_line = "('  Packing:|0 ',tr60,'100%|' )",&
                                  prog2_line = "('   Moving:|0 ',tr60,'100%|' )"

end module ahestetic
