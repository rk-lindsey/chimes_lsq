      program radf
      implicit none
c
c     *radf* computes a radial distribution function also an
c      intermolecular rdf
c
      real*8 zero,one,two,three,four
      parameter (zero = 0.d0,
     +           one = 1.d0,
     +           two = 2.d0,
     +           three = 3.d0,
     +           four = 4.d0)
c
      character out1*1,out2*5,out3*2,out4*6,out5*3,out6*7,out7*4,out8*8
      integer   ir,nat,ncon,i,n,k,j,m,jj,ii,atm,typ,mol,
     +          inc1,dd,xmove,ymove,zmove,kk,inc
      real*8    x,y,z,xjunk,yjunk,zjunk,
     +          b2,b2a,b2b,b2c,b2_2,b2_2b,b2_2c,a,b,c,xlo,xhi,ylo,
     +          yhi,zlo,zhi, qq
c
      open(15,file='dump.coor',status='unknown')
c     open(15,file='moo',status='unknown')
c
      inc = 0
      do jj = 1, 20500
        do ii = 1, 3
          read(15,*,end=101)
        enddo
        inc = inc + 1
        if(inc.le.9) then
          write(out1,'(i1)') inc
          out2 = out1//'.pdb'
          write(6,*) out2
          open(93,file=out2)
        endif
        if(inc.ge.10.and.inc.le.99) then
          write(out3,'(i2)') inc
          out4 = out3//'.pdb'
          write(6,*) out4
          open(93,file=out4)
        endif
        if(inc.ge.100.and.inc.le.999) then
          write(out5,'(i3)') inc
          out6 = out5//'.pdb'
          write(6,*) out6
          open(93,file=out6)
        endif
        if(inc.ge.999.and.inc.le.9999) then
          write(out7,'(i4)') inc
          out8 = out7//'.pdb'
          write(6,*) out8
          open(93,file=out8)
        endif
        write(93,801) 'HEADER 500 MUST_PBC'
        read(15,*) nat
        ncon = nat-1
        write(6,*)'Number of atoms=',nat
        read(15,*)
        read(15,*)xlo,xhi
        a = xhi-xlo
        read(15,*)ylo,yhi
        b = yhi-ylo
        read(15,*)zlo,zhi
        c = zhi-zlo
        read(15,*)
        write(93,802) 'CRYST1',a,b,c,90.00,90.00,90.00,' P 1'
        write(6,*)'Box=',a,b,c
        b2a = a/two
        b2b = b/two
        b2c = c/two
        b2 = b2a
        if(b2b.lt.b2) then
          b2 = b2b
        elseif(b2c.lt.b2) then
          b2 = b2c
        endif
        b2_2 = b2a*b2a
        b2_2b = b2b*b2b
        b2_2c = b2c*b2c
        if(b2_2b.lt.b2_2) then
          b2_2 = b2_2b
        elseif(b2_2c.lt.b2_2) then
          b2_2 = b2_2c
        endif
        do j = 1, nat
c id mol type xs ys zs ix iy iz
          read(15,*) atm,typ,qq,xjunk,yjunk,zjunk
          xmove = 0
          ymove = 0
          zmove = 0
          k = atm
          x = a*(xjunk + dble(xmove))
          y = b*(yjunk + dble(ymove))
          z = c*(zjunk + dble(zmove))
c         x = a*xjunk
c         y = b*yjunk
c         z = c*zjunk
c1 12.0000
c2 1.0080
c3 15.9990
c4 14.0000
c         if(atm.le.13920) then
            if(typ.eq.1)write(93,803)'ATOM',k,'C',1,x,y,z
            if(typ.eq.2)write(93,803)'ATOM',k,'H',2,x,y,z
            if(typ.eq.3)write(93,803)'ATOM',k,'O',3,x,y,z
            if(typ.eq.4)write(93,803)'ATOM',k,'N',4,x,y,z
c         else
c           write(93,803)'ATOM',k,'P',1,x,y,z
c         endif
        enddo
        write(93,804) 'TER'
        write(93,804) 'END'
        close(93)
      enddo
c
 101  inc1 = inc
      write(6,*) inc1
c
 801  format(A19)
 802  format(A6,3f9.3,3f7.2,A4)
 803  format(A4,i7,2x,A1,5x,i7,5x,3f10.3)
 804  format(A3)
 902  format(8f10.5)
 903  format(1A)
 905  format(f12.5,4f11.4)
 906  format(f12.5,3f11.4)
 908  format(4f16.3)
      end
