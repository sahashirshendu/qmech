      program sgseq
      integer, allocatable :: comp(:, :)
      integer atn, up, down, n, i, j
      real rand

      print *, "SGX - Press 1"
      print *, "SGY - Press 2"
      print *, "SGZ - Press 3"
      print *,'Number of SG devices ='
      read *, n
      allocate(comp(n, 3))
      do i = 1, n
        print *, 'Choice of device =',i
        read *, rand
        comp(i, 1) = rand
      enddo
      print *, ''
      print *, 'Number of atoms :'
      read *, atn

      do i=1,n
        if(i .eq. 1) then
          up = 0
          do j=1,atn
            call random_number(rand)
            if (rand .gt. 0.5) then
              up = up + 1
            end if
          end do
          down = atn - up
          comp(i, 2) = up
          comp(i, 3) = down
        else
          atn = up
          if(comp(i,1) .ne. comp(i-1,1)) then
            up = 0
            do j = 1, atn
              call random_number(rand)
              if(rand .gt. 0.5) then
                up = up + 1
              end if
            end do
          else
          endif
          down = atn - up
          comp(i, 2) = up
          comp(i, 3) = down
        endif

        print *, i, (comp(i,j),j=1,3)
      enddo
      end
