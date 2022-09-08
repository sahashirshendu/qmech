c """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
c Description: Sequential Stern Gerlach Experiment
c Author: Shirshendu Saha
c CRN: 10
c Registration No.: 20209110010
c """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      program sg_seq
        implicit none
        integer i, j, ch, chpr, n, ne
        integer xpcnt, zpcnt, xmcnt, zmcnt
        real rand
        print '(A)', "Choices :"
10      format (2X, A)
        print 10, "SGZ - Magnetic Field oriented along Z axis : Enter 1"
        print 10, "SGX - Magnetic Field oriented along X axis : Enter 2"
        print '(/1X, A)', "The number of Stern-Gerlach devices:"
        read *, n
        print '(/1X, A)', "The number of electrons:"
        read *, ne

        xpcnt = 0
        xmcnt = 0
        zpcnt = 0
        zmcnt = 0

        chpr = 0
        do i = 1, n
          print '(/, A, 1X, I3, 1X, A)', "SG Device -", i, ":"
          print '(A)', "Enter SG Device Type :"
          read *, ch
          if (ch .ne. chpr) then
            if (ch .eq. 1) then
              xpcnt = 0
              xmcnt = 0
              do j = 1, ne
                call random_number(rand)
                if (rand .gt. 0.5) then
                  zpcnt = zpcnt + 1
                else
                  zmcnt = zmcnt + 1
                end if
              end do
              ne = zpcnt
            else if (ch .eq. 2) then
              zpcnt = 0
              zmcnt = 0
              do j = 1, ne
                call random_number(rand)
                if (rand .gt. 0.5) then
                  xpcnt = xpcnt + 1
                else
                  xmcnt = xmcnt + 1
                end if
              end do
              ne = xpcnt
            else
              print '(A)', "Enter a valid choice!"
            end if
          else if (ch .eq. chpr) then
            if (ch .eq. 1) then
              ne = zpcnt
              zmcnt = 0
            else if (ch .eq. 2) then
              ne = xpcnt
              xmcnt = 0
            else
              print '(A)', "Enter a valid choice!"
            end if
          end if
          print '(/, A, 1X, I10)', "SGZ(+) :", zpcnt
          print '(A, 1X, I10)', "SGZ(-) :", zmcnt
          print '(A, 1X, I10)', "SGX(+) :", xpcnt
          print '(A, 1X, I10)', "SGX(-) :", xmcnt
          chpr = ch
        end do
      end
