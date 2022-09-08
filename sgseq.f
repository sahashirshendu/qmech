      program sg_seq
        integer i, j, ch, chpr, n, ne
        integer xpl, zpl, xmi, zmi
        print *, "Choices: SGZ -> Press 1, SGZ -> Press 2"
        print *, ""
        print *, "The number of Stern-Gerlach devices:"
        read *, n
        print *, "The number of electrons:"
        read *, ne

        xpl = 0
        xmi = 0
        zpl = 0
        zmi = 0

        chpr = 0
        do i = 1, n
          print *, ""
          print *, "Device", i
          print *, "Device Choice ->"
          read *, ch
          if (ch .ne. chpr) then
            if (ch .eq. 1) then
              xpl = 0
              xmi = 0
              do j = 1, ne
                call random_number(r)
                if (r .gt. 0.5) then
                  zpl = zpl + 1
                else
                  zmi = zmi + 1
                end if
              end do
              ne = zpl
            else if (ch .eq. 2) then
              zpl = 0
              zmi = 0
              do j = 1, ne
                call random_number(r)
                if (r .gt. 0.5) then
                  xpl = xpl + 1
                else
                  xmi = xmi + 1
                end if
              end do
              ne = xpl
            end if
          else if (ch .eq. chpr) then
            if (ch .eq. 1) then
              ne = zpl
              zmi = 0
            else if (ch .eq. 2) then
              ne = xpl
              xmi = 0
            end if
          end if
          print *, "Z+ =", zpl
          print *, "Z- =", zmi
          print *, "X+ =", xpl
          print *, "X- =", xmi
          chpr = ch
        end do
      end
