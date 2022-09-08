c Sequential Stern-Gerlach Experiment Simulation
      program sgseq
        integer, allocatable :: component(:, :)
        integer init_num, up_spin, down_spin, n, i, j
        real rand
        character(len=3)::opc(2)

        print*,'Enter the number of Stern-Gerlach devices:'
        read*,n
        allocate(component(n, 3))
        do i=1,n
          print *, 'Index of non-zero B-field component of device -',i
          read(*,*)rand
          if(rand.ne.1.and.rand.ne.2.and.rand.ne.3) then 
            print*,'Wrong index input !!!!'
            print*,'Indices can either be 1, 2 or 3.'
            stop
          endif
          component(i, 1) = rand
        enddo
        print *, ''
        print *, 'Number of efflux atoms from the oven = ?'
        read(*,*)init_num

        print *, 'There are total',n,'		Stern-Gerlach devices.'
        print *, 'Number of efflux atoms from the oven =',init_num
        print *, ''
        print *, 'Stern-Gerlach_device_no.   Spin_component', 
     1	  '   ',
     2	  'No. of outcoming beams', 
     3	  '   ',
     4	  '(+)',
     3	  '   ',
     4	  '(-)'

        do i=1,n
          if(component(i,1)==1) then
            opc(1) = 'S_x'
          else if(component(i,1)==2) then
            opc(1) = 'S_y'
          else
            opc(1) = 'S_z'
          endif

          if(i==1) then  ! For the first device
            opc(2) = '2'
c           print*,i,init_num
            up_spin = 0
            do j=1,init_num
              call random_number(rand)
              ! random variable < 0.5 means spin-up electrons
              if (rand.lt.0.5) then
                up_spin = up_spin + 1
              endif
            enddo
            down_spin=init_num-up_spin
            component(i, 2) = up_spin
            component(i, 3) = down_spin

          else ! For the second and the following devices
            init_num = up_spin
c           print*,i,init_num
            if(component(i,1).ne.component(i-1,1)) then
              opc(2) = '2'
              up_spin = 0
              do j=1,init_num
                call random_number(rand)
                if(rand.lt.0.5) then  ! random variable < 0.5 means spin-up electrons
                  up_spin = up_spin + 1
                endif
              enddo
            else
              opc(2) = '1'
            endif
            down_spin = init_num-up_spin
            component(i, 2) = up_spin
            component(i, 3) = down_spin
          endif

          print *, i, opc(1), " ", opc(2), component(i,2),component(i,3)
        enddo
      end
