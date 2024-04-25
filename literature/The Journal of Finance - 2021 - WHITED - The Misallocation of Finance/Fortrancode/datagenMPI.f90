program datagen
     use datatype
     use sizes
     use globals
     use pickmoments
     use omp_lib
     implicit none

     include 'mpif.h'

     integer :: inode, num_nodes, ierror, tag, status(MPI_STATUS_SIZE)      !MPI stuff
     integer :: master_id = 0
     !This is an error code that gets passed along
     integer :: erc=0

     ! Now we specify the parameters that will attempt to estimate later on.
     real(dp) :: params(nop)
     real(dp), allocatable :: vguess(:,:,:,:)
     real(dp), allocatable :: momvector(:), draws(:), drawmx(:,:)


     ! These variables are for the comparative statics loop.
     integer :: param, ii, jj, iiseed(2), pnum, nout, ndraws = 550
     real(dp) :: glop,diter1,diter2,pposition,ppick,gamma
     real(dp) :: maxmin(nop,2),allmaxmin(nop,2)
     character (len = 35) momdfil

     character (len = 50) readfil

     character (len = 4) string
     integer :: programme

     integer :: date_time(8),time1(8)
     character (len = 12) real_clock(3)


     intrinsic random_seed, random_number
     call MPI_INIT(ierror)
     call MPI_COMM_RANK(MPI_COMM_WORLD, inode, ierror)
     call MPI_COMM_SIZE(MPI_COMM_WORLD, num_nodes, ierror)


    write( string, '(i4)' ) inode
    momdfil = "Outputcomp/datagen"//trim(string)//".txt"
    open(lout,file=trim(momdfil))

    !=============Set the clock running======================================
    call date_and_time (real_clock(1), real_clock(2), real_clock(3), date_time)
    call date_and_time(values=time1)
    diter1=time1(5)*3600+time1(6)*60 +time1(7)+time1(8)*0.001_dp

    !=========Open all of the files=========================================

    !=============Set the bounds for the comparative statics exercises.=========
    allmaxmin = 0.0_dp
    allmaxmin(1,:)  = (/ 0.7_dp, 0.9_dp    /) !   returns to scale
    allmaxmin(2,:)  = (/ 0.3_dp, 0.7_dp    /) !   labor share
    allmaxmin(3,:)  = (/ 0.4_dp, 0.9_dp    /) !   serial correlation
    allmaxmin(4,:)  = (/ 0.05_dp, 0.4_dp   /) !   standard deviation
    allmaxmin(5,:)  = (/ 0.05_dp, 0.2_dp   /) !   k depreciation
    allmaxmin(6,:)  = (/ 0.05_dp, 0.2_dp   /) !   n depreciation
    allmaxmin(7,:)  = (/ 0.001_dp, 0.5_dp  /) !   k adjustment costs
    allmaxmin(8,:)  = (/ 0.001_dp, 0.5_dp  /) !   n adjustment costs
    allmaxmin(9,:)  = (/ 0.001_dp, 1.0_dp  /) !   collateral
    allmaxmin(10,:) = (/ 0.0001_dp, 0.15_dp /) !   equity issuance cost
    allmaxmin(11,:) = (/ 0.0001_dp, 0.15_dp /) !   debt issuance cost
    allmaxmin(12,:) = (/ 0.0001_dp, 0.15_dp /) !   dividend smoothing
    allmaxmin(13,:) = (/ 0.8_dp, 0.8_dp    /) !   wage

    maxmin = allmaxmin

    !=======================Do the comparative statics loops====================
    allocate(vguess(nz,nk,nn,nb))
    allocate(momvector(nmom))
    allocate(draws(ndraws*nop))
    vguess = 0.0_dp


    readfil = "Output/estfil.txt"
    open(1215,file=trim(readfil))
    do ii = 1,nop
      read(1215,*) params(ii)
    enddo
    close(1215)

    !write (lout,*) "---------------------------------------------------------------------------"
    !write (lout,*) "                Basic Parameter Settings                                   "
    !write (lout,*) "---------------------------------------------------------------------------"
    !
    !do ii=1,nop
    !    write (lout,"(a11,1x,f10.6)") pname(ii),params(ii)
    !enddo

    iiseed  = (/ 234567890, 987654320 /)

    call random_seed(PUT=iiseed)

    do jj = 1,ndraws*nop

       call random_number(harvest=ppick)
       draws(jj) = ppick

    enddo
    allocate(drawmx(ndraws,nop))
    drawmx = reshape(draws,(/ndraws, nop/))

    do jj = 1,ndraws

       if (mod(real(jj-1),real(num_nodes))==real(inode)) then
             write(*,"('I am on trial ',i3)") jj

              open(1215,file=trim(readfil))
              do ii = 1,nop
                read(1215,*) params(ii)
              enddo
              close(1215)

              do ii = 1,nop-3
                 pposition = drawmx(jj,ii)
                 params(ii) = (1.0_dp - pposition)*maxmin(ii,1) + pposition*maxmin(ii,2)
              enddo
              !write(*,"(<nop>f7.4)") params
              !=============================================================
              !                  SOLVE MODEL
              !=============================================================
              momvector = 0.0_dp
              call momentgen(params,momvector,vguess)

              gamma = momvector(16)
              nout = nop + 1
              if (gamma > 0.0_dp .and. gamma < 5.0_dp) then  ! Get rid of outliers

                write(lout,"(<nout>f26.16)") gamma,params
                write(*,"('foo  ',<nout>f26.16)") gamma,params

              endif
       endif
    enddo

    deallocate(momvector)
    !===========================================================================


    write (*,"('Start date: ',i4,'-',i2,'-',i2)") date_time(1:3)
    write (*,"('Start time: ',i2,':',i2,':',i2)") date_time(5:7)
    !write (lout,"('Start date: ',i4,'-',i2,'-',i2)") date_time(1:3)
    !write (lout,"('Start time: ',i2,':',i2,':',i2)") date_time(5:7)
    call date_and_time (real_clock(1), real_clock(2), real_clock(3), date_time)

    call date_and_time(values=time1)
    diter2=time1(5)*3600+time1(6)*60 +time1(7)+time1(8)*0.001_dp
    write (*,"('End date: ',i4,'-',i2,'-',i2)") date_time(1:3)
    write (*,"('End time: ',i2,':',i2,':',i2)") date_time(5:7)
    !write (lout,"('End date: ',i4,'-',i2,'-',i2)") date_time(1:3)
    !write (lout,"('End time: ',i2,':',i2,':',i2)") date_time(5:7)
    write (*,"('Elapse of time = ',f12.5)") diter2-diter1

    close(lout)

    call MPI_FINALIZE(ierror)

end program datagen
