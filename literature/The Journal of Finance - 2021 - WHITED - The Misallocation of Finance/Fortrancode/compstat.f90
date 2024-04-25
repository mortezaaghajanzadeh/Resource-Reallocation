program compstat
     use datatype
     use sizes
     use globals
     use pickmoments
     use omp_lib
     implicit none

    !This is an error code that gets passed along
    integer :: erc=0

    ! Now we specify the parameters that will attempt to estimate later on.
    real(dp) :: params(nop)
    real(dp), allocatable :: vguess(:,:,:,:)
    real(dp), allocatable :: momvector(:)


    ! These variables are for the comparative statics loop.
    integer :: param, vary, startploop, endploop, startvloop, endvloop, ii, jj, kk, nvary
    real(dp) :: glop,diter1,diter2
    real(dp) :: maxmin(nop,2),allmaxmin(nop,2)
    real(dp), allocatable :: allmoments(:,:,:),allparams(:,:)
    character (len = 11), allocatable :: momname(:)
    character (len = 35) momdfil

    character (len = 50) readfil

    integer :: date_time(8),time1(8)
    character (len = 12) real_clock(3)

    !=============Set the clock running======================================
    call date_and_time (real_clock(1), real_clock(2), real_clock(3), date_time)
    call date_and_time(values=time1)
    diter1=time1(5)*3600+time1(6)*60 +time1(7)+time1(8)*0.001_dp

    !=========Open all of the files=========================================

    !=============Set the bounds for the comparative statics exercises.=========
    allmaxmin(1,:)  = (/ 0.7_dp, 0.95_dp   /) !   returns to scale
    allmaxmin(2,:)  = (/ 0.3_dp, 0.7_dp    /) !   intangible share
    allmaxmin(3,:)  = (/ 0.6_dp, 0.9_dp    /) !   serial correlation
    allmaxmin(4,:)  = (/ 0.1_dp, 0.5_dp   /) !   standard deviation
    allmaxmin(5,:)  = (/ 0.05_dp, 0.2_dp   /) !   k depreciation
    allmaxmin(6,:)  = (/ 0.05_dp, 0.2_dp   /) !   n depreciation
    allmaxmin(7,:)  = (/ 0.001_dp, 0.5_dp  /) !   k adjustment costs
    allmaxmin(8,:)  = (/ 0.001_dp, 0.5_dp  /) !   n adjustment costs
    allmaxmin(9,:)  = (/ 0.001_dp, 1.0_dp  /) !   collateral
    allmaxmin(10,:) = (/ 0.00_dp, 0.15_dp /) !   equity issuance cost
    allmaxmin(11,:) = (/ 0.0_dp, 0.1_dp    /) !   debt issuance cost
    allmaxmin(12,:) = (/ 0.0_dp, 1.0_dp    /) !   dividend smoothing parameter  (0,0.01)    / or financial distress parameter
    allmaxmin(13,:) = (/ 0.4_dp, 0.9_dp    /) !   wage
    maxmin = allmaxmin
    nvary = 20
    !=======================Do the comparative statics loops====================
    allocate(allmoments(nmom,nvary,nop))
    allocate(allparams(nvary,nop))
    allocate(vguess(nz,nk,nn,nb))
    allocate(momname(nmom))

    allmoments = 0.0_dp
    momname = amomname(pickout)
    startploop = 10
    endploop   = 10!nop
    startvloop = 1
    endvloop   = nvary

    vguess = 0.0_dp

    allocate(momvector(nmom))
    allparams = 0.0_dp

    readfil = "Output/estfil.txt"
    open(1215,file=trim(readfil))
    do ii = 1,nop
      read(1215,*) params(ii)
    enddo
    close(1215)
    if (endploop == nop .and. endvloop == nvary) then
        momdfil = "Outputcomp/compstat.txt"
    else
        momdfil = "Outputcomp/teststat.txt"
    endif

    open(lout,file=trim(momdfil))

    write (lout,*) "---------------------------------------------------------------------------"
    write (lout,*) "                Basic Parameter Settings                                   "
    write (lout,*) "---------------------------------------------------------------------------"
    readfil = "Output/estfil.txt"
    open (unit = 66, file = trim(readfil))
    do ii = 1, nop
       read (66, '(F25.16)') params(ii)
    end do
    close (unit = 66)

    do ii=1,nop
        write (lout,"(a11,1x,f10.6)") pname(ii),params(ii)
    enddo
    paramloop: do param = startploop,endploop
        varyloop: do vary = startvloop,endvloop

              write(*,"('I am on trial ',2i3)") param,vary

              open(1215,file=trim(readfil))
              do ii = 1,nop
                read(1215,*) params(ii)
              enddo
              close(1215)

              !==================Change one parameter at a time=========================

              glop = (maxmin(param,2)-maxmin(param,1))*real(vary)/(real(nvary)-1.0_dp) &
                   - (maxmin(param,2)-maxmin(param,1))/(real(nvary)-1.0_dp) + maxmin(param,1)

              if (endploop > 1 .and. endvloop > 1)  params(param) = glop


              !=============================================================
              !                  SOLVE MODEL
              !=============================================================
              momvector = 0.0_dp
              call momentgen(params,momvector,vguess)



              !=============================================================
              !                 Update huge moment matrix
              !=============================================================

              do ii=1,nmom
                   allmoments(ii,vary,param) = momvector(ii)
              enddo

              !=============================================================
              !                 Update parameter values

              !=============================================================

              allparams(vary,param) = glop

        enddo varyloop
    enddo paramloop

    deallocate(momvector)
    !===========================================================================


    write (lout,*) "==========Parameter Values=============="
    do jj=1,nvary
      write(lout,"(30(' & ',f15.6))") allparams(jj,:)
    enddo



    write (lout,*) "==========Moments======================="


    do ii = 1,nmom
         write(lout,"(1x,a11)") momname(ii)
         write(lout,*) "     "
         do jj=1,nvary
           write(lout,"(30(' & ',f15.6))") allmoments(ii,jj,:)
         enddo
         write(lout,*) "     "
         write(lout,*) "     "
         write(lout,*) "     "
         write(lout,*) "     "
    enddo

    write (*,"('Start date: ',i4,'-',i2,'-',i2)") date_time(1:3)
    write (*,"('Start time: ',i2,':',i2,':',i2)") date_time(5:7)
    write (lout,"('Start date: ',i4,'-',i2,'-',i2)") date_time(1:3)
    write (lout,"('Start time: ',i2,':',i2,':',i2)") date_time(5:7)
    call date_and_time (real_clock(1), real_clock(2), real_clock(3), date_time)

    call date_and_time(values=time1)
    diter2=time1(5)*3600+time1(6)*60 +time1(7)+time1(8)*0.001_dp
    write (*,"('End date: ',i4,'-',i2,'-',i2)") date_time(1:3)
    write (*,"('End time: ',i2,':',i2,':',i2)") date_time(5:7)
    write (lout,"('End date: ',i4,'-',i2,'-',i2)") date_time(1:3)
    write (lout,"('End time: ',i2,':',i2,':',i2)") date_time(5:7)
    write (*,"('Elapse of time = ',f12.5)") diter2-diter1

    close(lout)

    deallocate(allmoments)

end program compstat
