subroutine valfun(params,vg)
      use datatype
      use sizes
      use globals
      use omp_lib
      implicit none

      real(dp), intent(in) :: params(nop)
      real(dp), intent(inout) :: vg(nz,nk,nn,nb)


      real(dp), parameter :: toler = 0.000001_dp, disttol = 1.0e-12_dp                    !tolerances
      integer,  parameter :: capt = 2000, distmax = 1000                                  !maximum iterations

      real(dp) :: theta, beta, rho, sigma, delta, inlta, &                                ! parameters to estimate
                  adjK, adjN, xi, lambda, wage, dishcost, divsmooth

      real(dp) :: fixcost, aye, bee                                                       ! other parameters

      real(dp) :: invest, innova, ncf                                                     ! investment, intangible investment, net cash flow

      real(dp) :: beeta                                                                   ! discount rate parameters

      real(dp) :: const, alpha, skale                                                     ! annoying functions of parameters for making profits

      real(dp) :: diff = 100.0_dp                                                         ! misc variables and constants

      real(dp) :: zmax, zmin, kmax, kmin, nmax, nmin, bmax, bmin, emax, emin              ! constants for setting up the grids

      real(dp) :: kmaxstar, kminstar, kstar                                               ! grid anchors

      real(dp) :: adjfactor, mess, str, ve, disterr, start, xtemp                         ! temporary variables

      integer :: n, iz, ie, ik, in, ib, iik, iin, iib, ize, ii, jj, &
                 counter, error, id, n1, n2, &                                            ! I hate initializing integers
                 jk, jn, jb, u_k, u_n, u_b, l_k, l_n, l_b

      integer :: pdiff, inarow(5), suminarow                                              ! Stuff for policy function convergence

      integer :: outofbounds, iter, maxiter = 10                                          ! Integers for adjusting the bounds

      real(dp) :: stretch_khi, stretch_klo                                                ! Stuff for fixing the initial grid bounds

      real(dp), allocatable :: zg(:), eg(:), tmat(:,:)                                    ! stuff for the transition matrix

      real(dp), allocatable :: v(:,:,:,:), vnew(:,:,:,:), &                               ! old and new value function
                               queue(:,:,:,:),   &                                        ! expected value
                               queuelong(:,:,:,:), &                                      ! interpolated expected value
                               vdiff(:,:,:,:), &                                          ! value difference
                               profit(:,:,:,:,:,:,:), &                                   ! profits
                               dish(:,:,:,:), &                                           ! debt issuance
                               avsmall(:,:,:)
      real(dp), allocatable :: mew_new(:,:,:,:), mew(:,:,:,:)                             ! stationary distribution


      real(dp), allocatable :: kg(:), ng(:), bg(:), kpg(:), npg(:), bpg(:)                ! grids  "p" is for policy

      logical :: exitloop

      type(grids), dimension(gmax)  :: sgrid                                              ! state grids
      type(grids), dimension(egmax) :: pgrid                                              ! policy grids
      type(policies), dimension(2)  :: endpoints                                          ! grid endpoints
      type(policies), allocatable   :: pol(:,:,:,:)                                       ! policy functions
      type(weights), allocatable    :: interp_wgt(:,:,:,:,:,:,:)                          ! Interpolation weights in the right places
      type(ipolicies), allocatable  :: gidx(:,:,:,:), sidx(:,:,:,:), oidx(:,:,:,:)        ! policy indices: "g" is the main one, "s" is small, "o" is old
      type(ipolicies), allocatable  :: coarse_idx(:,:,:,:)                                ! policy indices on the coarse grid (one down)
      type(gridinterps) :: k_wgt(npk)                                                     ! Structures that hold
      type(gridinterps) :: n_wgt(npn)                                                     ! and corrsponding grid points
      type(gridinterps) :: b_wgt(npb)                                                     ! interpolation weights

      call OMP_SET_NUM_THREADS(71)

      start = omp_get_wtime()
!=========================Allocate all of the arrays=================================
      allocate(zg(nz))
      allocate(kg(nk))
      allocate(ng(nn))
      allocate(bg(nb))
      allocate(kpg(npk))
      allocate(npg(npn))
      allocate(bpg(npb))
      allocate(tmat(nz,nz))
      allocate(mew(nz,nk,nn,nb))
      allocate(mew_new(nz,nk,nn,nb))
      allocate(dish(nk,nb,npk,npb))
      allocate(profit(nz,nk,nn,nb,npk,npn,npb))
      allocate(avsmall(2*pad+1,2*pad+1,2*pad+1))
      allocate(v(nz,nk,nn,nb))
      allocate(vnew(nz,nk,nn,nb))
      allocate(vdiff(nz,nk,nn,nb))
      allocate(queue(nz,nk,nn,nb))
      allocate(queuelong(nz,npk,npn,npb))

      allocate(gidx(nz,nk,nn,nb))
      allocate(sidx(nz,nk,nn,nb))
      allocate(oidx(nz,nk,nn,nb))

      allocate(pol(nz,nk,nn,nb))
      !profit = 0.0_dp
      !allval = 0.0_dp
      !===========================================================================================================
      !===========================================================================================================
      !===========================================================================================================
      !         This section reads in the parameters and sets up the state spaces
      !===========================================================================================================
      !===========================================================================================================
      !===========================================================================================================

      theta      = params(1)
      beta       = params(2)
      rho        = params(3)
      sigma      = params(4)
      delta      = params(5)
      inlta      = params(6)
      adjK       = params(7)
      adjN       = params(8)
      xi         = params(9)
      lambda     = params(10)
      dishcost   = params(11)
      divsmooth  = params(12)
      wage       = params(13)

      open (23,file="Output/otherparam.txt")
      read(23,*) fixcost
      close(23)

      beeta = 1.0_dp/(1.0_dp+rf)

      !===========================================================================================================
      !       This part sets up the shock grid and the transition matrix
      !===========================================================================================================

      n1   = nz
      str  = 3.0_dp
      tmat = 0.0_dp
      zg   = 0.0_dp
      ve   = sigma**2.0_dp
      call markov(rho,ve,n1,str,tmat,zg)
      zg = exp(zg)

      zmin = minval(zg)
      zmax = maxval(zg)

      !===========================================================================================================
      !       This part sets up the grid endpoints
      !===========================================================================================================
      const  = ((1.0_dp-share)*theta/wage)**(1.0_dp/(1.0_dp-(1.0_dp-share)*theta))
      alpha  = share*theta/(1.0_dp-(1.0_dp-share)*theta)
      skale  = const**((1.0_dp-share)*theta) - const*wage
      aye    = alpha*(1.0_dp - beta)
      bee    = alpha*beta
      zmin   = skale*zmin
      zmax   = skale*zmax

      kmaxstar  = (bee*log(inlta+rf) - log(zmax) - log(aye)-  bee*log(delta+rf)+ bee*log(aye) - bee*log(bee) + log(delta+rf))/(aye + bee - 1.0_dp)
      kminstar  = (bee*log(inlta+rf) - log(zmin) - log(aye)-  bee*log(delta+rf)+ bee*log(aye) - bee*log(bee) + log(delta+rf))/(aye + bee - 1.0_dp)
      kstar     = (bee*log(delta+rf) - 0.0_dp    - log(aye)-  bee*log(delta+rf)+ bee*log(aye) - bee*log(bee) + log(delta+rf))/(aye + bee - 1.0_dp)

      stretch_khi = 0.75_dp!0.5_dp
      stretch_klo = 1.0_dp!10.0_dp
      kmax = kmaxstar+log(stretch_khi)
      kmin = kminstar+log(stretch_klo)

      if (kmax < kmin) then
        kmax = kmaxstar
        kmin = kminstar
      endif

      nmax = 1.4_dp
      nmin = 0.6_dp

      if (nb>3) then
         bmax = xi*(1.0_dp-delta)     ! Do these grids differently because there is no obvious set of technological bounds
         bmin = -xi
      else
         bmax = 0.00001_dp
         bmin = -0.00001_dp
      endif

      endpoints(1).k = exp(kmin)
      endpoints(2).k = exp(kmax)
      endpoints(1).n = nmin
      endpoints(2).n = nmax
      endpoints(1).b = bmin
      endpoints(2).b = bmax


      stretchloop: do iter=1,maxiter
           outofbounds = 1               ! Note: this ONLY gets switched to 0 (good) if there's a solution and the solution is in bounds.

           if (iter > 1) then
              kmin = log(endpoints(1).k)
              kmax = log(endpoints(2).k)
              nmin = endpoints(1).n
              nmax = endpoints(2).n
              bmin = endpoints(1).b
              bmax = endpoints(2).b
              !write(*,*)kmin,kmax,nmin,nmax,bmin,bmax
           endif
                 !if (iter > 1) write(*,"('max, min',6f12.3)")  kmax, kmin, nmax, nmin, bmax, bmin
           n = nk
           call makegrid(kmax,kmin,n,kg)

           n = npk
           call makegrid(kmax,kmin,n,kpg)

           kg = exp(kg)
           kpg = exp(kpg)

           n = nn
           call makegrid(nmax,nmin,n,ng)

           n = npn
           call makegrid(nmax,nmin,n,npg)

           bg  = 0.0_dp
           bpg = 0.0_dp
           if (nb > 1) then
              n = nb
              call makegrid(bmax,bmin,n,bg)
              n = npb
              call makegrid(bmax,bmin,n,bpg)
           endif

           ! Put all of the grids in a structure so that you can move them around more easily

           sgrid(1:nz).z = zg;             pgrid(1:nz).z  = zg
           sgrid(1:nk).k = kg;             pgrid(1:npk).k = kpg
           sgrid(1:nn).n = ng;             pgrid(1:npn).n = npg
           sgrid(1:nb).b = bg;             pgrid(1:npb).b = bpg

           ! Make the interpolation weights

           k_wgt.lo = 0      ;         n_wgt.lo = 0      ;         b_wgt.lo = 0
           k_wgt.hi = 0      ;         n_wgt.hi = 0      ;         b_wgt.hi = 0
           k_wgt.w  = 0.0_dp ;         n_wgt.w  = 0.0_dp ;         b_wgt.w  = 0.0_dp

           call inbetween(sgrid,pgrid,k_wgt,n_wgt,b_wgt)
          ! do ik=1,npk
          !    write(*,*) kpg(ik),k_wgt(ik).lo,k_wgt(ik).hi,k_wgt(ik).w
          ! enddo

           if (compstat == 1) then
              open (unit=53,file="Outputcomp/"//trim(period)//trim(siz)//trim(estype)//"statespaces.txt")
              write(53,*)"======================== profit shock grid =============================="
              write(53,"(f20.14)")zg
              write(53,*)"======================== transition matrix==============================="
              do ii=1,(nz)
                write(53,"(50f10.6)")tmat(:,ii)
              enddo
              write(53,*) "  "
              write(53,*) " capital state space"
              write(53,"(1x,f20.4)") kg
              write(53,*) " "
              write(53,*) " capital policy space"
              write(53,"(1x,f20.4)") kpg
              write(53,*) "  "
              write(53,*) " intangible state space"
              write(53,"(1x,f20.4)") ng
              write(53,*) " "
              write(53,*) " intangible policy space"
              write(53,"(1x,f20.4)") npg
              write(53,*) " "
              write(53,*) " debt state space"
              write(53,"(1x,f20.4)") bg
              write(53,*) " "
              write(53,*) " debt policy space"
              write(53,"(1x,f20.4)") bpg
              write(53,*) " "

           endif


           if (verbose == 1) write(*,*) "Done with setup of grids at      ",omp_get_wtime()-start," secs."
           !===========================================================================================
           ! Now make the profit flow and initialize everything for the value function iteration
           !===========================================================================================



           dish = 0.0_dp

           if (dishcost > 0.0_dp) then
              !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iik,ik,iib,ib)
              !$OMP DO COLLAPSE(4)
              do iib=1,npb
                 do iik=1,npk
                    do ib=1,nb
                       do ik=1,nk
                          if (bpg(iib) > 0.0_dp) then
                             if (bg(ib) < 0.0_dp) then
                                dish(ik,ib,iik,iib) = bpg(iib)*kpg(iik)
                             else
                                dish(ik,ib,iik,iib) = bpg(iib)*kpg(iik) - bg(ib)*kg(ik)
                             endif
                          endif
                       enddo
                    enddo
                 enddo
              enddo
              !$OMP END DO
              !$OMP END PARALLEL
           endif

           !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,in,ib,iik,iin,iib,iz,ie,ize,invest,innova,ncf)
           !$OMP DO COLLAPSE(6)
           do iib=1,npb
             do iin=1,npn
               do iik=1,npk
                 do ib=1,nb
                   do in=1,nn
                     do ik=1,nk
                       do iz=1,nz
                         !This is slow but not that much slower.
                         invest                                        = (kpg(iik)          - kg(ik)*(1.0_dp-delta))
                         innova                                        = (npg(iin)*kpg(iik) - ng(in)*kg(ik)*(1.0_dp-inlta))
                         !if (invest < 0.0_dp) invest = invest*adjK! = invest + adjK*kg(ik)
                         !if (innova < 0.0_dp) invest = innova*adjN! = innova + adjN*kg(ik)*ng(in)


                         profit(iz,ik,in,ib,iik,iin,iib) = skale*zg(iz)*(kg(ik)**aye)*((ng(in)*kg(ik))**bee)            &
                                                         - bg(ib)*kg(ik) + bpg(iib)*kpg(iik)/(1.0_dp + rf*(1.0_dp-tau)) &
                                                         - dish(ik,ib,iik,iib)*dishcost  &
                                                         - invest   &
                                                         - innova   &
                                                         - 0.5_dp*adjK*(invest**2.0_dp)/kg(ik)  &
                                                         - 0.5_dp*adjN*(innova**2.0_dp)/(kg(ik)*ng(in))  &
                                                         - fixcost*skale
                         if (profit(iz,ik,in,ib,iik,iin,iib) < 0.0_dp) then
                             profit(iz,ik,in,ib,iik,iin,iib) = profit(iz,ik,in,ib,iik,iin,iib)*(1.0_dp + lambda)     ! - kpg(iik)*(lambda)
                         endif
                       enddo
                     enddo
                   enddo
                 enddo
               enddo
             enddo
           enddo
           !$OMP END DO
           !$OMP END PARALLEL

           !$OMP PARALLEL
           !$OMP WORKSHARE
           where (profit .ne. 0.0_dp)
             profit = profit - divsmooth*profit**2.0_dp
           end where
           !$OMP END WORKSHARE
           !$OMP END PARALLEL

           if (verbose == 1) write(*,*) "Done with setup of profit flow   ",omp_get_wtime()-start," secs."

           v    = profit(1:nz,1:nk,1:nn,1:nb,1,1,1)
           vnew = 0.0_dp
          ! gidx.k = 1
          ! gidx.n = 1
          ! gidx.b = 1
          ! sidx.k = 1
          ! sidx.n = 1
          ! sidx.b = 1
           oidx.k = 1
           oidx.n = 1
           oidx.b = 1
           errcode = 0
           diff  = 100.0_dp
           pdiff = 1000
           inarow = 100
           suminarow = 1000
           counter = 0
           !======================================value function iteration=============================
           convergeloop: do

                !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(in,ib)
                !$OMP DO COLLAPSE(2)
                do ib = 1,nb
                   do in = 1,nn
                      queue(:,:,in,ib) = matmul(tmat,v(:,:,in,ib))
                   enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL
                !if (counter > 0) then
                !write(*,*) queue(1,1,1,1)
                !write(*,*) " "
                !write(*,*) sum(tmat(:,1)*v(:,1,1,1))
                !write(*,*) sum(tmat(1,:)*v(:,1,1,1))
                !write(*,*) sum(tmat(1,:))
                !write(*,*) sum(tmat(:,1))
                !stop
                !endif
                !write(*,*) "Done with matrix multiplication ",omp_get_wtime()-start," secs."

                if (nk == npk .and. nn == npn .and. nb == npb) then
                    queuelong = queue
                else
                    queuelong = 0.0_dp
                    call interp(queue,queuelong,k_wgt,n_wgt,b_wgt)
                    !if (counter > 0) then
                    !  write(*,"(f12.3)") queue(1,1,1,:)
                    !  write(*,*) " "
                    !  write(*,"(f12.3)") queuelong(1,1,1,:)
                    !  write(*,*) " "
                    !  write(*,"(f12.3)") queue(1,1,:,1)
                    !  write(*,*) " "
                    !  write(*,"(f12.3)") queuelong(1,1,:,1)
                    !  write(*,*) " "
                    !  write(*,"(f12.3)") queue(nz,:,nn,nb)
                    !  write(*,*) " "
                    !  write(*,"(f12.3)") queuelong(nz,:,npn,npb)
                    !  stop
                    !endif
                endif
                !write(*,*) "Done with interpolation ",omp_get_wtime()-start," secs."


               ! Now we allocate the variables to the policy dimensions
               ! Make the Bellman equation

!               if (counter .le. 100000 .or. suminarow < 1) then        ! The first if is to get a good guess before you restrict the search space
               if (counter .le. 1 .or. suminarow < 1) then        ! The first if is to get a good guess before you restrict the search space
                                                                  ! The second is if Howard
                  if (suminarow > 0) then

                     if (verbose == 1)      write(*,*) "starting maximization: ", omp_get_wtime()-start
                     vnew = -huge(0.0)
                     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(xtemp,iib,iin,iik,ib,in,ik,iz)
                     do iib=1,npb
                        do iin=1,npn
                           do iik=1,npk
                              !$OMP DO COLLAPSE(4)
                              do ib=1,nb
                                 do in=1,nn
                                    do ik=1,nk
                                       do iz=1,nz
                                          xtemp = beeta*queuelong(iz,iik,iin,iib)+profit(iz,ik,in,ib,iik,iin,iib)
                                          if (xtemp > vnew(iz,ik,in,ib)) then
                                            vnew(iz,ik,in,ib) = xtemp
                                            gidx(iz,ik,in,ib).k = iik
                                            gidx(iz,ik,in,ib).n = iin
                                            gidx(iz,ik,in,ib).b = iib
                                          endif
                                       enddo
                                    enddo
                                 enddo
                              enddo
                              !OMP END DO
                           enddo
                        enddo
                     enddo
                     !$OMP END PARALLEL
                     if (verbose == 1)      write(*,*) "ending maximization: ", omp_get_wtime()-start

                  else

                     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,in,ib,iz,iik,iin,iib)
                     !$OMP DO COLLAPSE(4)
                     do ib=1,nb
                        do in=1,nn
                           do ik=1,nk
                              do iz=1,nz
                                 iik = gidx(iz,ik,in,ib).k
                                 iin = gidx(iz,ik,in,ib).n
                                 iib = gidx(iz,ik,in,ib).b
                                 vnew(iz,ik,in,ib) = beeta*queuelong(iz,iik,iin,iib) + profit(iz,ik,in,ib,iik,iin,iib)
                              enddo
                           enddo
                        enddo
                     enddo
                     !$OMP END DO
                     !$OMP END PARALLEL

                  endif
                else
                  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(avsmall, iz, ik, in, ib, jk, jn, jb, u_k, u_n, u_b, l_k, l_n, l_b)
                  !$OMP DO COLLAPSE(4)
                  do ib = 1,nb
                     do in = 1,nn
                        do ik = 1,nk
                           do iz = 1, nz
                               jk = gidx(iz,ik,in,ib).k
                               jn = gidx(iz,ik,in,ib).n
                               jb = gidx(iz,ik,in,ib).b

                               u_k = jk + pad + max(1+pad-jk,0) - max(jk+pad-npk, 0)
                               l_k = jk - pad + max(1+pad-jk,0) - max(jk+pad-npk, 0)
                               u_n = jn + pad + max(1+pad-jn,0) - max(jn+pad-npn, 0)
                               l_n = jn - pad + max(1+pad-jn,0) - max(jn+pad-npn, 0)
                               u_b = jb + pad + max(1+pad-jb,0) - max(jb+pad-npb, 0)
                               l_b = jb - pad + max(1+pad-jb,0) - max(jb+pad-npb, 0)

                               avsmall = beeta * queuelong(iz, l_k:u_k, l_n:u_n, l_b:u_b)  &
                                       + profit(iz,ik,in,ib,l_k:u_k, l_n:u_n, l_b:u_b)
                               iik = maxloc(maxval(maxval(avsmall,dim=3),dim=2),dim=1)
                               iin = maxloc(maxval(maxval(avsmall,dim=3),dim=1),dim=1)
                               iib = maxloc(maxval(maxval(avsmall,dim=2),dim=1),dim=1)
                               sidx(iz,ik,in,ib).k = iik
                               sidx(iz,ik,in,ib).n = iin
                               sidx(iz,ik,in,ib).b = iib
                               vnew(iz,ik,in,ib)   = avsmall(iik,iin,iib)
                           end do
                        end do
                     end do
                  end do
                  !$OMP END DO
                  !$OMP END PARALLEL
                  gidx.k = gidx.k + sidx.k - pad - 1 + max(pad+1-gidx.k,0) - max(gidx.k+pad-npk,0)
                  gidx.n = gidx.n + sidx.n - pad - 1 + max(pad+1-gidx.n,0) - max(gidx.n+pad-npn,0)
                  gidx.b = gidx.b + sidx.b - pad - 1 + max(pad+1-gidx.b,0) - max(gidx.b+pad-npb,0)

                endif
                !write(*,*) "Done with maximization ",omp_get_wtime()-start," secs."
                pdiff = sum(abs(oidx.k-gidx.k))+sum(abs(oidx.n-gidx.n))+sum(abs(oidx.b-gidx.b))
                inarow(2:size(inarow)) = inarow(1:size(inarow)-1)
                inarow(1) = pdiff
                suminarow = sum(inarow)
                vdiff = vnew-v
                adjfactor  = 0.5_dp*( (beeta/(1.0_dp-beeta)) * minval(vdiff) + (beeta/(1.0_dp-beeta)) * maxval(vdiff))

                !vdiff = abs(vdiff)
                diff = maxval(abs(vdiff))

                v=vnew+adjfactor
                oidx.k = gidx.k
                oidx.n = gidx.n
                oidx.b = gidx.b
                !write(*,*) "Done with update ",omp_get_wtime()-start," secs."
                counter = counter + 1
                if (verbose == 1) write(*,*)counter,diff, pdiff
                if (diff > 10.0e20 .and. counter > 1)   errcode = 2
                if (counter > capt) errcode = 5


               exitloop = diff < toler .and. pdiff < 1
               if (exitloop .or. errcode > 0 ) then
                 if( errcode == 0) then
                   vg = v
                 else
                   vg = 0.0_dp
                 endif
                 exit convergeloop
               endif
           enddo convergeloop

           if (mod(real(nfcnev),100.)==0.) then
              write(*,"('Error Code = ',i4,3x,'Iterations = ',i4,3x,'SA iteration = ',i8,3x,'MC iteration = ',i8)") errcode, counter, nfcnev, icount
           endif

           if (compstat == 1) then
               write(53,*) " Error code                            ", errcode
               write(53,*) " Number of iterations to convergence   ", counter
               close(53)
           endif

           if (errcode == 0) then
               if (verbose == 1) write(*,*) "Done with VFI                ",omp_get_wtime()-start," secs."

               !=============================================================================================================
               ! Construct the policy functions from the policy index functions.
               !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iz,ik,in,ib,iik,iin,iib)
               !$OMP DO COLLAPSE(4)
               do ib=1,nb
                  do in=1,nn
                     do ik=1,nk
                        do iz=1,nz
                            iik=gidx(iz,ik,in,ib).k
                            iin=gidx(iz,ik,in,ib).n
                            iib=gidx(iz,ik,in,ib).b
                            pol(iz,ik,in,ib).k = kpg(iik)
                            pol(iz,ik,in,ib).n = npg(iin)
                            pol(iz,ik,in,ib).b = bpg(iib)
                            pol(iz,ik,in,ib).d = profit(iz,ik,in,ib,iik,iin,iib)
                        enddo
                     enddo
                  enddo
               enddo
               !$OMP END DO
               !$OMP END PARALLEL

               !=============================================================================================================

               if (verbose == 1) write(*,*) "Done with policy functions        ",omp_get_wtime()-start," secs."

               !=======Update the output===================

               !$OMP PARALLEL
               !$OMP WORKSHARE
               vg = vnew
               !$OMP END WORKSHARE
               !$OMP END PARALLEL

               if (verbose == 1) write(*,*) "Done with final updates ",omp_get_wtime()-start," secs."



           else
               vg   = 0.0_dp
               vnew = 0.0_dp
               tmat = 0.0_dp
               mew  = 0.0_dp
               outofbounds = 0 ! if it fails, there is no point adjusting the bounds
           endif

           if (errcode == 0) then
               simdata = 0.0_dp
               call simpanel(params,v,pol,sgrid,tmat)
               call checkbounds(params,pgrid,endpoints,outofbounds,iter,maxiter)

           endif
           if (outofbounds == 0) then
              exit stretchloop
           endif
      enddo stretchloop

      if (errcode == 0) then

         mew = 1.0_dp/real(size(mew))
         if (mewave == 1) then
             allocate(interp_wgt(nz,nk,nn,nb,npk,npn,npb))
             allocate(coarse_idx(nz,nk,nn,nb))
             !=============================================================================================================
             ! This part makes the matrix of interpolation weights and coarse indices for the stationary distribution.
             !=============================================================================================================
             call makeweight(gidx,pol,sgrid,pgrid,k_wgt,n_wgt,b_wgt,interp_wgt,coarse_idx)
             if (verbose == 1) write(*,*) "Done with interpolation weights   ",omp_get_wtime()-start," secs."
             !=============================================================================================================
             ! This part makes the stationary distribution
             !=============================================================================================================

             mew_new = mew

             do id = 1,distmax
               call makedist(coarse_idx,gidx,interp_wgt,tmat,mew_new)
               disterr = maxval(abs(mew-mew_new))
               if (disterr < disttol) then
                  exit
               else
                  mew = mew_new
               endif
               !write(*,*)"id,disterr",id,disterr
             enddo
             !write(*,*) sum(mew); stop
             mew = mew/sum(mew) !get rid of rounding error
             if (verbose == 1) write(*,*) "Done with stationary distribution ",omp_get_wtime()-start," secs."
         endif
         !=============================================================================================================

         pol.l = const*((pol.k**beta) * ((pol.n*pol.k)**(1.0_dp - beta)))**alpha
         pol.y = skale*(pol.k**aye) * ((pol.n*pol.k)**bee)
         !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iz)
         !$OMP DO
         do iz = 1,nz
           pol(iz,:,:,:).l = zg(iz)*pol(iz,:,:,:).l
           pol(iz,:,:,:).y = zg(iz)*pol(iz,:,:,:).y
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         pol.y = pol.y + wage*pol.l

         !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik)
         !$OMP DO
         do ik = 1,nk
           pol(:,ik,:,:).i = pol(:,ik,:,:).k - (1.0_dp - delta)*kg(ik)
         enddo
         !$OMP END DO
         !$OMP END PARALLEL

         !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(in,ik)
         !$OMP DO COLLAPSE(2)
         do ik = 1,nk
            do in = 1,nn
              pol(:,ik,in,:).r = pol(:,ik,in,:).k*pol(:,ik,in,:).n - (1.0_dp - inlta)*kg(ik)*ng(in)
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL

         pol.e = pol.d
         !$OMP PARALLEL
         !$OMP WORKSHARE



         where (pol.e > 0.0_dp)
              pol.e = 0.0_dp
         end where

         !where (pol.e < 0.0_dp)    !only for fixed costs
         !     pol.e = 1.0_dp
         !end where




        ! where (pol.d < 0.0_dp)
        !      pol.d = 0.0_dp
        ! end where
         !$OMP END WORKSHARE
         !$OMP END PARALLEL

         !pol.y = pol.y + tau*rf*pol.b - divsmooth*pol.d**2.0_dp  ! + pol.e*lambda           ! deadweight costs in this model
         pol.y = pol.y + tau*rf*pol.b + pol.e*lambda              ! deadweight costs in this model    linear costs
         !pol.y = pol.y + tau*rf*pol.b - pol.e*pol.k*lambda       ! deadweight costs in this model


         if (mewave==1) then
            call aggregates(params,mew,pol,sgrid)
         endif

         if (compstat == 1) then
            call polplots(params,v,pol,mew,sgrid,pgrid)
         endif
         if (verbose == 1) write(*,*) "Done with plots                   ",omp_get_wtime()-start," secs."

      endif

      deallocate(profit)
      deallocate(avsmall)
      deallocate(queue)
      deallocate(queuelong)
      deallocate(vdiff)

      !=========================Deallocate all of the arrays=================================
      !deallocate(tmat)
      !deallocate(zg)
      !deallocate(pg)
      !deallocate(kg)
      !deallocate(bg)
      !deallocate(kpg)
      !deallocate(bpg)
      !deallocate(profit)
      !deallocate(invest)
      !deallocate(allval)
      !deallocate(v)
      !deallocate(vnew)
      !deallocate(queue)
      !deallocate(queuelong)
      !deallocate(pol)
      !deallocate(gidx)
      !deallocate(interp_wgt)
      !deallocate(coarse_idx)
      !deallocate(mew)
      !deallocate(mew_new)


end subroutine valfun



subroutine makegrid(kmax,kmin,n,k)
      use datatype
      use sizes
      use globals
      use omp_lib
      implicit none

      real(dp), intent(in) :: kmax,kmin
      integer,  intent(in) :: n
      real(dp), intent(out) :: k(n)

      integer :: ii

      k = 0.0_dp
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
      !$OMP DO SCHEDULE(dynamic)
      do ii=1,n
          k(ii) = kmin+real(ii-1)*(kmax-kmin)/real(n-1)
      end do
      !$OMP END DO
      !$OMP END PARALLEL

end subroutine makegrid


subroutine inbetween(sgrid,pgrid,k_wgt,n_wgt,b_wgt)    !  FIXED
      use datatype
      use sizes
      use globals
      use omp_lib
    implicit none

    type(grids), intent(in), dimension(gmax) :: sgrid
    type(grids), intent(in), dimension(egmax) :: pgrid
    type(gridinterps), intent(out) :: k_wgt(npk)
    type(gridinterps), intent(out) :: n_wgt(npn)
    type(gridinterps), intent(out) :: b_wgt(npb)

    integer :: kdown, kup, ndown, nup, bdown, bup
    real(dp) :: kfrac, nfrac, bfrac

    integer :: ik, in, ib

    ! I should generalize this, but I'm feeling lazy
    do ik = 1,npk
       kdown = floor(real(nk-1)*(real(ik-1)/real(npk-1))) + 1
       if (ik == npk) then
         kup = kdown
         kfrac = 1.0_dp
       else
         kup = kdown + 1
         kfrac = (log(pgrid(ik).k)-log(sgrid(kdown).k)) /(log(sgrid(kup).k)-log(sgrid(kdown).k))
       endif
       k_wgt(ik).lo = kdown
       k_wgt(ik).hi = kup
       k_wgt(ik).w  = kfrac

    enddo

    do in = 1,npn
       ndown = floor(real(nn-1)*(real(in-1)/real(npn-1))) + 1
       if (in == npn) then
         nup = ndown
         nfrac = 1.0_dp
       else
         nup = ndown + 1
         nfrac = (pgrid(in).n-sgrid(ndown).n) /(sgrid(nup).n-sgrid(ndown).n)
       endif
       n_wgt(in).lo = ndown
       n_wgt(in).hi = nup
       n_wgt(in).w  = nfrac

    enddo


    do ib = 1,npb
       bdown = floor(real(nb-1)*(real(ib-1)/real(npb-1))) + 1
       if (ib == npb) then
         bup = bdown
         bfrac = 1.0_dp
       else
         bup = bdown + 1
         bfrac = (pgrid(ib).b-sgrid(bdown).b) /(sgrid(bup).b-sgrid(bdown).b)
       endif
       b_wgt(ib).lo = bdown
       b_wgt(ib).hi = bup
       b_wgt(ib).w  = bfrac

    enddo

end subroutine inbetween

subroutine interp(obj,objlong,k_wgt,n_wgt,b_wgt)         ! FIXED
      use datatype
      use sizes
      use globals
      use omp_lib
    implicit none

    real(dp), intent(in) :: obj(nz,nk,nn,nb)
    real(dp), intent(out) :: objlong(nz,npk,npn,npb)
    type(gridinterps), intent(in) :: k_wgt(npk)
    type(gridinterps), intent(in) :: n_wgt(npn)
    type(gridinterps), intent(in) :: b_wgt(npb)

    real(dp) :: obb1, obb2, obb11, obb22
    integer :: iz, ik, in, ib


    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iz, ik, in, ib, obb1, obb2, obb11, obb22)
    !$OMP DO COLLAPSE(4)
    do ib = 1,npb
       do in = 1,npn
          do ik = 1,npk
             do iz = 1,nz
                obb1              = b_wgt(ib).w*obj(iz,k_wgt(ik).hi, n_wgt(in).hi, b_wgt(ib).hi)  + (1.0_dp - b_wgt(ib).w)*obj(iz,k_wgt(ik).hi, n_wgt(in).hi, b_wgt(ib).lo)
                obb2              = b_wgt(ib).w*obj(iz,k_wgt(ik).hi, n_wgt(in).lo, b_wgt(ib).hi)  + (1.0_dp - b_wgt(ib).w)*obj(iz,k_wgt(ik).hi, n_wgt(in).lo, b_wgt(ib).lo)

                obb11             = n_wgt(in).w*obb1  + (1.0_dp - n_wgt(in).w)*obb2

                obb1              = b_wgt(ib).w*obj(iz,k_wgt(ik).lo, n_wgt(in).hi, b_wgt(ib).hi)  + (1.0_dp - b_wgt(ib).w)*obj(iz,k_wgt(ik).lo, n_wgt(in).hi, b_wgt(ib).lo)
                obb2              = b_wgt(ib).w*obj(iz,k_wgt(ik).lo, n_wgt(in).lo, b_wgt(ib).hi)  + (1.0_dp - b_wgt(ib).w)*obj(iz,k_wgt(ik).lo, n_wgt(in).lo, b_wgt(ib).lo)

                obb22             = n_wgt(in).w*obb1  + (1.0_dp - n_wgt(in).w)*obb2

                objlong(iz,ik,in,ib) = k_wgt(ik).w*obb11 + (1.0_dp - k_wgt(ik).w)*obb22
             enddo
          enddo
       enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

end subroutine interp


subroutine markov(a,s2,nn,m,trans,z)
   implicit none

   real(8), intent(in) :: a, s2, m
   integer, intent(in) :: nn
   real(8), intent(out), dimension(1:nn) :: z
   real(8), intent(out), dimension(1:nn,1:nn) :: trans

   real(8) :: sy, s, w, f1(1:nn,1:nn), f2(1:nn,1:nn), minif
   integer :: ii, jj

   f1 = 0.0
   f2 = 1.0

   z = 0.0
   trans = 0.0

   sy = sqrt(s2/(1.0-a**2))
   do ii=1,nn
     z(ii) = -m*sy + ((2.0*m*sy)/(real(nn)-1.0))*(real(ii)-1.0)
   enddo

   s = sqrt(s2)
   w = z(2) - z(1)

   do jj=1,nn
      do ii=2,nn
        minif = (z(ii)-a*z(jj)-w/2.0)/s
        f1(ii,jj)   = 0.5*(1.0 + derf(minif/sqrt(2.0)))
        f2(ii-1,jj) = 0.5*(1.0 + derf(minif/sqrt(2.0)))
      enddo
  !    do ii=1,(nn-1)
  !      minif = (z(ii)-a*z(jj)-w/2.0)/s
  !      f2(ii,jj) = 0.5*(1.0 + derf(minif/sqrt(2.0)))
  !    enddo
   enddo

   f1 = transpose(f1)
   f2 = transpose(f2)

   !do ii=1,nn
   ! write (21,"(15(f8.3))") f1(ii,:)
   ! write (*,"(15(f8.3))") f1(ii,:)
   !enddo
   !do ii=1,nn
   ! write (21,"(15(f8.3))") f2(ii,:)
   ! write (*,"(15(f8.3))") f2(ii,:)
   !enddo

   trans = f2 - f1

return
end subroutine markov

subroutine makedist(coarse_idx,gidx,interp_wgt,tmat,mew)      !!!! FIXED
   use datatype
   use sizes
   use globals
   implicit none

   type(weights), intent(in), dimension(nz,nk,nn,nb,npk,npn,npb) :: interp_wgt   ! Interpolation weights
   type(ipolicies), intent(in), dimension(nz,nk,nn,nb) :: coarse_idx             ! coarse policy index array
   type(ipolicies), intent(in), dimension(nz,nk,nn,nb) :: gidx                   ! actual policy index array
   real(dp), intent(in) :: tmat(nz,nz)                                        ! Transition array
   real(dp), intent(inout) :: mew(nz,nk,nn,nb)                                   ! stationary distribution

   real(dp) :: mew_old(nz,nk,nn,nb), k_weight, n_weight, b_weight
   integer :: k_increm, n_increm, b_increm, kidx, sidx, bidx
   integer :: iz, ik, in, ib, ize


   mew_old = mew
   mew = 0.0_dp  !initialize the new distribution
   !x$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iz, ik, in, ib, ize, k_weight, n_weight, b_weight, k_increm, n_increm, b_increm, kidx, sidx, bidx)
   !x$OMP DO COLLAPSE(4)
   do iz = 1,nz
      do ib = 1,nb
         do in = 1,nn
            do ik = 1,nk
               if (mew_old(iz,ik,in,ib)>0.0_dp) then
                  k_weight = interp_wgt(iz,ik,in,ib,gidx(iz,ik,in,ib).k,gidx(iz,ik,in,ib).n,gidx(iz,ik,in,ib).b).k
                  n_weight = interp_wgt(iz,ik,in,ib,gidx(iz,ik,in,ib).k,gidx(iz,ik,in,ib).n,gidx(iz,ik,in,ib).b).n
                  b_weight = interp_wgt(iz,ik,in,ib,gidx(iz,ik,in,ib).k,gidx(iz,ik,in,ib).n,gidx(iz,ik,in,ib).b).b
                  k_increm = min(nk-coarse_idx(iz,ik,in,ib).k,1)
                  n_increm = min(nn-coarse_idx(iz,ik,in,ib).n,1)
                  b_increm = min(nb-coarse_idx(iz,ik,in,ib).b,1) !The min thing is for this when the policy is always glued up against the constraint.
                  do ize = 1,nz
                     if (tmat(ize,iz)>0.0_dp) then     ! I don't feel like making 2^n automatic with alphabetic indices.

                        kidx=coarse_idx(iz,ik,in,ib).k
                        sidx=coarse_idx(iz,ik,in,ib).n
                        bidx=coarse_idx(iz,ik,in,ib).b

                        mew(ize,kidx,sidx,bidx) = &
                        mew(ize,kidx,sidx,bidx) +                               tmat(ize,iz)*mew_old(iz,ik,in,ib) * k_weight*n_weight*b_weight
                        if (npk > nk .or. npn > nn .or. npb > nb ) then
                           mew(ize,kidx+k_increm,sidx,bidx) = &
                           mew(ize,kidx+k_increm,sidx,bidx) +                   tmat(ize,iz)*mew_old(iz,ik,in,ib) * (1.0_dp-k_weight)*n_weight*b_weight

                           mew(ize,kidx,sidx,bidx+b_increm) = &
                           mew(ize,kidx,sidx,bidx+b_increm) +                   tmat(ize,iz)*mew_old(iz,ik,in,ib) * k_weight*n_weight*(1.0_dp-b_weight)

                           mew(ize,kidx,sidx+n_increm,bidx) = &
                           mew(ize,kidx,sidx+n_increm,bidx) +                   tmat(ize,iz)*mew_old(iz,ik,in,ib) * k_weight*(1.0_dp-n_weight)*b_weight

                           mew(ize,kidx+k_increm,sidx,bidx+b_increm) = &
                           mew(ize,kidx+k_increm,sidx,bidx+b_increm) +          tmat(ize,iz)*mew_old(iz,ik,in,ib) * (1.0_dp-k_weight)*n_weight*(1.0_dp-b_weight)

                           mew(ize,kidx+k_increm,sidx+n_increm,bidx) = &
                           mew(ize,kidx+k_increm,sidx+n_increm,bidx) +          tmat(ize,iz)*mew_old(iz,ik,in,ib) * (1.0_dp-k_weight)*(1.0_dp-n_weight)*b_weight

                           mew(ize,kidx,sidx+n_increm,bidx+b_increm) = &
                           mew(ize,kidx,sidx+n_increm,bidx+b_increm) +          tmat(ize,iz)*mew_old(iz,ik,in,ib) * k_weight*(1.0_dp-n_weight)*(1.0_dp-b_weight)

                           mew(ize,kidx+k_increm,sidx+n_increm,bidx+b_increm) = &
                           mew(ize,kidx+k_increm,sidx+n_increm,bidx+b_increm) + tmat(ize,iz)*mew_old(iz,ik,in,ib) * (1.0_dp-k_weight)*(1.0_dp-n_weight)*(1.0_dp-b_weight)
                        endif
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo
   enddo
   !x$OMP END DO
   !x$OMP END PARALLEL

end subroutine makedist

subroutine makeweight(gidx,pol,sgrid,pgrid,k_wgt,n_wgt,b_wgt,interp_wgt,coarse_idx)  ! FIXED!!!
use datatype
use sizes
use globals
use omp_lib
 type(ipolicies), intent(in), dimension(nz,nk,nn,nb) :: gidx
 type(policies), intent(in), dimension(nz,nk,nn,nb) :: pol
 type(grids), intent(in), dimension(gmax) :: sgrid
 type(grids), intent(in), dimension(egmax) :: pgrid
 type(gridinterps), intent(in) :: k_wgt(npk)
 type(gridinterps), intent(in) :: n_wgt(npn)
 type(gridinterps), intent(in) :: b_wgt(npb)

 type(weights), intent(out), dimension(nz,nk,nn,nb,npk,npn,npb) :: interp_wgt
 type(ipolicies), intent(out), dimension(nz,nk,nn,nb) :: coarse_idx

 integer :: ib, ik, in, iz, iik, iin, iib, klo, nlo, blo

    if (nk == npk .and. nb == npb.and. nn == npn) then

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iz,ik,in,ib,iik,iin,iib,klo,blo)
       !$OMP DO COLLAPSE(4)
       do ib = 1,nb      ! Turn this into a matrix of interpolation weights
          do in = 1,nn
             do ik = 1,nk
                do iz = 1,nz
                   ! This bit just assigns a 1 to the big array when the policy intersects with the grid value.
                   coarse_idx(iz,ik,in,ib).k = gidx(iz,ik,in,ib).k
                   coarse_idx(iz,ik,in,ib).n = gidx(iz,ik,in,ib).n
                   coarse_idx(iz,ik,in,ib).b = gidx(iz,ik,in,ib).b
                   do iib=1,npb
                   do iin=1,npn
                   do iik=1,npk
                      if (pol(iz,ik,in,ib).k == pgrid(iik).k) then
                         interp_wgt(iz,ik,in,ib,iik,iin,iib).k = 1.0_dp
                      endif
                      if (pol(iz,ik,in,ib).n == pgrid(iin).n) then
                         interp_wgt(iz,ik,in,ib,iik,iin,iib).n = 1.0_dp
                      endif
                      if (pol(iz,ik,in,ib).b == pgrid(iib).b) then
                         interp_wgt(iz,ik,in,ib,iik,iin,iib).b = 1.0_dp
                      endif
                   enddo
                   enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
       !$OMP END DO
       !$OMP END PARALLEL

    else


       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iz,ik,in,ib)
       !$OMP DO COLLAPSE(4)
       do ib = 1,nb      ! Turn this into a matrix of interpolation weights
          do in = 1,nn
             do ik = 1,nk
                do iz = 1,nz
                   !This bit does the same except instead of ones and zeros you have
                   !the fractional distance between the optimal point and the lower
                   !boundary grid point.
                   !klo = 1
                   !call hunt(sgrid(1:nk).k,nk,pol(iz,ik,in,ib).k,klo)
                   !write(*,*) klo  , "one"


                   !klo = k_wgt(gidx(iz,ik,in,ib).k).lo  ! The index is the entry in the STATE grid
                   !coarse_idx(iz,ik,in,ib).k = klo      ! that is just below the optimal policy
                   !klo = k_wgt(gidx(iz,ik,in,ib).k).lo                   ! The index is the entry in the STATE grid
                   coarse_idx(iz,ik,in,ib).k = k_wgt(gidx(iz,ik,in,ib).k).lo ! that is just below the optimal policy
                   !write(*,*) klo, "two";
                   !stop
                   !blo = b_wgt(gidx(iz,ik,in,ib).b).lo  ! The index is the entry in the STATE grid
                   !coarse_idx(iz,ik,in,ib).b = blo      ! that is just below the optimal policy
                   !blo = b_wgt(gidx(iz,ik,in,ib).b).lo                   ! The index is the entry in the STATE grid
                   coarse_idx(iz,ik,in,ib).n = n_wgt(gidx(iz,ik,in,ib).n).lo ! that is just below the optimal policy
                   coarse_idx(iz,ik,in,ib).b = b_wgt(gidx(iz,ik,in,ib).b).lo ! that is just below the optimal policy
                enddo
             enddo
          enddo
       enddo
       !$OMP END DO
       !$OMP END PARALLEL



       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iz,ik,in,ib,iik,iin,iib)
       !$OMP DO COLLAPSE(7)
       do ib = 1,nb      ! Turn this into a matrix of interpolation weights
          do in = 1,nn
             do ik = 1,nk
                do iz = 1,nz
                   do iib = 1,npb
                      do iin = 1,npn
                         do iik = 1,npk
                            if (pol(iz,ik,in,ib).k == pgrid(iik).k) then
                               interp_wgt(iz,ik,in,ib,iik,iin,iib).k = k_wgt(iik).w
                               ! NOTE: this is the distance to the LOW point as a fraction of the total
                            endif
                            if (pol(iz,ik,in,ib).n == pgrid(iin).n) then
                               interp_wgt(iz,ik,in,ib,iik,iin,iib).n = n_wgt(iin).w
                               ! NOTE: this is the distance to the LOW point as a fraction of the total
                            endif
                            if (pol(iz,ik,in,ib).b == pgrid(iib).b) then
                               interp_wgt(iz,ik,in,ib,iik,iin,iib).b = b_wgt(iib).w
                               ! NOTE: this is the distance to the LOW point as a fraction of the total
                            endif
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
       !$OMP END DO
       !$OMP END PARALLEL

    endif
    !write(*,*) sum(interp_wgt,dim=1)
    !open(2,file='interp_wgt.asc')
    !  do ic=1,nc
    !     do iz=1,nz
    !        write(2,"(i4,1x,<ncp>f6.2)")iz,interp_wgt(:,ic,iz)
    !        write(2,*) " "
    !     enddo
    !  enddo
    !close(2)

end subroutine makeweight

subroutine polplots(params,v,pol,mew,sgrid,pgrid)
      use datatype
      use sizes
      use globals
      implicit none

      real(dp), intent(in) :: params(nop), v(nz,nk,nn,nb), mew(nz,nk,nn,nb)
      type(grids), intent(in), dimension(gmax) :: sgrid
      type(grids), intent(in), dimension(egmax) :: pgrid
      type(policies), intent(in), dimension(nz,nk,nn,nb) :: pol

      real(dp) :: alpha, beta, rho, sigma, delta, inlta, lambda, adjK, adjN, xi, wage, &       ! parameters to estimate
                  divsmooth, dishcost, kstar, mew_k(nz,nk), mew_n(nz,nn), mew_b(nz,nb)
      integer ::  ik, in, ib, iz, ie, seqa(nz)


      alpha      = params(1)
      beta       = params(2)
      rho        = params(3)
      sigma      = params(4)
      delta      = params(5)
      inlta      = params(6)
      adjK       = params(7)
      adjN       = params(8)
      xi         = params(9)
      lambda     = params(10)
      dishcost   = params(11)
      divsmooth  = params(12)
      wage       = params(13)

      open(5,file='Outputcomp/kgrid.asc')     ! grids
      open(6,file='Outputcomp/ngrid.asc')
      open(7,file='Outputcomp/bgrid.asc')
      open(9,file='Outputcomp/zgrid.asc')
      do ik=1,nk
               write(5,"(f26.12)") sgrid(ik).k
      enddo
      do in=1,nn
               write(6,"(f26.12)") sgrid(in).n
      enddo
      do ib=1,nb
               write(7,"(f26.12)") sgrid(ib).b
      enddo
      do iz=1,nz
               write(9,"(f26.12)") sgrid(iz).z
      enddo
      close(5)
      close(6)
      close(7)
      close(9)

      open(4,file='Outputcomp/v.asc')     ! value function
      open(5,file='Outputcomp/k.asc')     ! optimal capital policy
      open(6,file='Outputcomp/n.asc')     ! optimal labor policy
      open(7,file='Outputcomp/b.asc')     ! optimal debt policy
      open(9,file='Outputcomp/d.asc')     ! optimal dividend policy
      do ik=1,nk
         do in=1,nn
            do ib=1,nb
               write(4,"(1x,3f14.3,1x,<nz>f18.4)") sgrid(ik).k,sgrid(in).n,sgrid(ib).b,  v(:,ik,in,ib)
               write(5,"(1x,3f14.3,1x,<nz>f18.4)") sgrid(ik).k,sgrid(in).n,sgrid(ib).b,pol(:,ik,in,ib).k
               write(6,"(1x,3f14.3,1x,<nz>f18.4)") sgrid(ik).k,sgrid(in).n,sgrid(ib).b,pol(:,ik,in,ib).n
               write(7,"(1x,3f14.3,1x,<nz>f18.4)") sgrid(ik).k,sgrid(in).n,sgrid(ib).b,pol(:,ik,in,ib).b
               write(9,"(1x,3f14.3,1x,<nz>f18.4)") sgrid(ik).k,sgrid(in).n,sgrid(ib).b,pol(:,ik,in,ib).d
            enddo
         enddo
      enddo
      close(4)
      close(5)
      close(6)
      close(7)
      close(9)


      open (20,file='Outputcomp/mew.txt')
         do ib=1,nb
         do in=1,nn
           write(20,*) "  "
           do iz=1,nz
              write(20,"(<nk>f15.10)") mew(iz,:,in,ib)
           enddo
        enddo
        enddo
      close(20)

      mew_k = sum(sum(mew,dim=4),dim=3)
      mew_b = sum(sum(mew,dim=3),dim=2)
      mew_n = sum(sum(mew,dim=4),dim=2)

      open (20,file='Outputcomp/mew_k.txt')
           do iz=1,nz
              write(20,"(<nk>f15.10)") mew_k(iz,:)
           enddo
      close(20)
      open (20,file='Outputcomp/mew_b.txt')
           do iz=1,nz
              write(20,"(<nb>f15.10)") mew_b(iz,:)
           enddo
      close(20)
      open (20,file='Outputcomp/mew_n.txt')
           do iz=1,nz
              write(20,"(<nb>f15.10)") mew_n(iz,:)
           enddo
      close(20)
end subroutine polplots

subroutine checkbounds(params,pgrid,endpoints,outofbounds,iter,maxiter)
      use datatype
      use sizes
      use globals
      implicit none

      type(policies), intent(inout) :: endpoints(2)
      type(grids),    intent(in)    :: pgrid(egmax)

      real(dp), intent(in) ::  params(allnop)
      integer, intent(inout) :: outofbounds
      integer, intent(in) :: iter, maxiter

      type(policies) :: simpoints(2)
      real(dp) :: kpg(npk), npg(npn), bpg(npb)
      real(dp) :: kmax,kmin,nmax,nmin,bmax,bmin
      real(dp) :: maxsimk, minsimk, maxsimn, minsimn, maxsimb, minsimb
      real(dp) :: stretchfac, shrinkfac
      logical  :: toobig_k, toobig_n, hithi_k, hitlo_k, hithi_n, hitlo_n, skewlo_k, skewhi_k, skewlo_n, skewhi_n, bskew, bsmall, notbig_k, notbig_n
      real(dp), parameter :: eps=1.0e-4
      character (len=1) :: what


      outofbounds = 0

      maxsimk = maxval(simdata(2,burn_in+1:nyears,:))
      minsimk = minval(simdata(2,burn_in+1:nyears,:))
      maxsimn = maxval(simdata(3,burn_in+1:nyears,:))
      minsimn = minval(simdata(3,burn_in+1:nyears,:))
      maxsimb = maxval(simdata(4,burn_in+1:nyears,:))
      minsimb = minval(simdata(4,burn_in+1:nyears,:))

      simpoints(1).k = minsimk
      simpoints(2).k = maxsimk
      simpoints(1).n = minsimn
      simpoints(2).n = maxsimn
      simpoints(1).b = minsimb
      simpoints(2).b = maxsimb

      kpg = pgrid(1:npk).k
      npg = pgrid(1:npn).n
      bpg = pgrid(1:npb).b

      if (verbose == 1) write(*,"('      ',2a25  )")  "simulated","current endpoint"
      if (verbose == 1) write(*,"('k max ',2f25.4)")  maxsimk,kpg(npk)
      if (verbose == 1) write(*,"('k min ',2f25.4)")  minsimk,kpg(1)
      if (verbose == 1) write(*,"('n max ',2f25.4)")  maxsimn,npg(npn)
      if (verbose == 1) write(*,"('n min ',2f25.4)")  minsimn,npg(1)

      if (verbose == 1) write(*,"('b max ',2f25.4)")  maxsimb,bpg(npb)
      if (verbose == 1) write(*,"('b min ',2f25.4)")  minsimb,bpg(1)


      !===================================================If the grids are either too narrow or too asymmetric, then reshape everything.

      toobig_k = ((maxsimk - minsimk) < 0.1_dp*(kpg(npk)-kpg(1)))
      toobig_n = ((maxsimn - minsimn) < 0.1_dp*(npg(npn)-npg(1)))

      hithi_k  = (maxsimk > (kpg(npk)-eps))
      hitlo_k  = (minsimk < (kpg(1  )+eps))
      hithi_n  = (maxsimn > (npg(npn)-eps))
      hitlo_n  = (minsimn < (npg(1  )+eps))

      skewlo_k = (maxsimk < kpg(floor(npk*0.5_dp)))
      skewhi_k = (minsimk > kpg(ceiling(npk*0.5_dp)))
      skewlo_n = (maxsimn < npg(floor(npn*0.5_dp)))
      skewhi_n = (minsimn > npg(ceiling(npn*0.5_dp)))

      ! Do debt differently
      !bskew    = (minsimb > bpg(ceiling(npb*0.5_dp)) .and. (maxsimb >(minsimb+0.05_dp)))     ! This is so that you don't compress the state space against the collateral constraint.
      bsmall    = (minsimb < (bpg(1  ) + eps))

      if (verbose == 1) write(*,*)"toobig_k ", toobig_k
      if (verbose == 1) write(*,*)"toobig_n ", toobig_n
      if (verbose == 1) write(*,*)"hithi_k  ", hithi_k
      if (verbose == 1) write(*,*)"hitlo_k  ", hitlo_k
      if (verbose == 1) write(*,*)"hithi_n  ", hithi_n
      if (verbose == 1) write(*,*)"hitlo_n  ", hitlo_n
      if (verbose == 1) write(*,*)"skewlo_k ", skewlo_k
      if (verbose == 1) write(*,*)"skewhi_k ", skewhi_k
      if (verbose == 1) write(*,*)"skewlo_n ", skewlo_n
      if (verbose == 1) write(*,*)"skewhi_n ", skewhi_n
      !if (verbose == 1) write(*,*)"bskew    ", bskew
      if (verbose == 1) write(*,*)"bsmall   ", bsmall

      notbig_k  = .not. toobig_k
      notbig_n  = .not. toobig_n

      !========This part streches k==========

      if (toobig_k) then
          stretchfac = 1.25_dp
          shrinkfac  = 0.75_dp
          what = "k"
          call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
          outofbounds = 1
      elseif (notbig_k .and. hithi_k) then
          stretchfac = 1.25_dp
          shrinkfac  = 1.00_dp
          what = "k"
          call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
          outofbounds = 1
      elseif (notbig_k .and. hitlo_k) then
          stretchfac = 1.00_dp
          shrinkfac  = 0.75_dp
          what = "k"
          call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
          outofbounds = 1
      elseif (notbig_k .and. skewhi_k) then
          stretchfac = 1.00_dp
          shrinkfac  = 0.75_dp
          what = "k"
          call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
          outofbounds = 1
      elseif (notbig_k .and. skewlo_k) then
          stretchfac = 1.25_dp
          shrinkfac  = 1.00_dp
          what = "k"
          call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
          outofbounds = 1
      endif

      !========This part streches n==========

      if (toobig_n) then
          if (maxsimn>(minsimn+eps)) then      ! This is different from k because sometimes intangible ratios don't move
             stretchfac = 1.05_dp
             shrinkfac  = 0.95_dp
          else
             stretchfac = 1.25_dp
             shrinkfac  = 0.75_dp
          endif
          what = "n"
          call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
          outofbounds = 1
      elseif (notbig_n .and. hithi_n) then
          stretchfac = 1.25_dp
          shrinkfac  = 1.00_dp
          what = "n"
          call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
          outofbounds = 1
      elseif (notbig_n .and. hitlo_n) then
          stretchfac = 1.00_dp
          shrinkfac  = 0.75_dp
          what = "n"
          call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
          outofbounds = 1
      elseif (notbig_n .and. skewhi_n) then
          stretchfac = 1.00_dp
          shrinkfac  = 0.75_dp
          what = "n"
          call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
          outofbounds = 1
      elseif (notbig_n .and. skewlo_n) then
          stretchfac = 1.25_dp
          shrinkfac  = 1.00_dp
          what = "n"
          call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
          outofbounds = 1
      endif


      !===================Just do two rounds to make the endpoints tight. Better solution.

      if (iter==1 .and. toobig_k == .false. .and. hithi_k == .false. .and. hitlo_k == .false. .and. skewhi_k == .false. .and. skewlo_k == .false. ) then
          what = "k"
          stretchfac = 1.1_dp
          shrinkfac  = 0.9_dp
          call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
      endif
      if (iter==1 .and. toobig_n == .false. .and. hithi_n == .false. .and. hitlo_n == .false. .and. skewhi_n == .false. .and. skewlo_n == .false. ) then
          what = "n"
          stretchfac = 1.1_dp
          shrinkfac  = 0.9_dp
          call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
          outofbounds = 1
      endif

      !=======This part stretches b==========  !!! Order is important. You don't want to do this before the previous block.

    !  if (bskew) then
    !    what = "b"
    !    stretchfac = 1.00_dp
    !    shrinkfac  = 1.00_dp
    !    call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
    !    outofbounds = 1
      if (bsmall) then
        what = "b"
        stretchfac = 1.00_dp
        shrinkfac  = 1.00_dp
        call stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
        outofbounds = 1
      endif


      if (iter >= maxiter) then
         outofbounds = 0
      endif

end subroutine checkbounds
subroutine stretchbounds(params,endpoints,simpoints,stretchfac,shrinkfac,bskew,bsmall,what)
      use datatype
      use sizes
      use globals
      implicit none

         real(dp), intent(in) :: stretchfac, shrinkfac, params(nop)
         type(policies), intent(inout) :: simpoints(2), endpoints(2)
         logical, intent(in) :: bskew, bsmall
         character(1), intent(in) :: what
         real(dp) :: xi, delta

         xi = params(9)
         delta = params(5)

         if (what == "k") then
            endpoints(1).k = (simpoints(1).k*shrinkfac)
            endpoints(2).k = (simpoints(2).k*stretchfac)
         elseif (what == "n") then
            endpoints(1).n = (simpoints(1).n*shrinkfac)
            endpoints(2).n = (simpoints(2).n*stretchfac)
         elseif (what == "b") then
             if (bsmall) then
               endpoints(1).b = simpoints(1).b - abs(endpoints(2).b-endpoints(1).b)/real(nb)*2.0_dp
             endif
             if (bskew) then
               endpoints(1).b = simpoints(1).b - abs(endpoints(2).b-endpoints(1).b)/real(nb)
             endif
             endpoints(2).b = xi*(1.0_dp - delta)
        endif
end subroutine stretchbounds
