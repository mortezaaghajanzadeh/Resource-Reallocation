          subroutine simpanel(params,v,pol,sgrid,tmat)
          use datatype
          use sizes
          use globals
          implicit none

          type(policies), intent(in), dimension(nz,nk,nn,nb) :: pol
          type(grids),    intent(in), dimension(gmax)        :: sgrid
          real(dp), intent(in) :: tmat(nz,nz)
          real(dp), intent(in) :: v(nz,nk,nn,nb)
          real(dp), intent(in) :: params(nop)

          real(dp), allocatable :: cdf_wgt(:), phatcdf(:,:), nums(:), allv(:,:), allk(:,:), allz(:,:), allb(:,:), alln(:,:), diff(:), &
                                   draws(:,:), alld(:,:), alll(:,:)

          real(dp), allocatable :: k(:,:,:,:), b(:,:,:,:), d(:,:,:,:), n(:,:,:,:)


          integer,  allocatable ::  ls(:,:), zls(:,:), els(:,:)

          real(dp) :: kg(nk), zg(nz), bg(nb), ng(nn)
          real(dp) ::  theonek, theoneb, theonen, theonez, shocktemp, vold, kold, bold, zold, dold, nold, &
                       vprime, kprime, bprime, zprime, dprime, nprime, tw, eprime, eold

          integer :: iz, iiseed(1), ifi, iyi, iti, pickk, pickb, pickn, pickz, enn, erc

          intrinsic random_seed, random_number

          iiseed      = 123456789

          if (gendata == 1) then
             call random_seed(PUT=ooseed)
          else
             call random_seed(PUT=iiseed)
          endif

          !-----------------------------------Now simulate the model.-------------------------
          !
          simdata = 0.0_dp

          allocate (phatcdf(nz,nz))
          allocate (cdf_wgt(nz))
          allocate (draws(nYears+1,nfirms))
          allocate (diff(nz))
          allocate (ls(nyears+1,nfirms))
          allocate (zls(nyears+1,nfirms))
          allocate (k(nz,nk,nn,nb))
          allocate (b(nz,nk,nn,nb))
          allocate (d(nz,nk,nn,nb))
          allocate (n(nz,nk,nn,nb))

          k = pol.k
          n = pol.n
          b = pol.b
          d = pol.d
          kg = sgrid(1:nk).k
          ng = sgrid(1:nn).n
          bg = sgrid(1:nb).b
          zg = sgrid(1:nz).z

          ! Unconditional cdf of the shocks

          cdf_wgt = 1.0_dp/dble(nz)
          do iz = 2,nz
              cdf_wgt(iz) = cdf_wgt(iz) + cdf_wgt(iz-1)
          enddo

          ! Conditional shock cdf
          phatcdf = transpose(tmat)

          do iz=2,nz
            phatcdf(iz,:) = phatcdf(iz,:) + phatcdf(iz-1,:)
          enddo

          allocate(allv(nyears,nfirms))                                                                !
          allocate(allk(nyears,nfirms))                                                                !
          allocate(allz(nyears,nfirms))                                                                !
          allocate(allb(nyears,nfirms))                                                                !
          allocate(alln(nyears,nfirms))                                                                !
          allocate(alll(nyears,nfirms))                                                                !
          allocate(alld(nyears,nfirms))                                                                !
          allv          = 0.0_dp
          allk          = 0.0_dp
          allz          = 0.0_dp
          alln          = 0.0_dp
          alll          = 0.0_dp
          allb          = 0.0_dp
          alld          = 0.0_dp
          if (isnan(phatcdf(1,1))) then
               erc = 4
          else

               !$OMP PARALLEL PRIVATE(shocktemp,iyi,ifi)
               !$OMP DO COLLAPSE(2)
               do iyi = 1,nYears+1
                   do ifi = 1,nfirms
                       call random_number(harvest=shocktemp)
                       draws(iyi,ifi) = shocktemp
                   enddo
               enddo
               !$OMP END DO
               !$OMP END PARALLEL

! This just randomly assigns the starting locations from the unconditional distribution
               do ifi = 1,nfirms
                   ls(1,ifi) = nz
                   diff = draws(1,ifi) - cdf_wgt
                   oneloop: do iz = 1,nz
                     if (diff(iz) < 0.0_dp) then
                        ls(1,ifi) = iz
                        exit oneloop
                     endif
                   enddo oneloop
               enddo

               crosssection: do ifi=1,nfirms


                  ! This is the "time-zero" starting point.
                  call random_number(harvest=theonek)
                  call random_number(harvest=theonen)
                  call random_number(harvest=theoneb)

                  pickk = int(nk*theonek)+1
                  pickn = int(nn*theonen)+1
                  pickb = int(nb*theoneb)+1
                  pickz = ls(1,ifi)

                  vold = v(pickz,pickk,pickn,pickb)
                  kold = k(pickz,pickk,pickn,pickb)
                  bold = b(pickz,pickk,pickn,pickb)
                  dold = d(pickz,pickk,pickn,pickb)
                  nold = n(pickz,pickk,pickn,pickb)

                    timeseries: do iti=1,nyears
                       zls(iti,ifi) = ls(iti,ifi)

                       zold = zg(zls(iti,ifi))
                       !write(*,*)zls(iti,ifi),zold,eold
                       !This business pulls up next period's shock using the transition matrix
                       diff = draws(iti+1, ifi) - phatcdf(:, ls(iti, ifi))
                       ls(iti+1,ifi) = nz   !Initialize it at the end in case the loop never gets there
                       pickloop: do iz=1,nz
                         if (diff(iz) < 0.0_dp) then
                            ls(iti+1,ifi) = iz
                            exit pickloop
                         endif
                       enddo pickloop
                       ! We do the interpolation after selecting out the matrix for the relevant state only.


                       call interpol(ls(iti,ifi),kold,nold,bold,kg,ng,bg,v,vprime)
                       call interpol(ls(iti,ifi),kold,nold,bold,kg,ng,bg,k,kprime)
                       call interpol(ls(iti,ifi),kold,nold,bold,kg,ng,bg,b,bprime)
                       call interpol(ls(iti,ifi),kold,nold,bold,kg,ng,bg,d,dprime)
                       call interpol(ls(iti,ifi),kold,nold,bold,kg,ng,bg,n,nprime)

                       zls(iti+1,ifi) = ls(iti+1,ifi)
                       zprime = zg(zls(iti+1,ifi))

                       allv(iti,ifi) =vprime
                       allk(iti,ifi) =kprime
                       allb(iti,ifi) =bprime
                       allz(iti,ifi) =zprime
                       alld(iti,ifi) =dprime
                       alln(iti,ifi) =nprime
                       vold=vprime
                       kold=kprime
                       bold=bprime
                       dold=dprime
                       nold=nprime


                    enddo timeseries

               enddo  crosssection
          endif
         ! write(*,*) real(isdef/real(nyears))
          simdata(1,:,:) = allv(:,:)
          simdata(2,:,:) = allk(:,:)
          simdata(3,:,:) = alln(:,:)
          simdata(4,:,:) = allb(:,:)
          simdata(5,:,:) = allz(:,:)
          simdata(6,:,:) = alld(:,:)

          call labor(allk,alln,allz,params,alll)   !this computes the wage bill for a production
          simdata(7,:,:) = alll(:,:) ! Note that this is the total wage bill.

         if (erc .ne. 4) then
             if (compstat == 1) then 
                open (unit=262,file="Outputcomp/simdata.txt")
                do ifi=1,nfirms
                   do iti=burn_in,nyears
                        write(262,"(2i6,<numv>f18.3)") ifi,iti,simdata(:,iti,ifi)
                   enddo
                enddo
                close(262)
             endif
         else
             simdata = 0.0_dp
         endif

        end subroutine simpanel


     subroutine interpol(idz,kold,nold,bold,kg,ng,bg,v,vprime)
        use datatype
        use sizes
        implicit none
        integer,  intent(in) :: idz
        real(dp), intent(in) :: kold, nold, bold, kg(nk), ng(nn), bg(nb), v(nz,nk,nn,nb)
        real(dp), intent(out) :: vprime
        real(dp) :: kfrac, bfrac, nfrac, v1, v2, v11, v22
        integer :: idkhi, idklo, idnhi, idnlo, idbhi, idblo, enn


        enn = nk;  call pick(kold,kg,idklo,idkhi,enn,kfrac)
        enn = nn;  call pick(nold,ng,idnlo,idnhi,enn,nfrac)
        enn = nb;  call pick(bold,bg,idblo,idbhi,enn,bfrac)

    !    if (idklo == 0) then
    !       call wtfisgoingon(kold,nold,bold,kg,ng,bg,v)
    !    endif


        !==========================================================================

        v1 = bfrac*v(idz,idkhi,idnhi,idbhi) + (1.0_dp - bfrac)*v(idz,idkhi,idnhi,idblo)  !  ;write(*,"(4f12.3,2i2)") v1, v(idzhi,idkhi,idnhi,idblo), v(idzhi,idkhi,idnhi,idbhi), bfrac, idblo, idbhi
        v2 = bfrac*v(idz,idkhi,idnlo,idbhi) + (1.0_dp - bfrac)*v(idz,idkhi,idnlo,idblo)


        v11 = nfrac*v1 + (1.0_dp - nfrac)*v2


        v1 = bfrac*v(idz,idklo,idnhi,idbhi) + (1.0_dp - bfrac)*v(idz,idklo,idnhi,idblo)
        v2 = bfrac*v(idz,idklo,idnlo,idbhi) + (1.0_dp - bfrac)*v(idz,idklo,idnlo,idblo)

        v22 = nfrac*v1 + (1.0_dp - nfrac)*v2

        vprime = kfrac*v11 + (1.0_dp - kfrac)*v22


     end subroutine interpol


     subroutine pick(kold,kg,idklo,idkhi,enn,kfrac)
     use sizes
     use datatype
     use globals
     implicit none
     integer, intent(in)   :: enn
     real(dp), intent(in)  :: kold, kg(enn)
     real(dp), intent(out) :: kfrac
     integer,  intent(out) :: idklo, idkhi
     real(dp) ::  klo, khi
     real(dp), parameter :: eps = 1.0e-6
     integer :: ii

     idklo = 1
     do ii=enn,1,-1
       if (abs(kold-kg(ii))<eps) then
          idklo = ii
          exit
       endif
       if (kold > kg(ii)) then
          idklo = ii
          exit
       endif
     enddo
     if (idklo < enn) then
       idkhi = idklo + 1
       klo = kg(idklo)
       khi = kg(idkhi)
       kfrac = abs(kold-klo)/abs(khi-klo)
     else
        idkhi = enn
        kfrac = 1.0_dp
     endif
     end subroutine pick


  subroutine labor(k,n,z,params,l)   !this computes the wage bill for a production
     use datatype                    !function that is z^nu(k^alpha l^(1-alpha))^theta
     use sizes                       !nu = 1- (1-alpha)theta
     use globals
     implicit none
     real(dp), intent(in) :: k(nyears,nfirms), n(nyears,nfirms), z(nyears,nfirms), params(nop)
     real(dp), intent(out) :: l(nyears,nfirms)
     real(dp) :: aggk(nyears,nfirms)
     real(dp) :: theta, beta, wage, power1, power2, power3

     theta  = params(1)
     beta   = params(2)
     wage   = params(13)
     aggk   = (k**beta) * (n**(1.0_dp - beta))
     power1 = ((1.0_dp-share)*theta/wage)**(1.0_dp/(1.0_dp-(1.0_dp-share)*theta))
     power2 = (power1**((1.0_dp-share)*theta)-power1*wage)     !This is irrelevant here but whatever
     power3 = share*theta/(1.0_dp-(1.0_dp-share)*theta)

     l = power1*z*aggk**power3
     l = l*wage

  end subroutine labor
