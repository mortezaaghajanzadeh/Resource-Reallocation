subroutine aggregates(params,mew,pol,sgrid)
      use datatype
      use sizes
      use globals
      implicit none

      type(grids), intent(in), dimension(gmax) :: sgrid
      type(policies), intent(in), dimension(nz,nk,nn,nb) :: pol
      real(dp), intent(in) :: params(nop)
      real(dp), intent(in) :: mew(nz,nk,nn,nb)

      real(dp) :: theta, beta, rho, sigma, delta, inlta, adjK, adjN, xi, &
                  lambda, wage, dishcost, divsmooth     ! parameters to estimate

      real(dp) :: aye, bee, const, alpha, skale, fixcost
      real(dp) :: adjcosts(nz,nk,nn,nb)
      integer :: ik, in

      !=================================================================================================
      ! Extract the parameters
      !=================================================================================================
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


      const  = ((1.0_dp-share)*theta/wage)**(1.0_dp/(1.0_dp-(1.0_dp-share)*theta))
      skale  = const**((1.0_dp-share)*theta) - const*wage

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(in,ik)
      !$OMP DO COLLAPSE(2)
      do ik = 1,nk
         do in = 1,nn
           adjcosts(:,ik,in,:) = 0.5_dp*adjK*pol(:,ik,in,:).i**2.0_dp/sgrid(ik).k   +  0.5_dp*adjN*pol(:,ik,in,:).r**2.0_dp/(sgrid(ik).k*sgrid(in).n)
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      !=================================================================================================
      ! Make the actual moments
      !=================================================================================================

      aggY        = sum(pol.y*mew)  - skale*fixcost - sum(adjcosts*mew)
      aggI        = sum(pol.i*mew)
      aggL        = sum(pol.l*mew)
      aggN        = sum((pol.n*pol.k)*mew)
      aggKay      = sum(pol.k*mew)
      aggR        = sum(pol.r*mew)

      TFP         = log(aggY) - beta*share*theta*log(aggKay) - (1.0_dp-beta)*share*theta*log(aggN) - (1.0_dp-share)*theta*log(aggL)
      TFP         = exp(TFP)

end subroutine aggregates
