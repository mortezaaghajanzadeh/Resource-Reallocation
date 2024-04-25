subroutine makeaggregates(params)
      use datatype
      use sizes
      use globals
      implicit none

      real(dp), intent(in) :: params(nop)

      real(dp), allocatable :: v(:,:), k(:,:), n(:,:), b(:,:), d(:,:), z(:,:), w(:,:)

      real(dp), allocatable :: invest(:,:), sales(:,:), labor(:,:), intan(:,:), divs(:,:), adjcosts(:,:), glue_k(:,:)

      real(dp) :: theta, beta, rho, sigma, delta, inlta, adjK, adjN, xi, lambda, wage, dishcost, divsmooth     ! parameters to estimate

      real(dp) :: aye, bee, const, alpha, skale, fixcost

      integer :: nobs, ifi, iti
      !=================================================================================================
      ! Allocate the raw variable arrays. Note they are + 1 for lagging purposes.
      !=================================================================================================
      allocate (v(nyears-burn_in+1,nfirms))
      allocate (k(nyears-burn_in+1,nfirms))
      allocate (n(nyears-burn_in+1,nfirms))
      allocate (b(nyears-burn_in+1,nfirms))
      allocate (z(nyears-burn_in+1,nfirms))
      allocate (d(nyears-burn_in+1,nfirms))
      allocate (w(nyears-burn_in+1,nfirms))

      !=================================================================================================
      ! Extract the raw variables
      !=================================================================================================
      v(:,:) = simdata(1,burn_in:nyears,:)
      k(:,:) = simdata(2,burn_in:nyears,:)
      n(:,:) = simdata(3,burn_in:nyears,:)
      b(:,:) = simdata(4,burn_in:nyears,:)
      z(:,:) = simdata(5,burn_in:nyears,:)
      d(:,:) = simdata(6,burn_in:nyears,:)
      w(:,:) = simdata(7,burn_in:nyears,:)

      n = n*k
      b = b*k

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
      alpha  = share*theta/(1.0_dp-(1.0_dp-share)*theta)
      skale  = const**((1.0_dp-share)*theta) - const*wage
      aye    = alpha*(1.0_dp - beta)
      bee    = alpha*beta

      !=================================================================================================
      ! Allocate the constructed variables
      !=================================================================================================
      allocate (invest(nyears-burn_in,nfirms))
      allocate (sales(nyears-burn_in,nfirms))
      allocate (labor(nyears-burn_in,nfirms))
      allocate (intan(nyears-burn_in,nfirms))
      allocate (glue_k(nyears-burn_in,nfirms))
      allocate (adjcosts(nyears-burn_in,nfirms))


      !=================================================================================================
      ! Make the constructed variables
      !=================================================================================================

      invest = (k(2:nyears-burn_in+1,:) - (1.0_dp - delta)*k(1:nyears-burn_in,:))
      intan  = (n(2:nyears-burn_in+1,:) - (1.0_dp - inlta)*n(1:nyears-burn_in,:))
      labor  = w(2:nyears-burn_in+1,:)/wage

      adjcosts = 0.5_dp*adjK*invest**2.0_dp/k(1:nyears-burn_in,:) + 0.5_dp*adjN*intan**2.0_dp/n(1:nyears-burn_in,:)
      !k = (k**beta) * (n**(1.0_dp-beta))
      glue_k = (k(2:nyears-burn_in+1,:)**(1.0_dp-beta))*(n(2:nyears-burn_in+1,:)**beta)
      sales  = (z(2:nyears-burn_in+1,:)**(1.0_dp-theta+theta*share))*((glue_k**share)*(labor**(1.0_dp-share)))**theta - fixcost*skale
      !write(*,*) sum(sales)/real(size(sales))
      !sales  = (skale*z(2:nyears-burn_in+1,:))*(k(2:nyears-burn_in+1,:)**aye)*((n(2:nyears-burn_in+1,:))**bee) + w - fixcost*skale
      !write(*,*) sum(sales)/real(size(sales))
      sales = sales - adjcosts


      !=================================================================================================
      ! Make the actual moments
      !=================================================================================================

      nobs        = size(sales)
      aggY        = sum(sales)/nobs
      aggI        = sum(invest)/nobs
      aggL        = sum(labor)/nobs
      aggR        = sum(intan)/nobs
      aggKay      = sum(k(2:nyears-burn_in+1,:))/nobs
      aggN        = sum(n(2:nyears-burn_in+1,:))/nobs
      sales       = log(sales) - beta*share*theta*log(k(2:nyears-burn_in+1,:)) - (1.0_dp-beta)*share*theta*log(n(2:nyears-burn_in+1,:)) - (1.0_dp-share)*theta*log(labor)
      TFP         = sum(exp(sales))/nobs

      TFP         = log(aggY) - beta*share*theta*log(aggKay) - (1.0_dp-beta)*share*theta*log(aggN) - (1.0_dp-share)*theta*log(aggL)
      TFP         = exp(TFP)
      write(*,*) "simulated agg k and n", aggKay, aggN
end subroutine makeaggregates
