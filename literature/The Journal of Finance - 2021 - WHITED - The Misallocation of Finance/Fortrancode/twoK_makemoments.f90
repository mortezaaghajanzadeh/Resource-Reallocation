subroutine makemoments(params,simmoms)
      use datatype
      use sizes
      use globals
      use pickmoments
      implicit none

      real(dp), intent(in) :: params(nop)
      real(dp), intent(inout) :: simmoms(nmom)

      real(dp), allocatable :: v(:,:), k(:,:), n(:,:), b(:,:), d(:,:), z(:,:), w(:,:)


      real(dp), allocatable :: invest(:,:), debt(:,:), div(:,:), eg(:,:), eish(:,:), &
                               opinc(:,:), opinclag(:,:), mb(:,:), output(:,:)
      real(dp) :: allsimmoms(numinfl)

      real(dp) :: theta, beta, rho, sigma, delta, inlta, adjK, adjN, lambda, xi, wage, dishcost, divsmooth
      real(dp) :: gamma, ra_gain, infi_gain, hk_gain, r2
      real(dp) :: aye, bee, const, alpha, skale, fixcost ! permutations
      real(dp) :: mu_i, v_i, mu_lev, v_lev, mu_div, v_div, mu_eg, v_eg, gamma_est, mu_output, v_output, &
                  mu_inc, v_inc, rho_inc, mu_mb, v_mb, mu_linc, v_linc, gain, sludge, moresludge

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

      !b = b + (1-delta*xi)*k  this is a test this is only a test


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
      allocate (debt(nyears-burn_in,nfirms))
      allocate (div(nyears-burn_in,nfirms))
      allocate (eg(nyears-burn_in,nfirms))
      allocate (eish(nyears-burn_in,nfirms))
      allocate (opinc(nyears-burn_in,nfirms))
      allocate (opinclag(nyears-burn_in,nfirms))
      allocate (mb(nyears-burn_in,nfirms))
      allocate (output(nyears-burn_in,nfirms))


      !=================================================================================================
      ! Make the constructed variables
      !=================================================================================================

      invest     =  (k(2:nyears-burn_in+1,:) - (1.0_dp - delta)*k(1:nyears-burn_in,:))/k(1:nyears-burn_in,:)
      debt       =  b(2:nyears-burn_in+1,:)/(k(2:nyears-burn_in+1,:)+n(2:nyears-burn_in+1,:))
      div        =  d(2:nyears-burn_in+1,:)/(k(2:nyears-burn_in+1,:)+n(2:nyears-burn_in+1,:))
      mb         =  (v(2:nyears-burn_in+1,:)+b(2:nyears-burn_in+1,:))/(k(2:nyears-burn_in+1,:)+n(2:nyears-burn_in+1,:))
      eg         =  (n(2:nyears-burn_in+1,:) - (1.0_dp - inlta)*n(1:nyears-burn_in,:))/n(1:nyears-burn_in,:)
      opinc      =  (skale*z(2:nyears-burn_in+1,:)*(k(2:nyears-burn_in+1,:)**aye)*(n(2:nyears-burn_in+1,:)**bee) - fixcost*skale)/(k(2:nyears-burn_in+1,:)+n(2:nyears-burn_in+1,:))
      opinclag   =  (skale*z(1:nyears-burn_in+0,:)*(k(1:nyears-burn_in+0,:)**aye)*(n(1:nyears-burn_in+0,:)**bee) - fixcost*skale)/(k(1:nyears-burn_in,:)+n(1:nyears-burn_in,:))
      output     =  (skale*z(2:nyears-burn_in+1,:)*(k(2:nyears-burn_in+1,:)**aye)*(n(2:nyears-burn_in+1,:)**bee)+w(2:nyears-burn_in+1,:) - fixcost*skale)/(k(2:nyears-burn_in+1,:)+n(2:nyears-burn_in+1,:))

      eish = 0.0_dp
      where (div < 0.0_dp)
         eish = -div
      end where

      where (div < 0.0_dp)
         div = 0.0_dp
      end where

      if (compstat==1) then
          open (unit=262,file="Outputcomp/simvariables.txt")
          do ifi=1,1
             do iti=1,nyears-burn_in
                  write(262,"(12f15.5)") invest(iti,ifi),eg(iti,ifi),div(iti,ifi),eish(iti,ifi),debt(iti,ifi),mb(iti,ifi),opinc(iti,ifi)
             enddo
          enddo
          close(262)
      endif

      nobs = real(size(invest))
      !=================================================================================================
      ! Make the actual moments
      !=================================================================================================

      mu_i        = sum(invest)/nobs
      v_i         = sum(invest**2.0_dp)/nobs - mu_i**2.0_dp

      mu_lev      = sum(debt)/nobs
      v_lev       = sum(debt**2.0_dp)/nobs - mu_lev**2.0_dp

      mu_div      = sum(div)/nobs
      v_div       = sum(div**2.0_dp)/nobs - mu_div**2.0_dp

      mu_eg       = sum(eg)/nobs
      v_eg        = sum(eg**2.0_dp)/nobs - mu_eg**2.0_dp

      mu_inc      = sum(opinc)/nobs
      v_inc       = sum(opinc**2.0_dp)/nobs - mu_inc**2.0_dp

      mu_mb       = sum(mb)/nobs
      v_mb        = sum(mb**2.0_dp)/nobs - mu_mb**2.0_dp

      mu_output   = sum(output)/nobs
      v_output    = sum(log(output)**2.0_dp)/nobs -(sum(log(output))/nobs)**2.0_dp

      mu_linc     = sum(opinclag)/nobs
      v_linc      = sum(opinclag**2.0_dp)/nobs - mu_linc**2.0_dp
      rho_inc     = sum( (opinc-mu_inc)*(opinclag-mu_linc))/nobs/(sqrt(v_inc*v_linc))

      !==============This makes the gain======================

      deallocate (invest)
      deallocate (debt)
      deallocate (div)
      deallocate (eg)
      deallocate (eish)
      deallocate (opinc)
      deallocate (opinclag)
      deallocate (mb)
      deallocate (v)


      allocate (invest(nyears-burn_in+1,nfirms))
      allocate (eg(nyears-burn_in+1,nfirms))
      invest = simdata(2,burn_in:nyears,:) - (1.0_dp - delta)*simdata(2,burn_in-1:nyears-1,:)
      eg     = simdata(2,burn_in:nyears,:)*simdata(3,burn_in:nyears,:) - (1.0_dp - delta)*simdata(2,burn_in-1:nyears-1,:)*simdata(3,burn_in-1:nyears-1,:)


      gamma    = 1.519_dp
      call wz(k,n,b,w,d,z,gamma,invest,eg,params,ra_gain,gamma_est,r2)

      gamma    = 1000.0_dp
      call wz(k,n,b,w,d,z,gamma,invest,eg,params,infi_gain,sludge,moresludge)

      gamma = 1.519_dp
      call hk(k,n,b,w,z,invest,eg,params,hk_gain)
      !call hk3(k,n,b,w,z,params,gamma,hk_gain)

      allsimmoms(1) = mu_i         !  x  'mu_invest          ', &
      allsimmoms(2) = mu_eg        !  x  'mu_intang          ', &
      allsimmoms(3) = mu_lev       !  x  'mu_lev             ', &
      allsimmoms(4) = mu_inc       !  x  'mu_opinc           ', &
      allsimmoms(5) = mu_div       !  x  'mu_div             ', &
      allsimmoms(6) = mu_mb        !  x  'mu_mb              ', &
      allsimmoms(7) = v_i          !  x  'v_invest           ', &
      allsimmoms(8) = v_eg         !  x  'v_intang           ', &
      allsimmoms(9) = v_lev        !  x  'v_lev              ', &
      allsimmoms(10) = v_inc       !  x  'v_opinc            ', &
      allsimmoms(11) = v_div       !  x  'v_div              ', &
      allsimmoms(12) = v_mb        !  x  'v_mb               ', &
      allsimmoms(13) = rho_inc     !  x  'rho_opinc          ', &
      allsimmoms(14) = ra_gain     !  x  'gain               ', &
      allsimmoms(15) = infi_gain   !  x  'type gain          ', &
      allsimmoms(16) = gamma_est   !  x  'eslaticity of sub  ', &
      allsimmoms(17) = r2          !  x  'eslaticity of sub  ', &
      allsimmoms(18) = hk_gain     !  x  'eslaticity of sub  ', &
      allsimmoms(19) = aggY        !  x  'eslaticity of sub  ', &
      allsimmoms(20) = aggC        !  x  'eslaticity of sub  ', &
      allsimmoms(21) = wage
      allsimmoms(22) = mu_output
      allsimmoms(23) = v_output
      allsimmoms(24) = TFP
      simmoms = allsimmoms(pickout)

end subroutine makemoments


subroutine wz(k,n,b,w,e,z,gamma,invest,eg,params,ra_gain,gamma_est,r2)
      use datatype
      use sizes
      use globals
      use pickmoments
      implicit none

      real(dp), intent(in) :: k(nyears-burn_in+1,nfirms), n(nyears-burn_in+1,nfirms), z(nyears-burn_in+1,nfirms),   &
                              b(nyears-burn_in+1,nfirms), w(nyears-burn_in+1,nfirms), e(nyears-burn_in+1,nfirms),   &
                              invest(nyears-burn_in+1,nfirms), eg(nyears-burn_in+1,nfirms),                       &
                              params(nop), gamma

      real(dp), intent(out) :: ra_gain, gamma_est, r2

      real(dp) :: theta, beta, rho, sigma, delta, inlta, adjK, adjN, lambda, xi, wage, dishcost, divsmooth                                     ! parameters to estimate
      real(dp) :: cygma, r_s, lambda_s
      real(dp) :: aye, bee, const, alpha, skale, fixcost ! permutations

      real(dp) :: D_s, E_s, alpha_s, A_s_temp, F_s, F_s_eff, PF, theta_s, F_temp, F_eff_temp, F, F_eff


      real(dp), allocatable :: E_si(:,:), D_si(:,:), PF_si(:,:), A_si(:,:), &
                               D_si_eff(:,:), E_si_eff(:,:), F_si(:,:), F_si_eff(:,:), &
                               r_si(:,:), lambda_si(:,:), adj_si(:,:)

      real(dp), allocatable :: lD(:), lE(:), lF(:), lDE2(:), lhs(:), rhs(:,:), bhat(:)   ! stuff for the Kmenta regression

      integer :: len, nreg, panlen, ii


      panlen = nyears-burn_in+1
      allocate(E_si(panlen,nfirms))
      allocate(D_si(panlen,nfirms))
      allocate(PF_si(panlen,nfirms))
      allocate(A_si(panlen,nfirms))
      allocate(D_si_eff(panlen,nfirms))
      allocate(E_si_eff(panlen,nfirms))
      allocate(F_si(panlen,nfirms))
      allocate(F_si_eff(panlen,nfirms))
      allocate(r_si(panlen,nfirms))
      allocate(lambda_si(panlen,nfirms))
      allocate(adj_si(panlen,nfirms))

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


      cygma    = 1.77_dp
      r_s      = 0.03_dp
      lambda_s = 0.03_dp

      D_si = b
      E_si = k - b

      where (D_si < 0.0_dp)
         D_si = 0.0_dp
      end where

      PF_si      = (z**(1.0_dp-theta+theta*share)) *(((k**beta)*(n**(1.0_dp-beta)))**share * (w/wage)**(1.0_dp-share))**theta - fixcost*skale
      PF_si      = skale*z*(k**bee)*(n**aye) + w - fixcost*skale
      adj_si = 0.5_dp*adjK*invest**2.0_dp/k + 0.5_dp*adjN*eg**2.0_dp/n
      PF_si  = PF_si - adj_si


      D_s = sum(D_si)
      E_s = sum(E_si)


      alpha_s    = r_s*D_s**(1.0_dp/gamma)/(r_s*D_s**(1.0_dp/gamma)+lambda_s*E_s**(1.0_dp/gamma))  !r and lambda are from a previous version of the code, but they cancel out anyway.

      A_si       = PF_si**(cygma/(cygma - 1.0_dp))/(alpha_s*D_si**((gamma - 1.0_dp)/gamma)+(1.0_dp - alpha_s)*E_si**((gamma - 1.0_dp)/gamma))**(gamma/(gamma - 1.0_dp))

      A_s_temp   = sum(A_si**(cygma - 1.0_dp))
      D_si_eff   = D_s*A_si**(cygma - 1.0_dp)/A_s_temp
      E_si_eff   = E_s*A_si**(cygma - 1.0_dp)/A_s_temp

      F_si       = A_si*(alpha_s*D_si**((gamma - 1.0_dp)/gamma)+(1.0_dp - alpha_s)*E_si**((gamma - 1.0_dp)/gamma))**(gamma/(gamma - 1.0_dp))
      F_si_eff   = A_si*(alpha_s*D_si_eff**((gamma - 1.0_dp)/gamma)+(1.0_dp - alpha_s)*E_si_eff**((gamma - 1.0_dp)/gamma))**(gamma/(gamma - 1.0_dp))

      F_s        = sum(F_si**((cygma - 1.0_dp)/cygma))
      F_s        = F_s**(cygma/(cygma - 1.0_dp))
      F_s_eff    = sum(F_si_eff**((cygma - 1.0_dp)/cygma))
      F_s_eff    = F_s_eff**(cygma/(cygma - 1.0_dp))

      PF         = sum(PF_si)
      theta_s    = sum(PF_si)
      theta_s    = theta_s/PF

      F_temp     = F_s**theta_s
      F_eff_temp = F_s_eff**theta_s

      r_si       = alpha_s*(cygma - 1.0_dp)/cygma*PF_si/(alpha_s*D_si+(1.0_dp - alpha_s)*E_si**((gamma - 1.0_dp)/gamma)*D_si**(1.0_dp/gamma))
      lambda_si  = (1.0_dp - alpha_s)*(cygma - 1.0_dp)/cygma*PF_si/(alpha_s*D_si**((gamma - 1.0_dp)/gamma)*E_si**(1.0_dp/gamma)+(1.0_dp - alpha_s)*E_si)

      F          = (log(F_temp))
      F_eff      = (log(F_eff_temp))

      F          = exp(F)
      F_eff      = exp(F_eff)

      ra_gain    = F_s/F_s_eff

      !===========================This part does the KMENTA regression=======================================


      deallocate(F_si_eff)
      deallocate(D_si_eff)
      deallocate(E_si_eff)
      deallocate(A_si)
      deallocate(r_si)
      deallocate(lambda_si)

      panlen = (nyears-burn_in+1)*nfirms
      len = panlen
      nreg = 3
      allocate(lD(panlen))
      allocate(lE(panlen))
      allocate(lF(panlen))
      allocate(lDE2(panlen))
      lD = (reshape(D_si, (/ len /) ))
      lE = (reshape(E_si, (/ len /) ))
      lF = (reshape(PF_si, (/ len /) ))
      deallocate(F_si)
      deallocate(D_si)
      deallocate(E_si)
      !lD = asinh(lD)
      !lE = asinh(lE)
      !lF = asinh(lF)

      !lD = log(1.0_dp+lD)
      !lE = log(1.0_dp+lE)
      !lF = log(1.0_dp+lF)
      !lDE2 = (lD - lE)**2.0_dp
      !allocate(lhs(panlen))
      !allocate(rhs(panlen,nreg))
      !lhs = lF
      !rhs(:,1) = lD
      !rhs(:,2) = lE
      !rhs(:,3) = lDE2

      !deallocate(lD)
      !deallocate(lE)
      !deallocate(lF)
      !deallocate(lDE2)
      !
      !allocate(bhat(nreg+1))
      !call dools(len,nreg,lhs,rhs,bhat)
      !gamma_est = 1.0_dp + 2.0_dp*bhat(4)/(bhat(2)*bhat(3))
      len = count(lD > 0.0_dp)
      !lD = log(lD) This is the screening variable. Don't take the log
      lE = log(lE)
      lF = log(lF)
      lDE2 = (log(lD) - lE)**2.0_dp
      if (len > (nreg+1)) then
          allocate(lhs(len))
          allocate(rhs(len,nreg))
          rhs = 0.0_dp
          lhs = 0.0_dp
          lhs      = pack(lF,mask=lD>0.0_dp)
          rhs(:,1) = pack(log(lD),mask=lD>0.0_dp)
          rhs(:,2) = pack(lE,mask=lD>0.0_dp)
          rhs(:,3) = pack(lDE2,mask=lD>0.0_dp)

          !do ii=1,len;   write(20,"(4f20.12)") lhs(ii),rhs(ii,:);enddo;
          deallocate(lD)
          deallocate(lE)
          deallocate(lF)
          deallocate(lDE2)

          allocate(bhat(nreg+1))
          call dools(len,nreg,lhs,rhs,bhat,r2)
          gamma_est = 1.0_dp + 2.0_dp*bhat(4)/(bhat(2)*bhat(3))
     else
          gamma_est = -100.0_dp
     endif

end subroutine wz

subroutine hk(k,n,b,w,z,invest,eg,params,ra_gain)
      use datatype
      use sizes
      use globals
      use pickmoments
      implicit none

      real(dp), intent(in) :: k(nyears-burn_in+1,nfirms),n(nyears-burn_in+1,nfirms),z(nyears-burn_in+1,nfirms),   &
                              invest(nyears-burn_in+1,nfirms),eg(nyears-burn_in+1,nfirms),                      &
                              b(nyears-burn_in+1,nfirms),w(nyears-burn_in+1,nfirms),params(nop)
      real(dp), intent(out) :: ra_gain

      real(dp) :: theta, beta, rho, sigma, delta, inlta, adjK, adjN, lambda, xi, wage, dishcost, divsmooth                                     ! parameters to estimate
      real(dp) :: cygma, i_s, w_s
      real(dp) :: aye, bee, const, alpha, skale, fixcost ! permutations

      real(dp) :: L_s, K_s, alpha_s, A_s_temp, Y_s, Y_s_eff, &
                  PY, theta_s, Y_temp, Y_eff_temp, Y, Y_eff


      real(dp), allocatable :: K_si(:,:), L_si(:,:), PY_si(:,:), A_si(:,:), &
                               L_si_eff(:,:), K_si_eff(:,:), Y_si(:,:), Y_si_eff(:,:), &
                               i_si(:,:), w_si(:,:), adj_si(:,:)


      integer :: len, nreg, panlen


      panlen = nyears-burn_in+1
      allocate(K_si(panlen,nfirms))
      allocate(L_si(panlen,nfirms))
      allocate(PY_si(panlen,nfirms))
      allocate(A_si(panlen,nfirms))
      allocate(K_si_eff(panlen,nfirms))
      allocate(L_si_eff(panlen,nfirms))
      allocate(Y_si(panlen,nfirms))
      allocate(Y_si_eff(panlen,nfirms))
      allocate(i_si(panlen,nfirms))
      allocate(w_si(panlen,nfirms))
      allocate(adj_si(panlen,nfirms))

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
      aye = alpha*(1.0_dp - beta)
      bee = alpha*beta


      cygma = 1.77
      i_s = 1.0!0.1146   !used to be 5
      w_s = 1.0!0.4034

      K_si   = k!(k**beta)*(n**(1.0_dp - beta))        n is unobservable
      PY_si  = (z**(1.0_dp-theta+theta*share)) *(((k**beta)*(n**(1.0_dp-beta)))**share * (w/wage)**(1.0_dp-share))**theta - fixcost*skale
      PY_si  = skale*z*(k**bee)*(n**aye) + w - fixcost*skale
      L_si   = w/wage

      adj_si = 0.5_dp*adjK*invest**2.0_dp/k + 0.5_dp*adjN*eg**2.0_dp/n
      PY_si  = PY_si - adj_si



      K_s = sum(K_si)
      L_s = sum(L_si)
      alpha_s = i_s*K_s/(i_s*K_s+w_s*L_s)

      A_si = PY_si**(cygma/(cygma-1.0_dp))/(K_si**alpha_s*L_si**(1.0_dp-alpha_s))

      A_s_temp = sum(A_si**(cygma-1.0_dp))
      K_si_eff = K_s*A_si**(cygma-1.0_dp)/A_s_temp
      L_si_eff = L_s*A_si**(cygma-1.0_dp)/A_s_temp

      Y_si = A_si*K_si**alpha_s*L_si**(1.0_dp-alpha_s)
      Y_si_eff = A_si*K_si_eff**alpha_s*L_si_eff**(1.0_dp-alpha_s)

      Y_s = sum(Y_si**((cygma-1.0_dp)/cygma))
      Y_s = Y_s**(cygma/(cygma-1.0_dp))
      Y_s_eff = sum(Y_si_eff**((cygma-1.0_dp)/cygma))
      Y_s_eff = Y_s_eff**(cygma/(cygma-1.0_dp))

      PY = sum(PY_si)
      theta_s = sum(PY_si)
      theta_s = theta_s/PY

      Y_temp = Y_s**theta_s
      Y_eff_temp = Y_s_eff**theta_s

      i_si = alpha_s*(cygma-1.0_dp)/cygma*PY_si/K_si
      w_si = (1.0_dp-alpha_s)*(cygma-1.0_dp)/cygma*PY_si/L_si

      Y = log(Y_temp)
      Y_eff = log(Y_eff_temp)

      Y = exp(Y)
      Y_eff = exp(Y_eff)

      ra_gain = Y/Y_eff

end subroutine hk



subroutine hk3(k,n,b,w,z,params,gamma,ra_gain)
      use datatype
      use sizes
      use globals
      use pickmoments
      implicit none

      real(dp), intent(in) :: k(nyears-burn_in+1,nfirms),n(nyears-burn_in+1,nfirms),z(nyears-burn_in+1,nfirms),   &
                              b(nyears-burn_in+1,nfirms),w(nyears-burn_in+1,nfirms),params(nop)
      real(dp), intent(in)  :: gamma
      real(dp), intent(out) :: ra_gain

      real(dp) :: theta, beta, rho, sigma, delta, inlta, adjK, adjN, lambda, xi, wage, dishcost, divsmooth                                     ! parameters to estimate
      real(dp) :: cygma, i_s, r_s, w_s
      real(dp) :: aye, bee, const, alpha,skale ! permutations

      real(dp) :: L_s, K_s, N_s, alpha_s, beta_s, A_s_temp, Y_s, Y_s_eff, &
                  PY, theta_s, Y_temp, Y_eff_temp, Y, Y_eff


      real(dp), allocatable :: K_si(:,:), N_si(:,:), L_si(:,:), PY_si(:,:), A_si(:,:), &
                               L_si_eff(:,:), K_si_eff(:,:), N_si_eff(:,:), Y_si(:,:), Y_si_eff(:,:), &
                               i_si(:,:), w_si(:,:), r_si(:,:)


      integer :: len, nreg, panlen


      panlen = nyears-burn_in+1
      allocate(K_si(panlen,nfirms))
      allocate(N_si(panlen,nfirms))
      allocate(L_si(panlen,nfirms))
      allocate(PY_si(panlen,nfirms))
      allocate(A_si(panlen,nfirms))
      allocate(K_si_eff(panlen,nfirms))
      allocate(L_si_eff(panlen,nfirms))
      allocate(N_si_eff(panlen,nfirms))
      allocate(Y_si(panlen,nfirms))
      allocate(Y_si_eff(panlen,nfirms))
      allocate(i_si(panlen,nfirms))
      allocate(r_si(panlen,nfirms))
      allocate(w_si(panlen,nfirms))

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

      const  = ((1.0_dp-share)*theta/wage)**(1.0_dp/(1.0_dp-(1.0_dp-share)*theta))
      alpha  = share*theta/(1.0_dp-(1.0_dp-share)*theta)
      skale  = const**((1.0_dp-share)*theta) - const*wage
      aye = alpha*(1.0_dp - beta)
      bee = alpha*beta


      cygma = 1.77
      i_s = 1.0! 0.1146   !used to be 5     !1.0!
      r_s = 1.0! 0.4034                     !1.0!
      w_s = 1.0!0.01497                     !1.0!

      K_si   = k
      N_si   = n
      PY_si  = skale*z*(k**bee)*(n**aye) + w
      !PY_si  = (z**(1.0_dp-theta+theta*share)) *(((k**beta)*(n**(1.0_dp-beta)))**share * (w/wage)**(1.0_dp-share))**theta! the same
      L_si   = w/wage


      K_s = sum(K_si)
      N_s = sum(N_si)
      L_s = sum(L_si)

      alpha_s = i_s*K_s**(1.0_dp/gamma)/(i_s*K_s**(1.0_dp/gamma)+r_s*N_s**(1.0_dp/gamma))
      beta_s = (i_s*K_s+r_s*N_s)/(i_s*K_s+r_s*N_s+w_s*L_s)

      A_si = PY_si**(cygma/(cygma-1.0_dp))/((alpha_s*K_si**((gamma-1.0_dp)/gamma)+(1.0_dp-alpha_s)*N_si**((gamma-1.0_dp)/gamma))**(beta_s*gamma/(gamma-1.0_dp))*L_si**(1.0_dp-beta_s))

      A_s_temp = sum(A_si**(cygma-1.0_dp))
      K_si_eff = K_s*A_si**(cygma-1.0_dp)/A_s_temp
      N_si_eff = N_s*A_si**(cygma-1.0_dp)/A_s_temp
      L_si_eff = L_s*A_si**(cygma-1.0_dp)/A_s_temp

      Y_si = A_si*(alpha_s*K_si**((gamma-1.0_dp)/gamma)+(1.0_dp-alpha_s)*N_si**((gamma-1.0_dp)/gamma))**(beta_s*gamma/(gamma-1.0_dp))*L_si**(1.0_dp-beta_s)
      Y_si_eff = A_si*(alpha_s*K_si_eff**((gamma-1.0_dp)/gamma)+(1.0_dp-alpha_s)*N_si_eff**((gamma-1.0_dp)/gamma))**(beta_s*gamma/(gamma-1.0_dp))*L_si_eff**(1.0_dp-beta_s)

      Y_s = sum(Y_si**((cygma-1.0_dp)/cygma))
      Y_s = Y_s**(cygma/(cygma-1.0_dp))
      Y_s_eff = sum(Y_si_eff**((cygma-1.0_dp)/cygma))
      Y_s_eff = Y_s_eff**(cygma/(cygma-1.0_dp))

      PY = sum(PY_si)
      theta_s = sum(PY_si)
      theta_s = theta_s/PY

      Y_temp = Y_s**theta_s
      Y_eff_temp = Y_s_eff**theta_s

      i_si = alpha_s*beta_s*(sigma-1.0_dp)/sigma*PY_si/(alpha_s*K_si+(1.0_dp-alpha_s)*N_si**((gamma-1.0_dp)/gamma)*K_si**(1.0_dp/gamma))
      r_si = (1.0_dp-alpha_s)*beta_s*(sigma-1.0_dp)/sigma*PY_si/(alpha_s*K_si**((gamma-1.0_dp)/gamma)*N_si**(1.0_dp/gamma)+(1.0_dp-alpha_s)*N_si)
      w_si =(1.0_dp-beta_s)*(sigma-1.0_dp)/sigma*PY_si/L_si


      Y = log(Y_temp)
      Y_eff = log(Y_eff_temp)

      Y = exp(Y)
      Y_eff = exp(Y_eff)

      ra_gain = Y/Y_eff

end subroutine hk3
