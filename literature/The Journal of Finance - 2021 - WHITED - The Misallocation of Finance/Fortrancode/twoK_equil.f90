module packmeup
use datatype
use sizes
use globals
implicit none


contains



real(dp) function gap(pee,vg)
    use datatype
    use sizes
    use globals
    implicit none
    real(dp), intent(in) :: pee(nop)
    real(dp), intent(in) :: vg(nz,nk,nn,nb)
    real(dp) :: phi, gamma
    phi = 2.5_dp
    gamma = 1.0_dp

    aggC = (pee(13)/phi)**(1.0_dp/gamma)

    call valfun(pee,vg)
    if (mewave == 0) call makeaggregates(pee)
    gap = aggY - aggI - aggC - aggR

    write(*,*) "output     ",aggY
    write(*,*) "investment ",aggI
    write(*,*) "consumption",aggC
    write(*,*) "intangibles",aggR
    write(*,*) "wage       ",pee(13)

end function gap



end module packmeup
subroutine momentgen(params,momvector,vg)
      use datatype
      use sizes
      use globals
      use omp_lib
      use packmeup
      implicit none

      real(dp), intent(in) :: params(nop)
      real(dp), intent(out) :: momvector(nmom)
      real(dp), intent(inout) :: vg(nz,nk,nn,nb)

      real(dp) :: pee(nop), peetemp(nop)          ! expanded parameter vector
      real(dp) :: oldgap, newgap, wage, w_old, w_new, f_old, f_new, testgap, huntdirection=-1.0_dp
      real(dp), parameter :: gaptol=0.01_dp, huntincrem=0.2_dp
      integer  :: ip, ii, is, maxiter
      character (len = 50) :: readfil
      logical :: c1, c2, c3, c4, disaster
      maxiter = 40                                     ! number of iterations to calculate equilibrium

      allocate(simdata(numv,nYears,nFirms))

      do ip = 1,nop
         pee(ip) = params(ip)
      enddo

      oldgap = gap(pee,vg)
      if     (abs(oldgap/aggY) <= gaptol) then                               !If all is well, make moments and exit
         write(*,"('consumption, investment, output, labor, intang, wage', 2x, 6f12.6)") aggC, aggI, aggY, aggL, aggR, wage
         call makemoments(pee,momvector)
      elseif (abs(oldgap/aggY) > gaptol) then                                !Only keep going if there is no convergence
         peetemp = pee 
         peetemp(13) = pee(13)*(1.0_dp + huntincrem*huntdirection)           !Try a new point
         newgap = gap(peetemp,vg)
         write(*,*) "---------------------------------"
         write(*,*) "new gap = ", newgap
         write(*,*) "old gap = ", oldgap
         write(*,*) "---------------------------------"


         if (abs(newgap/aggY) <= gaptol) then                             !If all is well, make moments and exit
            write(*,"('consumption, investment, output, labor, intang, wage', 2x, 6f12.6)") aggC, aggI, aggY, aggL, aggR, wage
            call makemoments(pee,momvector)
         elseif (abs(newgap/aggY) > gaptol) then                          !Only keep going if there is no convergence

            c1 = (newgap>0.0_dp) .and. (oldgap<0.0_dp)                    !Different signs?
            c2 = (newgap<0.0_dp) .and. (oldgap>0.0_dp)                    !Different signs?

            c3 = (newgap>0.0_dp) .and. (oldgap>0.0_dp)                    !Same signs?
            c4 = (newgap<0.0_dp) .and. (oldgap<0.0_dp)                    !Same signs?


            if (c1 .or. c2) then


               w_old = pee(13)
               w_new = peetemp(13)
               f_old = oldgap
               f_new = newgap
               wage = 0.5_dp*(w_old + w_new)

               call bisectme(pee,f_old,f_new,w_old,w_new,wage,gaptol,maxiter,momvector,vg)


            elseif (c3 .or. c4) then                                     !If the same sign, find the hunting direction
               if (abs(oldgap)<abs(newgap)) then                     !If you move further away, change the sign
                   huntdirection = -1.0_dp*huntdirection
               endif

               peetemp = pee
               huntloop: do ii=1,maxiter

                  is = floor(dble(ii+1)/2.0_dp)
                  peetemp(13) = pee(13)*(1.0_dp + dble(is)*huntdirection*huntincrem)

                  newgap = gap(peetemp,vg)

                  disaster = wage < 0.0_dp .or. abs(newgap)>20.0_dp .or. abs(oldgap)>20.0_dp 
                  if (disaster) exit huntloop

                  if (verbose == 1) then 
                     write(*,*) "---------------------------------"
                     write(*,*) "new gap = ", newgap
                     write(*,*) "old gap = ", oldgap
                     write(*,*) "---------------------------------"
                  endif

                  c1 = (newgap>0.0_dp) .and. (oldgap<0.0_dp)
                  c2 = (newgap<0.0_dp) .and. (oldgap>0.0_dp)
                  if (verbose == 1) write(*,*) ii, oldgap,newgap
                  if (c1 .or. c2) exit huntloop
               enddo huntloop

               w_old = pee(13)
               w_new = peetemp(13)
               f_old = oldgap
               f_new = newgap
               wage = 0.5_dp*(w_old + w_new)

               if (disaster) then 
                  momvector = 100.0_dp 
               else
                  call bisectme(pee,f_old,f_new,w_old,w_new,wage,gaptol,maxiter,momvector,vg)
               endif 
           endif
         endif
      endif
      deallocate(simdata)


end subroutine momentgen

subroutine bisectme(pee,f_old,f_new,w_old,w_new,wage,gaptol,maxiter,momvector,vg)
      use datatype
      use sizes
      use globals
      use omp_lib
      use packmeup
      implicit none

      real(dp), intent(inout) :: pee(nop), f_old, f_new, w_old, w_new, wage
      real(dp), intent(in) :: gaptol
      integer, intent(in) :: maxiter
      real(dp), intent(out) :: momvector(nmom)
      real(dp), intent(inout) :: vg(nz,nk,nn,nb)
      logical :: c1, c2
      integer :: ii
      real(dp) :: testgap, peetemp(nop) 
      if (verbose == 1)  write(*,*) "starting bisection"
      bisectloop: do ii=1,maxiter
         peetemp = pee 
         peetemp(13) = wage
         testgap = gap(peetemp,vg)
         if (verbose == 1) write(*,*) "wage endpoints ", w_old, w_new
         if (verbose == 1) write(*,*) "gap endpoints  ", f_old, f_new
         if (abs(testgap/aggY) < gaptol) exit bisectloop

         c1 = (testgap*f_old) < 0.0_dp
         c2 = (testgap*f_new) < 0.0_dp
         write(*,*) "positive negative? ", c1, c2
         if (c1)  then
            write(*,*) "c1 ", f_old,f_new,testgap
            write(*,*) "c1 ", w_old,w_new,wage
            f_new = testgap
            w_new = wage
         elseif (c2)  then
            write(*,*) "c2 ", f_old,f_new,testgap
            write(*,*) "c2 ", w_old,w_new,wage
            f_old = testgap
            w_old = wage
         endif
         wage = 0.5_dp*(w_old + w_new)

      enddo bisectloop
      write(*,"('consumption, investment, output, labor, intang, wage', 2x, 6f12.6)") aggC, aggI, aggY, aggL, aggR, wage
      pee(13) = wage
      call makemoments(pee,momvector)
end subroutine bisectme
