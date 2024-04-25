subroutine momentgen(params,momvector,vg)
      use datatype
      use sizes
      use globals
      use omp_lib
      implicit none

      real(dp), intent(in) :: params(nop)
      real(dp), intent(out) :: momvector(nmom)
      real(dp), intent(inout) :: vg(nz,nk,nn,nb)

      allocate(simdata(numv,nYears,nFirms))

      call valfun(params,vg)
      momvector = -10.0_dp

      call makemoments(params,momvector)
      deallocate(simdata)


end subroutine momentgen

