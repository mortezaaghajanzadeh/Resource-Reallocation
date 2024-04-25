
subroutine findinv(matrix, inverse, n, errorflag)
        implicit none
        !declarations
        integer, intent(in) :: n
        integer, intent(out) :: errorflag  !return error status. -1 for error, 0 for normal
        double precision, intent(in), dimension(n,n) :: matrix  !input matrix
        double precision, intent(out), dimension(n,n) :: inverse !inverted matrix

        logical :: flag = .true.
        integer :: i, j, k, l
        double precision :: m
        double precision, dimension(n,2*n) :: augmatrix !augmented matrix

        !augment input matrix with an identity matrix
        do i = 1, n
                do j = 1, 2*n
                        if (j <= n ) then
                                augmatrix(i,j) = matrix(i,j)
                        else if ((i+n) == j) then
                                augmatrix(i,j) = 1
                        else
                                augmatrix(i,j) = 0
                        endif
                end do
        end do

        !reduce augmented matrix to upper traingular form
        do k =1, n-1
                if (augmatrix(k,k) == 0) then
                        flag = .false.
                        do i = k+1, n
                                if (augmatrix(i,k) /= 0) then
                                        do j = 1,2*n
                                                augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                                        end do
                                        flag = .true.
                                        exit
                                endif
                                if (flag .eqv. .false.) then
                                        print*, "matrix is non - invertible"
                                        inverse = 0
                                        errorflag = -1
                                        return
                                endif
                        end do
                endif
                do j = k+1, n
                        m = augmatrix(j,k)/augmatrix(k,k)
                        do i = k, 2*n
                                augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
                        end do
                end do
        end do

        !test for invertibility
        do i = 1, n
                if (augmatrix(i,i) == 0) then
                        print*, "matrix is non - invertible"
                        inverse = 0
                        errorflag = -1
                        return
                endif
        end do

        !make diagonal elements as 1
        do i = 1 , n
                m = augmatrix(i,i)
                do j = i , (2 * n)
                           augmatrix(i,j) = (augmatrix(i,j) / m)
                end do
        end do

        !reduced right side half of augmented matrix to identity matrix
        do k = n-1, 1, -1
                do i =1, k
                m = augmatrix(i,k+1)
                        do j = k, (2*n)
                                augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
                        end do
                end do
        end do

        !store answer
        do i =1, n
                do j = 1, n
                        inverse(i,j) = augmatrix(i,j+n)
                end do
        end do
        errorflag = 0
    end subroutine findinv
subroutine unique(n,example,res,k)
implicit none
integer, intent(in) :: n
real(8), intent(in) :: example(n)         ! The input
real(8), intent(out) :: res(n)  ! The output
integer, intent(out) :: k ! number of unique elements
integer :: i, j
res = 0.0
k = 1
res(1) = example(1)
outer: do i=2,size(example)
    do j=1,k
       if (res(j) == example(i)) then           ! Found a match so start looking again
         cycle outer
       end if
    end do     ! No match found so add it to the output
    k = k + 1
    res(k) = example(i)
end do outer
end subroutine unique


subroutine normp ( z, p, pdf )
    ! original form: normp ( z, p, q, pdf )

    !*****************************************************************************80
    !
    !! NORMP computes the cumulative density of the standard normal distribution.
    !
    !  Discussion:
    !
    !    This is algorithm 5666 from Hart, et al.
    !
    !  Modified:
    !
    !    13 January 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Alan Miller.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    !    Charles Mesztenyi, John Rice, Henry Thacher,
    !    Christoph Witzgall,
    !    Computer Approximations,
    !    Wiley, 1968,
    !    LC: QA297.C64.
    !
    !  Parameters:
    !
    !    Input, real (8) Z, divides the real line into two
    !    semi-infinite intervals, over each of which the standard normal
    !    distribution is to be integrated.
    !
    !    Output, real (8) P, Q, the integrals of the standard normal
    !    distribution over the intervals ( - Infinity, Z] and
    !    [Z, + Infinity ), respectively.
    !
    !    Output, real (8) PDF, the value of the standard normal
    !    distribution at Z.
    !
    implicit none

    real(8) :: cutoff = 7.071D+00
    real(8) expntl
    real(8) p
    real(8) :: p0 = 220.2068679123761D+00
    real(8) :: p1 = 221.2135961699311D+00
    real(8) :: p2 = 112.0792914978709D+00
    real(8) :: p3 = 33.91286607838300D+00
    real(8) :: p4 = 6.373962203531650D+00
    real(8) :: p5 = 0.7003830644436881D+00
    real(8) :: p6 = 0.03526249659989109D+00
    real(8) pdf
    real(8) q
    real(8) :: q0 = 440.4137358247522D+00
    real(8) :: q1 = 793.8265125199484D+00
    real(8) :: q2 = 637.3336333788311D+00
    real(8) :: q3 = 296.5642487796737D+00
    real(8) :: q4 = 86.78073220294608D+00
    real(8) :: q5 = 16.06417757920695D+00
    real(8) :: q6 = 1.755667163182642D+00
    real(8) :: q7 = 0.08838834764831844D+00
    real(8) :: root2pi = 2.506628274631001D+00
    real(8) z
    real(8) zabs

    zabs = abs ( z )
    !
    !  37 < |Z|.
    !
    if ( 37.0D+00 < zabs ) then
        pdf = 0.0D+00
        p = 0.0D+00
        !
        !  |Z| <= 37.
        !
    else

        expntl = exp ( - 0.5D+00 * zabs * zabs )
        pdf = expntl / root2pi
        !
        !  |Z| < CUTOFF = 10 / sqrt(2).
        !
        if ( zabs < cutoff ) then
            p = expntl * (((((( &
            p6   * zabs &
            + p5 ) * zabs &
            + p4 ) * zabs &
            + p3 ) * zabs &
            + p2 ) * zabs &
            + p1 ) * zabs &
            + p0 ) / ((((((( &
            q7   * zabs &
            + q6 ) * zabs &
            + q5 ) * zabs &
            + q4 ) * zabs &
            + q3 ) * zabs &
            + q2 ) * zabs &
            + q1 ) * zabs &
            + q0 )
            !
            !  CUTOFF <= |Z|.
            !
        else
            p = pdf / ( &
            zabs + 1.0D+00 / ( &
            zabs + 2.0D+00 / ( &
            zabs + 3.0D+00 / ( &
            zabs + 4.0D+00 / ( &
            zabs + 0.65D+00 )))))
        end if

    end if

    if ( z < 0.0D+00 ) then
        q = 1.0D+00 - p
    else
        q = p
        p = 1.0D+00 - q
    end if

    return

    end subroutine normp



!        subroutine legendre(lo, hi, num_nodes, nodes, weights)
!
!              use sizes
!
!              implicit none
!              real(8), intent(in) :: lo, hi
!              integer, intent(in) :: num_nodes
!              real(8), intent(out) :: nodes(num_nodes), weights(num_nodes)
!
!              real(8) :: tol, xm, xl
!              real(8) :: z, z1, p1, p2, p3, pp
!              integer :: m, i, j
!
!              tol = eps
!              m = floor((num_nodes+1)/2.0)
!              xm = (hi+lo)/2.0
!              xl = (hi-lo)/2.0
!
!              do i = 1,m
!                     z = cos(pi*(real(i)-0.25)/(real(num_nodes)+0.5))
!                     ! initial guess at z
!                     z1 = z+0.1
!                     do while (abs(z-z1) > tol)
!                            p1 = 1.0
!                            p2 = 0.0
!                            do j = 1,num_nodes
!                                   p3 = p2
!                                   p2 = p1
!                                   p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
!                            enddo
!                            pp = num_nodes*(z*p1-p2)/(z*z-1.0)
!                            z1 = z
!                            z = z1-p1/pp
!                     enddo
!                     nodes(i) = xm-xl*z
!                     nodes(num_nodes+1-i) = xm+xl*z
!                     weights(i) = 2.0*xl/((1.0-z*z)*pp*pp)
!                     weights(num_nodes+1-i) = weights(i)
!              enddo
!
!       end subroutine legendre

       subroutine spline (x, y, b, c, d, n)
       !======================================================================
       !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
       !  for cubic spline interpolation
       !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
       !  for  x(i) <= x <= x(i+1)
       !  Alex G: January 2010
       !----------------------------------------------------------------------
       !  input..
       !  x = the arrays of data abscissas (in strictly increasing order)
       !  y = the arrays of data ordinates
       !  n = size of the arrays xi() and yi() (n>=2)
       !  output..
       !  b, c, d  = arrays of spline coefficients
       !  comments ...
       !  spline.f90 program is based on fortran version of program spline.f
       !  the accompanying function fspline can be used for interpolation
       !======================================================================

       implicit none
       integer :: n
       real(8) :: x(n), y(n), b(n), c(n), d(n)
       integer :: i, j, gap
       real(8) h

       gap = n-1
       ! check input
       if ( n < 2 ) return
       if ( n < 3 ) then
         b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
         c(1) = 0.
         d(1) = 0.
         b(2) = b(1)
         c(2) = 0.
         d(2) = 0.
         return
       end if
       !
       ! step 1: preparation
       !
       d(1) = x(2) - x(1)
       c(2) = (y(2) - y(1))/d(1)
       do i = 2, gap
         d(i) = x(i+1) - x(i)
         b(i) = 2.0*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
       end do
       !
       ! step 2: end conditions
       !
       b(1) = -d(1)
       b(n) = -d(n-1)
       c(1) = 0.0
       c(n) = 0.0
       if(n /= 3) then
         c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
         c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
         c(1) = c(1)*d(1)**2.0/(x(4)-x(1))
         c(n) = -c(n)*d(n-1)**2.0/(x(n)-x(n-3))
       end if
       !
       ! step 3: forward elimination
       !
       do i = 2, n
         h = d(i-1)/b(i-1)
         b(i) = b(i) - h*d(i-1)
         c(i) = c(i) - h*c(i-1)
       end do
       !
       ! step 4: back substitution
       !
       c(n) = c(n)/b(n)
       do j = 1, gap
         i = n-j
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
       end do
       !
       ! step 5: compute spline coefficients
       !
       b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
       do i = 1, gap
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.0*c(i)
       end do
       c(n) = 3.0*c(n)
       d(n) = d(n-1)
       end subroutine spline


      subroutine indnv(val,vec,n,whereis)

         implicit none

         integer, intent(in) :: n
         real(8), intent(in) :: val, vec(n)
         integer, intent(out) :: whereis

         integer :: ii, temp

         temp = 0
         ii = 1
         pickloop: do

            if (ii > n-1) then
              temp = n
              exit pickloop
            else
              if (val > vec(ii) .and. val <= vec(ii + 1)) then
                 temp = ii
                 exit pickloop
              endif
            endif
            ii = ii + 1
         enddo pickloop

         whereis = temp

      end subroutine indnv

subroutine dools(len,nreg,lhs,desp,bhat,r2)
         implicit none
         integer, intent(in) :: len, nreg
         real(8), intent(in)  :: lhs(len),desp(len,nreg)
         real(8), intent(out) :: bhat(nreg+1),r2

         real(8) :: des(len,nreg+1),tdes(nreg+1,len),dd(nreg+1,nreg+1),idd(nreg+1,nreg+1),uhat(len),seb(nreg+1), my, sigu, vy
         integer :: ii,errorflag
         des = 1.0
         des(:,2:(nreg+1)) = desp

         tdes     = transpose(des)
         dd = matmul(tdes,des)
         call findinv(dd, idd, nreg+1, errorflag)

         bhat = matmul(tdes,lhs)
         bhat = matmul(idd,bhat)
         uhat = lhs - matmul(des,bhat)

         sigu = sum(uhat**2.0)/real(size(uhat))
         my = sum(lhs)/real(size(lhs))
         vy = sum((lhs - my)**2.0)/real(size(lhs))
         r2 = 1.0- sigu/vy
         seb  = 0.0
         do ii=1,nreg+1
            seb(ii) = sqrt(sum(uhat**2)/real(size(uhat))*idd(ii,ii))
         enddo

end subroutine dools
