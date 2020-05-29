! =================================================
!  parametrization for chiral-odd transversity x h1^dval (x;Q0) at Q0 = 1 GeV  
!  from  Phys. Rev. Lett. 120 (2018) 192001, arXiv:1802.05212 [hep-ph]
!
!  input: x, array of 5 parameters
!  output: x h1^dval(x;Q0)
! =================================================

double precision function xh1down(x,norm,xpow,ceb1,ceb2,ceb3)

	implicit none
	double precision, intent(in):: x, norm, xpow, ceb1, ceb2, ceb3
	double precision:: xsoffb
	
	xsoffb = dabs(                                                  &
    &     -0.012d0 * (1.d0-x)**10 * x**0.164 * (1.d0+98.94d0*x) +   &
    &     ( (0.59964d0 * (1.d0-x)**8.8801 *                         &
    &        (1.d0-2.9012d0*dsqrt(x)+16.865d0*x) ) / x**0.16276     &
    &      -( 0.10302d0 * (1.d0-x)**13.242 *                        &
    &         (1.d0-2.9012d0*dsqrt(x)+16.865d0*x) ) / x**0.16276    &
    &      +17.8826d0 * (1.d0-x)**10.8801 *x**1.876 *               &
    &       (1.d0+8.4703d0*x-36.507d0*x**2) ) / 4.d0                &
	&   ) / 2.d0  +                                                 &
    &   dabs(                                                       &
    &     5.0903d0 * (1.d0-x)**5.1244 * x**0.71978 *                &
    &      (1.d0-4.3654d0*dsqrt(x)+7.473d0*x) +                     &
    &     0.012d0 * (1.d0-x)**10 * x**0.164 * (1.d0+98.94d0*x) -    &
    &     0.015d0 * (1.d0-x)**3.89 * x**0.164 *                     &
    &      (1.d0+22.4d0*dsqrt(x)+98.94d0*x) +                       &
    &     ( ( 0.59964d0 * (1.d0-x)**8.8801 *                        &
    &         (1.d0-2.9012d0*dsqrt(x)+16.865d0*x) ) / x**0.16276    &
    &      -( 0.10302d0 * (1.d0-x)**13.242 *                        &
    &         (1.d0-2.9012d0*dsqrt(x)+16.865d0*x) ) / x**0.16276    &
    &      +17.8826d0 * (1.d0-x)**10.8801 * x**1.876 *              &
    &        (1.d0+8.4703d0*x-36.507d0*x**2) ) / 4.d0               &
	&   ) / 2.d0
	
	xh1down = xsoffb * norm * x**xpow *       &
	&     ( 1.d0 + ceb1*x + ceb2*(2.d0*x*x-1.d0) + ceb3*(4.d0*x*x*x-3.d0*x) ) 
	
end function xh1down

