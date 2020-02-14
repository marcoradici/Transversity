! =================================================
!  parametrization for chiral-odd transversity x h1^uval (x;Q0) at Q0 = 1 GeV  
!  from  Phys. Rev. Lett. 120 (2018) 192001, arXiv:1802.05212 [hep-ph]
!
!  input: x, array of 5 parameters
!  output: x h1^uval(x;Q0)
! =================================================

double precision function xh1up(x,norm,xpow,ceb1,ceb2,ceb3)

	implicit none
	double precision, intent(in):: x, norm, xpow, ceb1, ceb2, ceb3
	double precision:: xsoffb
	
	xsoffb = dabs(                                               &
    &    0.295d0 * (1.d0-8.42d0*x) * (1.d0-x)**10. * x**0.692 +  &
    &    ( ( 0.59964d0 * (1.d0-x)**8.8801 *                      &
    &        (1.d0-2.9012d0*dsqrt(x)+16.865d0*x) ) / x**0.16276  &
    &     -( 0.10302d0 * (1 - x)**13.242 *                       &
    &        (1.d0-2.9012d0*dsqrt(x)+16.865d0*x) ) / x**0.16276  &
    &     -17.8826d0 * (1-x)**10.8801 * x**1.876 *               &
    &      (1.d0+8.4703d0*x-36.507d0*x**2) ) / 4.d0              &
	&    ) / 2.d0  +                                             &    
    &   dabs(                                                    &
    &     -0.295d0 * (1.d0-8.42d0*x) * (1.d0-x)**10. * x**0.692 + &
    &      1.4335d0 * (1.d0-x)**3.0409 * x**0.45232 *            &
    &       (1.d0-2.3737d0*dsqrt(x)+8.9924d0*x) +                &
    &      0.677d0 * (1.d0-x)**3.34 * x**0.692 *                 &
    &       (1.d0-2.18d0*dsqrt(x)+15.87d0*x) +                   &
    &    ( ( 0.59964d0 * (1.d0-x)**8.8801 *                      &
    &        (1.d0-2.9012d0*dsqrt(x)+16.865d0*x) ) / x**0.16276  &
    &     -( 0.10302d0 * (1.d0-x)**13.242 *                      &
    &        (1.d0-2.9012d0*dsqrt(x)+16.865d0*x) ) / x**0.16276  &
    &     -17.8826d0 * (1.d0-x)**10.8801 * x**1.876 *            &
    &     (1.d0+8.4703d0*x-36.507d0*x**2) ) / 4.d0               &
	&    ) / 2.d0
	
	xh1up = xsoffb * norm * x**xpow *       &
	&     ( 1.d0 + ceb1*x + ceb2*(2.d0*x*x-1.d0) + ceb3*(4.d0*x*x*x-3.d0*x) ) 
	
end function xh1up
