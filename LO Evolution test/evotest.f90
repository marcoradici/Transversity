! =================================================
!  test chiral-odd evolution APFEL++ using parametrization for 
!  transversity x h1^qval (x;Q0) at Q0 = 1 GeV  
!  from  Phys. Rev. Lett. 120 (2018) 192001, arXiv:1802.05212 [hep-ph]
!
!  input parameters from fit_params_g=u.dat
!  compare with xh1g=u_q0(/q4/q10).txt
! =================================================

program Evotest

	implicit none
	
	integer, parameter:: nrepl=600
	integer:: i, irepl, ix
    integer:: PerturbativeOrder
    double precision:: mu0, AlphaQCDRef, MuAlphaQCDRef, mc, mb, p(10)
	double precision:: x, q(3), tmp, hmax, xh1up, xh1down, xf(-6:6)
	character*49::  input
	character*100:: param
	
! Scales
	q(1) = 1.d0
	q(2) = dsqrt(4.d0)
	q(3) = dsqrt(1.d1)
	
! Input directory and file
    input='../output/SIDIS+STAR06/a13939/new_param/gTfinite/'
    param=input//'fit_params_g=u.dat'
    open(unit=1,file=param,status='old')
		
! Output
    open(unit=11,file='test_q0.txt')
    open(unit=12,file='test_q4.txt')
    open(unit=13,file='test_q10.txt')
	
! Evolution constants
    PerturbativeOrder = 0
    mu0 = 1d0
    AlphaQCDRef = 0.13939d0
    MuAlphaQCDRef = 91.1876d0
    mc = 1.4d0
    mb = 4.5d0

	! Loop on replicas
	do i=1,nrepl
! Parameters of transversity at Q0=1 GeV
! up
	  read(1,*) irepl,tmp,p(7),p(8),p(9),p(10),hmax
	  p(6) = tmp/hmax
! down
      read(1,*) irepl,tmp,p(2),p(3),p(4),p(5),hmax
	  p(1) = tmp/hmax

! Initialize APFEL++
	  call InitialiseEvolution(PerturbativeOrder, mu0, AlphaQCDRef, MuAlphaQCDRef, mc, mb, p)

! Evolve PDF at various scales and in equidistant grid; then, write output
	  do ix=1,100
		x = 1.d-3 + (ix-1) * 1.d-2
		write(11,*) i,x,q(1)*q(1),xh1up(x,p(6),p(7),p(8),p(9),p(10)),    &
	    &           xh1down(x,p(1),p(2),p(3),p(4),p(5))
	  end do

      do ix=-6,6
	    xf(ix) = 0.d0
	  end do
	  	
	  do ix=1,100
		x = 1.d-3 + (ix-1) * 1.d-2
		call EvolveTransversities(x, q(2), xf)
		write(12,*) i,x,q(2)*q(2),xf(2),xf(1)
	  end do

      do ix=-6,6
	    xf(ix) = 0.d0
	  end do

	  do ix=1,100
		x = 1.d-3 + (ix-1) * 1.d-2
		call EvolveTransversities(x, q(3), xf)
		write(13,*) i,x,q(3)*q(3),xf(2),xf(1)
	  end do

	end do
	
	close(unit=1)
	close(unit=11)
	close(unit=12)
	close(unit=13)
	
end program Evotest
	
	
	