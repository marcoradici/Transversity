module evo_grid   ! calculate grid in 1< Q< 100 GeV and 0.0001< x< 1 for evolution

	implicit none
	
	integer, parameter:: nQ=50, nx=100
	double precision, parameter:: qmin=1.d0, qmax=1.d2  ! consistent with MSTW08
	double precision, parameter:: xmin=1.d-4 			!   "         "     "
	double precision:: Qvec(nQ), xvec(nx)

	save
	contains
	
	subroutine compute_grids
	
		integer:: i
	    
	    do i=1,nQ
		  Qvec(i) = 0.d0
	    end do
	    do i=1,nx
		  xvec(i) = 0.d0
	    end do
	
	    do i=1,nQ
		  Qvec(i) = qmin * (qmax/qmin)**(dble(i-1)/dble(nQ))
	    end do
	
	    do i=1,nx-20           ! for 0.0001< x< 0.1  80 points log-distributed
		  xvec(i) = dexp( dlog(1.d-4)*dble(106-i+1)/dble(106) )
	    end do
	    do i=1,20			   ! for 0.1< x< 1  20 equidistant points 
		  xvec(nx-20+i) = 1.d-1 + (1.d0-1.d-1) * dble(i-1)/dble(20)
	    end do
		
	end subroutine compute_grids
	
end module evo_grid
	
