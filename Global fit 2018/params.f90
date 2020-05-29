module params   ! read input parameters for xh1(x,Q0) from  gTfinite/fit_params_g=u.dat

	implicit none
	
	integer, parameter:: nreplica=600
	double precision, dimension(2,nreplica):: norm, xpow, ceb1, ceb2, ceb3
	
	save
	contains
	
	subroutine read_h1_params
		integer:: i, dummy
		double precision:: tmp, hmax
		character(len=23):: input
	    
		input = 'grid/fit_params_g=u.dat'
	    open(unit=1,file=input,status='old')
		
	    do i=1,nreplica
! up
		  read(1,*) dummy,tmp,xpow(2,i),ceb1(2,i),ceb2(2,i),ceb3(2,i),hmax
		  norm(2,i) = tmp / hmax
! down
		  read(1,*) dummy,tmp,xpow(1,i),ceb1(1,i),ceb2(1,i),ceb3(1,i),hmax
		  norm(1,i) = tmp / hmax
	  	end do
		
    	close(unit=1)
	end subroutine read_h1_params

end module params
	
