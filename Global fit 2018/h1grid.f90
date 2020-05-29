! Code to produce a grid for transversity from global fit 2018
!
! Reference: PRL 120 (2018) 192001, arXiv:1802.05212v2 [hep-ph] 
! 
!  output format:   xh1_PV18_XX_REP.dat  where XX = LO / NLO  using MSTW08 parametrization and
!                                        and REP indicates the replica: REP = 001, 002, ... , 600  
!                                       (200 replicas for each of 3 scenarios for gluons, organized as 
!                  #1 repl=1 for D1g(Q0)=0 ,  #2 repl=1 for D1g(Q0)= D1u(Q0)/4 , #3 repl=1 for D1g(Q0)=D1u(Q0)
!                  #4 repl=2 for   "       ,  #5 repl=2 for      "    "        , #6 repl=2 for    "    "      ...)
!
! for each file the output is:   grid in Q ,  grid in x  , then xh1q for > q =  t-bar b-bar  c-bar  s-bar  u-bar  d-bar  g  d  u  s  c  b  t   
! where 1< Q< 100 GeV, 0.0001< x< 1. 
! u_sea = d_sea = 0 at Q0 = 1 and non-singlet evolution  ==> 1) at NLO, sea quark generated by evolution
!                                                            2) at LO,  u = uval and d = dval at any Q 
!  compilation: see Makefile 

program h1grid  

	use params
	use evo_grid
	
	implicit none
	
	integer:: i, j, k, l
	double precision:: p(10), xh1(-6:6), xh1up, xh1down
	
	character(len=27):: fname   ! LO
!	character(len=29):: fname   ! NLO


!    read parameters to build xh1(x,Q0) for all replicas
	call read_h1_params

!    prepare grids in Q , x
	call compute_grids

!   start loop on replicas
	do i=1,nreplica

!    input parameters for replica i of xh1(x,Q0) 
! down
      p(1) = norm(1,i)
      p(2) = xpow(1,i)
      p(3) = ceb1(1,i)
      p(4) = ceb2(1,i)
      p(5) = ceb3(1,i)
! up
      p(6) = norm(2,i)
      p(7) = xpow(2,i)
      p(8) = ceb1(2,i)
      p(9) = ceb2(2,i)
      p(10) = ceb3(2,i)
	  
! 	Initialize APFEL++
      call InitialiseEvolution(p)

! 	prepare output in xh1_PV18_XX_REP.dat with  XX=LO/NLO and REP=replica
      write(fname,'(a,i3.3,a)') 'grid/LO/xh1_PV18_LO_',i,'.dat'      
!      write(fname,'(a,i3.3,a)') 'grid/NLO/xh1_PV18_NLO_',i,'.dat'      
      open(unit=1,file=fname,status='unknown') 
      write(1,*) ' Q ' 
      write(1,'(5(e11.5,1x))') (Qvec(j),j=1,nQ)
      write(1,*)
      write(1,*) ' x '
      write(1,'(5(e11.5,1x))') (xvec(j),j=1,nx)
      write(1,*)
      write(1,*)
      write(1,*) ' anti-t       anti-b       anti-c       anti-s      ',&
      &          ' anti-u       anti-d       g            d           ',&
      &          ' u            s            c            b            t'

!    evolve xh1 at predefined grid in (Q,x) and write output
      do j=1,nQ
        do k=1,nx
           
		  do l=-6,6
            xh1(l) = 0.d0
          end do
           
		  if(j == 1) then
            xh1(1) = xh1down(xvec(k),p(1),p(2),p(3),p(4),p(5))
            xh1(2) = xh1up(xvec(k),p(6),p(7),p(8),p(9),p(10))
          else	 
            call EvolveTransversities(xvec(k), Qvec(j), xh1)
          end if
           
		  write(1,'(13(1x,e12.5))') (xh1(l),l=-6,6)
        
		end do
     end do
     
	 close(unit=1)
   
   end do

   stop

end program h1grid