program Driver
  implicit none
  integer i
  integer PerturbativeOrder
  double precision mu0
  double precision AlphaQCDRef
  double precision MuAlphaQCDRef
  double precision mc
  double precision mb
  double precision p(10)

  double precision x
  double precision mu
  double precision xf(-6:6)

  PerturbativeOrder = 0
  mu0 = 1d0
  AlphaQCDRef = 0.13939d0
  MuAlphaQCDRef = 91.1876d0
  mc = 1.4d0
  mb = 4.5d0
  ! Parameters. The first 5 are for d_valence the second 5 for
  ! u_valence
  do i = 1, 10
     p(i) = 10
  enddo
  call InitialiseEvolution(PerturbativeOrder, mu0, AlphaQCDRef, MuAlphaQCDRef, mc, mb, p)

  x = 0.1d0
  mu = 10d0
  call EvolveTransversities(x, mu, xf);
  do i = -6, 6
     write(6,*) i,xf(i)
  enddo
end program
