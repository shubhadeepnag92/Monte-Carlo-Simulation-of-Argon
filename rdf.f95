!	FORTRAN PROGRAM TO CALCULATE RDF OF LIQUID Ar

	program mc

	implicit none
	integer isd,i,j,m,l,k,n,natom,ncyc,nx,bin,count,nbin(100000)
	real boxl,d_m,dE,ef,Efa,ei,Eia,ene,En,eo,eps,p_r,rsqr,rfhex,rfsqr
	real rftwelve,rhex,random,rihex,risqr,T,RT,R,rnsqr,rntwelve
	real sgmahex,sgmatwlv,rnhex,rtwelve,ritwelve,yt,unitcelllength
	real a(10000,3),a_p(10000,3),x(10000),yi(10000),OE,b,c,rdfav(100000)
	real g(10000),h(10000),xm(10000),E(10000),xf(10000),yn(10000),Kb
	character(2) s(5000)

!	OPEN FILES
	open(unit=1,file='inputparameter.dat')
	open(unit=2,file='in.dat')
	open(unit=3,file='finalconfi.xyz')
	open(unit=6,file='rdf.dat')
	open(unit=15,file='positionofneighbor.dat')

!	
	b=0;c=0
	boxl=5*4.81


	read(1,*)
        read(1,*) nx,unitcelllength

        ! READ THE INPUT DATA 
!
        read(2,*)
        read(2,*) isd,d_m,ncyc,eps,sgmahex,sgmatwlv,R,T

        boxl=nx*unitcelllength
        natom=nx*nx*nx*4


! 	 READ THE INITIAL CONFIGURATION FILE 
	read(3,*)
	do i=1,natom
		read (3,*) s(i),a(i,1),a(i,2),a(i,3)
	enddo !i
	
	
! 	CALCULATE THE rdf 

	do i=1,natom-1
		do j=i+1,natom                     !natom is the no of atoms in the box
				do m=1,3
					x(m)=a(i,m)-a(j,m) !x() is the distance between ith and jth particle
					x(m)=x(m)-boxl*anint(x(m)/boxl) !Minimum image convention
				enddo !m
			rsqr=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)

			bin=sqrt(rsqr)/0.001d0
			nbin(bin)=nbin(bin)+1
			write(15,*) i,j,sqrt(rsqr),bin,sqrt(rsqr)/3.41,nbin(bin)
		enddo !j 
	enddo !i	!OE is the total initial energy
	
	count=(5*4.81)/0.001d0
	do i=1,count
		rdfav(i)=nbin(i)*1.0d0/(natom*2)
		write(6,*) i*0.001d0,nbin(i),rdfav(i)
	enddo

	end program rdf
