!	MC PROGRAM TO SIMULATE LIQUID ARGON

	program mc

	implicit none
	integer isd,i,j,m,l,k,n,natom,ncyc,nx,bin,count,nbin(1000000),neq
	real*4 boxl,d_m,dE,ef,Efa,ei,Eia,ene,En,eo,eps,p_r,rsqr,rfhex,rfsqr
	real*4 rftwelve,rhex,random,rihex,risqr,T,RT,R,rnsqr,rntwelve,rupper,rlower
	real*4 sgmahex,sgmatwlv,rnhex,rtwelve,ritwelve,yt,unitcelllength,const,xr(3),xrsq
	real*4 a(100000,3),a_p(100000,3),x(3),yi(3),OE,b,c,rdfav(1000000),nr
	real*4 g(100000),h(100000),xm(3),E(100000),xf(3),yn(3),Kb,nid,gr,rn,rf,ri
	character(2) s(50000)

!	OPEN FILES
	open(unit=1,file='inputparameter.dat')
	open(unit=2,file='in_mcargon.dat')
	open(unit=3,file='simu.xyz')
	open(unit=4,file='ene.dat')
	open(unit=5,file='output.dat')
	open(unit=10,file='accpprob.dat')
	open(unit=7,file='finalconfi.xyz')

!	
	call random_seed(isd)
	b=0;c=0
	boxl=5*4.81


	read(1,*)
        read(1,*) nx,unitcelllength

        ! READ THE INPUT DATA 
!
        read(2,*)
        read(2,*) isd,d_m,ncyc,neq,eps,sgmahex,sgmatwlv,R,T

        boxl=nx*unitcelllength
        natom=nx*nx*nx*4

! 	WRITE IMPORTANT PARAMTERS OF THESE MC SIMULATION 	
	
	write(5,*) 'IMPORTANT PARAMTERS OF MC SIMULATION'
	write(5,*) 'iseed	dr_max(nanometer)	no of cycles	noofatoms'	
	write(5,*) isd,d_m,ncyc,natom
	write(5,*) 'Parameters of L-J Potential'
	write(5,*) 'eps	sigmahex	sigmatwlv  (sigma is in nanometer)'
	write(5,*) eps,sgmahex,sgmatwlv
	write(5,*) 'Boltzmann constant	Temperature	avogadro number'
	write(5,*) R,T
	

! 	 READ THE INITIAL CONFIGURATION FILE 
	read(3,*)
	read(3,*)
	do i=1,natom
		read (3,*) s(i),a(i,1),a(i,2),a(i,3)
	enddo !i
	

	write (5,*) '------------------------------------'

	
! 	CALCULATE THE INITIAL ENERGY 

!	lj energy	!
	OE=0.0 !Set initial energy at zero
	do i=1,natom-1
		do j=i+1,natom                     !natom is the no of atoms in the box
				eo=0.0d0
				do m=1,3
					x(m)=a(i,m)-a(j,m) !x() is the distance between ith and jth particle
					x(m)=x(m)-boxl*anint(x(m)/boxl) !Minimum image convention
				enddo !m
			rsqr=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)

			rhex=rsqr*rsqr*rsqr
			rtwelve=rhex*rhex
			eo=4*eps*((sgmatwlv/rtwelve)-(sgmahex/rhex))	! eo is the interaction energy between ith and jth particle
			OE=OE+eo
20		enddo !j 
	enddo !i	!OE is the total initial energy	


		
! 	WRITE INITIAL ENERGY IN OUTPUT 
	
	write(5,*) 'Initial energy of the system in Kj/mole'
	write(5,*) OE
	write(5,*) 'Initial energy of each particle in Kj/mole'
	write(5,*) OE/natom
	
	write(5,*) '----------------------------------------'


	
			! ---------------- MC simulation ------------------!
	


	do n=1,ncyc   ! ncycle is the no of simulation
	Efa=0.0 ! Set final energy at zero at the begining of each loop


	do i=1,natom	
	Eia=0.0;En=0.0 			! Set initial and new energy at zero at the begining of each loop
					! Ei and En are the initial and new energy between ith and all other particles 
		

	! CALCULATE INITIAL ENERGY(Ei) OF iTH PARTICLE WITH OTHERS  !
			
			do j=1,natom
				if (i.ne.j) then
					ei=0.0d0
					do m=1,3
						yi(m)=a(i,m)-a(j,m)
						yi(m)=yi(m)-boxl*anint(yi(m)/boxl) 
					enddo !m
				risqr=yi(1)*yi(1)+yi(2)*yi(2)+yi(3)*yi(3)	
				
			!	if(sqrt(risqr).le.2.5*3.41) then
				rihex=risqr*risqr*risqr
				ritwelve=rihex*rihex                  ! 'i' stands for initial
				ei=4*eps*((sgmatwlv/ritwelve)-(sgmahex/rihex))
			!ei is the initial energy between ith and jth particle
				Eia=Eia+ei
			!	endif
				endif 
			enddo !j


		! MOVE iTH PARTICLE RANDOMLY   !
			
			do j=1,3
				call random_number(yt)
				a_p(i,j)=a(i,j)+(2*yt-1)*d_m  ! move ith particle randomly
				if(a_p(i,j).gt.boxl) a_p(i,j)=a_p(i,j)-boxl*1.0d0
				if(a_p(i,j).lt.0.0d0) a_p(i,j)=a_p(i,j)+boxl*1.0d0
			enddo !j
		!


		! CALCULATE NEW ENERGY(En) OF iTH PARTICLE WITH OTHERS  !
			do j=1,natom
				ene=0.0d0
				if (i.ne.j) then
					do m=1,3
						yn(m)=a_p(i,m)-a(j,m)
						yn(m)=yn(m)-boxl*anint(yn(m)/boxl)
					enddo !m 
				rnsqr=yn(1)*yn(1)+yn(2)*yn(2)+yn(3)*yn(3)
			!	if(sqrt(rnsqr).le.3.41*2.5) then
				rnhex=rnsqr*rnsqr*rnsqr
				rntwelve=rnhex*rnhex			! 'n' stands for new
				ene=4*eps*((sgmatwlv/rntwelve)-(sgmahex/rnhex))
		!En is the new energy between ith and jth particle 
				En=En+ene
				!write(20,*) n,a_p(i,1),a_p(i,2),a_p(i,3),rnsqr,sgmatwlv/rntwelve,sgmahex/rnhex,En
			!	endif		
				endif
			enddo !j
		
		!write (*,*) n,'new energy',En


		! SAMPLING  !
			dE=En-Eia   ! dE is the difference between new and initial energy
			if (dE.gt.0) then
				call random_number(yt)
				p_r=exp(-dE/(R*T)) !Kb is the Boltzmann factor and T is the temperature
			endif 
			if (dE.lt.0.or.yt.lt.p_r) then
				do j=1,3
					a(i,j)=a_p(i,j)	!accept the new configuration
				enddo !j
				b=b+1
			else
			c=c+1
			endif
		
		enddo !i
	


	! CALCULATE THE FINAL ENERGY AT nTH CYCLE  !
		do i=1,natom-1
			do j=i+1,natom
				ef=0.0d0
			!	if (i.eq.j) goto 30
					do m=1,3
						xf(m)=a(i,m)-a(j,m)
						xf(m)=xf(m)-boxl*anint(xf(m)/boxl)
					enddo !m
				rfsqr=xf(1)*xf(1)+xf(2)*xf(2)+xf(3)*xf(3)
			!	if(sqrt(rfsqr).le.2.5*3.41) then
					rfhex=rfsqr*rfsqr*rfsqr
					rftwelve=rfhex*rfhex		! 'f' stands for final
					ef=4*eps*((sgmatwlv/rftwelve)-(sgmahex/rfhex))
					Efa=Efa+ef
			!	endif
30			enddo !j

		enddo !i

		if(n.ge.neq) then
		do i=1,natom
			write(40,*) n,s(i),a(i,1),a(i,2),a(i,3)
		enddo
		endif

!---- write no of nth cycle,Energy at nth cycle,Energy of each particle after nth cycle  in 'energy.dat' file to plot---!
		write(10,*) n,b/(b+c),c/(b+c)	
		write(*,*) n,Efa/natom    ! !Efa is the total final energy of nth cycle
	enddo !n
!


! --------- WRITE THE FINAL CONFIGURATION FILE ------------!
	write(7,*) natom
	write(7,*)
	do i=1,natom ! natom is the no of atom in simulation box
		write(7,*) s(i),a(i,1),a(i,2),a(i,3)
	enddo !i


	! -----------------WRITE FINAL ENERGY IN OUTPUT ----------------!
	write(5,*) 'Final energy of the system in Kj/mole',Efa
	write(5,*) 'Final energy of each particle in Kj/mole'
	write(5,*) Efa/natom

	end program mc
