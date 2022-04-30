
!**********************************************************************
!    This program solves 1D Poisson equation by Hyperbolic approach
!	
!	 Written by Sainath Ch, Email: s.chamarthi@al.t.u-tokyo.ac.jp
!**********************************************************************

!**********************************************************************
!    Roe, Rusanov, or HLL methods can be used. In 1D both are same
!	
!	 Fluxes can be either interpolated or reconstructed using appropriate 5th order compact formulas
	
	! 1. subroutine FX
	! 2. subroutine reconstruction3
	
	! Limiters are not used for the flux computation. I would say they should not be used

	! 3. Exact solution in subroutine exact

	! 4. Computation stops when the L1 norm is reduced by 10 orders of magnitude
	! 3rd order TVD-Runge Kutta is used for time integration. CFL number can be increased by using higher order time integration 

	! 5. subroutine boundarycondition
	! 5th order boundary conditions are used. The values in the ghost cells are extrapolated by Lagrange interpolation
!**********************************************************************



program hyperpoi3
		implicit none


		double precision, parameter 		:: pi=acos(-1.0d0)

		integer, parameter					:: NTMAX=300000
		integer 							:: i, j, z
		integer, parameter 					:: NX =24*24*4, n_eqn = 2, ghostp = 5						! No.of cells
		double precision, parameter			:: l1_target = 1.0d-12
		
		double precision 					:: dx, X(1:NX),xmin, xmax
		double precision 					:: cons(n_eqn,-5:NX+6),cons_old(n_eqn,-5:NX+6)
		double precision					:: cons_exact(n_eqn,1:NX)
										
		double precision					:: source(n_eqn,1:NX)
		double precision 					:: residual(n_eqn,1:NX)
				
		double precision 					:: dt, CFL						
		
		double precision 					:: error_inf, error_l1, error_l2

   		double precision					:: res_norm_initial(n_eqn,3),res_norm(n_eqn,3)
	
		integer								:: L1 = 1, L2 = 2, Linf = 3, i_Linf, i_Linf_initial

		character(len=3) 					:: gridsize
		character(len=5) 					:: saver
   		character(len=1024) 				:: fname
		
		integer 							:: time_ini,time_end,time_calc
	    double precision 					:: start, finish
	    integer 							:: file_number, k



	    call cpu_time(start)
	    write(*,*) 'Program start...'

	    call system_clock(count = time_ini)


	    xmin= 0.0d0; xmax =1.0d0
	    dx = (xmax-xmin)/float(NX)
	    

	    	do i=1,NX
	        	X(i)=(i-0.5d0)*dx
	    	enddo

	    call timestep		  (dx,dt,CFL,NX)				! Steady state
	    

	    ! Initial conditions for Poisson equation
	    

		call initial(NX,cons,n_eqn)

		! Compute the exact solution
		call exact(NX,X,cons_exact,n_eqn)


		


	!***********************************************************************
	!*****         			 TVD- Runge Kutta 				       	   *****
	!***********************************************************************
		
	iteration : do j=1,NTMAX

			cons_old = cons
			

			call boundarycondition(NX,cons,n_eqn,ghostp)
			call FX (NX,dx,x,cons,residual,n_eqn,ghostp)
			call sourceterm(NX,x,residual,cons,n_eqn)

!	 Runge Kutta step-1
		do z =  1,n_eqn

			do i =	1,NX
						cons(z,i) 	=	cons_old(z,i) + dt * (residual(z,i))
			enddo 

		enddo
			

			call boundarycondition(NX,cons,n_eqn,ghostp)
			call FX (NX,dx,x,cons,residual,n_eqn,ghostp)
			call sourceterm(NX,x,residual,cons,n_eqn)

! Runge Kutta step-2
		do z =  1,n_eqn

			do i =	1,NX
					cons(z,i) = (3.0d0/4.0d0) * cons_old(z,i) + (1.0d0/4.0d0)* cons(z,i) + (1.0d0/4.0d0)*dt*(residual(z,i))
			enddo

		enddo	

			call boundarycondition(NX,cons,n_eqn,ghostp)
			call FX (NX,dx,x,cons,residual,n_eqn,ghostp)
			call sourceterm(NX,x,residual,cons,n_eqn)

! Runge Kutta step-3
		do z =  1,n_eqn

			do i =	1,NX
				cons(z,i) = (1.0d0/3.0d0)*cons_old(z,i) + (2.0d0/3.0d0)* cons(z,i) + (2.0d0/3.0d0)*dt*(residual(z,i))
			enddo

		enddo	

	
	!***************************************************************************************************
	!*****                      ! Monitor residual norms to check steady convergence.			   *****
	!***** 	 I.e., see how well the numerical solution satisfies the target discrete equaton.      *****            
	!***************************************************************************************************
		

		 res_norm(:,L1  ) =  0.0d0
		 res_norm(:,L2  ) =  0.0d0
		 res_norm(:,Linf) = -1.0d0

		 do z = 1, NX
		  do k = 1, n_eqn ! Two equations
			   
			   res_norm(k,L1) = res_norm(k,L1) + abs(residual(k,z))
			   res_norm(k,L2) = res_norm(k,L2) + abs(residual(k,z))**2
			   
			   if (abs(residual(k,i-1)) >  res_norm(k,Linf) ) then
			                   i_Linf = z
			         res_norm(k,Linf) = abs(residual(k,z))
			   endif
		  end do
		 end do

		  res_norm(:,L1) =       res_norm(:,L1)/real(NX)
		  res_norm(:,L2) = sqrt( res_norm(:,L2)/real(NX) )

		 if (j==1) then 

			  res_norm_initial = res_norm

			  write(*,*) " Initial residual (after 1 iteration):"
			  write(*,*) "  Residual norm (L1) = ", res_norm(:,L1)
			  write(*,*) "  Residual norm (L2) = ", res_norm(:,L2)
			  write(*,*) "  Residual norm (Li) = ", res_norm(:,Linf), " i = ", i_Linf
			  write(*,*)
			  i_Linf_initial = i_Linf

		 endif

		 if (mod(j,20000)==0) then
		  write(*,*) " j = ", j
		  write(*,*) "  Residual norm (L1) = ", res_norm(:,L1  )
		  write(*,*) "  Residual norm (L2) = ", res_norm(:,L2  )
		  write(*,*) "  Residual norm (Li) = ", res_norm(:,Linf)
		  write(*,*)
		 endif

		 !Stop when the residual norm is reduced by 10 orders of magnitude in the L1 norm:
		 if ( minval( res_norm(:,L1)/res_norm_initial(:,L1) )  < l1_target) then
		  exit iteration
		 endif


	enddo iteration
		

	

		write(*,*)
		write(*,*) " Initial residual norm (after 1 iteration):"
		write(*,'(a,2es15.5)') "  Residual norm (L1) = ", res_norm_initial(:,L1)
		write(*,'(a,2es15.5)') "  Residual norm (L2) = ", res_norm_initial(:,L2)
		write(*,'(a,2es15.5,a,i6)') "  Residual norm (Li) = ", res_norm_initial(:,Linf), " i = ", i_Linf_initial
		write(*,*)

		write(*,*)
		write(*,*) " Final residual norm:"
		write(*,'(a,2es15.5)') "  Residual norm (L1) = ", res_norm(:,L1)
		write(*,'(a,2es15.5)') "  Residual norm (L2) = ", res_norm(:,L2)
		write(*,'(a,2es15.5,a,i6)') "  Residual norm (Li) = ", res_norm(:,Linf), " i = ", i_Linf
		write(*,*)

		
	!***********************************************************************
	!*****                Error w.r.t exact solution             	   *****
	!***********************************************************************
	    error_inf	= 0.0d0
	    error_l1	= 0.0d0
	    error_l2	= 0.0d0
	    
	    do i=1,NX

	        error_inf	=	MAX(error_inf,ABS(cons(1,i)-cons_exact(1,i)))
	        error_l1	=   ABS(cons(1,i)-cons_exact(1,i))+error_l1
	        error_l2	=   (cons(1,i)-cons_exact(1,i))**2+error_l2
	    
	    enddo
	    
	    error_l1	= error_l1/float(NX-1)
	    error_l2	= (error_l2/float(NX-1))**0.5
    
        write(*,*) '1D Poisson equation by hyperbolic approach, Roe - 5th order'
        write(*,*) '*********************************************'
        write(*,*) 'ø LInfinity Normal is:', error_inf
        write(*,*) 'ø L1        Normal is:', error_l1
        write(*,*) 'ø L2        Normal is:', error_l2

	    error_inf	= 0.0d0
	    error_l1	= 0.0d0
	    error_l2	= 0.0d0
	    
	    do i=1,NX

	        error_inf	=	MAX(error_inf,ABS(cons(2,i)-cons_exact(2,i)))
	        error_l1	=	ABS(cons(2,i)-cons_exact(2,i))+error_l1
	        error_l2	=	(cons(2,i)-cons_exact(2,i))**2+error_l2
	    
	    enddo
	    error_l1		=	error_l1/float(NX-1)
	    error_l2		=	(error_l2/float(NX-1))**0.5

        write(*,*) 'u LInfinity Normal is:', error_inf
        write(*,*) 'u L1        Normal is:', error_l1
        write(*,*) 'u L2        Normal is:', error_l2
        print*,'****************************************'

        call output(NX,x,cons,cons_exact,n_eqn,ghostp)


    
		call system_clock(count = time_end)
    
    	time_calc = time_end-time_ini
    
    	write(*,'(A20,I10,A)')   'Calculation time ',time_calc,' [CPU ticks]'

    	call cpu_time(finish)
    	write(*,*) " Total CPU time to solution = ", finish-start, " seconds"


		if (j .gt. NTMAX) then
		write(*,*) " Not converged... :( ..."
		write(*,*) "   Maximum iterations reached... NTMAX=", NTMAX
		write(*,*) "   Increase NTMAX, and try again."
		endif

		write(*,*) " Converged."
		write(*,*) " Final iteration      =", j-1
		write(*,*)
		write(*,*) "Finished Possion solver... See you later alligator!"
		write(*,*) 'Program ends...'


	end program hyperpoi3


	!***********************************************************************
	!*****                       Compute time step                	   *****
	!***********************************************************************

 	subroutine timestep(dx,dt,CFL,NX)
 	implicit none


 	integer 				:: i, j,NX
	double precision 		:: dx

	double precision 		:: dt, CFL

	double precision 		::  Lr, Tr, nu

	double precision, parameter 		:: pi=acos(-1.0d0)


		  nu 	= 	1.0d0
  		  Lr 	= 	1.0d0/(2.0d0*pi)
  		  Lr 	=   2.0d0/(pi*((pi/NX)+4))
  		  Tr 	= 	Lr**2
  		  CFL 	=	0.48d0
  		  dt 	= 	CFL*dx/sqrt(nu/Tr)

 	end subroutine timestep

	!***********************************************************************
	!*****                       Initial conditions                    *****
	!***********************************************************************


	subroutine initial(NX,cons,n_eqn)
		implicit none

		integer 							:: NX, n_eqn
		double precision 					:: cons(n_eqn,-5:NX+6)


		cons(1,:) =  000.0d0
		cons(2,:) =  000.0d0	

	end subroutine initial


	!***********************************************************************
	!*****                       Exact Solution	                       *****
	!***********************************************************************


	subroutine exact(NX,X,cons_exact,n_eqn)

	implicit none

		double precision, parameter 		:: pi=acos(-1.0d0)

		integer 		 					:: i, NX, n_eqn
		double precision 					:: X(1:NX)
		double precision					:: cons_exact(n_eqn,1:NX)
		double precision 					:: phianode, phicathode
		double precision					:: A, L, M
		

			A=10.0d0
			L=1.0d0
			M=3.0d0

		phianode		=	2.0d0
		phicathode		=	1.0d0
			
			do i=1,NX
				cons_exact(1,i)	=(((L/(2*pi*M))**2)*(-A)*cos(2*pi*M*x(i)/L)) &
							     +((x(i) * (phicathode- phianode))/L) + phianode  + ((L/(2*pi*M))**2)*A

				cons_exact(2,i) = (L/(2*pi*M)) * A * sin(2*pi*M*x(i)/L)   + ((phicathode- phianode)/L)

			enddo 

	end subroutine exact

	!***********************************************************************
	!*****                       Boundary conditions                   *****
	!***********************************************************************


	subroutine boundarycondition(NX,cons,n_eqn,ghostp)
	implicit none

	integer 							:: NX, n_eqn, ghostp
	double precision 					:: cons(n_eqn,-5:NX+6)
	double precision 					:: phi_anode,phi_cathode
		
		phi_anode		= 2.0d0
		phi_cathode		= 1.0d0

! Set left boundary condition
  

		cons(1,0)   = (128.0d0/35.0d0)*(phi_anode - (35.0d0/32.0d0)*cons(1,1) + (35.0d0/64.0d0)*cons(1,2)-&
					  (7.0d0/32.0d0)*cons(1,3) +(5.0d0/128.0d0)*cons(1,4))
		cons(2,0)    = 5.0d0 * cons(2,1) - 10.0d0*cons(2,2) + 10.0d0*cons(2,3) -5.0d0*cons(2,4) + cons(2,5)
		cons(:,-1)   = 5.0d0 * cons(:,0) - 10.0d0*cons(:,1) + 10.0d0*cons(:,2) -5.0d0*cons(:,3) + cons(:,4)
		cons(:,-2)   = 5.0d0 * cons(:,-1)-10.0d0*cons(:,0)  + 10.0d0*cons(:,1) -5.0d0*cons(:,2) + cons(:,3)
		cons(:,-3)   = 5.0d0 * cons(:,-2)-10.0d0*cons(:,-1) + 10.0d0*cons(:,0) -5.0d0*cons(:,1) + cons(:,2)
		! cons(:,-4)   = 5.0d0 * cons(:,-3)-10.0d0*cons(:,-2) + 10.0d0*cons(:,-1) -5.0d0*cons(:,0) + cons(:,1)
		! cons(:,-5)   = 5.0d0 * cons(:,-4)-10.0d0*cons(:,-3) + 10.0d0*cons(:,-2) -5.0d0*cons(:,-1) + cons(:,0)


! Set right boundary condition
       
		cons(1,NX+1) = (128.0d0/35.0d0)*(phi_cathode +(5.0d0/128.0d0)*cons(1,NX-3) - (7.0d0/32.0d0)*cons(1,NX-2) +&
					   (35.0d0/64.0d0)*cons(1,NX-1) -(35.0d0/32.d0)*cons(1,NX))
		cons(2,NX+1)  = cons(2,NX-4) - 5.0d0*cons(2,NX-3) + 10.0d0*cons(2,NX-2) -10.0d0*cons(2,NX-1)+5.0d0*cons(2,NX)
		cons(:,NX+2)  = cons(:,NX-3) - 5.0d0*cons(:,NX-2) + 10.0d0*cons(:,NX-1) -10.0d0*cons(:,NX)  +5.0d0*cons(:,NX+1)
		cons(:,NX+3)  = cons(:,NX-2) - 5.0d0*cons(:,NX-1) + 10.0d0*cons(:,NX)   -10.0d0*cons(:,NX+1)+5.0d0*cons(:,NX+2)
		cons(:,NX+4)  = cons(:,NX-1) - 5.0d0*cons(:,NX)   + 10.0d0*cons(:,NX+1) -10.0d0*cons(:,NX+2)+5.0d0*cons(:,NX+3)
		! cons(:,NX+5)  = cons(:,NX) - 5.0d0*cons(:,NX+1)   + 10.0d0*cons(:,NX+2) -10.0d0*cons(:,NX+3)+5.0d0*cons(:,NX+4)



			cons(1,0 )   = (8.0d0*phi_anode -6.0d0*cons(1,1) + cons(1,2)) * (1.0d0/3.0d0)
        
        cons(2,0 )   = cons(2,3)-3.0d0*cons(2,2)+3.0d0*cons(2,1)
        cons(:,-1)   = cons(:,2)-3.0d0*cons(:,1)+3.0d0*cons(:,0)
        cons(:,-2)   = cons(:,1)-3.0d0*cons(:,0)+3.0d0*cons(:,-1)
        cons(:,-3)   = cons(:,0)-3.0d0*cons(:,-1)+3.0d0*cons(:,-2)
        cons(:,-4)   = cons(:,-1)-3.0d0*cons(:,-2)+3.0d0*cons(:,-3)
       

! ! Set right boundary condition
       
! 		cons(1,NX+1) = (8.0d0*phi_cathode +1.0d0*cons(1,NX-1) -6.0d0* cons(1,NX)) * (1.0d0/3.0d0)
        
        cons(2,NX+1) = 3.0d0*cons(2,NX)   - 3.0d0*cons(2,NX-1)+cons(2,NX-2)
        cons(:,NX+2) = 3.0d0*cons(:,NX+1) - 3.0d0*cons(:,NX)  +cons(:,NX-1)
        cons(:,NX+3) = 3.0d0*cons(:,NX+2) - 3.0d0*cons(:,NX+1)  +cons(:,NX)
        cons(:,NX+4) = 3.0d0*cons(:,NX+3) - 3.0d0*cons(:,NX+2)  +cons(:,NX+1)
        cons(:,NX+5) = 3.0d0*cons(:,NX+4) - 3.0d0*cons(:,NX+3)  +cons(:,NX+2)

    end subroutine boundarycondition

 	

 	!***********************************************************************
	!*****                      Flux in X-direction                    *****
	!***********************************************************************


	subroutine FX (NX,dx,x,cons,residual,n_eqn,ghostp)
	implicit none

	integer								:: i,NX
	integer 							:: n_eqn, ghostp
	double precision, parameter 		:: pi=acos(-1.0d0)
	double precision 					:: dx, X(1:NX)

	double precision 					:: cons(n_eqn,-5:NX+6),residual(n_eqn,1:NX)

	double precision   					:: flux_x_half(n_eqn,-5:NX+6),flux_node_x(n_eqn,-5:NX+6)

	real*8 :: nu,Lr,Tr

		nu = 1.0d0
		Lr = 1.0d0/(2.0d0*acos(-1.0d0))
		Lr =2.0d0/(pi*((pi/NX)+4))

		Tr = Lr**2.0d0/nu


			do i = -ghostp,NX+1+ghostp

				flux_node_x(1,i) = -cons(2,i)
				flux_node_x(2,i) = -cons(1,i)

			enddo
	     
	call reconstruction3(NX,cons,flux_x_half,n_eqn,ghostp)


							
	    
	  	do i=1,NX
			
			!Finite volume or Reconstruction approach 
	        ! residual(:,i) = -((flux_x_half(:,i)-flux_x_half(:,i-1))/dx)

	         
	        
	        !Finite difference or interpolation approach 

! 4th order difference formula. Midpoint to node 	        

	        ! residual(:,i) = -((4.0/3.0d0)*(flux_x_half(:,i  )-flux_x_half(:,i-1)) +&
					    !      (-1.0/6.0d0)*(flux_node_x(:,i+1)-flux_node_x(:,i-1)))/dx

! 6th order difference formula. Midpoint to node 

			residual(:,i) =	-((3.0d0/2.0d0)*(flux_x_half(:,i)-flux_x_half(:,i-1))/dx + &
							(1.0d0/30.0d0)*(flux_x_half(:,i+1)-flux_x_half(:,i-2))/dx &
							+(-3.0d0/10.0d0)*(flux_node_x(:,i+1)-flux_node_x(:,i-1))/dx)

			! residual(:,i) = -((75.0d0/64.0d0)*(flux_x_half(:,i)-flux_x_half(:,i-1))/dx - &
			! 					(25.0d0/384.0d0)*(flux_x_half(:,i+1)-flux_x_half(:,i-2))/dx &
			! 	               + (3.0d0/640.0d0)*(flux_x_half(:,i+2)-flux_x_half(:,i-3))/dx)

   		

	    enddo

   

	end subroutine FX

	!***************************************************************************************************************
	!***** 					Roe, HLL or Rusanov + 5th order upwind reconstruction or interplation        	   *****
	!***************************************************************************************************************


	subroutine reconstruction3 (NX,cons,flux, n_eqn,ghostp)
		implicit none

		integer								:: i,NX,k
		integer								:: n_eqn, ghostp, total
		double precision, parameter 		:: pi=acos(-1.0d0)

		double precision 					:: cons(n_eqn,-5:NX+6)
		double precision 					:: flux(n_eqn,-5:NX+6)
		double precision 					:: consl(n_eqn,-5:NX+6), consr(n_eqn,-5:NX+6)

		double precision   					:: fleft(n_eqn,-5:NX+6), fright(n_eqn,-5:NX+6)

		double precision					:: righteigen(2,2), lefteigen(2,2),dissipation(n_eqn,-3:NX+3)
		double precision 					:: lambda1(-5:NX+6), lambda2(-5:NX+6)
		double precision  					:: delphi(-5:NX+6), delu(-5:NX+6)
		double precision 					:: alpha1, alpha2
		double precision 					:: SLm(-5:NX+6), SRp(-5:NX+6),laxmax,jacobian(2,2)

		double precision					:: a(0:NX), b(0:NX), c(0:NX), r(0:NX),temp(0:NX)


		real*8 :: nu,Lr,Tr

		nu = 1.0d0
		Lr = 1.0d0/(2.0d0*acos(-1.0d0))
		Lr =2.0d0/(pi*((pi/NX)+4))

		Tr = Lr**2.0d0/nu

		do k=1,2
        	do i =-3,NX+3

        	! Reconstruction formula
			! consl(k,i) = (2.0d0/60.0d0)*(cons(k,i-2)) + (-13.0d0/60.0d0)*(cons(k,i-1)) + (47.0d0/60.0d0)*(cons(k,i)) + &
			! 		  (27.0d0/60.0d0)*(cons(k,i+1)) + (-3.0d0/60.0d0)*(cons(k,i+2))

			! consr(k,i) = (-3.0d0/60.0d0)*(cons(k,i-1)) + (27.0d0/60.0d0)*(cons(k,i)) + (47.0d0/60.0d0)*(cons(k,i+1)) + &
			! 	      (-13.0d0/60.0d0)*(cons(k,i+2)) + (2.0d0/60.0d0)*(cons(k,i+3))

  	       ! Interpolation formula		
			consl(k,i) = (3.0d0/128.0d0)*(cons(k,i-2)) - (5.0d0/32.0d0)*(cons(k,i-1)) + (45.0d0/64.0d0)*(cons(k,i)) + &
					     (15.0d0/32.0d0)*(cons(k,i+1)) - (5.0d0/128.0d0)*(cons(k,i+2))
			
			consr(k,i) = (-5.0d0/128.0d0)*(cons(k,i-1)) + (15.0d0/32.0d0)*(cons(k,i)) + (45.0d0/64.0d0)*(cons(k,i+1)) - &
					     (5.0d0/32.0d0)*(cons(k,i+2)) + (3.0d0/128.0d0)*(cons(k,i+3))

			  ! consl(k,i) = -1.0d0/8.0d0*cons(k,i-1)+6.0d0/8.0d0*cons(k,i  )+3.0d0/8.0d0*cons(k,i+1)
     !          consr(k,i) =  3.0d0/8.0d0*cons(k,i  )+6.0d0/8.0d0*cons(k,i+1)-1.0d0/8.0d0*cons(k,i+2)

        
         	enddo
        enddo

	    do k = 1,n_eqn
	    	

	    do i = 0, NX
		      a(i) = 0.5d0
		      b(i) = 1.0d0
		      c(i) = 0.1d0
	   		 end do
	    
	    i = 0  ; a(i) = 0.0d+00
	    i = NX  ; c(i) = 0.0d+00
	    total = NX+1
	    
	    r(0)		= 0.1d0*cons(k, -1)+1.0d0*cons(k,0)+0.5d0*cons(k,  1)-0.5d0*consl(k, -1)
	    r(NX)		= 0.1d0*cons(k,NX-1)+1.0d0*cons(k,NX)+0.5d0*cons(k,NX+1)-0.1d0*consl(k,NX+1)
	    
	    do i=1,NX-1
	    r(i)		= 0.1d0*cons(k, i-1)+1.0d0*cons(k,i)+0.5d0*cons(k,  i+1)
	    enddo
	    
	    call tridiag(a,b,c,r,temp,total)
	        
	    consl(k,0:NX)=temp(0:NX) 

		    do i = 0, NX
		      a(i) = 0.1d0
		      b(i) = 1.0d0
		      c(i) = 0.5d0
		    end do

	    i = 0  ; a(i) = 0.0d+00
	    i = NX  ; c(i) = 0.0d+00

		r(0)		= 0.5d0*cons(k,0)+1.0d0*cons(k,1  )+0.1d0*cons(k,  2)-0.1d0*consr(k, -1)
	    r(NX)		= 0.5d0*cons(k,NX)+1.0d0*cons(k,NX+1)+0.1d0*cons(k,NX+2)-0.5d0*consr(k,NX+1)
	   
	    do i=1,NX-1
	    r(i)		= 0.5d0*cons(k, i)+1.0d0*cons(k,i+1)+0.1d0*cons(k,  i+2)
	    enddo
	    
	    call tridiag(a,b,c,r,temp,total)
	        
	    consr(k,0:NX)=temp(0:NX)

	    enddo


		do i=-3,NX+3

			fleft(1,i)	   = -consl(2,i)
			fleft(2,i)	   = -consl(1,i)

			fright(1,i)	   = -consr(2,i)
			fright(2,i)	   = -consr(1,i)


			delphi(i)			= consr(1,i)-consl(1,i)
			delu(i)				= consr(2,i)-consl(2,i)

	!***************************************
	!***** 		Roe solver	       	   *****
	!***************************************
		! Eigen vectors . We can also use the absolute Jacobian matrix
		righteigen(1,1) = 1.0d0
		righteigen(1,2) = 1.0d0
		righteigen(2,1) = 1.0d0
		righteigen(2,2) = -1.0d0


		lefteigen(1,1) = 0.50d0
		lefteigen(1,2) = 0.50d0
		lefteigen(2,1) = 0.50d0
		lefteigen(2,2) = -0.50d0

		alpha1 = lefteigen(1,1)*delphi(i) + lefteigen(1,2)*delu(i)
		alpha2 = lefteigen(2,1)*delphi(i) + lefteigen(2,2)*delu(i)

		lambda1(i) = abs(-1.0d0)
		lambda2(i) = abs(1.0d0)

		dissipation (1,i) = ((lambda1(i)*alpha1*righteigen(1,1)) + (lambda2(i)*alpha2*righteigen(1,2)))

		dissipation (2,i) = ((lambda1(i)*alpha1*righteigen(2,1)) + (lambda2(i)*alpha2*righteigen(2,2)))




		jacobian(1,1) = sqrt(1.0d0/Tr)
		jacobian(1,2) = 0.0d0
		jacobian(2,1) = 0.0d0
		jacobian(2,2) = sqrt(1.0d0/Tr)*Tr

		dissipation (1,i) = (delphi(i)*jacobian(1,1) + delu(i) *jacobian(1,2))

		dissipation (2,i) = (delphi(i)*jacobian(2,1) + delu(i) *jacobian(2,2))
		
		flux(1,i) =  0.5d0*((fright(1,i) +fleft(1,i) -  dissipation(1,i)))	
        flux(2,i) =  0.5d0*((fright(2,i) +fleft(2,i) -  dissipation(2,i)))	

    !***************************************
	!***** 		Rusanov solver	       *****
	!***************************************

        ! flux(:,i) =  0.5d0*((fright(:,i) +fleft(:,i) - (consr(:,i)-consl(:,i))))
    
    !***************************************
	!***** 		HLL solver	       *****
	!***************************************

   !      	SLm(i) 		   = -1.0d0
			! SRp(i) 		   =  1.0d0

			! SLm(i) 			= MIN(-1.0d0 , 1.0d0)
			! SRp(i) 			= MAX(-1.0d0 , 1.0d0)

			! laxmax         = MAX(-1,1)

			! if(0.0d0 .le. SLm(i)) then

				
			! 	flux(:,i) = fleft(:,i)
			
			! elseif(SLm(i) .le. 0.0d0 .and. 0.0d0 .le. SRp(i))	then
				
			! 	! flux(:,i) =  0.5d0*((fright(:,i) +fleft(:,i) + laxmax*(-consr(:,i)+consl(:,i))))

			! 	flux(:,i) = ( SRp(i)*fleft(:,i) - SLm(i)*fright(:,i) + SLm(i)*SRp(i)*(consr(:,i) - consl(:,i)) )/(SRp(i)-SLm(i))

			! elseif(0.0d0 .ge. SRp(i))then

			! 	flux(:,i) = fright(:,i)
			! endif
        

	enddo

	end subroutine reconstruction3

	!***********************************************************************
	!*****   Tridiagonal solver for compact scheme			           *****
	!***********************************************************************

	subroutine tridiag(l, d, u, r, x, n)

    
	! local variables are not implicit by default
	!
	    implicit none

	! input/output arguments
	!
	    integer                   , intent(in)  :: n
	    double precision, dimension(n), intent(in)  :: l, d, u, r
	    double precision, dimension(n), intent(out) :: x

	! local variables
	!
	    integer                    :: i, j
	    double precision               :: t
	    double precision, dimension(n) :: g
	!
	!-------------------------------------------------------------------------------
	!
	! decomposition and forward substitution
	!
	    t    = d(1)
	    x(1) = r(1) / t
	    do i = 2, n
	      j  = i - 1
	      g(i) = u(j) / t
	      t    = d(i) - l(i) * g(i)
	      if (t == 0.0d+00) then
	        write (*,*)"algebra::tridiag", "solution failed!"
	        stop
	      end if
	      x(i) = (r(i) - l(i) * x(j)) / t
	    end do

	! backsubstitution
	!
	    do i = n - 1, 1, -1
	      j    = i + 1
	      x(i) = x(i) - g(j) * x(j)
	    end do

	!-------------------------------------------------------------------------------
	!
  end subroutine tridiag


	!***********************************************************************
	!*****                       Source term 			               *****
	!***********************************************************************

	subroutine sourceterm(NX,X,residual,cons,n_eqn)

	implicit none
	double precision, parameter 		:: pi=acos(-1.0d0)
	double precision 					:: X(1:NX)
	
	integer								:: i, NX, n_eqn
	
	double precision					:: source(n_eqn,1:NX)	
	double precision 					:: cons(n_eqn,-5:NX+6),residual(n_eqn,1:NX)
	double precision					:: A, L, M

	real*8 :: nu,Lr,Tr

		nu = 1.0d0
		Lr = 1.0d0/(2.0d0*acos(-1.0d0))
		Lr = 2.0d0/(pi*((pi/NX)+4))

		Tr = Lr**2.0d0/nu



			A=10.0d0
			L=1.0d0
			M=3.0d0

			do i=1,NX


			source(1,i)		= (-A)*cos(2*pi*M*x(i)/L)
						! source(1,i) = -x(i)*x(i)*sin(1.0d0/(x(i)))
						! source(1,i) = -sqrt(x(i))
						! source(1,i) = -1.0d0/(x(i)+0.0)
						! source(1,i) = -sin((1/x(i)))
						! source(1,i) = (-A)*cos(2*pi*M*x(i)/L)*sin(x(i))
			source(2,i)		= -cons(2,i)
			residual(1,i)	= residual(1,i) + source(1,i)
			residual(2,i)	= (residual(2,i) + source(2,i))/Tr

			enddo

	end subroutine sourceterm

	

	!***********************************************************************
	!*****                       Output 			                   *****
	!***********************************************************************

	subroutine output(NX,x,cons,cons_exact,n_eqn,ghostp)

	integer 				:: i, j

	integer					:: NX, n_eqn, ghostp

	double precision 		:: X(1:NX),cons(n_eqn,-5:NX+6)
	double precision		:: cons_exact(n_eqn,1:NX)

	double precision		:: potential(1:NX), u(1:NX)
	double precision		:: potential_exact(1:NX), u_exact(1:NX)

	do i=1,NX
		
		potential(i)			= cons(1,i)
		u(i)					= cons(2,i)

		potential_exact(i)		= cons_exact(1,i)
		u_exact(i)				= cons_exact(2,i)

	enddo



	open(unit=9,file='Result-1b.plt',form='formatted',status='replace')
	
	! Write the data    
	    do i=1,NX
        	write(9,'(5es25.10)') x(i),potential(i), potential_exact(i), u_exact(i), u(i)
    	enddo
 	close(9)




 	end subroutine output


	!*******************************************************************************************************************************
	!*****                       We can also use reconstruction formula but not considered here 			                   *****
	!*******************************************************************************************************************************


!  	subroutine interpolations

! Compact reconstruction
!  	! do k = 1,2
    	

!     ! 	do i = 0, N
! 	   !    a(i) = 1.0d0/2.0d0
! 	   !    b(i) = 1.0d0
! 	   !    c(i) = 1.0d0/6.0d0
!    	! 	 end do
    
!     ! i = 0  ; a(i) = 0.0d+00
!     ! i = N  ; c(i) = 0.0d+00
!     ! total = N+1
    
!     ! r(0)		= (1.0d0/18.0d0)*un(k, -1)+(19.0d0/18.0d0)*un(k,0)+(10.0d0/18.0d0)*un(k,  1)-(1.0d0/2.0d0)*ul(k, -1)
!     ! r(N)		= (1.0d0/18.0d0)*un(k,N-1)+(19.0d0/18.0d0)*un(k,N)+(10.0d0/18.0d0)*un(k,N+1)-(1.0d0/6.0d0)*ul(k,N+1)
    
!     ! do i=1,N-1
!     ! r(i)		= (1.0d0/18.0d0)*un(k, i-1)+(19.0d0/18.0d0)*un(k,i)+(5.0d0/9.0d0)*un(k,  i+1)
!     ! enddo
    
!     ! call tridiag(a,b,c,r,temp,total)
        
!     ! ul(k,0:N)=temp(0:N) 

! 	   !  do i = 0, N
! 	   !    a(i) = 1.0d0/6.0d0
! 	   !    b(i) = 1.0d0
! 	   !    c(i) = 1.0d0/2.0d0
! 	   !  end do

!     ! i = 0  ; a(i) = 0.0d+00
!     ! i = N  ; c(i) = 0.0d+00

!     ! r(0)		= (10.0d0/18.0d0)*un(k,0)+(19.0d0/18.0d0)*un(k,  1)+(1.0d0/18.0d0)*un(k,  2)-(1.0d0/6.0d0)*ur(k, -1)
!     ! r(N)		= (10.0d0/18.0d0)*un(k,N)+(19.0d0/18.0d0)*un(k,N+1)+(1.0d0/18.0d0)*un(k,N+2)-(1.0d0/2.0d0)*ur(k,N+1)
   
!     ! do i=1,N-1
!     ! r(i)		= (10.0d0/18.0d0)*un(k,i)+(19.0d0/18.0d0)*un(k,i+1)+(1.0d0/18.0d0)*un(k,i+2)
!     ! enddo
    
!     ! call tridiag(a,b,c,r,temp,total)
        
!     ! ur(k,0:N)=temp(0:N)

!     ! enddo

! Compact reconstruction low dissipation
!     do k = 1,2
    	

!     	do i = 0, N
! 	      a(i) = 5.0d0/20.0d0
! 	      b(i) = 12.0d0/20.0d0
! 	      c(i) = 3.0d0/20.0d0
!    		 end do
    
!     i = 0  ; a(i) = 0.0d+00
!     i = N  ; c(i) = 0.0d+00
!     total = N+1
    
!     r(0)		= (3.0d0/120.0d0)*un(k, -1)+(67.0d0/120.0d0)*un(k,0)+(49.0d0/120.0d0)*un(k,  1)+(1.0d0/120.0d0)*un(k,  2)-(5.0d0/20.0d0)*ul(k, -1)
!     r(N)		= (3.0d0/120.0d0)*un(k,N-1)+(67.0d0/120.0d0)*un(k,N)+(49.0d0/120.0d0)*un(k,N+1)+(1.0d0/120.0d0)*un(k,  N+2)-(3.0d0/20.0d0)*ul(k,N+1)
    
!     do i=1,N-1
!     r(i)		= (3.0d0/120.0d0)*un(k,i-1)+(67.0d0/120.0d0)*un(k,i)+(49.0d0/120.0d0)*un(k,i+1)+(1.0d0/120.0d0)*un(k,  i+2)
!     enddo
    
!     call tridiag(a,b,c,r,temp,total)
        
!     ul(k,0:N)=temp(0:N) 

! 	    do i = 0, N
! 	      a(i) = 3.0d0/20.0d0
! 	      b(i) = 12.0d0/20.0d0
! 	      c(i) = 5.0d0/20.0d0
! 	    end do

!     i = 0  ; a(i) = 0.0d+00
!     i = N  ; c(i) = 0.0d+00

!     !  r(0)		= (1.0d0/120.0d0)*un(k, 0)+(49.0d0/120.0d0)*un(k,1)+(67.0d0/120.0d0)*un(k,  2)+(3.0d0/120.0d0)*un(k,  3)-(3.0d0/20.0d0)*ur(k, -1)
!     ! r(N)		= (1.0d0/120.0d0)*un(k,N)+(49.0d0/120.0d0)*un(k,N+1)+(67.0d0/120.0d0)*un(k,N+2)+(3.0d0/120.0d0)*un(k,  N+3)-(5.0d0/20.0d0)*ur(k,N+1)

!     r(0)		= (1.0d0/120.0d0)*un(k, -1)+(49.0d0/120.0d0)*un(k,0)+(67.0d0/120.0d0)*un(k,  1)+(3.0d0/120.0d0)*un(k,  2)-(3.0d0/20.0d0)*ur(k, -1)
!     r(N)		= (1.0d0/120.0d0)*un(k,N-1)+(49.0d0/120.0d0)*un(k,N)+(67.0d0/120.0d0)*un(k,N+1)+(3.0d0/120.0d0)*un(k,  N+2)-(5.0d0/20.0d0)*ur(k,N+1)
   
!     do i=1,N-1
!     r(i)		= (1.0d0/120.0d0)*un(k,i-1)+(49.0d0/120.0d0)*un(k,i)+(67.0d0/120.0d0)*un(k,i+1)+(3.0d0/120.0d0)*un(k,  i+2)
!     enddo
    
!     call tridiag(a,b,c,r,temp,total)
        
!     ur(k,0:N)=temp(0:N)

!     enddo


! ! Seventh order reconstruction

! !  		d1	=	(-1.0d0/204.0d0)
! ! 		d2	=	(+3.0d0/34.0d0)
! ! 		d3	=	(35.0d0/34.0d0)
! ! 		d4	=	(65.0d0/102.0d0)
! ! 		d5	=	(1.0d0/68.0d0)

! ! ! ! seventh order

! !      do k = 1,2
    	

! !     	do i = 0, N
! ! 	      a(i) = 9.0d0/17.0d0
! ! 	      b(i) = 1.0d0
! ! 	      c(i) = 4.0d0/17.0d0
! !    		 end do
    
! !     i = 0  ; a(i) = 0.0d+00
! !     i = N  ; c(i) = 0.0d+00
! !     total = N+1
    
! !     r(0)		= d1*un(k,-2)+d2*un(k, -1)+d3*un(k,0)+d4*un(k,  1)+d5*un(k,2)-(9.0d0/17.0d0)*ul(k, -1)
! !     r(N)		= d1*un(k,N-2)+d2*un(k, N-1)+d3*un(k,N)+d4*un(k,  N+1)+d5*un(k,N+2)-(4.0d0/17.0d0)*ul(k, N+1)
    
! !     do i=1,N-1
! !     r(i)		= d1*un(k,i-2)+d2*un(k, i-1)+d3*un(k,i)+d4*un(k,  i+1)+d5*un(k,i+2)
! !     enddo
    
! !     call tridiag(a,b,c,r,temp,total)
        
! !     ul(k,0:N)=temp(0:N) 

! ! 	    do i = 0, N
! ! 	      a(i) = 4.0d0/17.0d0
! ! 	      b(i) = 1.0d0
! ! 	      c(i) = 9.0d0/17.0d0
! ! 	    end do

! !     i = 0  ; a(i) = 0.0d+00
! !     i = N  ; c(i) = 0.0d+00

! ! 	! r(0)		= 0.5d0*un(k,0)+1.0d0*un(k,1  )+0.1d0*un(k,  2)-0.1d0*ur(k, -1)
! !  !    r(N)		= 0.5d0*un(k,N)+1.0d0*un(k,N+1)+0.1d0*un(k,N+2)-0.5d0*ur(k,N+1)

! !     r(0)		= d5*un(k,-1) +d4*un(k, 0)+d3*un(k,1)  +d2*un(k,  2)  +d1*un(k,3)-(4.0d0/17.0d0)*ur(k, -1)
! !     r(N)		= d5*un(k,N-1)+d4*un(k, N)+d3*un(k,N+1)+d2*un(k,  N+2)+d1*un(k,N+3)-(9.0d0/17.0d0)*ur(k, N+1)
   
! !     do i=1,N-1
! !     r(i)		= d5*un(k,i-1)+d4*un(k, i)+d3*un(k,i+1)+d2*un(k,  i+2)+d1*un(k,i+3)

! !     enddo
    
! !     call tridiag(a,b,c,r,temp,total)
        
! !     ur(k,0:N)=temp(0:N)



! !     enddo

! 	end subroutine interpolations


