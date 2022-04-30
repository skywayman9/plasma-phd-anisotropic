
!**********************************************************************
!    This program solves 2D Laplace equation by Hyperbolic approach
!	
!	 Written by Sainath Ch, Email: s.chamarthi@al.t.u-tokyo.ac.jp


	program hyperpoi2d
	
		implicit none


		double precision, parameter 		:: pi=acos(-1.0d0)

		integer, parameter					:: NTMAX=100000,l=1
		integer 							:: i, j, z, k, NT
		
		integer, parameter					:: NX = 16*l,NY = 16*l, n_eqn = 3, ghostp = 10 	

		double precision, parameter			:: l1_target = 1.0d-13, epsilon = 1.0d0-14
		double precision 					:: CFL, dt

		double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
		double precision 					:: dx, dy, xmin, xmax, ymin, ymax

		double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),cons_old(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)


		double precision 					:: cons_exact(n_eqn,1:NX,1:NY)

		double precision 					:: residual(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision 					:: flux_x(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),flux_y(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision					:: res_norm_initial(n_eqn,3), res_norm(n_eqn,3)

		integer 							:: L1 = 1, L2 = 2, Linf = 3, i_Linf, i_Linf_initial

		double precision					:: dnorm(n_eqn),iter_diff,denominator

		double precision 					:: error_inf, error_l1, error_l2

		integer 							:: time_ini, time_end, time_calc

		double precision					:: start, finish


		common /domain/ xmin, xmax, ymin,ymax
		common /grid/ dx, dy

		call cpu_time(start)
		write(*,*) 'Program start...'

		call system_clock(count = time_ini)


		! Domain. Keep it like this so that for other problems (4 contacts or R-T) it will be easier to modify the domain

		xmin = 0.0d0
		xmax = 1.0d0

		ymin =  0.0d0
		ymax =  1.0d0

		! Generate simple grid

	
	    dx = (xmax - xmin)/(NX)
		dy = (ymax - ymin)/(NY)


			do i=-ghostp,NX+ghostp
	        	X(i)=(i-0.5d0)*dx+xmin
	        	! write(*,*) x(i),i
	    	enddo
! 	    	write(*,*) dx

! pause
	    	do j=-ghostp,NY+ghostp
	        	y(j)=(j-0.5d0)*dy+ymin
	    	enddo


	!***********************************************************************
	!*****                       call initial conditions               *****
	!***********************************************************************


		call initialconditions(NX,NY,ghostp,cons,n_eqn)

		call exactsolution(NX,NY,ghostp,x,y,cons_exact,n_eqn)

		call timestep(NX,NY,dt,dx,dy,CFL,x,y,ghostp,cons,n_eqn)


	!***********************************************************************
	!*****         			 TVD- Runge Kutta 				       	   *****
	!***********************************************************************
		
	iteration : do NT=1,NTMAX

		call timestep(NX,NY,dt,dx,dy,CFL,x,y,ghostp,cons,n_eqn)
		cons_old = cons
		
	!***********************************************************************
	!*****                       Time step -1	                       *****
	!***********************************************************************

			call boundarycondition(NX,NY,cons,n_eqn,ghostp,x,y)
			
			call FX(NX,NY,dx,x,y,cons,flux_x,n_eqn,ghostp)
			
			call GY(NX,NY,dy,x,y,cons,flux_y,n_eqn,ghostp)

			call sourceterm(NX,NY,x,y,cons,flux_x,flux_y,residual,n_eqn,ghostp)
						   


			do k = 1,n_eqn

				do j = 1, NY

					do i = 1, NX

						cons(k,i,j) = cons_old(k,i,j) + dt* residual(k,i,j)
					
					enddo
				
				enddo

			enddo

	!***********************************************************************
	!*****                       Time step - 2	                       *****
	!***********************************************************************	

			call boundarycondition(NX,NY,cons,n_eqn,ghostp,x,y)
			
			call FX(NX,NY,dx,x,y,cons,flux_x,n_eqn,ghostp)
			
			call GY(NX,NY,dy,x,y,cons,flux_y,n_eqn,ghostp)

			call sourceterm(NX,NY,x,y,cons,flux_x,flux_y,residual,n_eqn,ghostp)


			do k = 1,n_eqn

				do j = 1, NY

					do i = 1, NX

						cons(k,i,j) = (3.0d0/4.0d0)*cons_old(k,i,j) + (1.0d0/4.0d0)*( cons(k,i,j) + dt* residual(k,i,j) )
					
					enddo
				
				enddo

			enddo

	!***********************************************************************
	!*****                       Time step -3	                       *****
	!***********************************************************************

			call boundarycondition(NX,NY,cons,n_eqn,ghostp,x,y)
			
			call FX(NX,NY,dx,x,y,cons,flux_x,n_eqn,ghostp)
			
			call GY(NX,NY,dy,x,y,cons,flux_y,n_eqn,ghostp)

			call sourceterm(NX,NY,x,y,cons,flux_x,flux_y,residual,n_eqn,ghostp)


			do k = 1,n_eqn

				do j = 1, NY

					do i = 1, NX

						cons(k,i,j) = (1.0d0/3.0d0)*cons_old(k,i,j) + (2.0d0/3.0d0)*( cons(k,i,j) + dt* residual(k,i,j) )
					
					enddo
				
				enddo

			enddo

		!***************************************************************************************************
		!*****                      ! Monitor residual norms to check steady convergence.			   *****
		!***** 	 I.e., see how well the numerical solution satisfies the target discrete equaton.      *****            
		!***************************************************************************************************

		 res_norm(:,L1  ) =  0.0d0
		 res_norm(:,L2  ) =  0.0d0
		 res_norm(:,Linf) =  -1.0d0

		do k = 1, n_eqn ! Three equations
			
			do j = 1,NY

				do i = 1, NX
		  
			   
			   	   res_norm(k,L1) = res_norm(k,L1) + abs(residual(k,i,j))
				   res_norm(k,L2) = res_norm(k,L2) + abs(residual(k,i,j))**2

				   	if (abs(residual(k,i-1,j-1)) >  res_norm(k,Linf) ) then
				                   i_Linf = i
				         res_norm(k,Linf) = abs(residual(k,i,j))
				   endif

		  		end do
			end do

		end do

		  res_norm(:,L1) =       res_norm(:,L1)/real(NX)/real(NY)
		  res_norm(:,L2) = sqrt( res_norm(:,L2)/real(NX)/real(NY) )

		error_inf	= 0.0d0
	    error_l1	= 0.0d0
	    error_l2	= 0.0d0
	    
	    do j =1, NY
	    	
	    	do i=1,NX

		        error_inf	=	MAX(error_inf,ABS(cons(1,i,j)-cons_exact(1,i,j)))
		        error_l1	=   ABS(cons(1,i,j)-cons_exact(1,i,j))+error_l1
		        error_l2	=   (cons(1,i,j)-cons_exact(1,i,j))**2+error_l2
	    
	    	end do

	    end do	
	    
	    error_l1	= error_l1/float(NX-1)/float(NY-1)
	    error_l2	= (error_l2/float(NX-1)/float(NY-1))**0.5

	    

		 !Stop when the residual norm is reduced by 12 orders of magnitude in the L1 norm:

		 if ( minval( res_norm(:,L1)/res_norm_initial(:,L1) )  < l1_target) then
		  	
		  	exit iteration

		 endif

		if (NT==1) then 

			  res_norm_initial = res_norm

	write(*,*) " Initial residual (after 1 iteration):"
	write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (L1) = ", res_norm(:,L1)
	write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (L2) = ", res_norm(:,L2)
	write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (Li) = ", res_norm(:,Linf)
			  
			  i_Linf_initial = i_Linf

		 endif

		 if (mod(NT,1000)==0) then
	write(*,*) NT
	write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (L1) = ", res_norm(:,L1  )
	write(*,'(a22,4es12.4,4es12.4,4es12.2)') "  Residual norm (L2) = ", res_norm(:,L2  )
	write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (Li) = ", res_norm(:,Linf)
	write(*,*)
  	write(*,*) error_l2, MINVAL(cons(1,1:NX,1:NY))
  	write(*,*)
		  call output(NX,NY,x,y,cons,cons_exact,n_eqn,ghostp)
		 endif


	!************************************************************
	!*****          Monitor normalized differences			*****
	!************************************************************


		 	dnorm=0.0d0 
			iter_diff=0.0d0
			denominator=0.0d0
			
		do j = 1,NY

			do i = 1, NX

				do k=1,n_eqn

					denominator=denominator + (cons_old(k,i,j)*cons_old(k,i,j))

					iter_diff= iter_diff + (cons(k,i,j)-cons_old(k,i,j))**2

					dnorm(k)=(iter_diff/(denominator+epsilon))

					dnorm(k)= (dnorm(k)/NX/NY)**0.5d0

				enddo

			enddo
		enddo	




		enddo iteration


		write(*,*)
		write(*,*) " Final residual norm:"
		write(*,'(a15,4es12.4,4es12.4,4es12.4)') "  Residual norm (L1) = ", res_norm(:,L1)
		write(*,'(a15,4es12.4,4es12.4,4es12.4)') "  Residual norm (L2) = ", res_norm(:,L2)
		write(*,'(a15,4es12.4,4es12.4,4es12.4)') "  Residual norm (Li) = ", res_norm(:,Linf)
		write(*,*)
	



	!***************************************************************************************************
	!*****                Error w.r.t exact Solution             	   							   *****
	!***************************************************************************************************


		error_inf	= 0.0d0
	    error_l1	= 0.0d0
	    error_l2	= 0.0d0
	    
	    do j =1, NY
	    	
	    	do i=1,NX

		        error_inf	=	MAX(error_inf,ABS(cons(1,i,j)-cons_exact(1,i,j)))
		        error_l1	=   ABS(cons(1,i,j)-cons_exact(1,i,j))+error_l1
		        error_l2	=   (cons(1,i,j)-cons_exact(1,i,j))**2+error_l2
	    
	    	end do

	    end do	
	    
	    error_l1	= error_l1/float(NX)/float(NY)
	    error_l2	= (error_l2/float(NX)/float(NY))**0.5
    	
    	
    	write(*,*) 'Non-linear diffusion equation'
    	write(*,*)
        write(*,*) 'Order of accuracy for Solution variable - U'
        write(*,*) '*********************************************'
        write(*,*) 'ø LInfinity Normal is:', error_inf
        write(*,*) 'ø L1        Normal is:', error_l1
        write(*,*) 'ø L2        Normal is:', error_l2
        write(*,*) '*********************************************'

       	error_inf	= 0.0d0
	    error_l1	= 0.0d0
	    error_l2	= 0.0d0
	    
	    do j =1, NY
	    	
	    	do i=1,NX

		        error_inf	=	MAX(error_inf,ABS(cons(2,i,j)-cons_exact(2,i,j)))
		        error_l1	=   ABS(cons(2,i,j)-cons_exact(2,i,j))+error_l1
		        error_l2	=   (cons(2,i,j)-cons_exact(2,i,j))**2+error_l2
	    
	    	enddo
	    enddo	
	    error_l1	= error_l1/float(NX)/float(NY)
	    error_l2	= (error_l2/float(NX)/float(NY))**0.5

	    write(*,*) 'Gradient variable - p'
        write(*,*) '*********************************************'
        write(*,*) 'ø LInfinity Normal is:', error_inf
        write(*,*) 'ø L1        Normal is:', error_l1
        write(*,*) 'ø L2        Normal is:', error_l2
        write(*,*) '*********************************************'

       	error_inf	= 0.0d0
	    error_l1	= 0.0d0
	    error_l2	= 0.0d0
	    
	    do j =1, NY
	    	
	    	do i=1,NX

		        error_inf	=	MAX(error_inf,ABS(cons(3,i,j)-cons_exact(3,i,j)))
		        error_l1	=   ABS(cons(3,i,j)-cons_exact(3,i,j))+error_l1
		        error_l2	=   (cons(3,i,j)-cons_exact(3,i,j))**2+error_l2
	    
	    	enddo
	    enddo	
	    error_l1	= error_l1/float(NX)/float(NY)
	    error_l2	= (error_l2/float(NX)/float(NY))**0.5

	    write(*,*) 'Gradient variable - q'
        write(*,*) '*********************************************'
        write(*,*) 'ø LInfinity Normal is:', error_inf
        write(*,*) 'ø L1        Normal is:', error_l1
        write(*,*) 'ø L2        Normal is:', error_l2
        write(*,*) '*********************************************'
	    


        call output(NX,NY,x,y,cons,cons_exact,n_eqn,ghostp)


		call system_clock(count = time_end)
    
    	time_calc = time_end-time_ini
    
    	write(*,'(A20,I10,A)')   'Calculation time ',time_calc,' [CPU ticks]'

    	call cpu_time(finish)
    	write(*,*) " Total CPU time to solution = ", finish-start, " seconds"


		if (j .gt. NTMAX) then
		write(*,*) "   Not converged... =                        ..."
		write(*,*) "   Maximum iterations reached... NTMAX=", NTMAX
		write(*,*) "   Increase NTMAX, and try again or write implicit code "
		else

		write(*,*) " Converged                     "
		write(*,*) " Final iteration        =", NT-1
		write(*,*)
		write(*,*) "Finished Possion solver... See you later alligator! "
		endif
		write(*,*) 'Program ends...'

		write(*,*)MINVAL(cons(1,1:NX,1:NY))




	end program hyperpoi2d


	!***********************************************************************
	!*****                       Initial conditions                    *****
	!***********************************************************************


	subroutine initialconditions(NX,NY,ghostp,cons,n_eqn)


	implicit none

	integer 							:: i, j, k
	integer								:: NX, NY, n_eqn, ghostp
	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision, parameter 		:: pi=acos(-1.0d0)


	do k =1,n_eqn

		do i = -ghostp, NX+ghostp
			
			do j =  -ghostp, NY+ghostp

				cons(k,i,j)		= 0.0d0

			enddo	
		enddo
	enddo

	end subroutine initialconditions


	!***********************************************************************
	!*****                       Exact Solution	                       *****
	!***********************************************************************


	subroutine exactsolution(NX,NY,ghostp,x,y,cons_exact,n_eqn)

	implicit none

	integer								:: i, j
	integer								:: NX,NY, n_eqn, ghostp

	double precision, parameter 		:: pi=acos(-1.0d0)

	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
	double precision 					:: cons_exact(n_eqn,1:NX,1:NY)
	double precision 					:: dx, dy, xmin, xmax, ymin, ymax, L


		L = 1.0d0

		cons_exact = 0.0d0


		do i = 1, NY
			
			do j = 1, NX

				cons_exact(1,i,j) = (sin(pi*(x(i)))*sin(pi*y(j)))
				cons_exact(2,i,j) = pi*(cos(pi*(x(i)))*sin(pi*y(j)))
				cons_exact(3,i,j) = pi*(sin(pi*(x(i)))*cos(pi*y(j)))
			
			enddo
		enddo

		
	end subroutine exactsolution


	!***********************************************************************
	!*****                       Compute time step                	   *****
	!***********************************************************************

 	subroutine timestep(NX,NY,dt,dx,dy,CFL,x,y,ghostp,cons,n_eqn)
 	implicit none

 	integer 							:: i, j,n_eqn
 	
 	integer								:: NX, NY,ghostp
	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

 	double precision 					:: dx, dy
	double precision 					:: dt, CFL, dtnew
	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision, parameter 		:: pi=acos(-1.0d0)

	double precision					:: x_velocity, y_velocity
	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx, d_xy,d_yx,d_yy

		dt = 1.0d10;
		
		do i = 1, NX

			do j = 1, NY

					d_xx = 3.0d0+cons(1,i,j)**(5.0d0/1.0d0)+3.0d0*cons(1,i,j)**2.0d0

				! d_xx = sin(3.14159265358979d0*x(i))**2*sin(3.14159265358979d0*y(j))**2 &
    !   + sin(3.14159265358979d0*x(i))*sin(3.14159265358979d0*y(j)) + 1

					Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
					Tr = (1.0d0*Lr**2.0d0)/(d_xx)

				  x_velocity =  DABS(sqrt(d_xx/Tr))
				  y_velocity =  DABS(sqrt(d_xx/Tr))

				dtnew = min(dx/x_velocity, dy/y_velocity)
				if(dtnew .lt. dt) dt = dtnew

				
			enddo
		enddo


  		  CFL 	=	0.99d0
  		  dt 	= 	CFL*dt


 	end subroutine timestep

 	!***********************************************************************
	!*****                       Boundary conditions                   *****
	!***********************************************************************

	subroutine boundarycondition(NX,NY,cons,n_eqn,ghostp,x,y)

	implicit none

	integer 							:: i, j, k
	double precision, parameter 		:: pi=acos(-1.0d0)
	integer								:: NX, NY, n_eqn, ghostp

	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	
	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)



    ! ! Left boundary conditions

	do j = 1,NY
		
		i=0     
	    
	    cons(1,i,j )   = (8.0d0*(dsin((x(i)+x(i+1))*0.5d0*pi)*dsin(y(j)*pi)) - 6.0d0*cons(1,i+1,j) + cons(1,i+2,j)) * (1.0d0/3.0d0)

	    cons(2,i,j )   = cons(2,i+3,j)-3.0d0*cons(2,i+2,j)+3.0d0*cons(2,i+1,j)
        
        cons(3,i,j )   = cons(3,i+3,j)-3.0d0*cons(3,i+2,j)+3.0d0*cons(3,i+1,j)
    enddo

    do j = 1,NY

    	do i = -1,-ghostp,-1
        
        	cons(:,i,j)   = cons(:,i+3,j)-3.0d0*cons(:,i+2,j)+3.0d0*cons(:,i+1,j)

        enddo
    enddo


	! Right boundary conditions


	do j = 1, NY
		
		i = NX

		cons(1,i+1,j) = (8.0d0*(dsin((x(i)+x(i+1))*0.5d0*pi)*dsin(y(j)*pi)) +1.0d0*cons(1,i-1,j) -6.0d0* cons(1,i,j)) * (1.0d0/3.0d0)
        
        cons(2,i+1,j) = 3.0d0*cons(2,i,j)   - 3.0d0*cons(2,i-1,j)+cons(2,i-2,j)
        cons(3,i+1,j) = 3.0d0*cons(3,i,j)   - 3.0d0*cons(3,i-1,j)+cons(3,i-2,j)
        
    enddo

    do j = 1, NY

    	do i= NX, NX+4

        cons(:,i+2,j) = 3.0d0*cons(:,i+1,j) - 3.0d0*cons(:,i,j)  +cons(:,i-1,j)

        enddo
    enddo


	! Bottom boundary conditions


	do i = 1,NX
		j=0     
	    cons(1,i,j )   = (8.0d0*((dsin(x(i)*pi)*dsin((y(j)+y(j+1))*0.5d0*pi))) - 6.0d0*cons(1,i,j+1) + cons(1,i,j+2)) * (1.0d0/3.0d0)

	    cons(2,i,j )   = cons(2,i,j+3)-3.0d0*cons(2,i,j+2)+3.0d0*cons(2,i,j+1)
        cons(3,i,j )   = cons(3,i,j+3)-3.0d0*cons(3,i,j+2)+3.0d0*cons(3,i,j+1)
    enddo

    do i = 1,NX  
    	do j = -1,-ghostp,-1
        
        	cons(:,i,j)   = cons(:,i,j+3)-3.0d0*cons(:,i,j+2)+3.0d0*cons(:,i,j+1)

        enddo
    enddo


	! Top boundary conditions
	do i = 1, NX
		 j = NY

		cons(1,i,j+1) = (8.0d0*((dsin(x(i)*pi)*dsin((y(j)+y(j+1))*0.5d0*pi))) +1.0d0*cons(1,i,j-1) -6.0d0* cons(1,i,j)) * (1.0d0/3.0d0)
        
        cons(2,i,j+1) = 3.0d0*cons(2,i,j)   - 3.0d0*cons(2,i,j-1)+cons(2,i,j-2)
        cons(3,i,j+1) = 3.0d0*cons(3,i,j)   - 3.0d0*cons(3,i,j-1)+cons(3,i,j-2)
        
    enddo

    do i = 1, NX

    	do j= NY, NY+4

        cons(:,i,j+2) = 3.0d0*cons(:,i,j+1) - 3.0d0*cons(:,i,j)  +cons(:,i,j-1)

        enddo
    enddo


	

	end subroutine boundarycondition



	!***********************************************************************
	!*****                      Flux in X-direction                    *****
	!***********************************************************************

	subroutine FX(NX,NY,dx,x,y,cons,flux_x,n_eqn,ghostp)

	implicit none

	double precision, parameter 		:: pi=acos(-1.0d0)
	integer								:: i, j ,k
	integer 							:: NX, NY, ghostp, n_eqn
	double precision					:: dx, x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision 					:: flux_x(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision					:: flux_node_x(n_eqn,-ghostp:NX+ghostp), flux_x_half(n_eqn,-ghostp:NX+ghostp)

	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),cons_local(n_eqn,-ghostp:NX+ghostp)

	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx(-ghostp:NX+ghostp), d_xy(-ghostp:NX+ghostp),d_yx(-ghostp:NX+ghostp),d_yy(-ghostp:NX+ghostp),beta




	do j =1, NY

		do i = -ghostp, NX+ghostp

			cons_local(1,i) = cons(1,i,j)
			cons_local(2,i) = cons(2,i,j)
			cons_local(3,i) = cons(3,i,j)

			d_xx(i) = 3.0d0+cons_local(1,i)**(5.0d0/1.0d0)+3.0d0*cons_local(1,i)**2.0d0

			! d_xx(i) = sin(3.14159265358979d0*x(i))**2*sin(3.14159265358979d0*y(j))**2 &
   !    + sin(3.14159265358979d0*x(i))*sin(3.14159265358979d0*y(j)) + 1

			flux_node_x(1,i) = -cons_local(2,i)*d_xx(i)
			flux_node_x(2,i) = -cons_local(1,i)
			flux_node_x(3,i) = 0.0d0
		end do

		call reconstruction5EX (NX,NY,flux_node_x,cons_local,flux_x_half, n_eqn,ghostp,x,y,j,d_xx,d_xy,d_yx,d_yy)

		do k = 1,n_eqn

			do i = 1, NX

				!Finite volume or Reconstruction approach 
	        	! flux_x(k,i,j) = -((flux_x_half(k,i)-flux_x_half(k,i-1))/dx)


				flux_x(k,i,j)  = -((9.0d0/8.0d0)*(flux_x_half(k,i)-flux_x_half(k,i-1))/dx &
					             - (1.0d0/24.0d0)*(flux_x_half(k,i+1)-flux_x_half(k,i-2))/dx)


			enddo
		enddo
	enddo

	end subroutine FX
			

	!***************************************************************************************************************
	!***** 					Roe, HLL or Rusanov + 5th order compact upwind reconstruction or interplation     *****
	!*****											in X-direction											   *****
	!***************************************************************************************************************

	subroutine reconstruction5EX (NX,NY,flux_node,cons_local,flux, n_eqn,ghostp,x,y,j,d_xx,d_xy,d_yx,d_yy)

	implicit none 

	integer, intent(in)					:: j
	integer								:: i, k	
	integer								:: NX,NY, n_eqn, ghostp, total
	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision					:: cons_local(n_eqn,-ghostp:NX+ghostp),flux(n_eqn,-ghostp:NX+ghostp),flux_node(n_eqn,-ghostp:NX+ghostp)

	double precision 					:: consl(n_eqn,-ghostp:NX+ghostp), consr(n_eqn,-ghostp:NX+ghostp)

	double precision   					:: fleft(n_eqn,-ghostp:NX+ghostp), fright(n_eqn,-ghostp:NX+ghostp)

	double precision					:: dissipation(n_eqn,-3:NX+3),jacobian(3,3)

	double precision  					:: delphi(-ghostp:NX+ghostp), delp(-ghostp:NX+ghostp), delq(-ghostp:NX+ghostp)


	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx(-ghostp:NX+ghostp), d_xy(-ghostp:NX+ghostp)
	double precision 					:: d_yy(-ghostp:NY+ghostp), d_yx(-ghostp:NY+ghostp)
	
	double precision 					:: d_xxl(-ghostp:NX+ghostp), d_xxr(-ghostp:NX+ghostp)
	double precision 					:: d_xyl(-ghostp:NX+ghostp), d_xyr(-ghostp:NX+ghostp)

	double precision 					:: d_yxl(-ghostp:NY+ghostp),d_yxr(-ghostp:NY+ghostp)
	double precision 					:: d_yyl(-ghostp:NY+ghostp),d_yyr(-ghostp:NY+ghostp)


		do k=1,n_eqn

        	do i =-2,NX+2

              consl(k,i) = -1.0d0/8.0d0*cons_local(k,i-1)+6.0d0/8.0d0*cons_local(k,i  )+3.0d0/8.0d0*cons_local(k,i+1)
              consr(k,i) =  3.0d0/8.0d0*cons_local(k,i  )+6.0d0/8.0d0*cons_local(k,i+1)-1.0d0/8.0d0*cons_local(k,i+2)

	         	! ! Reconstruction formula
           !    fleft(k,i) =  -1.0d0/8.0d0*flux_node(k,i-1)+6.0d0/8.0d0*flux_node(k,i  )+3.0d0/8.0d0*flux_node(k,i+1)
           !    fright(k,i) =  3.0d0/8.0d0*flux_node(k,i  )+6.0d0/8.0d0*flux_node(k,i+1)-1.0d0/8.0d0*flux_node(k,i+2)

            d_xxl(i) = (-1.0d0/8.0d0)*(d_xx(i-1)) + (6.0d0/8.0d00)*(d_xx(i)) + (3.0d0/8.0d0)*(d_xx(i+1))
			
			d_xxr(i) = (3.0d0/8.0d0)*(d_xx(i)) + (6.0d0/8.0d00)*(d_xx(i+1)) - (1.0d0/8.0d0)*(d_xx(i+2))
       
         	enddo
        enddo

    
		do i=-2,NX+2

				Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
				Tr = (1.0d0*Lr**2.0d0)/( ((d_xxl(i)+d_xxr(i))*0.5d0) )

			  fleft (1,i)  =  -consl(2,i)*d_xxl(i)
			  fright(1,i)  =  -consr(2,i)*d_xxr(i)

			  fleft (2,i)  =  -consl(1,i)
			  fright(2,i)  =  -consr(1,i)

			  fleft (3,i)  = 0.0d0
			  fright(3,i)  = 0.0d0




			delphi(i)			= consr(1,i)-consl(1,i)
			delp(i)				= consr(2,i)-consl(2,i)
			delq(i)				= consr(3,i)-consl(3,i)

	!***************************************
	!***** 		Roe solver	       	   *****
	!***************************************

		 !***** Roe solver          *****
		 !The upwind dissipation matrix is |A| = R*|Lambda|*L
		 jacobian(1,1) = abs(sqrt(((d_xxl(i)+d_xxr(i))*0.5d0)/Tr))
		 jacobian(1,2) = 0.0d0
		 jacobian(1,3) = 0.0d0
		 jacobian(2,1) = 0.0d0
		 jacobian(2,2) = Tr*abs(sqrt(((d_xxl(i)+d_xxr(i))*0.5d0)/Tr))
		 jacobian(2,3) = 0.0d0
		 jacobian(3,1) = 0.0d0
		 jacobian(3,2) = 0.0d0
		 jacobian(3,3) = 0.0d0

		dissipation (1,i) = (delphi(i)*jacobian(1,1) + delp(i) *jacobian(1,2)+delq(i)*jacobian(1,3))

		dissipation (2,i) = (delphi(i)*jacobian(2,1) + delp(i) *jacobian(2,2)+delq(i)*jacobian(2,3))

		dissipation (3,i) = (delphi(i)*jacobian(3,1) + delp(i) *jacobian(3,2)+delq(i)*jacobian(3,3))

		
		
		flux(1,i) =  0.5d0*((fright(1,i) +fleft(1,i) -  dissipation(1,i)))	
        flux(2,i) =  0.5d0*((fright(2,i) +fleft(2,i) -  dissipation(2,i)))	
        flux(3,i) =  0.5d0*((fright(3,i) +fleft(3,i) -  dissipation(3,i)))
        

	enddo



	end subroutine reconstruction5EX


 	!***********************************************************************
	!*****                      Flux in Y-direction                    *****
	!***********************************************************************

	subroutine GY(NX,NY,dy,x,y,cons,flux_y,n_eqn,ghostp)

	implicit none
	double precision, parameter 		:: pi=acos(-1.0d0)
	integer								:: i, j, k

	integer 							:: NX, NY, ghostp, n_eqn

	double precision					:: dy, x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision 					:: flux_y(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	
	double precision					:: flux_node_y(n_eqn,-ghostp:NY+ghostp), flux_y_half(n_eqn,-ghostp:NY+ghostp)

	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),cons_local(n_eqn,-ghostp:NY+ghostp)

	double precision 					:: nu,Lr,Tr
	 double precision 					:: d_xx(-ghostp:NX+ghostp), d_xy(-ghostp:NX+ghostp),d_yx(-ghostp:NX+ghostp),d_yy(-ghostp:NX+ghostp),beta




	
	do j =1, NX
		
		do i = -ghostp, NY+ghostp

			cons_local(1,i) = cons(1,j,i)
			cons_local(2,i) = cons(2,j,i)
			cons_local(3,i) = cons(3,j,i)


			d_xx(i) = 3.0d0+cons_local(1,i)**(5.0d0/1.0d0)+3.0d0*cons_local(1,i)**2.0d0

			! d_xx(i) = sin(3.14159265358979d0*x(j))**2*sin(3.14159265358979d0*y(i))**2 &
   !    + sin(3.14159265358979d0*x(j))*sin(3.14159265358979d0*y(i)) + 1


			flux_node_y(1,i) =  -cons_local(3,i)*d_xx(i)
			flux_node_y(2,i) =  0.0d0
			flux_node_y(3,i) =  -cons_local(1,i)
		end do

		call reconstruction5EY (NX,NY,flux_node_y,cons_local,flux_y_half, n_eqn,ghostp,x,y,j,d_xx,d_xy,d_yx,d_yy)

		do k = 1,n_eqn

			do i = 1, NY

				!Finite volume or Reconstruction approach 
	        	! flux_y(k,j,i) = -((flux_y_half(k,i)-flux_y_half(k,i-1))/dy)

	        	flux_y(k,j,i)  = -((9.0d0/8.0d0)*(flux_y_half(k,i)-flux_y_half(k,i-1))/dy &
				                - (1.0d0/24.0d0)*(flux_y_half(k,i+1)-flux_y_half(k,i-2))/dy)
			enddo
		enddo
	enddo
	





	end subroutine GY


	!***************************************************************************************************************
	!***** 					Roe, HLL or Rusanov + 5th order compact upwind reconstruction or interplation     *****
	!*****											in Y-direction											   *****
	!***************************************************************************************************************

	subroutine reconstruction5EY (NX,NY,flux_node,cons_local,flux, n_eqn,ghostp,x,y,j,d_xx,d_xy,d_yx,d_yy)

	implicit none 		

	integer, intent(in)					:: j
	integer								:: i,k	
	integer								:: NX,NY, n_eqn, ghostp, total
	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision					:: cons_local(n_eqn,-ghostp:NY+ghostp),flux(n_eqn,-ghostp:NY+ghostp),flux_node(n_eqn,-ghostp:NY+ghostp)

	double precision 					:: consl(n_eqn,-ghostp:NY+ghostp), consr(n_eqn,-ghostp:NY+ghostp)

	double precision   					:: fleft(n_eqn,-ghostp:NY+ghostp), fright(n_eqn,-ghostp:NY+ghostp)

	double precision					:: dissipation(n_eqn,-3:NY+3),jacobian(3,3)

	double precision  					:: delphi(-ghostp:NY+ghostp), delp(-ghostp:NY+ghostp), delq(-ghostp:NY+ghostp)


	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx(-ghostp:NX+ghostp), d_xy(-ghostp:NX+ghostp)
	double precision 					:: d_yy(-ghostp:NY+ghostp), d_yx(-ghostp:NY+ghostp)
	
	double precision 					:: d_xxl(-ghostp:NX+ghostp), d_xxr(-ghostp:NX+ghostp)
	double precision 					:: d_xyl(-ghostp:NX+ghostp), d_xyr(-ghostp:NX+ghostp)

	double precision 					:: d_yxl(-ghostp:NY+ghostp),d_yxr(-ghostp:NY+ghostp)
	double precision 					:: d_yyl(-ghostp:NY+ghostp),d_yyr(-ghostp:NY+ghostp)


		do k=1,n_eqn

        	do i =-2,NY+2
              consl(k,i) = -1.0d0/8.0d0*cons_local(k,i-1)+6.0d0/8.0d0*cons_local(k,i  )+3.0d0/8.0d0*cons_local(k,i+1)
              consr(k,i) =  3.0d0/8.0d0*cons_local(k,i  )+6.0d0/8.0d0*cons_local(k,i+1)-1.0d0/8.0d0*cons_local(k,i+2)

	         	! Reconstruction formula
              ! fleft(k,i) =  -1.0d0/8.0d0*flux_node(k,i-1)+6.0d0/8.0d0*flux_node(k,i  )+3.0d0/8.0d0*flux_node(k,i+1)
              ! fright(k,i) =  3.0d0/8.0d0*flux_node(k,i  )+6.0d0/8.0d0*flux_node(k,i+1)-1.0d0/8.0d0*flux_node(k,i+2)

            d_xxl(i) = (-1.0d0/8.0d0)*(d_xx(i-1)) + (6.0d0/8.0d00)*(d_xx(i+0)) + (3.0d0/8.0d0)*(d_xx(i+1))
			d_xxr(i) = (+3.0d0/8.0d0)*(d_xx(i+0)) + (6.0d0/8.0d00)*(d_xx(i+1)) - (1.0d0/8.0d0)*(d_xx(i+2))

         	enddo
        enddo

    
		do i=-2,NY+2

				Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
				Tr = (Lr**2.0d0)/( ((d_xxl(i)+d_xxr(i))*0.5d0))


			  fleft (1,i)  =  -consl(3,i)*d_xxl(i)
			  fright(1,i)  =  -consr(3,i)*d_xxr(i)

			  fleft (2,i)  =  0.0d0
			  fright(2,i)  =  0.0d0

			  fleft (3,i)  = -consl(1,i)
			  fright(3,i)  = -consr(1,i)




			delphi(i)			= consr(1,i)-consl(1,i)
			delp(i)				= consr(2,i)-consl(2,i)
			delq(i)				= consr(3,i)-consl(3,i)

	!***************************************
	!***** 		Roe solver	       	   *****
	!***************************************

		 !***** Roe solver          *****
		 !The upwind dissipation matrix is |A| = R*|Lambda|*L
		 jacobian(1,1) = abs(sqrt(((d_xxl(i)+d_xxr(i))*0.5d0)/Tr))
		 jacobian(1,2) = 0.0d0
		 jacobian(1,3) = 0.0d0
		 jacobian(2,1) = 0.0d0
		 jacobian(2,2) = 0.0d0
		 jacobian(2,3) = 0.0d0
		 jacobian(3,1) = 0.0d0
		 jacobian(3,2) = 0.0d0
		 jacobian(3,3) = Tr*abs(sqrt(((d_xxl(i)+d_xxr(i))*0.5d0)/Tr))


		dissipation (1,i) = (delphi(i)*jacobian(1,1) + delp(i) *jacobian(1,2)+delq(i)*jacobian(1,3))

		dissipation (2,i) = (delphi(i)*jacobian(2,1) + delp(i) *jacobian(2,2)+delq(i)*jacobian(2,3))

		dissipation (3,i) = (delphi(i)*jacobian(3,1) + delp(i) *jacobian(3,2)+delq(i)*jacobian(3,3))

		
		
		flux(1,i) =  0.5d0*((fright(1,i) +fleft(1,i) -  dissipation(1,i)))	
        flux(2,i) =  0.5d0*((fright(2,i) +fleft(2,i) -  dissipation(2,i)))	
        flux(3,i) =  0.5d0*((fright(3,i) +fleft(3,i) -  dissipation(3,i)))
        

	enddo




	end subroutine reconstruction5EY



 	!***********************************************************************
	!*****                       Source term 			               *****
	!***********************************************************************

 	subroutine sourceterm(NX,NY,x,y,cons,flux_x,flux_y,residual,n_eqn,ghostp)

 		implicit none


 		double precision, parameter 		:: pi=acos(-1.0d0)
 		integer								:: NX, NY, n_eqn, ghostp
		double precision 					:: X(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp), L
	
		
		integer								:: i,j
		
		
		double precision					:: sourcet(3,1:NX,1:NY)	
		double precision 					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision 					:: residual(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision 					:: flux_x(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),flux_y(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision              		:: d_xx, d_xy,d_yx,d_yy
		double precision              		:: nu,Lr,Tr, my_function
			

			do j = 1, NY

				do i = 1, NX

			  d_xx = 3.0d0+cons(1,i,j)**(5.0d0/1.0d0)+3.0d0*cons(1,i,j)**2.0d0

			  ! d_xx = sin(3.14159265358979d0*x(i))**2*sin(3.14159265358979d0*y(j))**2 &
     !  + sin(3.14159265358979d0*x(i))*sin(3.14159265358979d0*y(j)) + 1

			  Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
	  		  Tr = (Lr**2.0d0)/(d_xx)


				sourcet(1,i,j)      = my_function(x(i), y(j))

				sourcet(2,i,j)		= -cons(2,i,j)
				sourcet(3,i,j)		= -cons(3,i,j)

					residual(1,i,j)		=  flux_x(1,i,j) + flux_y(1,i,j) + sourcet(1,i,j)
					residual(2,i,j)		= (flux_x(2,i,j) + flux_y(2,i,j) + sourcet(2,i,j))/Tr
					residual(3,i,j)		= (flux_x(3,i,j) + flux_y(3,i,j) + sourcet(3,i,j))/Tr
				enddo

		enddo
	end subroutine sourceterm
	
double precision function my_function(x, y)
	implicit none
	double precision, intent(in) :: x
	double precision, intent(in) :: y

	double precision, parameter :: pi = 3.14159265358979d0
	my_function = pi**2*(12*sin(3.14159265358979d0*x)**5*sin( &
      3.14159265358979d0*y)**5 - 5*sin(3.14159265358979d0*x)**5*sin( &
      3.14159265358979d0*y)**3 - 5*sin(3.14159265358979d0*x)**3*sin( &
      3.14159265358979d0*y)**5 + 18*sin(3.14159265358979d0*x)**2*sin( &
      3.14159265358979d0*y)**2 - 6*sin(3.14159265358979d0*x)**2 - 6*sin &
      (3.14159265358979d0*y)**2 + 6)*sin(3.14159265358979d0*x)*sin( &
      3.14159265358979d0*y)


end function


	!***********************************************************************
	!*****                       Output 			                   *****
	!***********************************************************************

	subroutine output(NX,NY,x,y,cons,cons_exact,n_eqn,ghostp)

	integer 				:: i, j

	integer					:: NX,NY, n_eqn, ghostp

	double precision 		:: X(-ghostp:NX+ghostp),Y(-ghostp:NY+ghostp),cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: cons_exact(n_eqn,1:NX,1:NY)


	open(unit=1,file='solution_4x.plt',form='formatted',status='replace')
	
	write(1,*) 'TITLE=" "Anistropic"'
    write(1,*) 'VARIABLES = "x","y","Pote","vx","vy","exac"'
	write(1,*) "ZONE I=",NX," J=",NY," F=POINT"


    do j = 1, NY
        do i = 1, NX

          write(1,'(6F25.8)') x(i), y(j), cons(1,i,j), cons(2,i,j), cons(3,i,j), cons_exact(1,i,j)

	  enddo
    enddo
	
    close(1)



 	end subroutine output




