! 
!**********************************************************************
!    This program solves 2D Diffusion equation with Tensor
! Generalized and non-aligned tensor
!**********************************************************************


	program hyperpoi2d
	
		implicit none


		double precision, parameter 		:: pi=acos(-1.0d0)

		integer, parameter					:: NTMAX=40000000,l=2
		integer 							:: i, j, z, k, NT
		
		integer, parameter					:: NX = 20*l,NY = 20*l, n_eqn = 3, ghostp = 10 	

		double precision, parameter			:: l1_target = 1.0d-12, epsilon = 1.0d0-14
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
		ymax =  0.005d0

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

		call timestep(NX,NY,dt,dx,dy,CFL,x,y,ghostp)


	!***********************************************************************
	!*****         			 TVD- Runge Kutta 				       	   *****
	!***********************************************************************
		
	iteration : do NT=1,NTMAX
	call timestep(NX,NY,dt,dx,dy,CFL,x,y,ghostp)

		
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

		 if ( MAXVAL( res_norm(:,L1)/res_norm_initial(:,L1) )  < l1_target) then
		  	
		  	exit iteration

		 endif

		if (NT==1) then 

			  res_norm_initial = res_norm

	write(*,*) " Initial residual (after 1 iteration):"
	write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (L1) = ", res_norm_initial(:,L1)
	write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (L2) = ", res_norm_initial(:,L2)
	write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (Li) = ", res_norm_initial(:,Linf)
			  
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

		write(*,*) " Initial residual (after 1 iteration):"
		write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (L1) = ", res_norm_initial(:,L1)
		write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (L2) = ", res_norm_initial(:,L2)
		write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (Li) = ", res_norm_initial(:,Linf)

		write(*,*)
		write(*,*) " Final residual norm:"
		write(*,'(a15,4es12.4,4es12.4,4es12.4)') "  Residual norm (L1) = ", res_norm(:,L1)
		write(*,'(a15,4es12.4,4es12.4,4es12.4)') "  Residual norm (L2) = ", res_norm(:,L2)
		write(*,'(a15,4es12.4,4es12.4,4es12.4)') "  Residual norm (Li) = ", res_norm(:,Linf)
		write(*,*)
	



	!***************************************************************************************************
	!*****                Error w.r.t exact solution, doing it only for Potential             	   *****
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
    
        write(*,*) '2D Poisson equation by hyperbolic approach, Roe - 5th order Compact'
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

		write(*,*) MINVAL(cons(1,1:NX,1:NY))




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

				cons(k,i,j)		= 1.0d0

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
				

				cons_exact(1,i,j) = (sin(x(i)*pi)*sin(y(j)*pi*2000.0d0))
			enddo
		enddo

		
	end subroutine exactsolution


	!***********************************************************************
	!*****                       Compute time step                	   *****
	!***********************************************************************

 	subroutine timestep(NX,NY,dt,dx,dy,CFL,x,y,ghostp)
 	implicit none

 	integer 							:: i, j
 	
 	integer								:: NX, NY,ghostp

 	double precision 					:: dx, dy
	double precision 					:: dt, CFL, dtnew
	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision, parameter 		:: pi=acos(-1.0d0), beta =  atan(tand(5.0d0)), D_parallel = 1.0

	double precision					:: x_velocity, y_velocity
	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx, d_xy,d_yx,d_yy




	dt = 1.0d010
		
		do i = 1, NX

			do j = 1, NY

				d_xx = D_parallel
				d_xy =   0.0d0
				d_xy =   0.0d0
				d_yy =   1.0d0
				Lr = (0.005/sqrt(1.0d0+0.005d0**2))/(2.0d0*acos(-1.0d0))
				Tr = Lr**2.0d0

				  x_velocity =  DABS(sqrt(d_xx/Tr))
				  y_velocity =  DABS(sqrt(d_yy/Tr))
				  dtnew = min( (x(i+1)-x(i))/x_velocity, (y(j+1)-y(j))/y_velocity)
				    
				   if(dtnew .lt. dt) dt = dtnew

				
			enddo
		enddo

  		  CFL 	=	0.45d0
  		  dt 	= 	CFL*dt


 	end subroutine timestep

 	!***********************************************************************
	!*****                       Boundary conditions                   *****
	!***********************************************************************

	subroutine boundarycondition(NX,NY,cons,n_eqn,ghostp,x,y)

	implicit none

	integer 							:: i, j, k
	double precision, parameter 		:: pi=acos(-1.0d0), beta =  atan(tand(5.0d0)), D_parallel = 1.0
	integer								:: NX, NY, n_eqn, ghostp

	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)	
	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	! Primary variable for Dirichlet boundary condition. Conservative variable 1

	double precision 					:: phi_left,phi_right, phi_top, phi_bottom

! Gradient variables for Neumann boundary conditions. Conservative variables 2 and 3 for boundaries

	double precision					:: px_top, px_bottom, px_right, px_left

	double precision					:: py_top, py_bottom, py_left, py_right

! Kuzmin problem
 ! ! Left boundary conditions

	do j = 1,NY
		
		i=0     
	    

        cons(1,i,j)   = (128.0d0/35.0d0)*( 0.0d0 - (35.0d0/32.0d0)*cons(1,i+1,j) + (35.0d0/64.0d0)*cons(1,i+2,j)-&
					  (7.0d0/32.0d0)*cons(1,i+3,j) +(5.0d0/128.0d0)*cons(1,i+4,j))
		
		cons(2,i,j)    = 5.0d0 * cons(2,i+1,j) - 10.0d0*cons(2,i+2,j) + 10.0d0*cons(2,i+3,j) -5.0d0*cons(2,i+4,j) + cons(2,i+5,j)
		cons(3,i,j)    = 5.0d0 * cons(3,i+1,j) - 10.0d0*cons(3,i+2,j) + 10.0d0*cons(3,i+3,j) -5.0d0*cons(3,i+4,j) + cons(3,i+5,j)
    enddo

    do j = 1,NY

    	do i = -1,-ghostp,-1
        
        	cons(1:3,i,j)   = 5.0d0 * cons(1:3,i+1,j) - 10.0d0*cons(1:3,i+2,j) + 10.0d0*cons(1:3,i+3,j) -5.0d0*cons(1:3,i+4,j) + cons(1:3,i+5,j)

        enddo
    enddo

! ! Right boundary conditions


	do j = 1, NY
		
		i = NY

        cons(1,i+1,j) = (128.0d0/35.0d0)*(0.0d0+(5.0d0/128.0d0)*cons(1,i-3,j) - (7.0d0/32.0d0)*cons(1,i-2,j) +&
					   (35.0d0/64.0d0)*cons(1,i-1,j) -(35.0d0/32.d0)*cons(1,i,j))

		cons(2,i+1,j)  = cons(2,i-4,j) - 5.0d0*cons(2,i-3,j) + 10.0d0*cons(2,i-2,j) -10.0d0*cons(2,i-1,j)+5.0d0*cons(2,i,j)
		cons(3,i+1,j)  = cons(3,i-4,j) - 5.0d0*cons(3,i-3,j) + 10.0d0*cons(3,i-2,j) -10.0d0*cons(3,i-1,j)+5.0d0*cons(3,i,j)
        
    ! enddo

    ! do j = 1, NY

    	do i= NY, NY+6      

        cons(1:3,i+2,j)  = cons(1:3,i-3,j) - 5.0d0*cons(1:3,i-2,j) + 10.0d0*cons(1:3,i-1,j) -10.0d0*cons(1:3,i,j)  +5.0d0*cons(1:3,i+1,j)

        enddo
    enddo




	! Bottom boundary conditions

	do i = 1,NX
		j=0     

        cons(1,i,j)   = (128.0d0/35.0d0)*(0.0d0 - (35.0d0/32.0d0)*cons(1,i,j+1) + (35.0d0/64.0d0)*cons(1,i,j+2)-&
					    (7.0d0/32.0d0)*cons(1,i,j+3) +(5.0d0/128.0d0)*cons(1,i,j+4))
		
		cons(2,i,j)    = 5.0d0 * cons(2,i,j+1) - 10.0d0*cons(2,i,j+2) + 10.0d0*cons(2,i,j+3) -5.0d0*cons(2,i,j+4) + cons(2,i,j+5)
		cons(3,i,j)    = 5.0d0 * cons(3,i,j+1) - 10.0d0*cons(3,i,j+2) + 10.0d0*cons(3,i,j+3) -5.0d0*cons(3,i,j+4) + cons(3,i,j+5)
    enddo

    do i = 1,NX  
    	do j = -1,-ghostp,-1
        
        	
        	cons(1:3,i,j)   = 5.0d0 * cons(1:3,i,j+1) - 10.0d0*cons(1:3,i,j+2) + 10.0d0*cons(1:3,i,j+3) -5.0d0*cons(1:3,i,j+4) + cons(1:3,i,j+5)

        enddo
    enddo


	! Top boundary conditions
	do i = 1, NX
		 j = NX


        cons(1,i,j+1) = (128.0d0/35.0d0)*((sin(x(i)*pi)*sin((y(j)+y(j+1))*0.5d0*pi*2000.0d0)) +(5.0d0/128.0d0)*cons(1,i,j-3) - (7.0d0/32.0d0)*cons(1,i,j-2) +&
					   (35.0d0/64.0d0)*cons(1,i,j-1) -(35.0d0/32.d0)*cons(1,i,j))

		cons(2,i,j+1)  = cons(2,i,j-4) - 5.0d0*cons(2,i,j-3) + 10.0d0*cons(2,i,j-2) -10.0d0*cons(2,i,j-1)+5.0d0*cons(2,i,j)
		cons(3,i,j+1)  = cons(3,i,j-4) - 5.0d0*cons(3,i,j-3) + 10.0d0*cons(3,i,j-2) -10.0d0*cons(3,i,j-1)+5.0d0*cons(3,i,j)
        
    enddo

    do i = 1, NX

    	do j= NX, NX+6

        
        cons(1:3,i,j+2)  = cons(1:3,i,j-3) - 5.0d0*cons(1:3,i,j-2) + 10.0d0*cons(1:3,i,j-1) -10.0d0*cons(1:3,i,j)  +5.0d0*cons(1:3,i,j+1)
 

        enddo
    enddo

! phi_left = 0.0d0; phi_top = 0.0d0; 	phi_right = 0.0d0;  phi_bottom = 0.0d0;		

! 	! Bottom boundary conditions (actually bottom)

! 	do i = 1,NX
! 		j=0     
! 	    cons(1,i,j )   = (8.0d0*phi_bottom - 6.0d0*cons(1,i,j+1) + cons(1,i,j+2)) * (1.0d0/3.0d0)

! 	    cons(2,i,j )   = cons(2,i,j+3)-3.0d0*cons(2,i,j+2)+3.0d0*cons(2,i,j+1)
!         cons(3,i,j )   = cons(3,i,j+3)-3.0d0*cons(3,i,j+2)+3.0d0*cons(3,i,j+1)
!     enddo

!     do i = 1,NX  
!     	do j = -1,-ghostp,-1
        
!         	cons(:,i,j)   = cons(:,i,j+3)-3.0d0*cons(:,i,j+2)+3.0d0*cons(:,i,j+1)

!         enddo
!     enddo


! 	! Top boundary conditions (actually top)
! 	do i = 1, NX
! 		 j = NX

! 		cons(1,i,j+1) = (8.0d0*phi_top +1.0d0*cons(1,i,j-1) -6.0d0* cons(1,i,j)) * (1.0d0/3.0d0)
        
!         cons(2,i,j+1) = 3.0d0*cons(2,i,j)   - 3.0d0*cons(2,i,j-1)+cons(2,i,j-2)
!         cons(3,i,j+1) = 3.0d0*cons(3,i,j)   - 3.0d0*cons(3,i,j-1)+cons(3,i,j-2)
        
!     enddo

!     do i = 1, NX

!     	do j= NX, NX+2

!         cons(:,i,j+2) = 3.0d0*cons(:,i,j+1) - 3.0d0*cons(:,i,j)  +cons(:,i,j-1)

!         enddo
!     enddo


! 	! Left boundary conditions (actually left)

! 	do j = 1,NY
		
! 		i=0     
	    
! 	    cons(1,i,j )   = (8.0d0*phi_left - 6.0d0*cons(1,i+1,j) + cons(1,i+2,j)) * (1.0d0/3.0d0)

! 	    cons(2,i,j )   = cons(2,i+3,j)-3.0d0*cons(2,i+2,j)+3.0d0*cons(2,i+1,j)
        
!         cons(3,i,j )   = cons(3,i+3,j)-3.0d0*cons(3,i+2,j)+3.0d0*cons(3,i+1,j)
!     enddo

!     do j = 1,NY

!     	do i = -1,-ghostp,-1
        
!         	cons(:,i,j)   = cons(:,i+3,j)-3.0d0*cons(:,i+2,j)+3.0d0*cons(:,i+1,j)

!         enddo
!     enddo


! 	! Right boundary conditions (actually right)


! 	do j = 1, NY
		
! 		i = NY

! 		cons(1,i+1,j) = (8.0d0*phi_right +1.0d0*cons(1,i-1,j) -6.0d0* cons(1,i,j)) * (1.0d0/3.0d0)
        
!         cons(2,i+1,j) = 3.0d0*cons(2,i,j)   - 3.0d0*cons(2,i-1,j)+cons(2,i-2,j)
!         cons(3,i+1,j) = 3.0d0*cons(3,i,j)   - 3.0d0*cons(3,i-1,j)+cons(3,i-2,j)
        
!     enddo

!     do j = 1, NY

!     	do i= NY, NY+2

!         cons(:,i+2,j) = 3.0d0*cons(:,i+1,j) - 3.0d0*cons(:,i,j)  +cons(:,i-1,j)

!         enddo
!     enddo


	

	end subroutine boundarycondition



	!***********************************************************************
	!*****                      Flux in X-direction                    *****
	!***********************************************************************

	subroutine FX(NX,NY,dx,x,y,cons,flux_x,n_eqn,ghostp)

	implicit none

	integer								:: i, j ,k
	integer 							:: NX, NY, ghostp, n_eqn
	double precision					:: dx, x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision 					:: flux_x(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision					:: flux_node_x(n_eqn,-ghostp:NX+ghostp), flux_x_half(n_eqn,-ghostp:NX+ghostp)

	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),cons_local(n_eqn,-ghostp:NX+ghostp)

	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx, d_xy,d_yx,d_yy
	double precision, parameter 		:: pi=acos(-1.0d0), beta =  atan(tand(5.0d0)), D_parallel = 1.0



	do j =1, NY

		do i = -ghostp, NX+ghostp

  d_xx = D_parallel
  d_xy =   0.0d0
  d_xy =   0.0d0
  d_yy =   1.0d0

			cons_local(1,i) = cons(1,i,j)
			cons_local(2,i) = cons(2,i,j)
			cons_local(3,i) = cons(3,i,j)

				Lr = (0.005/sqrt(1.0d0+0.005d0**2))/(2.0d0*acos(-1.0d0))
				Tr = Lr**2.0d0

			flux_node_x(1,i) = -cons_local(2,i)*d_xx-cons_local(3,i)*d_xy
			flux_node_x(2,i) = -cons_local(1,i)/Tr
			flux_node_x(3,i) = 0.0d0
		end do

		call reconstruction5EX (NX,NY,cons_local,flux_x_half, n_eqn,ghostp,x,y,j)

		do k = 1,n_eqn

			do i = 1, NX

				!Finite volume or Reconstruction approach 
	        	flux_x(k,i,j) = -((flux_x_half(k,i)-flux_x_half(k,i-1))/dx)

	        	!Finite difference or interpolation approach
				! flux_x(k,i,j) = -((4.0/3.0d0)*(flux_x_half(k,i  )-flux_x_half(k,i-1)) +&
				!                  (-1.0/6.0d0)*(flux_node_x(k,i+1)-flux_node_x(k,i-1)))/dx

				flux_x(k,i,j) = -((9.0/8.0d0)*(flux_x_half(k,i  )-flux_x_half(k,i-1)) +&
				                 (-1.0/24.0d0)*(flux_x_half(k,i+1)-flux_x_half(k,i-2)))/dx

				! flux_x(k,i,j) = -((3.0d0/2.0d0)*(flux_x_half(k,i)-flux_x_half(k,i-1))/dx + &
				! 			    (1.0d0/30.0d0)*(flux_x_half(k,i+1)-flux_x_half(k,i-2))/dx &
				! 			    +(-3.0d0/10.0d0)*(flux_node_x(k,i+1)-flux_node_x(k,i-1))/dx)

				flux_x(k,i,j) = -((75.0d0/64.0d0)*(flux_x_half(k,i)-flux_x_half(k,i-1))/dx - (25.0d0/384.0d0)*(flux_x_half(k,i+1)-flux_x_half(k,i-2))/dx &
				              	+ (3.0d0/640.0d0)*(flux_x_half(k,i+2)-flux_x_half(k,i-3))/dx)



			enddo
		enddo
	enddo

	end subroutine FX
			

	!***************************************************************************************************************
	!***** 					Roe, HLL or Rusanov + 5th order compact upwind reconstruction or interplation     *****
	!*****											in X-direction											   *****
	!***************************************************************************************************************

	subroutine reconstruction5EX (NX,NY,cons_local,flux, n_eqn,ghostp,x,y,j)

	implicit none 

	integer, intent(in)					:: j
	integer								:: i, k	
	integer								:: NX,NY, n_eqn, ghostp, total
	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision					:: cons_local(n_eqn,-ghostp:NX+ghostp),flux(n_eqn,-ghostp:NX+ghostp)

	double precision 					:: consl(n_eqn,-ghostp:NX+ghostp), consr(n_eqn,-ghostp:NX+ghostp)

	double precision   					:: fleft(n_eqn,-ghostp:NX+ghostp), fright(n_eqn,-ghostp:NX+ghostp)

	double precision					:: righteigen(3,3), lefteigen(3,3)

	double precision					:: dissipation(n_eqn,-2:NX+2),jacobian(3,3)

	double precision 					:: lambda1(-ghostp:NX+ghostp), lambda2(-ghostp:NX+ghostp)
	double precision  					:: delphi(-ghostp:NX+ghostp), delp(-ghostp:NX+ghostp), delq(-ghostp:NX+ghostp)
	double precision 					:: alpha1, alpha2
	double precision 					:: SLm(-ghostp:NX+ghostp), SRp(-ghostp:NX+ghostp),laxmax

	double precision					:: a(0:NX+1), b(0:NX+1), c(0:NX+1), r(0:NX+1),temp(0:NX+1)

	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx, d_xy,d_yx,d_yy

	double precision					:: xf, yf

	double precision					:: xl, xr, yl, yr

	double precision, parameter			:: epsilon = 1.0d-23
    double precision 					:: recon_poly5_r(0:3), recon_poly5_l(0:3)
    double precision 					:: alpha_l(0:3),alpha_r(0:3)
    double precision 					:: beta_w(0:3), omega_l(0:3), omega_r(0:3)

    double precision, parameter 		:: pi=acos(-1.0d0), beta =  atan(tand(5.0d0)), D_parallel = 1.0
   
    double precision, parameter 		:: Constant=1.0d0,p =1.0d0






		do k=1,n_eqn

        	do i =-2,NX+2

        	   ! ! Interpolation formula		
            !   consl(k,i) = -1.0d0/8.0d0*cons_local(k,i-1)+6.0d0/8.0d0*cons_local(k,i  )+3.0d0/8.0d0*cons_local(k,i+1)
            !   consr(k,i) =  3.0d0/8.0d0*cons_local(k,i  )+6.0d0/8.0d0*cons_local(k,i+1)-1.0d0/8.0d0*cons_local(k,i+2)

        	! Reconstruction formula
   !      	consl(k,i) = (2.0d0/60.0d0)*(cons_local(k,i-2)) + (-13.0d0/60.0d0)*(cons_local(k,i-1)) + (47.0d0/60.0d0)*(cons_local(k,i)) + &
			! 		  (27.0d0/60.0d0)*(cons_local(k,i+1)) + (-3.0d0/60.0d0)*(cons_local(k,i+2))

			! consr(k,i) = (-3.0d0/60.0d0)*(cons_local(k,i-1)) + (27.0d0/60.0d0)*(cons_local(k,i)) + (47.0d0/60.0d0)*(cons_local(k,i+1)) + &
			! 	      (-13.0d0/60.0d0)*(cons_local(k,i+2)) + (2.0d0/60.0d0)*(cons_local(k,i+3))

   !      	xl = (2.0d0/60.0d0)*(x(i-2)) + (-13.0d0/60.0d0)*(x(i-1)) + (47.0d0/60.0d0)*(x(i)) + &
			! 		  (27.0d0/60.0d0)*(x(i+1)) + (-3.0d0/60.0d0)*(x(i+2))

			! xr = (-3.0d0/60.0d0)*(x(i-1)) + (27.0d0/60.0d0)*(x(i)) + (47.0d0/60.0d0)*(x(i+1)) + &
			! 	      (-13.0d0/60.0d0)*(x(i+2)) + (2.0d0/60.0d0)*(x(i+3))


  	       ! Interpolation formula		
			consl(k,i) = (3.0d0/128.0d0)*(cons_local(k,i-2)) - (5.0d0/32.0d0)*(cons_local(k,i-1)) + (45.0d0/64.0d0)*(cons_local(k,i)) + &
					     (15.0d0/32.0d0)*(cons_local(k,i+1)) - (5.0d0/128.0d0)*(cons_local(k,i+2))
			
			consr(k,i) = (-5.0d0/128.0d0)*(cons_local(k,i-1)) + (15.0d0/32.0d0)*(cons_local(k,i)) + (45.0d0/64.0d0)*(cons_local(k,i+1)) - &
					     (5.0d0/32.0d0)*(cons_local(k,i+2)) + (3.0d0/128.0d0)*(cons_local(k,i+3))
			! consl(k,i) = (1.0d0/60.0d0)*(cons_local(k,i-2)) - (8.0d0/60.0d0)*(cons_local(k,i-1)) + (37.0d0/60.0d0)*(cons_local(k,i)) +&
   !     		 		  (37.0d0/60.0d0)*(cons_local(k,i+1)) - (8.0d0/60.0d0)*cons_local(k,i+2)+ (1.0d0/60.0d0)*(cons_local(k,i+3))
	  !   	consr(k,i) = (1.0d0/60.0d0)*(cons_local(k,i-2)) - (8.0d0/60.0d0)*(cons_local(k,i-1)) + (37.0d0/60.0d0)*(cons_local(k,i)) +&
   !     			      (37.0d0/60.0d0)*(cons_local(k,i+1)) - (8.0d0/60.0d0)*cons_local(k,i+2)+ (1.0d0/60.0d0)*(cons_local(k,i+3))

		 !  	xl = (3.0d0/128.0d0)*(x(i-2)) - (5.0d0/32.0d0)*(x(i-1)) + (45.0d0/64.0d0)*(x(i)) + &
			! 		     (15.0d0/32.0d0)*(x(i+1)) - (5.0d0/128.0d0)*(x(i+2))
			
			! xr = (-5.0d0/128.0d0)*(x(i-1)) + (15.0d0/32.0d0)*(x(i)) + (45.0d0/64.0d0)*(x(i+1)) - &
			! 		     (5.0d0/32.0d0)*(x(i+2)) + (3.0d0/128.0d0)*(x(i+3))

			! WCNS interpolation

!  				beta_w(0) 	= (13.0d0/12.0d0)*(cons_local(k,i)-2.0d0*cons_local(k,i+1)+cons_local(k,i+2))**2+ &
! 	            				(1.0/4.0d0)*(3.0d0*cons_local(k,i)-4.0d0*cons_local(k,i+1)+cons_local(k,i+2))**2

! 	            beta_w(1) 	= (13.0d0/12.0d0)*(cons_local(k,i-1)-2.d0*cons_local(k,i)+cons_local(k,i+1))**2+ &
! 	            				(1.0/4.0d0)*(cons_local(k,i-1)-cons_local(k,i+1))**2

! 	            beta_w(2)		= (13.0d0/12.0d0)*(cons_local(k,i-2)-2.d0*cons_local(k,i-1)+cons_local(k,i))**2 &
! 	            				+(1.0/4.0d0)*(cons_local(k,i-2)-4.0d0*cons_local(k,i-1)+3.d0*cons_local(k,i))**2


! 	            ! recon_poly5_l(0) = 3.0d0/8.0d0*cons_local(k,i)+6.0d0/8.0d0*cons_local(k,i+1)-1.0d0/8.0d0*cons_local(k,i+2)
! 	            ! recon_poly5_l(1) = -1.0d0/8.0d0*cons_local(k,i-1)+6.0d0/8.0d0*cons_local(k,i)+3.0d0/8.0d0*cons_local(k,i+1)
! 	            ! recon_poly5_l(2) = 3.0d0/8.0d0*cons_local(k,i-2)-10.0d0/8.0d0*cons_local(k,i-1)+15.0d0/8.0d0*cons_local(k,i)

! 	           	recon_poly5_l(0) =1.0d0/8.0d0*( 3.0d0*cons_local(k,i  )+ 6.0d0*cons_local(k,i+1)-1.0d0 *cons_local(k,i+2))
! 	            recon_poly5_l(1) =1.0d0/8.0d0*(-1.0d0*cons_local(k,i-1)+ 6.0d0*cons_local(k,i  )+3.0d0 *cons_local(k,i+1))
! 	            recon_poly5_l(2) =1.0d0/8.0d0*( 3.0d0*cons_local(k,i-2)-10.0d0*cons_local(k,i-1)+15.0d0*cons_local(k,i  ))

! 	! Jiang-Shu linear weights
! 	            ! alpha_l(0) 		= (5.0/16.0d0)
! 	            ! alpha_l(1) 		= (10.0/16.0d0)
! 	            ! alpha_l(2) 		= (1.0/16.0d0)
! ! ! Borges
! 	        alpha_l(0) = 5.0d0/16.0d0  * (Constant + (abs(beta(0)-beta(2)) / (epsilon + beta(0)))**p)
!             alpha_l(1) = 10.0d0/16.0d0 * (Constant + (abs(beta(0)-beta(2)) / (epsilon + beta(1)))**p)
!             alpha_l(2) = 1.0d0/16.0d0  * (Constant + (abs(beta(0)-beta(2)) / (epsilon + beta(2)))**p)

! 	            omega_l(0)=alpha_l(0)/(alpha_l(0)+alpha_l(1)+alpha_l(2))
! 	            omega_l(1)=alpha_l(1)/(alpha_l(0)+alpha_l(1)+alpha_l(2))
! 	            omega_l(2)=alpha_l(2)/(alpha_l(0)+alpha_l(1)+alpha_l(2))

! 	            consl(k,i) = omega_l(0) * recon_poly5_l(0) + omega_l(1) * recon_poly5_l(1) + omega_l(2) * recon_poly5_l(2)



! 	! right of i+1/2
! 	            beta_w(0)	= (13.0d0/12.0d0)*(cons_local(k,i+1)-2.0d0*cons_local(k,i+2)+cons_local(k,i+3))**2+ &
! 	            						(1.0/4.0d0)*(3.0d0*cons_local(k,i+1)-4.0d0*cons_local(k,i+2)+cons_local(k,i+3))**2

! 	            beta_w(1)	= (13.0d0/12.0d0)*(cons_local(k,i)-2.d0*cons_local(k,i+1)+cons_local(k,i+2))**2+ &
! 	           							(1.0/4.0d0)*(cons_local(k,i)-cons_local(k,i+2))**2

! 	            beta_w(2) = (13.0d0/12.0d0)*(cons_local(k,i-1)-2.d0*cons_local(k,i)+cons_local(k,i+1))**2 &
! 	            						+(1.0/4.0d0)*(cons_local(k,i-1)-4.0d0*cons_local(k,i)+3.d0*cons_local(k,i+1))**2


! 				! recon_poly5_r(0) = 15.0d0/8.0d0*cons_local(k,i+1)-10.0d0/8.0d0*cons_local(k,i+2)+3.0d0/8.0d0*cons_local(k,i+3)
! 	   !          recon_poly5_r(1) = 3.0d0/8.0d0*cons_local(k,i)+6.0d0/8.0d0*cons_local(k,i+1)-1.0d0/8.0d0*cons_local(k,i+2)
! 	   !          recon_poly5_r(2) = -1.0d0/8.0d0*cons_local(k,i-1)+6.0d0/8.0d0*cons_local(k,i)+3.0d0/8.0d0*cons_local(k,i+1)

! 	   			recon_poly5_r(0) = 1.0d0/8.0d0*(15.0d0*cons_local(k,i+1)-10.0d0*cons_local(k,i+2)+3.0d0*cons_local(k,i+3))
! 	            recon_poly5_r(1) = 1.0d0/8.0d0*( 3.0d0*cons_local(k,i  )+ 6.0d0*cons_local(k,i+1)-1.0d0*cons_local(k,i+2))
! 	            recon_poly5_r(2) = 1.0d0/8.0d0*(-1.0d0*cons_local(k,i-1)+ 6.0d0*cons_local(k,i  )+3.0d0*cons_local(k,i+1))

! 	! Jiang-Shu linear weights
! 	            ! alpha_r(0) = (1.0/16.0d0)
! 	            ! alpha_r(1) = (10.0/16.0d0)
! 	            ! alpha_r(2) = (5.0/16.0d0)

! ! Borges
! 	        alpha_r(0) = 1.0d0/16.0d0  * (Constant + (abs(beta(0)-beta(2)) / (epsilon + beta(0)))**p)
!             alpha_r(1) = 10.0d0/16.0d0 * (Constant + (abs(beta(0)-beta(2)) / (epsilon + beta(1)))**p)
!             alpha_r(2) = 5.0d0/16.0d0  * (Constant + (abs(beta(0)-beta(2)) / (epsilon + beta(2)))**p)

! 	! Non-linear weights
! 	            omega_r(0)=alpha_r(0)/(alpha_r(0)+alpha_r(1)+alpha_r(2))
! 	            omega_r(1)=alpha_r(1)/(alpha_r(0)+alpha_r(1)+alpha_r(2))
! 	            omega_r(2)=alpha_r(2)/(alpha_r(0)+alpha_r(1)+alpha_r(2))

	            
! 	            consr(k,i) = omega_r(0) * recon_poly5_r(0) + omega_r(1) * recon_poly5_r(1) + omega_r(2) * recon_poly5_r(2) 

        
         	enddo
        enddo

    ! Upwind Comapct Reconstruction

		! do k = 1,n_eqn
		    	

		!     	do i = 0, NX
		! 	      a(i) = 1.0d0/2.0d0
		! 	      b(i) = 1.0d0
		! 	      c(i) = 1.0d0/6.0d0
		!    		 end do
		    
		!     i = 0  ; a(i) = 0.0d+00
		!     i = NX  ; c(i) = 0.0d+00
		!     total = NX+1
		    
		!     r(0)		= (1.0d0/18.0d0)*cons_local(k, -1)+(19.0d0/18.0d0)*cons_local(k,0)+(10.0d0/18.0d0)*cons_local(k,  1)-(1.0d0/2.0d0)*consl(k, -1)
		!     r(NX)		= (1.0d0/18.0d0)*cons_local(k,NX-1)+(19.0d0/18.0d0)*cons_local(k,NX)+(10.0d0/18.0d0)*cons_local(k,NX+1)-(1.0d0/6.0d0)*consl(k,NX+1)
		    
		!     do i=1,NX-1
		!     r(i)		= (1.0d0/18.0d0)*cons_local(k, i-1)+(19.0d0/18.0d0)*cons_local(k,i)+(5.0d0/9.0d0)*cons_local(k,  i+1)
		!     enddo
		    
		!     call tridiag(a,b,c,r,temp,total)
		        
		!     consl(k,0:NX)=temp(0:NX) 

		! 	    do i = 0, NX
		! 	      a(i) = 1.0d0/6.0d0
		! 	      b(i) = 1.0d0
		! 	      c(i) = 1.0d0/2.0d0
		! 	    end do

		!     i = 0  ; a(i) = 0.0d+00
		!     i = NX  ; c(i) = 0.0d+00

		!     r(0)		= (10.0d0/18.0d0)*cons_local(k,0)+(19.0d0/18.0d0)*cons_local(k,  1)+(1.0d0/18.0d0)*cons_local(k,  2)-(1.0d0/6.0d0)*consr(k, -1)
		!     r(NX)		= (10.0d0/18.0d0)*cons_local(k,NX)+(19.0d0/18.0d0)*cons_local(k,NX+1)+(1.0d0/18.0d0)*cons_local(k,NX+2)-(1.0d0/2.0d0)*consr(k,NX+1)
		   
		!     do i=1,NX-1
		!     r(i)		= (10.0d0/18.0d0)*cons_local(k,i)+(19.0d0/18.0d0)*cons_local(k,i+1)+(1.0d0/18.0d0)*cons_local(k,i+2)
		!     enddo
		    
		!     call tridiag(a,b,c,r,temp,total)
		        
		!     consr(k,0:NX)=temp(0:NX)

		! enddo

	 ! Compact reconstruction low dissipation

		    ! do k = 1,n_eqn
		    	

		    ! 	do i = 0, NX
			   !    a(i) = 5.0d0/20.0d0
			   !    b(i) = 12.0d0/20.0d0
			   !    c(i) = 3.0d0/20.0d0
		   	! 	 end do
		    
		    ! i = 0  ; a(i) = 0.0d+00
		    ! i = NX  ; c(i) = 0.0d+00
		    ! total = NX+1
		    
		    ! r(0)		= (3.0d0/120.0d0)*cons_local(k, -1)+(67.0d0/120.0d0)*cons_local(k,0)+(49.0d0/120.0d0)*cons_local(k,  1)+(1.0d0/120.0d0)*cons_local(k,  2)-(5.0d0/20.0d0)*consl(k, -1)
		    ! r(NX)		= (3.0d0/120.0d0)*cons_local(k,NX-1)+(67.0d0/120.0d0)*cons_local(k,NX)+(49.0d0/120.0d0)*cons_local(k,NX+1)+(1.0d0/120.0d0)*cons_local(k,  NX+2)-(3.0d0/20.0d0)*consl(k,NX+1)
		    
		    ! do i=1,NX-1
		    ! r(i)		= (3.0d0/120.0d0)*cons_local(k,i-1)+(67.0d0/120.0d0)*cons_local(k,i)+(49.0d0/120.0d0)*cons_local(k,i+1)+(1.0d0/120.0d0)*cons_local(k,  i+2)
		    ! enddo
		    
		    ! call tridiag(a,b,c,r,temp,total)
		        
		    ! consl(k,0:NX)=temp(0:NX) 

			   !  do i = 0, NX
			   !    a(i) = 3.0d0/20.0d0
			   !    b(i) = 12.0d0/20.0d0
			   !    c(i) = 5.0d0/20.0d0
			   !  end do

		    ! i = 0  ; a(i) = 0.0d+00
		    ! i = NX  ; c(i) = 0.0d+00

		    ! r(0)		= (1.0d0/120.0d0)*cons_local(k, -1)+(49.0d0/120.0d0)*cons_local(k,0)+(67.0d0/120.0d0)*cons_local(k,  1)+(3.0d0/120.0d0)*cons_local(k,  2)-(3.0d0/20.0d0)*consr(k, -1)
		    ! r(NX)		= (1.0d0/120.0d0)*cons_local(k,NX-1)+(49.0d0/120.0d0)*cons_local(k,NX)+(67.0d0/120.0d0)*cons_local(k,NX+1)+(3.0d0/120.0d0)*cons_local(k,  NX+2)-(5.0d0/20.0d0)*consr(k,NX+1)
		   
		    ! do i=1,NX-1
		    ! r(i)		= (1.0d0/120.0d0)*cons_local(k,i-1)+(49.0d0/120.0d0)*cons_local(k,i)+(67.0d0/120.0d0)*cons_local(k,i+1)+(3.0d0/120.0d0)*cons_local(k,  i+2)
		    ! enddo
		    
		    ! call tridiag(a,b,c,r,temp,total)
		        
		    ! consr(k,0:NX)=temp(0:NX)

		    ! enddo

	! Upwind interpolation
     do k = 1,n_eqn
	    	

	    	do i = 0, NX
		      a(i) = 0.5d0
		      b(i) = 1.0d0
		      c(i) = 0.1d0
	   		 end do
	    
	    i = 0  ; a(i) = 0.0d+00
	    i = NX  ; c(i) = 0.0d+00
	    total = NX+1
	    
	    r(0)		= 0.1d0*cons_local(k, -1)+1.0d0*cons_local(k,0)+0.5d0*cons_local(k,  1)-0.5d0*consl(k, -1)
	    r(NX)		= 0.1d0*cons_local(k,NX-1)+1.0d0*cons_local(k,NX)+0.5d0*cons_local(k,NX+1)-0.1d0*consl(k,NX+1)
	    
	    do i=1,NX-1
	    r(i)		= 0.1d0*cons_local(k, i-1)+1.0d0*cons_local(k,i)+0.5d0*cons_local(k,  i+1)
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

		r(0)		= 0.5d0*cons_local(k,0)+1.0d0*cons_local(k,1  )+0.1d0*cons_local(k,  2)-0.1d0*consr(k, -1)
	    r(NX)		= 0.5d0*cons_local(k,NX)+1.0d0*cons_local(k,NX+1)+0.1d0*cons_local(k,NX+2)-0.5d0*consr(k,NX+1)
	   
	    do i=1,NX-1
	    r(i)		= 0.5d0*cons_local(k, i)+1.0d0*cons_local(k,i+1)+0.1d0*cons_local(k,  i+2)
	    enddo
	    
	    call tridiag(a,b,c,r,temp,total)
	        
	    consr(k,0:NX)=temp(0:NX)

    enddo


		


		do i=-3,NX+3

  d_xx = D_parallel
  d_xy =   0.0d0
  d_xy =   0.0d0
  d_yy =   1.0d0
  	Lr = (0.005/sqrt(1.0d0+0.005d0**2))/(2.0d0*acos(-1.0d0))
	Tr = Lr**2.0d0



			fleft(1,i)	   = -consl(2,i)*d_xx-consl(3,i)*d_xy
			fleft(2,i)	   = -consl(1,i)/Tr
			fleft(3,i)	   = 0.0d0

			fright(1,i)	   = -consr(2,i)*d_xx-consr(3,i)*d_xy
			fright(2,i)	   = -consr(1,i)/Tr
			fright(3,i)	   = 0.0d0


			delphi(i)			= consr(1,i)-consl(1,i)
			delp(i)				= consr(2,i)-consl(2,i)
			delq(i)				= consr(3,i)-consl(3,i)

	!***************************************
	!***** 		Roe solver	       	   *****
	!***************************************

		 !***** Roe solver          *****
		 !The upwind dissipation matrix is |A| = R*|Lambda|*L
		 jacobian(1,1) = abs(sqrt(d_xx/Tr))
		 jacobian(1,2) = 0.0d0
		 jacobian(1,3) = 0.0d0
		 jacobian(2,1) = 0.0d0
		 jacobian(2,2) = abs(sqrt(d_xx/Tr))*d_xx/d_xx
		 jacobian(2,3) = abs(sqrt(d_xx/Tr))*d_xy/d_xx
		 jacobian(3,1) = 0.0d0
		 jacobian(3,2) = 0.0d0
		 jacobian(3,3) = 0.0d0

		dissipation (1,i) = (delphi(i)*jacobian(1,1) + delp(i) *jacobian(1,2)+delq(i)*jacobian(1,3))

		dissipation (2,i) = (delphi(i)*jacobian(2,1) + delp(i) *jacobian(2,2)+delq(i)*jacobian(2,3))

		dissipation (3,i) = (delphi(i)*jacobian(3,1) + delp(i) *jacobian(3,2)+delq(i)*jacobian(3,3))

		
		
		flux(1,i) =  0.5d0*((fright(1,i) +fleft(1,i) -  0.1d0*dissipation(1,i)))	
        flux(2,i) =  0.5d0*((fright(2,i) +fleft(2,i) -  0.1d0*dissipation(2,i)))	
        flux(3,i) =  0.5d0*((fright(3,i) +fleft(3,i) -  0.1d0*dissipation(3,i)))

    !***************************************
	!***** 		Rusanov solver	       *****
	!***************************************

        ! flux(:,i) =  0.5d0*((fright(:,i) +fleft(:,i) - (consr(:,i)-consl(:,i))))
    
    !***************************************
	!***** 		HLL solver	       *****
	!***************************************

   			! SLm(i) 		   = -1.0d0
			! SRp(i) 		   =  1.0d0

			! SLm(i) 		= MIN(-1.0d0 , 1.0d0)
			! SRp(i) 		= MAX(-1.0d0 , 1.0d0)

			! laxmax        	= MAX(-1,1)

			! if(0.0d0 .le. SLm(i)) then

				
			! 	flux(:,i) = fleft(:,i)
			
			! elseif(SLm(i) .le. 0.0d0 .and. 0.0d0 .le. SRp(i))	then
				
			! 	! flux(:,i) =  0.5d0*((fright(:,i) +fleft(:,i) + laxmax*(-consr(:,i)+consl(:,i))))

			! 	flux(:,i) = ( SRp(i)*fleft(:,i) - SLm(i)*fright(:,i) + SLm(i)*SRp(i)*(consr(:,i) - consl(:,i)) )/(SRp(i)-SLm(i))

			! elseif(0.0d0 .ge. SRp(i))then

			! 	flux(:,i) = fright(:,i)
			! endif
        

	enddo



	end subroutine reconstruction5EX


 	!***********************************************************************
	!*****                      Flux in Y-direction                    *****
	!***********************************************************************

	subroutine GY(NX,NY,dy,x,y,cons,flux_y,n_eqn,ghostp)

	implicit none

	integer								:: i, j, k

	integer 							:: NX, NY, ghostp, n_eqn

	double precision					:: dy, x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision 					:: flux_y(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	
	double precision					:: flux_node_y(n_eqn,-ghostp:NY+ghostp), flux_y_half(n_eqn,-ghostp:NY+ghostp)

	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),cons_local(n_eqn,-ghostp:NY+ghostp)

	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx, d_xy,d_yx,d_yy

	double precision, parameter 		:: pi=acos(-1.0d0), beta =  atan(tand(5.0d0)), D_parallel = 1.0



	
	do j =1, NX
		
		do i = -ghostp, NY+ghostp

  d_xx = D_parallel
  d_xy =   0.0d0
  d_xy =   0.0d0
  d_yy =   1.0d0
  	Lr = (0.005/sqrt(1.0d0+0.005d0**2))/(2.0d0*acos(-1.0d0))
	Tr = Lr**2.0d0

			cons_local(1,i) = cons(1,j,i)
			cons_local(2,i) = cons(2,j,i)
			cons_local(3,i) = cons(3,j,i)


			flux_node_y(1,i) =  -cons_local(3,i)*d_yy-cons_local(2,i)*d_yx
			flux_node_y(2,i) =  0.0d0
			flux_node_y(3,i) =  -cons_local(1,i)/Tr
		end do

		call reconstruction5EY (NX,NY,cons_local,flux_y_half, n_eqn,ghostp,x,y,j)

		do k = 1,n_eqn

			do i = 1, NY

				!Finite volume or Reconstruction approach 
	        	flux_y(k,j,i) = -((flux_y_half(k,i)-flux_y_half(k,i-1))/dy)

				! flux_y(k,j,i) = -((4.0/3.0d0)*(flux_y_half(k,i  )-flux_y_half(k,i-1)) +&
				!                 (-1.0/6.0d0)*(flux_node_y(k,i+1)-flux_node_y(k,i-1)))/dy

				flux_y(k,j,i) = -((9.0/8.0d0)*(flux_y_half(k,i  )-flux_y_half(k,i-1)) +&
				                (-1.0/24.0d0)*(flux_y_half(k,i+1)-flux_y_half(k,i-2)))/dy


				! flux_y(k,j,i) = -((3.0d0/2.0d0)*(flux_y_half(k,i)-flux_y_half(k,i-1))/dy + &
				! 			    (1.0d0/30.0d0)*(flux_y_half(k,i+1)-flux_y_half(k,i-2))/dy &
				! 			    +(-3.0d0/10.0d0)*(flux_node_y(k,i+1)-flux_node_y(k,i-1))/dy)

				flux_y(k,j,i) = -((75.0d0/64.0d0)*(flux_y_half(k,i)-flux_y_half(k,i-1))/dy - (25.0d0/384.0d0)*(flux_y_half(k,i+1)-flux_y_half(k,i-2))/dy &
				              	+ (3.0d0/640.0d0)*(flux_y_half(k,i+2)-flux_y_half(k,i-3))/dy)




			enddo
		enddo
	enddo
	





	end subroutine GY


	!***************************************************************************************************************
	!***** 					Roe, HLL or Rusanov + 5th order compact upwind reconstruction or interplation     *****
	!*****											in Y-direction											   *****
	!***************************************************************************************************************

	subroutine reconstruction5EY (NX,NY,cons_local,flux, n_eqn,ghostp,x,y,j)

	implicit none 		

	integer, intent(in)					:: j
	integer								:: i,k	
	integer								:: NX,NY, n_eqn, ghostp, total
	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision					:: cons_local(n_eqn,-ghostp:NY+ghostp),flux(n_eqn,-ghostp:NY+ghostp)

	double precision 					:: consl(n_eqn,-ghostp:NY+ghostp), consr(n_eqn,-ghostp:NY+ghostp)

	double precision   					:: fleft(n_eqn,-ghostp:NY+ghostp), fright(n_eqn,-ghostp:NY+ghostp)

	double precision					:: righteigen(3,3), lefteigen(3,3)

	double precision					:: dissipation(n_eqn,-2:NY+2),jacobian(3,3)

	double precision 					:: lambda1(-ghostp:NY+ghostp), lambda2(-ghostp:NY+ghostp)
	double precision  					:: delphi(-ghostp:NY+ghostp), delp(-ghostp:NY+ghostp), delq(-ghostp:NY+ghostp)
	double precision 					:: alpha1, alpha2
	double precision 					:: SLm(-ghostp:NY+ghostp), SRp(-ghostp:NY+ghostp),laxmax

	double precision					:: a(0:NY+1), b(0:NY+1), c(0:NY+1), r(0:NY+1),temp(0:NY+1)

	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx, d_xy,d_yx,d_yy

	double precision					:: xf, yf

	double precision					:: xl, xr, yl, yr

	double precision, parameter			:: epsilon = 1.0d-23
    double precision 					:: recon_poly5_r(0:3), recon_poly5_l(0:3)
    double precision 					:: alpha_l(0:3),alpha_r(0:3)
    double precision 					:: beta_w(0:3), omega_l(0:3), omega_r(0:3)
   
    double precision, parameter 			:: Constant=1.0d0,p =1.0d0

    double precision, parameter 		:: pi=acos(-1.0d0), beta =  atan(tand(5.0d0)), D_parallel = 1.0




		do k=1,n_eqn

        	do i =-3,NY+3

        	   ! Interpolation formula		
              ! consl(k,i) = -1.0d0/8.0d0*cons_local(k,i-1)+6.0d0/8.0d0*cons_local(k,i  )+3.0d0/8.0d0*cons_local(k,i+1)
              ! consr(k,i) =  3.0d0/8.0d0*cons_local(k,i  )+6.0d0/8.0d0*cons_local(k,i+1)-1.0d0/8.0d0*cons_local(k,i+2)


        	! Reconstruction formula
   !      	consl(k,i) = (2.0d0/60.0d0)*(cons_local(k,i-2)) + (-13.0d0/60.0d0)*(cons_local(k,i-1)) + (47.0d0/60.0d0)*(cons_local(k,i)) + &
			! 		  (27.0d0/60.0d0)*(cons_local(k,i+1)) + (-3.0d0/60.0d0)*(cons_local(k,i+2))

			! consr(k,i) = (-3.0d0/60.0d0)*(cons_local(k,i-1)) + (27.0d0/60.0d0)*(cons_local(k,i)) + (47.0d0/60.0d0)*(cons_local(k,i+1)) + &
			! 	      (-13.0d0/60.0d0)*(cons_local(k,i+2)) + (2.0d0/60.0d0)*(cons_local(k,i+3))

   !      	yl = (2.0d0/60.0d0)*(y(i-2)) + (-13.0d0/60.0d0)*(y(i-1)) + (47.0d0/60.0d0)*(y(i)) + &
			! 		  (27.0d0/60.0d0)*(y(i+1)) + (-3.0d0/60.0d0)*(y(i+2))

			! yr = (-3.0d0/60.0d0)*(y(i-1)) + (27.0d0/60.0d0)*(y(i)) + (47.0d0/60.0d0)*(y(i+1)) + &
			! 	      (-13.0d0/60.0d0)*(y(i+2)) + (2.0d0/60.0d0)*(y(i+3))

  	       ! Interpolation formula		
			consl(k,i) = (3.0d0/128.0d0)*(cons_local(k,i-2)) - (5.0d0/32.0d0)*(cons_local(k,i-1)) + (45.0d0/64.0d0)*(cons_local(k,i)) + &
					     (15.0d0/32.0d0)*(cons_local(k,i+1)) - (5.0d0/128.0d0)*(cons_local(k,i+2))
			
			consr(k,i) = (-5.0d0/128.0d0)*(cons_local(k,i-1)) + (15.0d0/32.0d0)*(cons_local(k,i)) + (45.0d0/64.0d0)*(cons_local(k,i+1)) - &
					     (5.0d0/32.0d0)*(cons_local(k,i+2)) + (3.0d0/128.0d0)*(cons_local(k,i+3))

			! consl(k,i) = (1.0d0/60.0d0)*(cons_local(k,i-2)) - (8.0d0/60.0d0)*(cons_local(k,i-1)) + (37.0d0/60.0d0)*(cons_local(k,i)) +&
   !     		 		  (37.0d0/60.0d0)*(cons_local(k,i+1)) - (8.0d0/60.0d0)*cons_local(k,i+2)+ (1.0d0/60.0d0)*(cons_local(k,i+3))
	  !   	consr(k,i) = (1.0d0/60.0d0)*(cons_local(k,i-2)) - (8.0d0/60.0d0)*(cons_local(k,i-1)) + (37.0d0/60.0d0)*(cons_local(k,i)) +&
   !     			      (37.0d0/60.0d0)*(cons_local(k,i+1)) - (8.0d0/60.0d0)*cons_local(k,i+2)+ (1.0d0/60.0d0)*(cons_local(k,i+3))

		 !  	yl = (3.0d0/128.0d0)*(y(i-2)) - (5.0d0/32.0d0)*(y(i-1)) + (45.0d0/64.0d0)*(y(i)) + &
			! 		     (15.0d0/32.0d0)*(y(i+1)) - (5.0d0/128.0d0)*(y(i+2))
			
			! yr = (-5.0d0/128.0d0)*(y(i-1)) + (15.0d0/32.0d0)*(y(i)) + (45.0d0/64.0d0)*(y(i+1)) - &
			! 		     (5.0d0/32.0d0)*(y(i+2)) + (3.0d0/128.0d0)*(y(i+3))

			! WCNS interpolation

!  				beta_w(0) 	= (13.0d0/12.0d0)*(cons_local(k,i)-2.0d0*cons_local(k,i+1)+cons_local(k,i+2))**2+ &
! 	            				(1.0/4.0d0)*(3.0d0*cons_local(k,i)-4.0d0*cons_local(k,i+1)+cons_local(k,i+2))**2

! 	            beta_w(1) 	= (13.0d0/12.0d0)*(cons_local(k,i-1)-2.d0*cons_local(k,i)+cons_local(k,i+1))**2+ &
! 	            				(1.0/4.0d0)*(cons_local(k,i-1)-cons_local(k,i+1))**2

! 	            beta_w(2)		= (13.0d0/12.0d0)*(cons_local(k,i-2)-2.d0*cons_local(k,i-1)+cons_local(k,i))**2 &
! 	            				+(1.0/4.0d0)*(cons_local(k,i-2)-4.0d0*cons_local(k,i-1)+3.d0*cons_local(k,i))**2


! 	            ! recon_poly5_l(0) = 3.0d0/8.0d0*cons_local(k,i)+6.0d0/8.0d0*cons_local(k,i+1)-1.0d0/8.0d0*cons_local(k,i+2)
! 	            ! recon_poly5_l(1) = -1.0d0/8.0d0*cons_local(k,i-1)+6.0d0/8.0d0*cons_local(k,i)+3.0d0/8.0d0*cons_local(k,i+1)
! 	            ! recon_poly5_l(2) = 3.0d0/8.0d0*cons_local(k,i-2)-10.0d0/8.0d0*cons_local(k,i-1)+15.0d0/8.0d0*cons_local(k,i)

! 	           	recon_poly5_l(0) =1.0d0/8.0d0*( 3.0d0*cons_local(k,i  )+ 6.0d0*cons_local(k,i+1)-1.0d0 *cons_local(k,i+2))
! 	            recon_poly5_l(1) =1.0d0/8.0d0*(-1.0d0*cons_local(k,i-1)+ 6.0d0*cons_local(k,i  )+3.0d0 *cons_local(k,i+1))
! 	            recon_poly5_l(2) =1.0d0/8.0d0*( 3.0d0*cons_local(k,i-2)-10.0d0*cons_local(k,i-1)+15.0d0*cons_local(k,i  ))

! 	! Jiang-Shu linear weights
! 	            ! alpha_l(0) 		= (5.0/16.0d0)
! 	            ! alpha_l(1) 		= (10.0/16.0d0)
! 	            ! alpha_l(2) 		= (1.0/16.0d0)
! ! ! Borges
! 	        alpha_l(0) = 5.0d0/16.0d0  * (Constant + (abs(beta(0)-beta(2)) / (epsilon + beta(0)))**p)
!             alpha_l(1) = 10.0d0/16.0d0 * (Constant + (abs(beta(0)-beta(2)) / (epsilon + beta(1)))**p)
!             alpha_l(2) = 1.0d0/16.0d0  * (Constant + (abs(beta(0)-beta(2)) / (epsilon + beta(2)))**p)

! 	            omega_l(0)=alpha_l(0)/(alpha_l(0)+alpha_l(1)+alpha_l(2))
! 	            omega_l(1)=alpha_l(1)/(alpha_l(0)+alpha_l(1)+alpha_l(2))
! 	            omega_l(2)=alpha_l(2)/(alpha_l(0)+alpha_l(1)+alpha_l(2))

! 	            consl(k,i) = omega_l(0) * recon_poly5_l(0) + omega_l(1) * recon_poly5_l(1) + omega_l(2) * recon_poly5_l(2)



! 	! right of i+1/2
! 	            beta_w(0)	= (13.0d0/12.0d0)*(cons_local(k,i+1)-2.0d0*cons_local(k,i+2)+cons_local(k,i+3))**2+ &
! 	            						(1.0/4.0d0)*(3.0d0*cons_local(k,i+1)-4.0d0*cons_local(k,i+2)+cons_local(k,i+3))**2

! 	            beta_w(1)	= (13.0d0/12.0d0)*(cons_local(k,i)-2.d0*cons_local(k,i+1)+cons_local(k,i+2))**2+ &
! 	           							(1.0/4.0d0)*(cons_local(k,i)-cons_local(k,i+2))**2

! 	            beta_w(2) = (13.0d0/12.0d0)*(cons_local(k,i-1)-2.d0*cons_local(k,i)+cons_local(k,i+1))**2 &
! 	            						+(1.0/4.0d0)*(cons_local(k,i-1)-4.0d0*cons_local(k,i)+3.d0*cons_local(k,i+1))**2


! 				! recon_poly5_r(0) = 15.0d0/8.0d0*cons_local(k,i+1)-10.0d0/8.0d0*cons_local(k,i+2)+3.0d0/8.0d0*cons_local(k,i+3)
! 	   !          recon_poly5_r(1) = 3.0d0/8.0d0*cons_local(k,i)+6.0d0/8.0d0*cons_local(k,i+1)-1.0d0/8.0d0*cons_local(k,i+2)
! 	   !          recon_poly5_r(2) = -1.0d0/8.0d0*cons_local(k,i-1)+6.0d0/8.0d0*cons_local(k,i)+3.0d0/8.0d0*cons_local(k,i+1)

! 	   			recon_poly5_r(0) = 1.0d0/8.0d0*(15.0d0*cons_local(k,i+1)-10.0d0*cons_local(k,i+2)+3.0d0*cons_local(k,i+3))
! 	            recon_poly5_r(1) = 1.0d0/8.0d0*( 3.0d0*cons_local(k,i  )+ 6.0d0*cons_local(k,i+1)-1.0d0*cons_local(k,i+2))
! 	            recon_poly5_r(2) = 1.0d0/8.0d0*(-1.0d0*cons_local(k,i-1)+ 6.0d0*cons_local(k,i  )+3.0d0*cons_local(k,i+1))

! 	! Jiang-Shu linear weights
! 	            ! alpha_r(0) = (1.0/16.0d0)
! 	            ! alpha_r(1) = (10.0/16.0d0)
! 	            ! alpha_r(2) = (5.0/16.0d0)

! ! Borges
! 	        alpha_r(0) = 1.0d0/16.0d0  * (Constant + (abs(beta(0)-beta(2)) / (epsilon + beta(0)))**p)
!             alpha_r(1) = 10.0d0/16.0d0 * (Constant + (abs(beta(0)-beta(2)) / (epsilon + beta(1)))**p)
!             alpha_r(2) = 5.0d0/16.0d0  * (Constant + (abs(beta(0)-beta(2)) / (epsilon + beta(2)))**p)

! 	! Non-linear weights
! 	            omega_r(0)=alpha_r(0)/(alpha_r(0)+alpha_r(1)+alpha_r(2))
! 	            omega_r(1)=alpha_r(1)/(alpha_r(0)+alpha_r(1)+alpha_r(2))
! 	            omega_r(2)=alpha_r(2)/(alpha_r(0)+alpha_r(1)+alpha_r(2))

	            
! 	            consr(k,i) = omega_r(0) * recon_poly5_r(0) + omega_r(1) * recon_poly5_r(1) + omega_r(2) * recon_poly5_r(2) 


        
         	enddo
        enddo

    ! Upwind Comapct Reconstruction
		! do k = 1,n_eqn
		    	

		!     	do i = 0, NY
		! 	      a(i) = 1.0d0/2.0d0
		! 	      b(i) = 1.0d0
		! 	      c(i) = 1.0d0/6.0d0
		!    		end do
		    
		!     i = 0  ; a(i) = 0.0d+00
		!     i = NY  ; c(i) = 0.0d+00
		!     total = NY+1
		    
		!     r(0)		= (1.0d0/18.0d0)*cons_local(k, -1)+(19.0d0/18.0d0)*cons_local(k,0)+(10.0d0/18.0d0)*cons_local(k,  1)-(1.0d0/2.0d0)*consl(k, -1)
		!     r(NY)		= (1.0d0/18.0d0)*cons_local(k,NY-1)+(19.0d0/18.0d0)*cons_local(k,NY)+(10.0d0/18.0d0)*cons_local(k,NY+1)-(1.0d0/6.0d0)*consl(k,NY+1)
		    
		!     do i=1,NY-1
		!     r(i)		= (1.0d0/18.0d0)*cons_local(k, i-1)+(19.0d0/18.0d0)*cons_local(k,i)+(5.0d0/9.0d0)*cons_local(k,  i+1)
		!     enddo
		    
		!     call tridiag(a,b,c,r,temp,total)
		        
		!     consl(k,0:NY)=temp(0:NY) 

		! 	    do i = 0, NY
		! 	      a(i) = 1.0d0/6.0d0
		! 	      b(i) = 1.0d0
		! 	      c(i) = 1.0d0/2.0d0
		! 	    end do

		!     i = 0  ; a(i) = 0.0d+00
		!     i = NY  ; c(i) = 0.0d+00

		!     r(0)		= (10.0d0/18.0d0)*cons_local(k,0)+(19.0d0/18.0d0)*cons_local(k,  1)+(1.0d0/18.0d0)*cons_local(k,  2)-(1.0d0/6.0d0)*consr(k, -1)
		!     r(NY)		= (10.0d0/18.0d0)*cons_local(k,NY)+(19.0d0/18.0d0)*cons_local(k,NY+1)+(1.0d0/18.0d0)*cons_local(k,NY+2)-(1.0d0/2.0d0)*consr(k,NY+1)
		   
		!     do i=1,NY-1
		!     r(i)		= (10.0d0/18.0d0)*cons_local(k,i)+(19.0d0/18.0d0)*cons_local(k,i+1)+(1.0d0/18.0d0)*cons_local(k,i+2)
		!     enddo
		    
		!     call tridiag(a,b,c,r,temp,total)
		        
		!     consr(k,0:NY)=temp(0:NY)

		! enddo


	 ! Compact reconstruction low dissipation
		    ! do k = 1,n_eqn
		    	

		    ! 	do i = 0, NY
			   !    a(i) = 5.0d0/20.0d0
			   !    b(i) = 12.0d0/20.0d0
			   !    c(i) = 3.0d0/20.0d0
		   	! 	 end do
		    
		    ! i = 0  ; a(i) = 0.0d+00
		    ! i = NY  ; c(i) = 0.0d+00
		    ! total = NY+1
		    
		    ! r(0)		= (3.0d0/120.0d0)*cons_local(k, -1)+(67.0d0/120.0d0)*cons_local(k,0)+(49.0d0/120.0d0)*cons_local(k,  1)+(1.0d0/120.0d0)*cons_local(k,  2)-(5.0d0/20.0d0)*consl(k, -1)
		    ! r(NY)		= (3.0d0/120.0d0)*cons_local(k,NY-1)+(67.0d0/120.0d0)*cons_local(k,NY)+(49.0d0/120.0d0)*cons_local(k,NY+1)+(1.0d0/120.0d0)*cons_local(k,  NY+2)-(3.0d0/20.0d0)*consl(k,NY+1)
		    
		    ! do i=1,NY-1
		    ! r(i)		= (3.0d0/120.0d0)*cons_local(k,i-1)+(67.0d0/120.0d0)*cons_local(k,i)+(49.0d0/120.0d0)*cons_local(k,i+1)+(1.0d0/120.0d0)*cons_local(k,  i+2)
		    ! enddo
		    
		    ! call tridiag(a,b,c,r,temp,total)
		        
		    ! consl(k,0:NY)=temp(0:NY) 

			   !  do i = 0, NY
			   !    a(i) = 3.0d0/20.0d0
			   !    b(i) = 12.0d0/20.0d0
			   !    c(i) = 5.0d0/20.0d0
			   !  end do

		    ! i = 0  ; a(i) = 0.0d+00
		    ! i = NY  ; c(i) = 0.0d+00

		    ! r(0)		= (1.0d0/120.0d0)*cons_local(k, -1)+(49.0d0/120.0d0)*cons_local(k,0)+(67.0d0/120.0d0)*cons_local(k,  1)+(3.0d0/120.0d0)*cons_local(k,  2)-(3.0d0/20.0d0)*consr(k, -1)
		    ! r(NY)		= (1.0d0/120.0d0)*cons_local(k,NY-1)+(49.0d0/120.0d0)*cons_local(k,NY)+(67.0d0/120.0d0)*cons_local(k,NY+1)+(3.0d0/120.0d0)*cons_local(k,  NY+2)-(5.0d0/20.0d0)*consr(k,NY+1)
		   
		    ! do i=1,NY-1
		    ! r(i)		= (1.0d0/120.0d0)*cons_local(k,i-1)+(49.0d0/120.0d0)*cons_local(k,i)+(67.0d0/120.0d0)*cons_local(k,i+1)+(3.0d0/120.0d0)*cons_local(k,  i+2)
		    ! enddo
		    
		    ! call tridiag(a,b,c,r,temp,total)
		        
		    ! consr(k,0:NY)=temp(0:NY)

		    ! enddo

! Upwind interpolation
	 do k = 1,n_eqn
	    	

	    	do i = 0, NY
		      a(i) = 0.5d0
		      b(i) = 1.0d0
		      c(i) = 0.1d0
	   		 end do
	    
	    i = 0  ; a(i) = 0.0d+00
	    i = NY  ; c(i) = 0.0d+00
	    total = NY+1
	    
	    r(0)		= 0.1d0*cons_local(k, -1)+1.0d0*cons_local(k,0)+0.5d0*cons_local(k,  1)-0.5d0*consl(k, -1)
	    r(NY)		= 0.1d0*cons_local(k,NY-1)+1.0d0*cons_local(k,NY)+0.5d0*cons_local(k,NY+1)-0.1d0*consl(k,NY+1)
	    
	    do i=1,NY-1
	    r(i)		= 0.1d0*cons_local(k, i-1)+1.0d0*cons_local(k,i)+0.5d0*cons_local(k,  i+1)
	    enddo
	    
	    call tridiag(a,b,c,r,temp,total)
	        
	    consl(k,0:NY)=temp(0:NY) 

		    do i = 0, NY
		      a(i) = 0.1d0
		      b(i) = 1.0d0
		      c(i) = 0.5d0
		    end do

	    i = 0  ; a(i) = 0.0d+00
	    i = NY  ; c(i) = 0.0d+00

		r(0)		= 0.5d0*cons_local(k,0)+1.0d0*cons_local(k,1  )+0.1d0*cons_local(k,  2)-0.1d0*consr(k, -1)
	    r(NY)		= 0.5d0*cons_local(k,NY)+1.0d0*cons_local(k,NY+1)+0.1d0*cons_local(k,NY+2)-0.5d0*consr(k,NY+1)
	   
	    do i=1,NY-1
	    r(i)		= 0.5d0*cons_local(k, i)+1.0d0*cons_local(k,i+1)+0.1d0*cons_local(k,  i+2)
	    enddo
	    
	    call tridiag(a,b,c,r,temp,total)
	        
	    consr(k,0:NY)=temp(0:NY)

	    enddo

		do i=-3,NY+3

  d_xx = D_parallel
  d_xy =   0.0d0
  d_xy =   0.0d0
  d_yy =   1.0d0
  	Lr = (0.005/sqrt(1.0d0+0.005d0**2))/(2.0d0*acos(-1.0d0))
	Tr = Lr**2.0d0

			fleft(1,i)	   = -(consl(3,i)*d_yy+consl(2,i)*d_yx)
			fleft(2,i)	   = 0.0d0
			fleft(3,i)	   = -consl(1,i)/Tr

			fright(1,i)	   = -(consr(3,i)*d_yy+consr(2,i)*d_yx)
			fright(2,i)	   = 0.0d0
			fright(3,i)	   = -consr(1,i)/Tr


			delphi(i)			= consr(1,i)-consl(1,i)
			delp(i)				= consr(2,i)-consl(2,i)
			delq(i)				= consr(3,i)-consl(3,i)

	!***************************************
	!***** 		Roe solver	       	   *****
	!***************************************

		 !***** Roe solver          *****
		 !The upwind dissipation matrix is |A| = R*|Lambda|*L
		 jacobian(1,1) = abs(sqrt(d_yy/Tr))
		 jacobian(1,2) = 0.0d0
		 jacobian(1,3) = 0.0d0
		 jacobian(2,1) = 0.0d0
		 jacobian(2,2) = 0.0d0
		 jacobian(2,3) = 0.0d0
		 jacobian(3,1) = 0.0d0
		 jacobian(3,2) = abs(sqrt(d_yy/Tr))*d_yx/d_yy
		 jacobian(3,3) = abs(sqrt(d_yy/Tr))*d_yy/d_yy

		dissipation (1,i) = (delphi(i)*jacobian(1,1) + delp(i) *jacobian(1,2)+delq(i)*jacobian(1,3))

		dissipation (2,i) = (delphi(i)*jacobian(2,1) + delp(i) *jacobian(2,2)+delq(i)*jacobian(2,3))

		dissipation (3,i) = (delphi(i)*jacobian(3,1) + delp(i) *jacobian(3,2)+delq(i)*jacobian(3,3))

		
		
		flux(1,i) =  0.5d0*((fright(1,i) +fleft(1,i) -  0.1d0*dissipation(1,i)))	
        flux(2,i) =  0.5d0*((fright(2,i) +fleft(2,i) -  0.1d0*dissipation(2,i)))	
        flux(3,i) =  0.5d0*((fright(3,i) +fleft(3,i) -  0.1d0*dissipation(3,i)))
        

	enddo




	end subroutine reconstruction5EY



 	!***********************************************************************
	!*****                       Source term 			               *****
	!***********************************************************************

 	subroutine sourceterm(NX,NY,x,y,cons,flux_x,flux_y,residual,n_eqn,ghostp)

 		implicit none

 		integer								:: NX, NY, n_eqn, ghostp
		double precision 					:: X(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp), L
	
		
		integer								:: i,j
		
		
		double precision,allocatable		:: sourcet(:,:,:)	
		double precision 					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision 					:: residual(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision 					:: flux_x(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),flux_y(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision              		:: d_xx, d_xy,d_yx,d_yy
		double precision              		:: nu,Lr,Tr,Tr_1,Tr_2

		double precision, parameter 		:: pi=acos(-1.0d0), beta =  atan(tand(5.0d0)), D_parallel = 1.0
			
		allocate (sourcet(3,1:NX,1:NY))
			do j = 1, NY

				do i = 1, NX
				  
				  d_xx = D_parallel
				  d_xy =   0.0d0
				  d_xy =   0.0d0
				  d_yy =   1.0d0
				  Lr = (0.005/sqrt(1.0d0+0.005d0**2))/(2.0d0*acos(-1.0d0))
				  Tr = Lr**2.0d0


				sourcet(1,i,j)		= sin(x(i)*pi)*sin(y(j)*pi*2000.0d0)*(pi**2)*(1.0d0+2000.0d0**2)
				sourcet(2,i,j)		= -cons(2,i,j)/Tr
				sourcet(3,i,j)		= -cons(3,i,j)/Tr

					residual(:,i,j)		= flux_x(:,i,j) + flux_y(:,i,j) + sourcet(:,i,j)
				enddo

		enddo
	end subroutine sourceterm


	!***********************************************************************
	!*****                       Output 			                   *****
	!***********************************************************************

	subroutine output(NX,NY,x,y,cons,cons_exact,n_eqn,ghostp)

	integer 				:: i, j

	integer					:: NX,NY, n_eqn, ghostp

	double precision 		:: X(-ghostp:NX+ghostp),Y(-ghostp:NY+ghostp),cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: cons_exact(n_eqn,1:NX,1:NY)


	open(unit=1,file='Result-5a.plt')
	
	write(1,*) 'TITLE=" "Normal"'
    write(1,*) 'VARIABLES = "x","y","Pote","vx","vy","exac"'
	write(1,*) "ZONE I=",NX," J=",NY," F=POINT"


    do j = 1, NY
        do i = 1, NX

          write(1,'(6F25.8)') x(i), y(j), cons(1,i,j), cons(2,i,j), cons(3,i,j), cons_exact(1,i,j)

	  enddo
    enddo
	
    close(1)



 	end subroutine output

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
	    integer                    	   :: i, j
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




