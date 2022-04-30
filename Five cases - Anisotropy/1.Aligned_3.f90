! 
!**********************************************************************
! This program solves 2D anisotropic Diffusion equation with Tensor
! Generalized and non-aligned tensor
!**********************************************************************


	program hyperpoi2d
	
		implicit none


		double precision, parameter 		:: pi=acos(-1.0d0)

		integer, parameter					:: NTMAX=40000000,l=2
		integer 							:: i, j, z, k, NT
		
		integer, parameter					:: NX = 8*l,NY = 8*l, n_eqn = 3, ghostp = 15	

		double precision, parameter			:: l1_target = 1.0d-14, epsilon = 1.0d0-14
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


		call initialconditions(NX,NY,ghostp,cons,n_eqn,x,y)

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


	subroutine initialconditions(NX,NY,ghostp,cons,n_eqn,x,y)


	implicit none

	integer 							:: i, j, k
	integer								:: NX, NY, n_eqn, ghostp

	double precision,intent(in)			::x(-ghostp:NX+ghostp),y(-ghostp:NY+ghostp)
	
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


		do j = 1, NY
			
			do i = 1, NX
				
				cons_exact(1,i,j) = (1.0d0/(2.0d0*pi**2))*sin(pi*(x(i)))*sin(pi*y(j))
			
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

	double precision, parameter 		:: pi=acos(-1.0d0), beta =  00.0d0, D_parallel = 1.0d6

	double precision					:: x_velocity, y_velocity
	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx, d_xy,d_yx,d_yy



	dt = 1.0d010
		
		do i = 1, NX

			do j = 1, NY

				  d_xx =   D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
				  d_xy =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)				  
				  d_yx =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yy =   D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)

  					Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
					Tr = (2.0d0*(Lr**2.0d0))/(d_xx**1.0d0+d_yy**1.0d0+2.0d0*d_xy)
				  
				  x_velocity =  DABS(sqrt(d_xx/Tr))
				  y_velocity =  DABS(sqrt(d_yy/Tr))

				  dtnew = min( (x(i+1)-x(i))/x_velocity, (y(j+1)-y(j))/y_velocity)
				    
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
	double precision, parameter 		:: pi=acos(-1.0d0), beta =  00.0d0, D_parallel = 1.0d6
	integer								:: NX, NY, n_eqn, ghostp

	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)	
	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)



 !    ! Left boundary conditions (actually left)

	do j = 1,NY
		
		i=0     
	    
	    cons(1,i,j )   = (8.0d0*((1.0d0/(2.0d0*pi**2))*sin(pi*0.5d0*(x(i)+x(i+1)))*sin(pi*y(j))) - 6.0d0*cons(1,i+1,j) + cons(1,i+2,j)) * (1.0d0/3.0d0)

	    cons(2,i,j )   = cons(2,i+3,j)-3.0d0*cons(2,i+2,j)+3.0d0*cons(2,i+1,j)
        
        cons(3,i,j )   = cons(3,i+3,j)-3.0d0*cons(3,i+2,j)+3.0d0*cons(3,i+1,j)
    enddo

    do j = 1,NY

    	do i = -1,-ghostp,-1
        
        	cons(:,i,j)   = cons(:,i+3,j)-3.0d0*cons(:,i+2,j)+3.0d0*cons(:,i+1,j)

        enddo
    enddo


	! Right boundary conditions (actually right)


	do j = 1, NY
		
		i = NY

		cons(1,i+1,j) = (8.0d0*((1.0d0/(2.0d0*pi**2))*sin(pi*0.5d0*(x(i)+x(i+1)))*sin(pi*y(j))) +1.0d0*cons(1,i-1,j) -6.0d0* cons(1,i,j)) * (1.0d0/3.0d0)
        
        cons(2,i+1,j) = 3.0d0*cons(2,i,j)   - 3.0d0*cons(2,i-1,j)+cons(2,i-2,j)
        cons(3,i+1,j) = 3.0d0*cons(3,i,j)   - 3.0d0*cons(3,i-1,j)+cons(3,i-2,j)
        
    enddo

    do j = 1, NY

    	do i= NY, NY+6

        cons(:,i+2,j) = 3.0d0*cons(:,i+1,j) - 3.0d0*cons(:,i,j)  +cons(:,i-1,j)

        enddo
    enddo


    ! Bottom boundary conditions (actually bottom)

	do i = 1,NX
		j=0     
	    cons(1,i,j )   = (8.0d0*((1.0d0/(2.0d0*pi**2))*sin(pi*0.5d0*(y(j)+y(j+1)))*sin(pi*x(i))) - 6.0d0*cons(1,i,j+1) + cons(1,i,j+2)) * (1.0d0/3.0d0)

	    cons(2,i,j )   = cons(2,i,j+3)-3.0d0*cons(2,i,j+2)+3.0d0*cons(2,i,j+1)
        cons(3,i,j )   = cons(3,i,j+3)-3.0d0*cons(3,i,j+2)+3.0d0*cons(3,i,j+1)
    enddo

    do i = 1,NX  
    	do j = -1,-ghostp,-1
        
        	cons(:,i,j)   = cons(:,i,j+3)-3.0d0*cons(:,i,j+2)+3.0d0*cons(:,i,j+1)

        enddo
    enddo


	! Top boundary conditions (actually top)
	do i = 1, NX
		 j = NX

		cons(1,i,j+1) = (8.0d0*((1.0d0/(2.0d0*pi**2))*sin(pi*0.5d0*(y(j)+y(j+1)))*sin(pi*x(i))) +1.0d0*cons(1,i,j-1) -6.0d0* cons(1,i,j)) * (1.0d0/3.0d0)
        
        cons(2,i,j+1) = 3.0d0*cons(2,i,j)   - 3.0d0*cons(2,i,j-1)+cons(2,i,j-2)
        cons(3,i,j+1) = 3.0d0*cons(3,i,j)   - 3.0d0*cons(3,i,j-1)+cons(3,i,j-2)
        
    enddo

    do i = 1, NX

    	do j= NX, NX+6

        cons(:,i,j+2) = 3.0d0*cons(:,i,j+1) - 3.0d0*cons(:,i,j)  +cons(:,i,j-1)

        enddo
    enddo



	

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
	double precision, parameter 		:: pi=acos(-1.0d0), beta =  00.0d0, D_parallel = 1.0d6

	do j =1, NY

		do i = -ghostp, NX+ghostp

				  d_xx =   D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
				  d_xy =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yx =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yy =   D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)

					Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))


					Tr = (2.0d0*(Lr**2.0d0))/(d_xx**1.0d0+d_yy**1.0d0+2.0d0*d_xy)

			cons_local(1,i) = cons(1,i,j)
			cons_local(2,i) = cons(2,i,j)
			cons_local(3,i) = cons(3,i,j)

			flux_node_x(1,i) = -cons_local(2,i)*d_xx-cons_local(3,i)*d_xy
			flux_node_x(2,i) = -cons_local(1,i)
			flux_node_x(3,i) = 0.0d0
		end do

		call reconstruction5EX (NX,NY,cons_local,flux_x_half, n_eqn,ghostp,x,y,j)

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

    double precision, parameter 		:: pi=acos(-1.0d0), beta =  00.0d0, D_parallel = 1.0d6
   
    double precision, parameter 		:: Constant=1.0d0,p =1.0d0

    double precision 					:: d1, d2 ,d3, d4, d5


		do k=1,n_eqn

        	do i =-4,NX+4

        	 ! Interpolation formula		
              consl(k,i) = -1.0d0/8.0d0*cons_local(k,i-1)+6.0d0/8.0d0*cons_local(k,i  )+3.0d0/8.0d0*cons_local(k,i+1)
              consr(k,i) =  3.0d0/8.0d0*cons_local(k,i  )+6.0d0/8.0d0*cons_local(k,i+1)-1.0d0/8.0d0*cons_local(k,i+2)
              

         	enddo
        enddo

   

		do i=-4,NX+4

				  d_xx =   D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
				  d_xy =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yx =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yy =   D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)


				  	Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
	
					Tr = (2.0d0*(Lr**2.0d0))/(d_xx**1.0d0+d_yy**1.0d0+2.0d0*d_xy)


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

	integer								:: i, j, k

	integer 							:: NX, NY, ghostp, n_eqn

	double precision					:: dy, x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision 					:: flux_y(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	
	double precision					:: flux_node_y(n_eqn,-ghostp:NY+ghostp), flux_y_half(n_eqn,-ghostp:NY+ghostp)

	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),cons_local(n_eqn,-ghostp:NY+ghostp)

	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx, d_xy,d_yx,d_yy

	double precision, parameter 		:: pi=acos(-1.0d0), beta =  00.0d0, D_parallel = 1.0d6



	
	do j =1, NX
		
		do i = -ghostp, NY+ghostp

				  d_xx =   D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
				  d_xy =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yx =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yy =   D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)

				  	Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
					Tr = (2.0d0*(Lr**2.0d0))/(d_xx**1.0d0+d_yy**1.0d0+2.0d0*d_xy)

			cons_local(1,i) = cons(1,j,i)
			cons_local(2,i) = cons(2,j,i)
			cons_local(3,i) = cons(3,j,i)


			flux_node_y(1,i) =  -cons_local(3,i)*d_yy-cons_local(2,i)*d_yx
			flux_node_y(2,i) =  0.0d0
			flux_node_y(3,i) =  -cons_local(1,i)
		end do

		call reconstruction5EY (NX,NY,cons_local,flux_y_half, n_eqn,ghostp,x,y,j)

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

    double precision, parameter 		:: pi=acos(-1.0d0), beta =  00.0d0, D_parallel = 1.0d6

    double precision 					:: d1, d2 ,d3, d4, d5





		do k=1,n_eqn

        	do i =-4,NY+4

        	 ! Interpolation formula		
              consl(k,i) = -1.0d0/8.0d0*cons_local(k,i-1)+6.0d0/8.0d0*cons_local(k,i  )+3.0d0/8.0d0*cons_local(k,i+1)
              consr(k,i) =  3.0d0/8.0d0*cons_local(k,i  )+6.0d0/8.0d0*cons_local(k,i+1)-1.0d0/8.0d0*cons_local(k,i+2)
        
         	enddo
        enddo





		do i=-4,NY+4

				  d_xx =   D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
				  d_xy =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yx =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yy =   D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)


				  	Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))

					Tr = (2.0d0*(Lr**2.0d0))/(d_xx**1.0d0+d_yy**1.0d0+2.0d0*d_xy)

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

 		integer								:: NX, NY, n_eqn, ghostp
		double precision 					:: X(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp), L
	
		
		integer								:: i,j
		
		
		double precision					:: sourcet(3,1:NX,1:NY)	
		double precision 					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision 					:: residual(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision 					:: flux_x(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),flux_y(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision              		:: d_xx, d_xy,d_yx,d_yy
		double precision              		:: nu,Lr,Tr,Tr_1,Tr_2

		double precision, parameter 		:: pi=acos(-1.0d0), beta =  00.0d0, D_parallel = 1.0d6

		double precision					:: A, B, C, D, E, F,frac, my_function

	
			

			do j = 1, NY

				do i = 1, NX


				  sourcet(1,i,j)      = my_function(D_parallel,x(i),y(j))
				  
				  d_xx =   D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
				  d_xy =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yx =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yy =   D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)

					Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
					Tr = (2.0d0*(Lr**2.0d0))/(d_xx**1.0d0+d_yy**1.0d0+2.0d0*d_xy)	


  				sourcet(2,i,j)		= -cons(2,i,j)/Tr
				sourcet(3,i,j)		= -cons(3,i,j)/Tr


					residual(1,i,j)		= flux_x(1,i,j) + flux_y(1,i,j) + sourcet(1,i,j)
					residual(2,i,j)		= (flux_x(2,i,j) + flux_y(2,i,j) + sourcet(2,i,j))
					residual(3,i,j)		= (flux_x(3,i,j) + flux_y(3,i,j) + sourcet(3,i,j))
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


	open(unit=1,file='solution_1.plt')
	
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

 	

!******************************************************************************
!*                      Code generated with sympy 1.1.1                       *
!*                                                                            *
!*              See http://www.sympy.org/ for more information.               *
!*                                                                            *
!*                       This file is part of 'project'                       *
!******************************************************************************

double precision function my_function(D_para, x, y)
implicit none
double precision, intent(in) :: D_para
double precision, intent(in) :: x
double precision, intent(in) :: y

my_function = (1.0d0/2.0d0)*D_para*sin(3.14159265358979d0*x)*sin( &
      3.14159265358979d0*y) + (1.0d0/2.0d0)*sin(3.14159265358979d0*x)* &
      sin(3.14159265358979d0*y)

end function


