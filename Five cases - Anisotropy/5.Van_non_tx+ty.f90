
!**********************************************************************
!    This program solves 2D Laplace equation by Hyperbolic approach
!	
!	 Written by Sainath Ch, Email: s.chamarthi@al.t.u-tokyo.ac.jp


	program hyperpoi2d
	
		implicit none


		double precision, parameter 		:: pi=acos(-1.0d0)

		integer, parameter					:: NTMAX=100000,l=1
		integer 							:: i, j, z, k, NT
		
		integer, parameter					:: NX = 16*l,NY = 16*l, n_eqn = 3, ghostp = 20 	

		double precision, parameter			:: l1_target = 1.0d-10, epsilon = 1.0d0-14
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
	open(unit=4, file="residual.txt",action="write",status="replace")

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
	write(4,'(i8,4es12.4,4es12.4,4es12.4)') NT,(res_norm(1,L2)),(res_norm(2,L2)),(res_norm(3,L2))
			  
			  i_Linf_initial = i_Linf

		 endif

		 if (mod(NT,1000)==0) then
	write(*,*) NT
	write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (L1) = ", res_norm(:,L1  )
	write(*,'(a22,4es12.4,4es12.4,4es12.2)') "  Residual norm (L2) = ", res_norm(:,L2  )
	write(*,'(a22,4es12.4,4es12.4,4es12.4)') "  Residual norm (Li) = ", res_norm(:,Linf)
	write(4,'(i8,4es12.4,4es12.4,4es12.4)') NT,(res_norm(1,L2)),(res_norm(2,L2)),(res_norm(3,L2))
	write(*,*)
  	write(*,*) error_l2, MINVAL(cons(1,1:NX,1:NY))
  	write(*,*)
		  call output(NX,NY,x,y,cons,cons_exact,n_eqn,ghostp)
		 endif







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
    	
    	
    	write(*,*) 'Non-linear Anisotropic diffusion equation'
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

	double precision, parameter 		:: pi=acos(-1.0d0), gamma = 10.0d0

	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
	double precision 					:: cons_exact(n_eqn,1:NX,1:NY)
	double precision 					:: dx, dy, xmin, xmax, ymin, ymax, L


		L = 1.0d0

		cons_exact = 0.0d0


		do i = 1, NY
			
			do j = 1, NX

				cons_exact(1,i,j) = (x(i)*y(j)*(sin(pi*x(i))*sin(pi*y(j)))**gamma)
				cons_exact(2,i,j) = pi*gamma*x(i)*y(j)*(sin(3.14159265358979d0*x(i))*sin( &
      								3.14159265358979d0*y(j)))**gamma*cos(3.14159265358979d0*x(i))/sin( &
      								3.14159265358979d0*x(i)) + y(j)*(sin(3.14159265358979d0*x(i))*sin( &
      								3.14159265358979d0*y(j)))**gamma
      			cons_exact(3,i,j) = pi*gamma*x(i)*y(j)*(sin(3.14159265358979d0*x(i))*sin( &
      								3.14159265358979d0*y(j)))**gamma*cos(3.14159265358979d0*y(j))/sin( &
      								3.14159265358979d0*y(j)) + x(i)*(sin(3.14159265358979d0*x(i))*sin( &
      								3.14159265358979d0*y(j)))**gamma
			
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

	double precision, parameter 		:: pi=acos(-1.0d0), beta =  30.0d0, D_parallel = 1.0d9

	double precision					:: x_velocity, y_velocity
	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx, d_xy,d_yx,d_yy,p,q

		dt = 1.0d10;

		do i = 1, NX

			do j = 1, NY



					d_xx = D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+(1.0d0+cons(2,i,j)+cons(3,i,j))*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
					d_xy = 0.5d0*(D_parallel-1.0d0-cons(2,i,j)-cons(3,i,j))*sin(2.0d0*beta*pi/180.0d0)
					d_yx = 0.5d0*(D_parallel-1.0d0-cons(2,i,j)-cons(3,i,j))*sin(2.0d0*beta*pi/180.0d0)
					d_yy = D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)*(1.0d0+cons(2,i,j)+cons(3,i,j))

					Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
					Tr = (2.0d0*Lr**2.0d0)/(d_xx+d_yy+2.0d0*d_xy)

					! write(*,*) d_yy,d_xx, DABS(sqrt(d_xx/Tr)),DABS(sqrt(d_yy/Tr))
					! pause
				  x_velocity =  DABS(sqrt(d_xx/Tr))
				  y_velocity =  DABS(sqrt(d_yy/Tr))

				dtnew = min(dx/x_velocity, dy/y_velocity)
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
	double precision, parameter 		:: pi=acos(-1.0d0), beta =  30.0d0, D_parallel = 1.0d9, gamma = 10.0d0
	integer								:: NX, NY, n_eqn, ghostp

	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)	
	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	do j = 1,NY
		
		i=0     
	    
	    cons(1,i,j )   = (8.0d0*((x(i)+x(i+1))*0.5d0*y(j)*(sin(pi*(x(i)+x(i+1))*0.5d0)*sin(pi*y(j)))**gamma) - 6.0d0*cons(1,i+1,j) + cons(1,i+2,j)) * (1.0d0/3.0d0)

	    cons(2,i,j )   = cons(2,i+3,j)-3.0d0*cons(2,i+2,j)+3.0d0*cons(2,i+1,j)
        
        cons(3,i,j )   = cons(3,i+3,j)-3.0d0*cons(3,i+2,j)+3.0d0*cons(3,i+1,j)
    enddo

    do j = 1,NY

    	do i = -1,-ghostp,-1
        
        	cons(:,i,j)   = cons(:,i+3,j)-3.0d0*cons(:,i+2,j)+3.0d0*cons(:,i+1,j)

        enddo
    enddo


	! ! Right boundary conditions (actually right)


	do j = 1, NY
		
		i = NX

		cons(1,i+1,j) = (8.0d0*((x(i)+x(i+1))*0.5d0*y(j)*(sin(pi*(x(i)+x(i+1))*0.5d0)*sin(pi*y(j)))**gamma) +1.0d0*cons(1,i-1,j) -6.0d0* cons(1,i,j)) * (1.0d0/3.0d0)
        
        cons(2,i+1,j) = 3.0d0*cons(2,i,j)   - 3.0d0*cons(2,i-1,j)+cons(2,i-2,j)
        cons(3,i+1,j) = 3.0d0*cons(3,i,j)   - 3.0d0*cons(3,i-1,j)+cons(3,i-2,j)
        
    enddo

    do j = 1, NY

    	do i= NX, NX+8

        cons(:,i+2,j) = 3.0d0*cons(:,i+1,j) - 3.0d0*cons(:,i,j)  +cons(:,i-1,j)

        enddo
    enddo


	! ! Top boundary conditions (actually top)
	do i = 1, NX
		 j = NY

		cons(1,i,j+1) = (8.0d0*(x(i)*(y(j)+y(j+1))*0.5d0*(sin(pi*x(i))*sin(pi*(y(j)+y(j+1))*0.5d0))**gamma) +1.0d0*cons(1,i,j-1) -6.0d0* cons(1,i,j)) * (1.0d0/3.0d0)
        
        cons(2,i,j+1) = 3.0d0*cons(2,i,j)   - 3.0d0*cons(2,i,j-1)+cons(2,i,j-2)
        cons(3,i,j+1) = 3.0d0*cons(3,i,j)   - 3.0d0*cons(3,i,j-1)+cons(3,i,j-2)
        
    enddo

    do i = 1, NX

    	do j= NY, NY+8

        cons(:,i,j+2) = 3.0d0*cons(:,i,j+1) - 3.0d0*cons(:,i,j)  +cons(:,i,j-1)

        enddo
    enddo

   ! Bottom boundary conditions (actually bottom)

	do i = 1,NX
		j=0     
	    cons(1,i,j )   = (8.0d0*(x(i)*(y(j)+y(j+1))*0.5d0*(sin(pi*x(i))*sin(pi*(y(j)+y(j+1))*0.5d0))**gamma) - 6.0d0*cons(1,i,j+1) + cons(1,i,j+2)) * (1.0d0/3.0d0)

	    cons(2,i,j )   = cons(2,i,j+3)-3.0d0*cons(2,i,j+2)+3.0d0*cons(2,i,j+1)
        cons(3,i,j )   = cons(3,i,j+3)-3.0d0*cons(3,i,j+2)+3.0d0*cons(3,i,j+1)
    enddo

    do i = 1,NX  
    	do j = -1,-ghostp,-1
        
        	cons(:,i,j)   = cons(:,i,j+3)-3.0d0*cons(:,i,j+2)+3.0d0*cons(:,i,j+1)

        enddo
    enddo

!     ! Tensor

! 	do j = 1,NY
		
! 		i=0     
	    

!         cons(1,i,j)   = (128.0d0/35.0d0)*( ((x(i)+x(i+1))*0.5d0*y(j)*(sin(pi*(x(i)+x(i+1))*0.5d0)*sin(pi*y(j)))**gamma)- (35.0d0/32.0d0)*cons(1,i+1,j) + (35.0d0/64.0d0)*cons(1,i+2,j)-&
! 					  (7.0d0/32.0d0)*cons(1,i+3,j) +(5.0d0/128.0d0)*cons(1,i+4,j))
		
! 		cons(2,i,j)    = 5.0d0 * cons(2,i+1,j) - 10.0d0*cons(2,i+2,j) + 10.0d0*cons(2,i+3,j) -5.0d0*cons(2,i+4,j) + cons(2,i+5,j)
! 		cons(3,i,j)    = 5.0d0 * cons(3,i+1,j) - 10.0d0*cons(3,i+2,j) + 10.0d0*cons(3,i+3,j) -5.0d0*cons(3,i+4,j) + cons(3,i+5,j)
!     enddo

!     do j = 1,NY

!     	do i = -1,-ghostp,-1
        
!         	! cons(1,i,j)   = 5.0d0 * cons(1,i+1,j) - 10.0d0*cons(1,i+2,j) + 10.0d0*cons(1,i+3,j) -5.0d0*cons(1,i+4,j) + cons(1,i+5,j)
!         	cons(2,i,j)   = 5.0d0 * cons(2,i+1,j) - 10.0d0*cons(2,i+2,j) + 10.0d0*cons(2,i+3,j) -5.0d0*cons(2,i+4,j) + cons(2,i+5,j)
!         	cons(3,i,j)   = 5.0d0 * cons(3,i+1,j) - 10.0d0*cons(3,i+2,j) + 10.0d0*cons(3,i+3,j) -5.0d0*cons(3,i+4,j) + cons(3,i+5,j)

!         enddo
!     enddo

! ! ! ! ! Right boundary conditions


! 	do j = 1, NY
		
! 		i = NX

!         ! cons(1,i+1,j) = (128.0d0/35.0d0)*(((x(i)+x(i+1))*0.5d0*y(j)*(sin(pi*(x(i)+x(i+1))*0.5d0)*sin(pi*y(j)))**gamma)+(5.0d0/128.0d0)*cons(1,i-3,j) - (7.0d0/32.0d0)*cons(1,i-2,j) +&
! 					   ! (35.0d0/64.0d0)*cons(1,i-1,j) -(35.0d0/32.d0)*cons(1,i,j))

! 		cons(2,i+1,j)  = cons(2,i-4,j) - 5.0d0*cons(2,i-3,j) + 10.0d0*cons(2,i-2,j) -10.0d0*cons(2,i-1,j)+5.0d0*cons(2,i,j)
! 		cons(3,i+1,j)  = cons(3,i-4,j) - 5.0d0*cons(3,i-3,j) + 10.0d0*cons(3,i-2,j) -10.0d0*cons(3,i-1,j)+5.0d0*cons(3,i,j)
        
!     ! enddo

!     ! do j = 1, NY

!     	do i= NX, NX+8      

!         ! cons(1,i+2,j)  = cons(1,i-3,j) - 5.0d0*cons(1,i-2,j) + 10.0d0*cons(1,i-1,j) -10.0d0*cons(1,i,j)  +5.0d0*cons(1,i+1,j)
!         cons(2,i+2,j)  = cons(2,i-3,j) - 5.0d0*cons(2,i-2,j) + 10.0d0*cons(2,i-1,j) -10.0d0*cons(2,i,j)  +5.0d0*cons(2,i+1,j)
!         cons(3,i+2,j)  = cons(3,i-3,j) - 5.0d0*cons(3,i-2,j) + 10.0d0*cons(3,i-1,j) -10.0d0*cons(3,i,j)  +5.0d0*cons(3,i+1,j)

!         enddo
!     enddo




! ! 	! Bottom boundary conditions

! 	do i = 1,NX
! 		j=0     

!         cons(1,i,j)   = (128.0d0/35.0d0)*((x(i)*(y(j)+y(j+1))*0.5d0*(sin(pi*x(i))*sin(pi*(y(j)+y(j+1))*0.5d0))**gamma) - (35.0d0/32.0d0)*cons(1,i,j+1) + (35.0d0/64.0d0)*cons(1,i,j+2)-&
! 					    (7.0d0/32.0d0)*cons(1,i,j+3) +(5.0d0/128.0d0)*cons(1,i,j+4))
		
! 		cons(2,i,j)    = 5.0d0 * cons(2,i,j+1) - 10.0d0*cons(2,i,j+2) + 10.0d0*cons(2,i,j+3) -5.0d0*cons(2,i,j+4) + cons(2,i,j+5)
! 		cons(3,i,j)    = 5.0d0 * cons(3,i,j+1) - 10.0d0*cons(3,i,j+2) + 10.0d0*cons(3,i,j+3) -5.0d0*cons(3,i,j+4) + cons(3,i,j+5)
!     enddo

!     do i = 1,NX  
!     	do j = -1,-ghostp,-1
        
        	
!         	cons(1,i,j)   = 5.0d0 * cons(1,i,j+1) - 10.0d0*cons(1,i,j+2) + 10.0d0*cons(1,i,j+3) -5.0d0*cons(1,i,j+4) + cons(1,i,j+5)
!         	cons(2,i,j)   = 5.0d0 * cons(2,i,j+1) - 10.0d0*cons(2,i,j+2) + 10.0d0*cons(2,i,j+3) -5.0d0*cons(2,i,j+4) + cons(2,i,j+5)
!         	cons(3,i,j)   = 5.0d0 * cons(3,i,j+1) - 10.0d0*cons(3,i,j+2) + 10.0d0*cons(3,i,j+3) -5.0d0*cons(3,i,j+4) + cons(3,i,j+5)

!         enddo
!     enddo


! ! ! ! 	! Top boundary conditions
! 	do i = 1, NX
! 		 j = NY


!         cons(1,i,j+1) = (128.0d0/35.0d0)*((x(i)*(y(j)+y(j+1))*0.5d0*(sin(pi*x(i))*sin(pi*(y(j)+y(j+1))*0.5d0))**gamma) +(5.0d0/128.0d0)*cons(1,i,j-3) - (7.0d0/32.0d0)*cons(1,i,j-2) +&
! 					   (35.0d0/64.0d0)*cons(1,i,j-1) -(35.0d0/32.d0)*cons(1,i,j))

! 		cons(2,i,j+1)  = cons(2,i,j-4) - 5.0d0*cons(2,i,j-3) + 10.0d0*cons(2,i,j-2) -10.0d0*cons(2,i,j-1)+5.0d0*cons(2,i,j)
! 		cons(3,i,j+1)  = cons(3,i,j-4) - 5.0d0*cons(3,i,j-3) + 10.0d0*cons(3,i,j-2) -10.0d0*cons(3,i,j-1)+5.0d0*cons(3,i,j)
        
!     enddo

!     do i = 1, NX

!     	do j= NY, NY+8

        
!         cons(1,i,j+2)  = cons(1,i,j-3) - 5.0d0*cons(1,i,j-2) + 10.0d0*cons(1,i,j-1) -10.0d0*cons(1,i,j)  +5.0d0*cons(1,i,j+1)
!         cons(2,i,j+2)  = cons(2,i,j-3) - 5.0d0*cons(2,i,j-2) + 10.0d0*cons(2,i,j-1) -10.0d0*cons(2,i,j)  +5.0d0*cons(2,i,j+1)
!         cons(3,i,j+2)  = cons(3,i,j-3) - 5.0d0*cons(3,i,j-2) + 10.0d0*cons(3,i,j-1) -10.0d0*cons(3,i,j)  +5.0d0*cons(3,i,j+1)

!         enddo
!     enddo


	

	end subroutine boundarycondition



	!***********************************************************************
	!*****                      Flux in X-direction                    *****
	!***********************************************************************

	subroutine FX(NX,NY,dx,x,y,cons,flux_x,n_eqn,ghostp)

	implicit none

	double precision, parameter 		:: pi=acos(-1.0d0), beta =  30.0d0, D_parallel = 1.0d9
	integer								:: i, j ,k
	integer 							:: NX, NY, ghostp, n_eqn
	double precision					:: dx, x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision 					:: flux_x(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision					:: flux_node_x(n_eqn,-ghostp:NX+ghostp), flux_x_half(n_eqn,-ghostp:NX+ghostp)

	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),cons_local(n_eqn,-ghostp:NX+ghostp)

	double precision 					:: nu,Lr,Tr,p,q
	double precision 					:: d_xx(-ghostp:NX+ghostp), d_xy(-ghostp:NX+ghostp),d_yx(-ghostp:NX+ghostp),d_yy(-ghostp:NX+ghostp)




	do j =1, NY

		do i = -ghostp, NX+ghostp

			cons_local(1,i) = cons(1,i,j)
			cons_local(2,i) = cons(2,i,j)
			cons_local(3,i) = cons(3,i,j)

					d_xx(i) = D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+(1.0d0+cons(2,i,j)+cons(3,i,j))*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
					d_xy(i) = 0.5d0*(D_parallel-1.0d0-cons(2,i,j)-cons(3,i,j))*sin(2.0d0*beta*pi/180.0d0)
					d_yx(i) = 0.5d0*(D_parallel-1.0d0-cons(2,i,j)-cons(3,i,j))*sin(2.0d0*beta*pi/180.0d0)
					d_yy(i) = D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)*(1.0d0+cons(2,i,j)+cons(3,i,j))

			flux_node_x(1,i) = -cons_local(2,i)*d_xx(i)-cons_local(3,i)*d_xy(i)
			flux_node_x(2,i) = -cons_local(1,i)
			flux_node_x(3,i) = 0.0d0
		end do

		call reconstruction5EX (NX,NY,flux_node_x,cons_local,flux_x_half, n_eqn,ghostp,x,y,j,d_xx,d_xy,d_yx,d_yy)

		do k = 1,n_eqn

			do i = 1, NX


				!Finite volume or Reconstruction approach 
	        	flux_x(k,i,j) = -((flux_x_half(k,i)-flux_x_half(k,i-1))/dx)


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

	double precision					:: a(0:NX+1), b(0:NX+1), c(0:NX+1), r(0:NX+1),temp(0:NX+1)

	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx(-ghostp:NX+ghostp), d_xy(-ghostp:NX+ghostp)
	double precision 					:: d_yy(-ghostp:NY+ghostp), d_yx(-ghostp:NY+ghostp)
	
	double precision 					:: d_xxl(-ghostp:NX+ghostp), d_xxr(-ghostp:NX+ghostp)
	double precision 					:: d_xyl(-ghostp:NX+ghostp), d_xyr(-ghostp:NX+ghostp)

	double precision 					:: d_yxl(-ghostp:NY+ghostp),d_yxr(-ghostp:NY+ghostp)
	double precision 					:: d_yyl(-ghostp:NY+ghostp),d_yyr(-ghostp:NY+ghostp)


		do k=1,n_eqn

        	do i =-2,NX+2

              consl(k,i) = (3.0d0/128.0d0)*(cons_local(k,i-2)) - (5.0d0/32.0d0)*(cons_local(k,i-1)) + (45.0d0/64.0d0)*(cons_local(k,i)) + &
					     (15.0d0/32.0d0)*(cons_local(k,i+1)) - (5.0d0/128.0d0)*(cons_local(k,i+2))
			
			consr(k,i) = (-5.0d0/128.0d0)*(cons_local(k,i-1)) + (15.0d0/32.0d0)*(cons_local(k,i)) + (45.0d0/64.0d0)*(cons_local(k,i+1)) - &
					     (5.0d0/32.0d0)*(cons_local(k,i+2)) + (3.0d0/128.0d0)*(cons_local(k,i+3))


  			fleft(k,i) = (3.0d0/128.0d0)*(flux_node(k,i-2)) - (5.0d0/32.0d0)*(flux_node(k,i-1)) + (45.0d0/64.0d0)*(flux_node(k,i)) + &
					     (15.0d0/32.0d0)*(flux_node(k,i+1)) - (5.0d0/128.0d0)*(flux_node(k,i+2))
			
			fright(k,i) = (-5.0d0/128.0d0)*(flux_node(k,i-1)) + (15.0d0/32.0d0)*(flux_node(k,i)) + (45.0d0/64.0d0)*(flux_node(k,i+1)) - &
					     (5.0d0/32.0d0)*(flux_node(k,i+2)) + (3.0d0/128.0d0)*(flux_node(k,i+3))

		  	
            
       
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
    
		do i=-2,NX+2

	

			d_xxl(i) = (3.0d0/128.0d0)*(d_xx(i-2)) - (5.0d0/32.0d0)*(d_xx(i-1)) + (45.0d0/64.0d0)*(d_xx(i)) + &
					     (15.0d0/32.0d0)*(d_xx(i+1)) - (5.0d0/128.0d0)*(d_xx(i+2))
			
			d_xxr(i) = (-5.0d0/128.0d0)*(d_xx(i-1)) + (15.0d0/32.0d0)*(d_xx(i)) + (45.0d0/64.0d0)*(d_xx(i+1)) - &
					     (5.0d0/32.0d0)*(d_xx(i+2)) + (3.0d0/128.0d0)*(d_xx(i+3))

		  	d_xyl(i) = (3.0d0/128.0d0)*(d_xy(i-2)) - (5.0d0/32.0d0)*(d_xy(i-1)) + (45.0d0/64.0d0)*(d_xy(i)) + &
					     (15.0d0/32.0d0)*(d_xy(i+1)) - (5.0d0/128.0d0)*(d_xy(i+2))
			
			d_xyr(i) = (-5.0d0/128.0d0)*(d_xy(i-1)) + (15.0d0/32.0d0)*(d_xy(i)) + (45.0d0/64.0d0)*(d_xy(i+1)) - &
					     (5.0d0/32.0d0)*(d_xy(i+2)) + (3.0d0/128.0d0)*(d_xy(i+3))

		  	d_yyl(i) = (3.0d0/128.0d0)*(d_yy(i-2)) - (5.0d0/32.0d0)*(d_yy(i-1)) + (45.0d0/64.0d0)*(d_yy(i)) + &
					     (15.0d0/32.0d0)*(d_yy(i+1)) - (5.0d0/128.0d0)*(d_yy(i+2))
			
			d_yyr(i) = (-5.0d0/128.0d0)*(d_yy(i-1)) + (15.0d0/32.0d0)*(d_yy(i)) + (45.0d0/64.0d0)*(d_yy(i+1)) - &
					     (5.0d0/32.0d0)*(d_yy(i+2)) + (3.0d0/128.0d0)*(d_yy(i+3))

		  	d_yxl(i) = (3.0d0/128.0d0)*(d_yx(i-2)) - (5.0d0/32.0d0)*(d_yx(i-1)) + (45.0d0/64.0d0)*(d_yx(i)) + &
					     (15.0d0/32.0d0)*(d_yx(i+1)) - (5.0d0/128.0d0)*(d_yx(i+2))
			
			d_yxr(i) = (-5.0d0/128.0d0)*(d_yx(i-1)) + (15.0d0/32.0d0)*(d_yx(i)) + (45.0d0/64.0d0)*(d_yx(i+1)) - &
					     (5.0d0/32.0d0)*(d_yx(i+2)) + (3.0d0/128.0d0)*(d_yx(i+3))


			d_xxl(i) = (-1.0d0/8.0d0)*(d_xx(i-1)) + (6.0d0/8.0d00)*(d_xx(i)) + (3.0d0/8.0d0)*(d_xx(i+1))
			
			d_xxr(i) = (3.0d0/8.0d0)*(d_xx(i)) + (6.0d0/8.0d00)*(d_xx(i+1)) - (1.0d0/8.0d0)*(d_xx(i+2))

		  	d_xyl(i) = (-1.0d0/8.0d0)*(d_xy(i-1)) + (6.0d0/8.0d00)*(d_xy(i)) + (3.0d0/8.0d0)*(d_xy(i+1))
			
			d_xyr(i) = (3.0d0/8.0d0)*(d_xy(i)) + (6.0d0/8.0d00)*(d_xy(i+1)) - (1.0d0/8.0d0)*(d_xy(i+2))

		  	d_yyl(i) = (-1.0d0/8.0d0)*(d_yy(i-1)) + (6.0d0/8.0d00)*(d_yy(i)) + (3.0d0/8.0d0)*(d_yy(i+1))
			
			d_yyr(i) = (3.0d0/8.0d0)*(d_yy(i)) + (6.0d0/8.0d00)*(d_yy(i+1)) - (1.0d0/8.0d0)*(d_yy(i+2))

		  	d_yxl(i) = (-1.0d0/8.0d0)*(d_yx(i-1)) + (6.0d0/8.0d00)*(d_yx(i)) + (3.0d0/8.0d0)*(d_yx(i+1))
			
			d_yxr(i) = (3.0d0/8.0d0)*(d_yx(i)) + (6.0d0/8.0d00)*(d_yx(i+1)) - (1.0d0/8.0d0)*(d_yx(i+2))


			! fleft(1,i) = -consl(2,i)*d_xxl(i)-consl(3,i)*d_xyl(i)
			! fleft(2,i) = -consl(1,i)
			! fleft(3,i) = 0.0d0

			! fright(1,i) = -consr(2,i)*d_xxr(i)-consr(3,i)*d_xyr(i)
			! fright(2,i) = -consr(1,i)
			! fright(3,i) = 0.0d0

	Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
	Tr = (2.0d0*Lr**2.0d0)/( ((d_xxl(i)+d_xxr(i))*0.5d0) + ((d_yyl(i)+d_yyr(i))*0.5d0) + 2.0d0*(((d_xyl(i)+d_xyr(i))*0.5d0)))


			delphi(i)			= consr(1,i)-consl(1,i)
			delp(i)				= consr(2,i)-consl(2,i)
			delq(i)				= consr(3,i)-consl(3,i)

	!***************************************
	!***** 		Roe solver	       	   *****
	!***************************************

		 !***** Roe solver          *****
		 !The upwind dissipation matrix is |A| = R*|Lambda|*L
		 jacobian(1,1) = abs((((d_xxl(i)+d_xxr(i))*0.5d0)/Lr))
		 jacobian(1,2) = 0.0d0
		 jacobian(1,3) = 0.0d0
		 jacobian(2,1) = 0.0d0
		 jacobian(2,2) = Tr*abs((((d_xxl(i)+d_xxr(i))*0.5d0)/Lr))
		 jacobian(2,3) = Tr*abs((((d_xxl(i)+d_xxr(i))*0.5d0)/Lr))*(d_xyl(i)+d_xyr(i))*0.5d0/((d_xxl(i)+d_xxr(i))*0.5d0)
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
	double precision, parameter 		:: pi=acos(-1.0d0), beta =  30.0d0, D_parallel = 1.0d9
	integer								:: i, j, k

	integer 							:: NX, NY, ghostp, n_eqn

	double precision					:: dy, x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision 					:: flux_y(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	
	double precision					:: flux_node_y(n_eqn,-ghostp:NY+ghostp), flux_y_half(n_eqn,-ghostp:NY+ghostp)

	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),cons_local(n_eqn,-ghostp:NY+ghostp)

	double precision 					:: nu,Lr,Tr,p,q
	 double precision 					:: d_xx(-ghostp:NX+ghostp), d_xy(-ghostp:NX+ghostp),d_yx(-ghostp:NX+ghostp),d_yy(-ghostp:NX+ghostp)




	
	do j =1, NX
		
		do i = -ghostp, NY+ghostp

			cons_local(1,i) = cons(1,j,i)
			cons_local(2,i) = cons(2,j,i)
			cons_local(3,i) = cons(3,j,i)

					d_xx(i) = D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+(1.0d0+cons(2,j,i)+cons(3,j,i))*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
					d_xy(i) = 0.5d0*(D_parallel-1.0d0-cons(2,j,i)-cons(3,j,i))*sin(2.0d0*beta*pi/180.0d0)
					d_yx(i) = 0.5d0*(D_parallel-1.0d0-cons(2,j,i)-cons(3,j,i))*sin(2.0d0*beta*pi/180.0d0)
					d_yy(i) = D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)*(1.0d0+cons(2,j,i)+cons(3,j,i))


			flux_node_y(1,i) =  -cons_local(3,i)*d_yy(i)-cons_local(2,i)*d_yx(i)
			flux_node_y(2,i) =  0.0d0
			flux_node_y(3,i) =  -cons_local(1,i)
		end do

		call reconstruction5EY (NX,NY,flux_node_y,cons_local,flux_y_half, n_eqn,ghostp,x,y,j,d_xx,d_xy,d_yx,d_yy)

		do k = 1,n_eqn

			do i = 1, NY


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

	double precision					:: a(0:NY+1), b(0:NY+1), c(0:NY+1), r(0:NY+1),temp(0:NY+1)
	double precision 					:: nu,Lr,Tr
	double precision 					:: d_xx(-ghostp:NX+ghostp), d_xy(-ghostp:NX+ghostp)
	double precision 					:: d_yy(-ghostp:NY+ghostp), d_yx(-ghostp:NY+ghostp)
	
	double precision 					:: d_xxl(-ghostp:NX+ghostp), d_xxr(-ghostp:NX+ghostp)
	double precision 					:: d_xyl(-ghostp:NX+ghostp), d_xyr(-ghostp:NX+ghostp)

	double precision 					:: d_yxl(-ghostp:NY+ghostp),d_yxr(-ghostp:NY+ghostp)
	double precision 					:: d_yyl(-ghostp:NY+ghostp),d_yyr(-ghostp:NY+ghostp)


		do k=1,n_eqn

        	do i =-2,NY+2
              consl(k,i) = (3.0d0/128.0d0)*(cons_local(k,i-2)) - (5.0d0/32.0d0)*(cons_local(k,i-1)) + (45.0d0/64.0d0)*(cons_local(k,i)) + &
					     (15.0d0/32.0d0)*(cons_local(k,i+1)) - (5.0d0/128.0d0)*(cons_local(k,i+2))
			
			consr(k,i) = (-5.0d0/128.0d0)*(cons_local(k,i-1)) + (15.0d0/32.0d0)*(cons_local(k,i)) + (45.0d0/64.0d0)*(cons_local(k,i+1)) - &
					     (5.0d0/32.0d0)*(cons_local(k,i+2)) + (3.0d0/128.0d0)*(cons_local(k,i+3))


  			fleft(k,i) = (3.0d0/128.0d0)*(flux_node(k,i-2)) - (5.0d0/32.0d0)*(flux_node(k,i-1)) + (45.0d0/64.0d0)*(flux_node(k,i)) + &
					     (15.0d0/32.0d0)*(flux_node(k,i+1)) - (5.0d0/128.0d0)*(flux_node(k,i+2))
			
			fright(k,i) = (-5.0d0/128.0d0)*(flux_node(k,i-1)) + (15.0d0/32.0d0)*(flux_node(k,i)) + (45.0d0/64.0d0)*(flux_node(k,i+1)) - &
					     (5.0d0/32.0d0)*(flux_node(k,i+2)) + (3.0d0/128.0d0)*(flux_node(k,i+3))

		  	

  
         	enddo
        enddo

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

    
		do i=-2,NY+2

			d_xxl(i) = (3.0d0/128.0d0)*(d_xx(i-2)) - (5.0d0/32.0d0)*(d_xx(i-1)) + (45.0d0/64.0d0)*(d_xx(i)) + &
					     (15.0d0/32.0d0)*(d_xx(i+1)) - (5.0d0/128.0d0)*(d_xx(i+2))
			
			d_xxr(i) = (-5.0d0/128.0d0)*(d_xx(i-1)) + (15.0d0/32.0d0)*(d_xx(i)) + (45.0d0/64.0d0)*(d_xx(i+1)) - &
					     (5.0d0/32.0d0)*(d_xx(i+2)) + (3.0d0/128.0d0)*(d_xx(i+3))

		  	d_xyl(i) = (3.0d0/128.0d0)*(d_xy(i-2)) - (5.0d0/32.0d0)*(d_xy(i-1)) + (45.0d0/64.0d0)*(d_xy(i)) + &
					     (15.0d0/32.0d0)*(d_xy(i+1)) - (5.0d0/128.0d0)*(d_xy(i+2))
			
			d_xyr(i) = (-5.0d0/128.0d0)*(d_xy(i-1)) + (15.0d0/32.0d0)*(d_xy(i)) + (45.0d0/64.0d0)*(d_xy(i+1)) - &
					     (5.0d0/32.0d0)*(d_xy(i+2)) + (3.0d0/128.0d0)*(d_xy(i+3))

		  	d_yyl(i) = (3.0d0/128.0d0)*(d_yy(i-2)) - (5.0d0/32.0d0)*(d_yy(i-1)) + (45.0d0/64.0d0)*(d_yy(i)) + &
					     (15.0d0/32.0d0)*(d_yy(i+1)) - (5.0d0/128.0d0)*(d_yy(i+2))
			
			d_yyr(i) = (-5.0d0/128.0d0)*(d_yy(i-1)) + (15.0d0/32.0d0)*(d_yy(i)) + (45.0d0/64.0d0)*(d_yy(i+1)) - &
					     (5.0d0/32.0d0)*(d_yy(i+2)) + (3.0d0/128.0d0)*(d_yy(i+3))

		  	d_yxl(i) = (3.0d0/128.0d0)*(d_yx(i-2)) - (5.0d0/32.0d0)*(d_yx(i-1)) + (45.0d0/64.0d0)*(d_yx(i)) + &
					     (15.0d0/32.0d0)*(d_yx(i+1)) - (5.0d0/128.0d0)*(d_yx(i+2))
			
			d_yxr(i) = (-5.0d0/128.0d0)*(d_yx(i-1)) + (15.0d0/32.0d0)*(d_yx(i)) + (45.0d0/64.0d0)*(d_yx(i+1)) - &
					     (5.0d0/32.0d0)*(d_yx(i+2)) + (3.0d0/128.0d0)*(d_yx(i+3))

			d_xxl(i) = (-1.0d0/8.0d0)*(d_xx(i-1)) + (6.0d0/8.0d00)*(d_xx(i)) + (3.0d0/8.0d0)*(d_xx(i+1))
			
			d_xxr(i) = (3.0d0/8.0d0)*(d_xx(i)) + (6.0d0/8.0d00)*(d_xx(i+1)) - (1.0d0/8.0d0)*(d_xx(i+2))

		  	d_xyl(i) = (-1.0d0/8.0d0)*(d_xy(i-1)) + (6.0d0/8.0d00)*(d_xy(i)) + (3.0d0/8.0d0)*(d_xy(i+1))
			
			d_xyr(i) = (3.0d0/8.0d0)*(d_xy(i)) + (6.0d0/8.0d00)*(d_xy(i+1)) - (1.0d0/8.0d0)*(d_xy(i+2))

		  	d_yyl(i) = (-1.0d0/8.0d0)*(d_yy(i-1)) + (6.0d0/8.0d00)*(d_yy(i)) + (3.0d0/8.0d0)*(d_yy(i+1))
			
			d_yyr(i) = (3.0d0/8.0d0)*(d_yy(i)) + (6.0d0/8.0d00)*(d_yy(i+1)) - (1.0d0/8.0d0)*(d_yy(i+2))

		  	d_yxl(i) = (-1.0d0/8.0d0)*(d_yx(i-1)) + (6.0d0/8.0d00)*(d_yx(i)) + (3.0d0/8.0d0)*(d_yx(i+1))
			
			d_yxr(i) = (3.0d0/8.0d0)*(d_yx(i)) + (6.0d0/8.0d00)*(d_yx(i+1)) - (1.0d0/8.0d0)*(d_yx(i+2))

			! fleft(1,i) =  -consl(3,i)*d_yyl(i)-consl(2,i)*d_yxl(i)
			! fleft(2,i) =  0.0d0
			! fleft(3,i) =  -consl(1,i)

			! fright(1,i) =  -consr(3,i)*d_yyr(i)-consr(2,i)*d_yxr(i)
			! fright(2,i) =  0.0d0
			! fright(3,i) =  -consr(1,i)



				Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
				Tr = (2.0d0*Lr**2.0d0)/( ((d_xxl(i)+d_xxr(i))*0.5d0) + ((d_yyl(i)+d_yyr(i))*0.5d0) + 2.0d0*(((d_xyl(i)+d_xyr(i))*0.5d0)))


			delphi(i)			= consr(1,i)-consl(1,i)
			delp(i)				= consr(2,i)-consl(2,i)
			delq(i)				= consr(3,i)-consl(3,i)

	!***************************************
	!***** 		Roe solver	       	   *****
	!***************************************

		 !***** Roe solver          *****
		 !The upwind dissipation matrix is |A| = R*|Lambda|*L
		 jacobian(1,1) = abs((((d_yyl(i)+d_yyr(i))*0.5d0)/Lr))
		 jacobian(1,2) = 0.0d0
		 jacobian(1,3) = 0.0d0
		 jacobian(2,1) = 0.0d0
		 jacobian(2,2) = 0.0d0
		 jacobian(2,3) = 0.0d0
		 jacobian(3,1) = 0.0d0
		 jacobian(3,2) = Tr*abs((((d_yyl(i)+d_yyr(i))*0.5d0)/Lr))*(d_xyl(i)+d_xyr(i))*0.5d0/((d_yyl(i)+d_yyr(i))*0.5d0)
		 jacobian(3,3) = Tr*abs((((d_yyl(i)+d_yyr(i))*0.5d0)/Lr))


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


 		double precision, parameter 		:: pi=acos(-1.0d0), beta =  30.0d0, D_parallel = 1.0d9, gamma = 10.0d0
 		integer								:: NX, NY, n_eqn, ghostp
		double precision 					:: X(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp), L
	
		
		integer								:: i,j
		
		
		double precision					:: sourcet(3,1:NX,1:NY)	
		double precision 					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision 					:: residual(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision 					:: flux_x(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),flux_y(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision             	    :: d_xx, d_xy,d_yx,d_yy
		double precision                    :: nu,Lr,Tr
		
		double precision					:: A, B, C, D, E, F,frac, my_function,p,q

	
			

			do j = 1, NY

				do i = 1, NX


				sourcet(1,i,j)		= my_function(D_parallel, beta*pi/180.0d0, gamma, x(i), y(j))

				p = cons(2,i,j)
				q = cons(3,i,j)


					d_xx = D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+(1.0d0+cons(2,i,j)+cons(3,i,j))*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
					d_xy = 0.5d0*(D_parallel-1.0d0-cons(2,i,j)-cons(3,i,j))*sin(2.0d0*beta*pi/180.0d0)
					d_yx = 0.5d0*(D_parallel-1.0d0-cons(2,i,j)-cons(3,i,j))*sin(2.0d0*beta*pi/180.0d0)
					d_yy = D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)*(1.0d0+cons(2,i,j)+cons(3,i,j))
 


			  Lr = 2.0d0/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))

	  		  Tr = (2.0d0*Lr**2.0d0)/(d_xx+d_yy+2.0d0*d_xy)


				sourcet(2,i,j)		= -cons(2,i,j)
				sourcet(3,i,j)		= -cons(3,i,j)

					residual(1,i,j)		= flux_x(1,i,j) + flux_y(1,i,j) + sourcet(1,i,j)
					residual(2,i,j)		= (flux_x(2,i,j) + flux_y(2,i,j) + sourcet(2,i,j))/Tr
					residual(3,i,j)		= (flux_x(3,i,j) + flux_y(3,i,j) + sourcet(3,i,j))/Tr
				enddo

		enddo
	end subroutine sourceterm

	double precision function my_function(D_para, beta, gamma, x, y)
implicit none
double precision, intent(in) :: D_para
double precision, intent(in) :: beta
double precision, intent(in) :: gamma
double precision, intent(in) :: x
double precision, intent(in) :: y

double precision, parameter :: pi = 3.14159265358979d0
my_function = (1.0d0/4.0d0)*(sin(3.14159265358979d0*x)*sin( &
      3.14159265358979d0*y))**gamma*(-2*pi*gamma*x*(D_para*sin(beta)**2 &
      *sin(3.14159265358979d0*x)*sin(3.14159265358979d0*y) + (pi*gamma* &
      x*y*(sin(3.14159265358979d0*x)*sin(3.14159265358979d0*y))**gamma* &
      sin(3.14159265358979d0*x + 3.14159265358979d0*y) + x*(sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y))**gamma*sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y) + y*(sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y))**gamma*sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y) + sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y))*cos(beta)**2)*( &
      pi*gamma*y*cos(6.28318530717959d0*y) + pi*gamma*y - 2*pi*y + 2* &
      sin(6.28318530717959d0*y))*sin(3.14159265358979d0*x)**2 - 2*pi* &
      gamma*y*(D_para*sin(3.14159265358979d0*x)*sin(3.14159265358979d0* &
      y)*cos(beta)**2 + (pi*gamma*x*y*(sin(3.14159265358979d0*x)*sin( &
      3.14159265358979d0*y))**gamma*sin(3.14159265358979d0*x + &
      3.14159265358979d0*y) + x*(sin(3.14159265358979d0*x)*sin( &
      3.14159265358979d0*y))**gamma*sin(3.14159265358979d0*x)*sin( &
      3.14159265358979d0*y) + y*(sin(3.14159265358979d0*x)*sin( &
      3.14159265358979d0*y))**gamma*sin(3.14159265358979d0*x)*sin( &
      3.14159265358979d0*y) + sin(3.14159265358979d0*x)*sin( &
      3.14159265358979d0*y))*sin(beta)**2)*(pi*gamma*x*cos( &
      6.28318530717959d0*x) + pi*gamma*x - 2*pi*x + 2*sin( &
      6.28318530717959d0*x))*sin(3.14159265358979d0*y)**2 - x*(sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y))**gamma*(pi*gamma &
      *y*cos(3.14159265358979d0*y) + sin(3.14159265358979d0*y))*(-pi**2 &
      *gamma**2*x*y*(sin(3.14159265358979d0*x - 6.28318530717959d0*y) - &
      sin(3.14159265358979d0*x + 6.28318530717959d0*y)) + 4*pi**2*gamma &
      *x*y*(gamma - 1)*sin(3.14159265358979d0*x)*cos(3.14159265358979d0 &
      *y)**2 + 4*pi*gamma*x*sin(3.14159265358979d0*y)**2*cos( &
      3.14159265358979d0*x) + pi*gamma*(2*x + y)*(cos( &
      3.14159265358979d0*x - 6.28318530717959d0*y) - cos( &
      3.14159265358979d0*x + 6.28318530717959d0*y)) + (-4*pi**2*gamma*x &
      *y + 4)*sin(3.14159265358979d0*x)*sin(3.14159265358979d0*y)**2)* &
      sin(3.14159265358979d0*x)**2*cos(beta)**2 + x*(sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y))**gamma*(pi*gamma &
      *y*cos(3.14159265358979d0*y) + sin(3.14159265358979d0*y))*(pi**2* &
      gamma**2*x*y*(sin(6.28318530717959d0*x - 3.14159265358979d0*y) + &
      sin(6.28318530717959d0*x + 3.14159265358979d0*y)) + 4*pi**2*gamma &
      *x*y*(gamma - 1)*sin(3.14159265358979d0*y)*cos(3.14159265358979d0 &
      *x)**2 + 4*pi*gamma*y*sin(3.14159265358979d0*x)**2*cos( &
      3.14159265358979d0*y) + pi*gamma*(x + 2*y)*(cos( &
      6.28318530717959d0*x - 3.14159265358979d0*y) - cos( &
      6.28318530717959d0*x + 3.14159265358979d0*y)) + (-4*pi**2*gamma*x &
      *y + 4)*sin(3.14159265358979d0*x)**2*sin(3.14159265358979d0*y))* &
      sin(beta)*sin(3.14159265358979d0*x)*sin(3.14159265358979d0*y)*cos &
      (beta) + y*(sin(3.14159265358979d0*x)*sin(3.14159265358979d0*y)) &
      **gamma*(pi*gamma*x*cos(3.14159265358979d0*x) + sin( &
      3.14159265358979d0*x))*(-pi**2*gamma**2*x*y*(sin( &
      3.14159265358979d0*x - 6.28318530717959d0*y) - sin( &
      3.14159265358979d0*x + 6.28318530717959d0*y)) + 4*pi**2*gamma*x*y &
      *(gamma - 1)*sin(3.14159265358979d0*x)*cos(3.14159265358979d0*y) &
      **2 + 4*pi*gamma*x*sin(3.14159265358979d0*y)**2*cos( &
      3.14159265358979d0*x) + pi*gamma*(2*x + y)*(cos( &
      3.14159265358979d0*x - 6.28318530717959d0*y) - cos( &
      3.14159265358979d0*x + 6.28318530717959d0*y)) + (-4*pi**2*gamma*x &
      *y + 4)*sin(3.14159265358979d0*x)*sin(3.14159265358979d0*y)**2)* &
      sin(beta)*sin(3.14159265358979d0*x)*sin(3.14159265358979d0*y)*cos &
      (beta) - y*(sin(3.14159265358979d0*x)*sin(3.14159265358979d0*y)) &
      **gamma*(pi*gamma*x*cos(3.14159265358979d0*x) + sin( &
      3.14159265358979d0*x))*(pi**2*gamma**2*x*y*(sin( &
      6.28318530717959d0*x - 3.14159265358979d0*y) + sin( &
      6.28318530717959d0*x + 3.14159265358979d0*y)) + 4*pi**2*gamma*x*y &
      *(gamma - 1)*sin(3.14159265358979d0*y)*cos(3.14159265358979d0*x) &
      **2 + 4*pi*gamma*y*sin(3.14159265358979d0*x)**2*cos( &
      3.14159265358979d0*y) + pi*gamma*(x + 2*y)*(cos( &
      6.28318530717959d0*x - 3.14159265358979d0*y) - cos( &
      6.28318530717959d0*x + 3.14159265358979d0*y)) + (-4*pi**2*gamma*x &
      *y + 4)*sin(3.14159265358979d0*x)**2*sin(3.14159265358979d0*y))* &
      sin(beta)**2*sin(3.14159265358979d0*y)**2 + 8*(pi**2*gamma**2*x*y &
      *cos(3.14159265358979d0*x)*cos(3.14159265358979d0*y) + pi*gamma*x &
      *sin(3.14159265358979d0*y)*cos(3.14159265358979d0*x) + pi*gamma*y &
      *sin(3.14159265358979d0*x)*cos(3.14159265358979d0*y) + sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y))*(-D_para*sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y) + pi*gamma*x*y*( &
      sin(3.14159265358979d0*x)*sin(3.14159265358979d0*y))**gamma*sin( &
      3.14159265358979d0*x + 3.14159265358979d0*y) + x*(sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y))**gamma*sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y) + y*(sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y))**gamma*sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y) + sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y))*sin(beta)*sin( &
      3.14159265358979d0*x)*sin(3.14159265358979d0*y)*cos(beta))/(sin( &
      3.14159265358979d0*x)**3*sin(3.14159265358979d0*y)**3)


end function


	!***********************************************************************
	!*****                       Output 			                   *****
	!***********************************************************************

	subroutine output(NX,NY,x,y,cons,cons_exact,n_eqn,ghostp)

	integer 				:: i, j

	integer					:: NX,NY, n_eqn, ghostp

	double precision 		:: X(-ghostp:NX+ghostp),Y(-ghostp:NY+ghostp),cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: cons_exact(n_eqn,1:NX,1:NY)


	open(unit=1,file='solution_5.plt',form='formatted',status='replace')
	
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




