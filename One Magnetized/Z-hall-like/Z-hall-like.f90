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
		
		integer, parameter					:: NX = 48*l,NY = 48*l, n_eqn = 3, ghostp = 10 	

		double precision, parameter			:: l1_target = 1.0d-10, epsilon = 1.0d0-14
		double precision 					:: CFL, dt

		double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
		double precision 					:: dx, dy, xmin, xmax, ymin, ymax

		double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),cons_old(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)


		double precision 					:: residual(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision 					:: flux_x(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp),flux_y(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		double precision					:: res_norm_initial(n_eqn,3), res_norm(n_eqn,3)

		integer 							:: L1 = 1, L2 = 2, Linf = 3, i_Linf, i_Linf_initial

		integer 							:: time_ini, time_end, time_calc

		double precision					:: start, finish, N_bar


		common /domain/ xmin, xmax, ymin,ymax
		common /grid/ dx, dy

		call cpu_time(start)
		write(*,*) 'Program start...'

		call system_clock(count = time_ini)


		! Domain. Keep it like this so that for other problems (4 contacts or R-T) it will be easier to modify the domain

		xmin = 0.0d0
		xmax = 0.020d0

		ymin =  0.01d0
		ymax =  0.020d0

		N_bar=1.0d0/sqrt((1.0d0/(xmax-xmin)**2.0d0)+(1.0d0/(ymax-ymin)**2.0d0))

		write(*,*) N_bar

		! stop

		! Generate simple grid

	
	    dx = (xmax - xmin)/(NX)
		dy = (ymax - ymin)/(NY)


			do i=-ghostp,NX+ghostp
	        	X(i)=(i-0.5d0)*dx+xmin
	    	enddo

	    	do j=-ghostp,NY+ghostp
	        	y(j)=(j-0.5d0)*dy+ymin
	    	enddo


	!***********************************************************************
	!*****                       call initial conditions               *****
	!***********************************************************************
		open(unit=4, file="residual.txt",action="write",status="replace")

		call initialconditions(NX,NY,ghostp,cons,n_eqn)

		call timestep(NX,NY,dt,dx,dy,CFL,x,y,ghostp)


	!***********************************************************************
	!*****         			 TVD- Runge Kutta 				       	   *****
	!***********************************************************************
		
	iteration : do NT=1,NTMAX

		
		cons_old = cons

	!***********************************************************************
	!*****                       Time step -1	                       *****
	!***********************************************************************

			call boundarycondition(NX,NY,cons,n_eqn,ghostp,x,y,dx,dy)
			
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

			call boundarycondition(NX,NY,cons,n_eqn,ghostp,x,y,dx,dy)
			
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

			call boundarycondition(NX,NY,cons,n_eqn,ghostp,x,y,dx,dy)
			
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
  	write(*,*) 'Minimum Potential: ', MINVAL(cons(1,1:NX,1:NY))
  	write(*,*)
		  call output(NX,NY,x,y,cons,n_eqn,ghostp)
		 endif






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
	




        call output(NX,NY,x,y,cons,n_eqn,ghostp)


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
	!*****                       Compute time step                	   *****
	!***********************************************************************

 	subroutine timestep(NX,NY,dt,dx,dy,CFL,x,y,ghostp)
 	implicit none

 	integer 							:: i, j
 	
 	integer								:: NX, NY,ghostp

 	double precision 					:: dx, dy
	double precision 					:: dt, CFL, dtnew
	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision, parameter 		:: pi=acos(-1.0d0), beta =  45.0d0, D_parallel = 1.0d3

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
					

					Lr = 1.0d0/(sqrt(4.0d0)*pi)
					Lr = (2.0d0*0.0089442d0*1.0d0)/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
					Tr = (Lr**2.0d0)/(d_xx+2.0d0*d_xy+d_yy)

  				  x_velocity =  DABS(sqrt(d_xx/Tr))	
				  
				

  					Lr = 1.0d0/(sqrt(4.0d0)*pi)
  					Lr = (2.0d0*0.0089442d0*1.0d0)/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NY)+4.0d0))
					Tr = (Lr**2.0d0)/(d_xx+2.0d0*d_xy+d_yy)
				  
				  y_velocity =  DABS(sqrt(d_yy/Tr))

				  dtnew = min( (x(i+1)-x(i))/x_velocity, (y(j+1)-y(j))/y_velocity)
				    
				   if(dtnew .lt. dt) dt = dtnew

				   ! write(*,*) d_xx, d_xy
				   ! pause

				
			enddo
		enddo

  		  CFL 	=	0.45d0
  		  dt 	= 	CFL*dt


 	end subroutine timestep

 	!***********************************************************************
	!*****                       Boundary conditions                   *****
	!***********************************************************************

	subroutine boundarycondition(NX,NY,cons,n_eqn,ghostp,x,y,dx,dy)

	implicit none

	integer 							:: i, j, k
	double precision, parameter 		:: pi=acos(-1.0d0), beta =  45.0d0, D_parallel = 1.0d3
	integer								:: NX, NY, n_eqn, ghostp

	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp),dx,dy	
	double precision					:: cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	! Primary variable for Dirichlet boundary condition. Conservative variable 1

	double precision 					:: phi_left,phi_right, phi_top, phi_bottom

	double precision					:: u0,u1,u2, mu

	double precision 					:: dbd(0:2),dkprdxk(0:2,0:2)

	double precision 					:: beltabd(0:2),alphabd(0:2),omigabd(0:2)

	double precision 					:: vch(3,3),vbd(3,0:2),ubd(3,0:2)

	phi_left = 1.0d0; phi_top = 0.0d0; 	phi_right = 0.0d0;  phi_bottom = 0.0d0;

! Left boundary conditions

	do j = 1,NY
		
		i=0     
	    
	    cons(1,i,j )   = (8.0d0*phi_left - 6.0d0*cons(1,i+1,j) + cons(1,i+2,j)) * (1.0d0/3.0d0)	   

	    cons(2,i,j )   = cons(2,i+3,j)-3.0d0*cons(2,i+2,j)+3.0d0*cons(2,i+1,j)        
        cons(3,i,j )   = cons(3,i+3,j)-3.0d0*cons(3,i+2,j)+3.0d0*cons(3,i+1,j)

    enddo

    do j = 1,NY

    	vch=cons(:,0:2,j)

	do k=1,3

        u0=vch(k,1);u1=vch(k,2);u2=vch(k,3)

        beltabd(0)=dy*dy
        beltabd(1)=(u0-u1)*(u0-u1)
        beltabd(2)=(61.d0*u0*u0-196.d0*u0*u1+74.d0*u0*u2+ &
                   160.d0*u1*u1-124.d0*u1*u2+25*u2*u2)/12.d0

        dbd(0)=dy**2; dbd(1)=dy**1; dbd(2)=1-(dbd(0)+dbd(1))

        do i=0,2
            alphabd(i)=dbd(i)/(1.d-6+beltabd(i))**3
        enddo

        do i=0,2
            omigabd(i)=alphabd(i)/(alphabd(0)+alphabd(1)+alphabd(2))
        enddo

        dkprdxk=0.d0
        dkprdxk(0,0)=u0
        dkprdxk(1,0)=1.5*u0 - 0.5*u1
        dkprdxk(2,0)=1.875*u0 - 1.25*u1 + 0.375*u2

        dkprdxk(1,1)=u0/dy-u1/dy
        dkprdxk(2,1)=2.0*u0/dy - 3.0*u1/dy + u2/dy
        dkprdxk(2,2)=u0/(dy*dy) - 2.0*u1/(dy*dy) + u2/(dy*dy)   

        vbd(k,:)=matmul(omigabd,dkprdxk)

      enddo

    	do i = -1,-ghostp,-1

    		 mu=(0.5d0+i)*dy
      		 cons(:,i,j)=vbd(:,0)+mu*vbd(:,1)+0.5*mu*mu*vbd(:,2) 
        
        	! cons(:,i,j)   = cons(:,i+3,j)-3.0d0*cons(:,i+2,j)+3.0d0*cons(:,i+1,j)
        	! cons(1,i,j)   = cons(1,i+3,j)-3.0d0*cons(1,i+2,j)+3.0d0*cons(1,i+1,j)

        enddo
    enddo


	! Right boundary conditions


	do j = 1, NY
		
		i = NX

		cons(1,i+1,j) = (8.0d0*phi_right +1.0d0*cons(1,i-1,j) -6.0d0* cons(1,i,j)) * (1.0d0/3.0d0)
        
        cons(2,i+1,j) = 3.0d0*cons(2,i,j)   - 3.0d0*cons(2,i-1,j)+cons(2,i-2,j)
        cons(3,i+1,j) = 3.0d0*cons(3,i,j)   - 3.0d0*cons(3,i-1,j)+cons(3,i-2,j)
        
    	enddo

    do j = 1, NY


    	vch=cons(:,NX-1:NX+1,j)

      do k=1,3

        u0=vch(k,3);u1=vch(k,2);u2=vch(k,1);

        beltabd(0)=dy*dy
        beltabd(1)=(u0-u1)*(u0-u1)
        beltabd(2)=(61.d0*u0*u0-196.d0*u0*u1+74.d0*u0*u2+&
                   160.d0*u1*u1-124.d0*u1*u2+25*u2*u2)/12.d0


        dbd(0)=dy**2; dbd(1)=dy**1; dbd(2)=1-(dbd(0)+dbd(1))

        do i=0,2
            alphabd(i)=dbd(i)/(1.e-6+beltabd(i))**3
        enddo

        do i=0,2
            omigabd(i)=alphabd(i)/(alphabd(0)+alphabd(1)&
                     +alphabd(2))
        enddo


        dkprdxk=0.d0
        dkprdxk(0,0)=u0
        dkprdxk(1,0)=1.5*u0 - 0.5*u1
        dkprdxk(2,0)=1.875*u0 - 1.25*u1 + 0.375*u2  

        dkprdxk(1,1)=u0/dy-u1/dy
        dkprdxk(2,1)=2.0*u0/dy - 3.0*u1/dy + u2/dy
        dkprdxk(2,2)=u0/(dy*dy) - 2.0*u1/(dy*dy) + u2/(dy*dy)


        vbd(k,:)=matmul(omigabd,dkprdxk)
      enddo

    	do i= 2, 4

    	 mu=(1.5d0-i)*dy
     	 cons(:,NX+i,j)=vbd(:,0)+mu*vbd(:,1)+0.5*mu*mu*vbd(:,2)

        ! cons(:,i+2,j) = 3.0d0*cons(:,i+1,j) - 3.0d0*cons(:,i,j)  +cons(:,i-1,j)

        enddo
    enddo

	! Bottom boundary conditions (actually bottom)

	do i = 1,NX
		j=0 


		cons(1,i,j )   = cons(1,i,j+3)-3.0d0*cons(1,i,j+2)+3.0d0*cons(1,i,j+1)
		cons(2,i,j )   = (8.0d0*phi_bottom - 6.0d0*cons(2,i,j+1) + cons(2,i,j+2)) * (1.0d0/3.0d0)
	    cons(3,i,j )   = (8.0d0*phi_bottom - 6.0d0*cons(3,i,j+1) + cons(3,i,j+2)) * (1.0d0/3.0d0)

    	do j = -1,-ghostp,-1
        
        	cons(1:3,i,j)   = cons(1:3,i,j+3)-3.0d0*cons(1:3,i,j+2)+3.0d0*cons(1:3,i,j+1)

        enddo
    enddo


	! Top boundary conditions (actually top)
	do i = 1, NX
		j = NY

		cons(1,i,j+1) = 3.0d0*cons(1,i,j)   - 3.0d0*cons(1,i,j-1)+cons(1,i,j-2)
		cons(2,i,j+1) = (8.0d0*phi_top +1.0d0*cons(2,i,j-1) -6.0d0* cons(2,i,j)) * (1.0d0/3.0d0)
		cons(3,i,j+1) = (8.0d0*phi_top +1.0d0*cons(3,i,j-1) -6.0d0* cons(3,i,j)) * (1.0d0/3.0d0)

    	do j= NY, NY+6

        cons(1:3,i,j+2) = 3.0d0*cons(1:3,i,j+1) - 3.0d0*cons(1:3,i,j)  +cons(1:3,i,j-1)
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

	double precision 					:: nu,Lr,Tr(-ghostp:NX+ghostp)

	double precision, parameter 		:: pi=acos(-1.0d0), beta =  45.0d0, D_parallel = 1.0d3
	double precision 					:: d_xx(-ghostp:NX+ghostp), d_xy(-ghostp:NX+ghostp),d_yx(-ghostp:NX+ghostp),d_yy(-ghostp:NX+ghostp)

	! write(*,*) beta
	! pause

	do j =1, NY

		do i = -ghostp, NX+ghostp

				  d_xx(i) =   D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
				  d_xy(i) =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yx(i) =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yy(i) =   D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)



			cons_local(1,i) = cons(1,i,j)
			cons_local(2,i) = cons(2,i,j)
			cons_local(3,i) = cons(3,i,j)

			flux_node_x(1,i) = -cons_local(2,i)
			flux_node_x(2,i) = -cons_local(1,i)*d_xx(i)
			flux_node_x(3,i) = -cons_local(1,i)*d_xy(i)
		end do

		call reconstruction5EX (NX,NY,cons_local,flux_x_half, n_eqn,ghostp,x,y,j,d_xx,d_xy,d_yx,d_yy)

		do k = 1,n_eqn

			do i = 1, NX

				!Finite volume or Reconstruction approach 
	        	flux_x(k,i,j) = -((flux_x_half(k,i)-flux_x_half(k,i-1))/dx)

	        	!Finite difference or interpolation approach
				! flux_x(k,i,j) = -((4.0/3.0d0)*(flux_x_half(k,i  )-flux_x_half(k,i-1)) +&
				!                  (-1.0/6.0d0)*(flux_node_x(k,i+1)-flux_node_x(k,i-1)))/dx

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

	subroutine reconstruction5EX (NX,NY,cons_local,flux, n_eqn,ghostp,x,y,j,d_xx,d_xy,d_yx,d_yy)

	implicit none 
	integer, intent(in)					:: j
	integer								:: i, k	
	integer								:: NX,NY, n_eqn, ghostp, total
	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision					:: cons_local(n_eqn,-ghostp:NX+ghostp),flux(n_eqn,-ghostp:NX+ghostp)

	double precision 					:: consl(n_eqn,-ghostp:NX+ghostp), consr(n_eqn,-ghostp:NX+ghostp)

	double precision   					:: fleft(n_eqn,-ghostp:NX+ghostp), fright(n_eqn,-ghostp:NX+ghostp)

	double precision					:: dissipation(n_eqn,-3:NX+3),jacobian(3,3)

	double precision  					:: delphi(-ghostp:NX+ghostp), delp(-ghostp:NX+ghostp), delq(-ghostp:NX+ghostp)


	double precision 					:: nu,Lr,Tr
	double precision					:: d_xx_1,d_xy_1,d_yy_1
	double precision 					:: d_xx(-ghostp:NX+ghostp), d_xy(-ghostp:NX+ghostp)
	double precision 					:: d_yy(-ghostp:NY+ghostp), d_yx(-ghostp:NY+ghostp)
	
	double precision 					:: d_xxl(-ghostp:NX+ghostp), d_xxr(-ghostp:NX+ghostp)
	double precision 					:: d_xyl(-ghostp:NX+ghostp), d_xyr(-ghostp:NX+ghostp)

	double precision 					:: d_yxl(-ghostp:NY+ghostp),d_yxr(-ghostp:NY+ghostp)
	double precision 					:: d_yyl(-ghostp:NY+ghostp),d_yyr(-ghostp:NY+ghostp)

	double precision, parameter			:: epsilon = 1.0d-23
    double precision 					:: recon_poly5_r(0:3), recon_poly5_l(0:3)
    double precision 					:: alpha_l(0:3),alpha_r(0:3)
    double precision 					:: beta_w(0:3), omega_l(0:3), omega_r(0:3)

    double precision, parameter 		:: pi=acos(-1.0d0), beta =  45.0d0, D_parallel = 1.0d3
   
    double precision, parameter 		:: Constant=1.0d0,p =2.0d0


		do k=1,n_eqn

        	do i =-4,NX+4

        	 	beta_w(0) 	= (13.0d0/12.0d0)*(cons_local(k,i)-2.0d0*cons_local(k,i+1)+cons_local(k,i+2))**2+ &
	            				(1.0/4.0d0)*(3.0d0*cons_local(k,i)-4.0d0*cons_local(k,i+1)+cons_local(k,i+2))**2

	            beta_w(1) 	= (13.0d0/12.0d0)*(cons_local(k,i-1)-2.d0*cons_local(k,i)+cons_local(k,i+1))**2+ &
	            				(1.0/4.0d0)*(cons_local(k,i-1)-cons_local(k,i+1))**2

	            beta_w(2)		= (13.0d0/12.0d0)*(cons_local(k,i-2)-2.d0*cons_local(k,i-1)+cons_local(k,i))**2 &
	            				+(1.0/4.0d0)*(cons_local(k,i-2)-4.0d0*cons_local(k,i-1)+3.d0*cons_local(k,i))**2


   	           	recon_poly5_l(0) =1.0d0/8.0d0*( 3.0d0*cons_local(k,i  )+ 6.0d0*cons_local(k,i+1)-1.0d0 *cons_local(k,i+2))
	            recon_poly5_l(1) =1.0d0/8.0d0*(-1.0d0*cons_local(k,i-1)+ 6.0d0*cons_local(k,i  )+3.0d0 *cons_local(k,i+1))
	            recon_poly5_l(2) =1.0d0/8.0d0*( 3.0d0*cons_local(k,i-2)-10.0d0*cons_local(k,i-1)+15.0d0*cons_local(k,i  ))

	! Jiang-Shu linear weights
	            ! alpha_l(0) 		= (5.0/16.0d0)
	            ! alpha_l(1) 		= (10.0/16.0d0)
	            ! alpha_l(2) 		= (1.0/16.0d0)
! ! Borges
	        alpha_l(0) = 5.0d0/16.0d0  * (Constant + (abs(beta_w(0)-beta_w(2)) / (epsilon + beta_w(0)))**p)
            alpha_l(1) = 10.0d0/16.0d0 * (Constant + (abs(beta_w(0)-beta_w(2)) / (epsilon + beta_w(1)))**p)
            alpha_l(2) = 1.0d0/16.0d0  * (Constant + (abs(beta_w(0)-beta_w(2)) / (epsilon + beta_w(2)))**p)
	        ! alpha_l(0) = 3.0d0/10.0d0  * (1.0d0 / (epsilon + beta(0))**2)
         !    alpha_l(1) = 6.0d0/10.0d0 * (1.0d0 / (epsilon + beta(1))**2)
         !    alpha_l(2) = 1.0d0/10.0d0  * (1.0d0 / (epsilon + beta(2))**2)

	            omega_l(0)=alpha_l(0)/(alpha_l(0)+alpha_l(1)+alpha_l(2))
	            omega_l(1)=alpha_l(1)/(alpha_l(0)+alpha_l(1)+alpha_l(2))
	            omega_l(2)=alpha_l(2)/(alpha_l(0)+alpha_l(1)+alpha_l(2))

	            consl(k,i) = omega_l(0) * recon_poly5_l(0) + omega_l(1) * recon_poly5_l(1) + omega_l(2) * recon_poly5_l(2)



	! right of i+1/2
	            beta_w(0)	= (13.0d0/12.0d0)*(cons_local(k,i+1)-2.0d0*cons_local(k,i+2)+cons_local(k,i+3))**2+ &
	            						(1.0/4.0d0)*(3.0d0*cons_local(k,i+1)-4.0d0*cons_local(k,i+2)+cons_local(k,i+3))**2

	            beta_w(1)	= (13.0d0/12.0d0)*(cons_local(k,i)-2.d0*cons_local(k,i+1)+cons_local(k,i+2))**2+ &
	           							(1.0/4.0d0)*(cons_local(k,i)-cons_local(k,i+2))**2

	            beta_w(2) = (13.0d0/12.0d0)*(cons_local(k,i-1)-2.d0*cons_local(k,i)+cons_local(k,i+1))**2 &
	            						+(1.0/4.0d0)*(cons_local(k,i-1)-4.0d0*cons_local(k,i)+3.d0*cons_local(k,i+1))**2


! 				

            	recon_poly5_r(0) = 1.0d0/8.0d0*(15.0d0*cons_local(k,i+1)-10.0d0*cons_local(k,i+2)+3.0d0*cons_local(k,i+3))
	            recon_poly5_r(1) = 1.0d0/8.0d0*( 3.0d0*cons_local(k,i  )+ 6.0d0*cons_local(k,i+1)-1.0d0*cons_local(k,i+2))
	            recon_poly5_r(2) = 1.0d0/8.0d0*(-1.0d0*cons_local(k,i-1)+ 6.0d0*cons_local(k,i  )+3.0d0*cons_local(k,i+1))

	        alpha_r(0) = 1.0d0/16.0d0  * (Constant + (abs(beta_w(0)-beta_w(2)) / (epsilon + beta_w(0)))**p)
            alpha_r(1) = 10.0d0/16.0d0 * (Constant + (abs(beta_w(0)-beta_w(2)) / (epsilon + beta_w(1)))**p)
            alpha_r(2) = 5.0d0/16.0d0  * (Constant + (abs(beta_w(0)-beta_w(2)) / (epsilon + beta_w(2)))**p)


	! Non-linear weights
	            omega_r(0)=alpha_r(0)/(alpha_r(0)+alpha_r(1)+alpha_r(2))
	            omega_r(1)=alpha_r(1)/(alpha_r(0)+alpha_r(1)+alpha_r(2))
	            omega_r(2)=alpha_r(2)/(alpha_r(0)+alpha_r(1)+alpha_r(2))

	            
	            consr(k,i) = omega_r(0) * recon_poly5_r(0) + omega_r(1) * recon_poly5_r(1) + omega_r(2) * recon_poly5_r(2)

        
         	enddo
        enddo

    

		do i=-4,NX+4

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

				d_xx_1 = (d_xxl(i)+d_xxr(i))*0.5d0

				d_xy_1 = (d_xyl(i)+d_xyr(i))*0.5d0

				d_yy_1 = (d_yyl(i)+d_yyr(i))*0.5d0

				  	Lr = 1.0d0/(sqrt(4.0d0)*pi)
				  	Lr = (2.0d0*0.0089442d0*1.0d0)/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
	

				Tr = (Lr**2.0d0)/(d_xx_1+2.0d0*d_xy_1+d_yy_1)


			fleft(1,i)	   = -consl(2,i)
			fleft(2,i)	   = -consl(1,i)*d_xxl(i)
			fleft(3,i)	   = -consl(1,i)*d_xyl(i)

			fright(1,i)	   = -consr(2,i)
			fright(2,i)	   = -consr(1,i)*d_xxr(i)
			fright(3,i)	   = -consr(1,i)*d_xyr(i)


			delphi(i)			= consr(1,i)-consl(1,i)
			delp(i)				= consr(2,i)-consl(2,i)
			delq(i)				= consr(3,i)-consl(3,i)

	!***************************************
	!***** 		Roe solver	       	   *****
	!***************************************

		 !***** Roe solver          *****
		 !The upwind dissipation matrix is |A| = R*|Lambda|*L
		 jacobian(1,1) = abs(sqrt(d_xx_1/Tr))
		 jacobian(1,2) = 0.0d0
		 jacobian(1,3) = 0.0d0
		 jacobian(2,1) = 0.0d0
		 jacobian(2,2) = Tr*abs(sqrt(d_xx_1/Tr))
		 jacobian(2,3) = 0.0d0
		 jacobian(3,1) = 0.0d0
		 jacobian(3,2) = Tr*abs(sqrt(d_xx_1/Tr))*d_xy_1/d_xx_1
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

	double precision 					:: nu,Lr,Tr(-ghostp:NY+ghostp)
	

	double precision, parameter 		:: pi=acos(-1.0d0), beta =  45.0d0, D_parallel = 1.0d3
	double precision 					:: d_xx(-ghostp:NX+ghostp), d_xy(-ghostp:NX+ghostp),d_yx(-ghostp:NX+ghostp),d_yy(-ghostp:NX+ghostp)


	
	do j =1, NX
		
		do i = -ghostp, NY+ghostp

				  d_xx(i) =   D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
				  d_xy(i) =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yx(i) =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yy(i) =   D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)



			cons_local(1,i) = cons(1,j,i)
			cons_local(2,i) = cons(2,j,i)
			cons_local(3,i) = cons(3,j,i)


			flux_node_y(1,i) =  -cons_local(3,i)
			flux_node_y(2,i) =  -cons_local(1,i)*d_yx(i)
			flux_node_y(3,i) =  -cons_local(1,i)*d_yy(i)
		end do

		call reconstruction5EY (NX,NY,cons_local,flux_y_half, n_eqn,ghostp,x,y,j,d_xx,d_xy,d_yx,d_yy)

		do k = 1,n_eqn

			do i = 1, NY

				!Finite volume or Reconstruction approach 
	        	flux_y(k,j,i) = -((flux_y_half(k,i)-flux_y_half(k,i-1))/dy)

				! flux_y(k,j,i) = -((4.0/3.0d0)*(flux_y_half(k,i  )-flux_y_half(k,i-1)) +&
				!                 (-1.0/6.0d0)*(flux_node_y(k,i+1)-flux_node_y(k,i-1)))/dy


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

	subroutine reconstruction5EY (NX,NY,cons_local,flux, n_eqn,ghostp,x,y,j,d_xx,d_xy,d_yx,d_yy)

	implicit none 		
		

	integer, intent(in)					:: j
	integer								:: i,k	
	integer								:: NX,NY, n_eqn, ghostp, total
	double precision					:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision					:: cons_local(n_eqn,-ghostp:NY+ghostp),flux(n_eqn,-ghostp:NY+ghostp)

	double precision 					:: consl(n_eqn,-ghostp:NY+ghostp), consr(n_eqn,-ghostp:NY+ghostp)

	double precision   					:: fleft(n_eqn,-ghostp:NY+ghostp), fright(n_eqn,-ghostp:NY+ghostp)

	double precision					:: dissipation(n_eqn,-3:NY+3),jacobian(3,3)


	double precision  					:: delphi(-ghostp:NY+ghostp), delp(-ghostp:NY+ghostp), delq(-ghostp:NY+ghostp)


	double precision 					:: nu,Lr,Tr
	double precision					:: d_xx_1,d_xy_1,d_yy_1,d_yx_1
	double precision 					:: d_xx(-ghostp:NX+ghostp), d_xy(-ghostp:NX+ghostp)
	double precision 					:: d_yy(-ghostp:NY+ghostp), d_yx(-ghostp:NY+ghostp)
	
	double precision 					:: d_xxl(-ghostp:NX+ghostp), d_xxr(-ghostp:NX+ghostp)
	double precision 					:: d_xyl(-ghostp:NX+ghostp), d_xyr(-ghostp:NX+ghostp)

	double precision 					:: d_yxl(-ghostp:NY+ghostp),d_yxr(-ghostp:NY+ghostp)
	double precision 					:: d_yyl(-ghostp:NY+ghostp),d_yyr(-ghostp:NY+ghostp)

	double precision, parameter			:: epsilon = 1.0d-23
    double precision 					:: recon_poly5_r(0:3), recon_poly5_l(0:3)
    double precision 					:: alpha_l(0:3),alpha_r(0:3)
    double precision 					:: beta_w(0:3), omega_l(0:3), omega_r(0:3)
   
    double precision, parameter 		:: Constant=1.0d0,p =1.0d0

    double precision, parameter 		:: pi=acos(-1.0d0), beta =  45.0d0, D_parallel = 1.0d3





		do k=1,n_eqn

        	do i =-4,NY+4

        	 
        	 	beta_w(0) 	= (13.0d0/12.0d0)*(cons_local(k,i)-2.0d0*cons_local(k,i+1)+cons_local(k,i+2))**2+ &
	            				(1.0/4.0d0)*(3.0d0*cons_local(k,i)-4.0d0*cons_local(k,i+1)+cons_local(k,i+2))**2

	            beta_w(1) 	= (13.0d0/12.0d0)*(cons_local(k,i-1)-2.d0*cons_local(k,i)+cons_local(k,i+1))**2+ &
	            				(1.0/4.0d0)*(cons_local(k,i-1)-cons_local(k,i+1))**2

	            beta_w(2)	= (13.0d0/12.0d0)*(cons_local(k,i-2)-2.d0*cons_local(k,i-1)+cons_local(k,i))**2 &
	            				+(1.0/4.0d0)*(cons_local(k,i-2)-4.0d0*cons_local(k,i-1)+3.d0*cons_local(k,i))**2


   	           	recon_poly5_l(0) =1.0d0/8.0d0*( 3.0d0*cons_local(k,i  )+ 6.0d0*cons_local(k,i+1)-1.0d0 *cons_local(k,i+2))
	            recon_poly5_l(1) =1.0d0/8.0d0*(-1.0d0*cons_local(k,i-1)+ 6.0d0*cons_local(k,i  )+3.0d0 *cons_local(k,i+1))
	            recon_poly5_l(2) =1.0d0/8.0d0*( 3.0d0*cons_local(k,i-2)-10.0d0*cons_local(k,i-1)+15.0d0*cons_local(k,i  ))


! ! Borges
	        alpha_l(0) = 5.0d0/16.0d0  * (Constant + (abs(beta_w(0)-beta_w(2)) / (epsilon + beta_w(0)))**p)
            alpha_l(1) = 10.0d0/16.0d0 * (Constant + (abs(beta_w(0)-beta_w(2)) / (epsilon + beta_w(1)))**p)
            alpha_l(2) = 1.0d0/16.0d0  * (Constant + (abs(beta_w(0)-beta_w(2)) / (epsilon + beta_w(2)))**p)


	            omega_l(0)=alpha_l(0)/(alpha_l(0)+alpha_l(1)+alpha_l(2))
	            omega_l(1)=alpha_l(1)/(alpha_l(0)+alpha_l(1)+alpha_l(2))
	            omega_l(2)=alpha_l(2)/(alpha_l(0)+alpha_l(1)+alpha_l(2))

	            consl(k,i) = omega_l(0) * recon_poly5_l(0) + omega_l(1) * recon_poly5_l(1) + omega_l(2) * recon_poly5_l(2)



	! right of i+1/2
	            beta_w(0)	= (13.0d0/12.0d0)*(cons_local(k,i+1)-2.0d0*cons_local(k,i+2)+cons_local(k,i+3))**2+ &
	            						(1.0/4.0d0)*(3.0d0*cons_local(k,i+1)-4.0d0*cons_local(k,i+2)+cons_local(k,i+3))**2

	            beta_w(1)	= (13.0d0/12.0d0)*(cons_local(k,i)-2.d0*cons_local(k,i+1)+cons_local(k,i+2))**2+ &
	           							(1.0/4.0d0)*(cons_local(k,i)-cons_local(k,i+2))**2

	            beta_w(2) = (13.0d0/12.0d0)*(cons_local(k,i-1)-2.d0*cons_local(k,i)+cons_local(k,i+1))**2 &
	            						+(1.0/4.0d0)*(cons_local(k,i-1)-4.0d0*cons_local(k,i)+3.d0*cons_local(k,i+1))**2


! 				

            	recon_poly5_r(0) = 1.0d0/8.0d0*(15.0d0*cons_local(k,i+1)-10.0d0*cons_local(k,i+2)+3.0d0*cons_local(k,i+3))
	            recon_poly5_r(1) = 1.0d0/8.0d0*( 3.0d0*cons_local(k,i  )+ 6.0d0*cons_local(k,i+1)-1.0d0*cons_local(k,i+2))
	            recon_poly5_r(2) = 1.0d0/8.0d0*(-1.0d0*cons_local(k,i-1)+ 6.0d0*cons_local(k,i  )+3.0d0*cons_local(k,i+1))

	        alpha_r(0) = 1.0d0/16.0d0  * (Constant + (abs(beta_w(0)-beta_w(2)) / (epsilon + beta_w(0)))**p)
            alpha_r(1) = 10.0d0/16.0d0 * (Constant + (abs(beta_w(0)-beta_w(2)) / (epsilon + beta_w(1)))**p)
            alpha_r(2) = 5.0d0/16.0d0  * (Constant + (abs(beta_w(0)-beta_w(2)) / (epsilon + beta_w(2)))**p)


	! Non-linear weights
	            omega_r(0)=alpha_r(0)/(alpha_r(0)+alpha_r(1)+alpha_r(2))
	            omega_r(1)=alpha_r(1)/(alpha_r(0)+alpha_r(1)+alpha_r(2))
	            omega_r(2)=alpha_r(2)/(alpha_r(0)+alpha_r(1)+alpha_r(2))

	            
	            consr(k,i) = omega_r(0) * recon_poly5_r(0) + omega_r(1) * recon_poly5_r(1) + omega_r(2) * recon_poly5_r(2)

        
         	enddo
        enddo

    

		do i=-4,NY+4

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

				d_xx_1 = (d_xxl(i)+d_xxr(i))*0.5d0

				d_xy_1 = (d_xyl(i)+d_xyr(i))*0.5d0

				d_yy_1 = (d_yyl(i)+d_yyr(i))*0.5d0

					
			  	Lr = (2.0d0*0.0089442d0*1.0d0)/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NY)+4.0d0))
				Tr = (Lr**2.0d0)/(d_xx_1+2.0d0*d_xy_1+d_yy_1)

			fleft(1,i)	   = -consl(3,i)
			fleft(2,i)	   = -consl(1,i)*d_xyl(i)
			fleft(3,i)	   = -consl(1,i)*d_yyl(i)

			fright(1,i)	   = -consr(3,i)
			fright(2,i)	   = -consr(1,i)*d_xyr(i)
			fright(3,i)	   = -consr(1,i)*d_yyr(i)


			delphi(i)			= consr(1,i)-consl(1,i)
			delp(i)				= consr(2,i)-consl(2,i)
			delq(i)				= consr(3,i)-consl(3,i)

	!***************************************
	!***** 		Roe solver	       	   *****
	!***************************************

		 !***** Roe solver          *****
		 !The upwind dissipation matrix is |A| = R*|Lambda|*L
		 jacobian(1,1) = abs(sqrt(d_yy_1/Tr))
		 jacobian(1,2) = 0.0d0
		 jacobian(1,3) = 0.0d0
		 jacobian(2,1) = 0.0d0
		 jacobian(2,2) = 0.0d0
		 jacobian(2,3) = Tr*abs(sqrt(d_yy_1/Tr))*d_xy_1/d_yy_1
		 jacobian(3,1) = 0.0d0
		 jacobian(3,2) = 0.0d0
		 jacobian(3,3) = Tr*abs(sqrt(d_yy_1/Tr))

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
		double precision              		:: nu,Lr,Tr

		double precision, parameter 		:: pi=acos(-1.0d0), beta =  45.0d0, D_parallel = 1.0d3



	
			

			do j = 1, NY

				do i = 1, NX

	
				sourcet(1,i,j)		= 0.0d0
				  
				  d_xx =   D_parallel*cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)+sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)
				  d_xy =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
	  			  d_yx =   0.5d0*sin(2.0d0*beta*pi/180.0d0)*(D_parallel-1.0d0)
				  d_yy =   D_parallel*sin(beta*pi/180.0d0)*sin(beta*pi/180.0d0)+cos(beta*pi/180.0d0)*cos(beta*pi/180.0d0)


					Lr = 1.0d0/(sqrt(4.0d0)*pi)
					Lr = (2.0d0*0.0089442d0*1.0d0)/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NX)+4.0d0))
					Tr = (Lr**2.0d0)/(d_xx+2.0d0*d_xy+d_yy)
  				  sourcet(2,i,j)		= -cons(2,i,j)
				  



					Lr = 1.0d0/(sqrt(4.0d0)*pi)
					Lr = (2.0d0*0.0089442d0*1.0d0)/(acos(-1.0d0)*(acos(-1.0d0)*(1.0d0/NY)+4.0d0))
					Tr = (Lr**2.0d0)/(d_xx+2.0d0*d_xy+d_yy)				

				sourcet(3,i,j)		= -cons(3,i,j)

					residual(1,i,j)		= flux_x(1,i,j) + flux_y(1,i,j) + sourcet(1,i,j)
					residual(2,i,j)		= (flux_x(2,i,j) + flux_y(2,i,j) + sourcet(2,i,j))/Tr
					residual(3,i,j)		= (flux_x(3,i,j) + flux_y(3,i,j) + sourcet(3,i,j))/Tr
				enddo

		enddo
	end subroutine sourceterm


	!***********************************************************************
	!*****                       Output 			                   *****
	!***********************************************************************

	subroutine output(NX,NY,x,y,cons,n_eqn,ghostp)

	integer 				:: i, j

	integer					:: NX,NY, n_eqn, ghostp

	double precision 		:: X(-ghostp:NX+ghostp),Y(-ghostp:NY+ghostp),cons(n_eqn,-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	


	open(unit=21, file="Potential.txt",action="write",status="replace")
	open(unit=22, file="U_x.txt",action="write",status="replace")
	open(unit=23, file="U_y.txt",action="write",status="replace")
	
	
	
		do i=1,NX
		    
		    write(21,'(1600F14.7)')(cons(1,i,j),j=1,NY)
		    write(22,'(1600F14.7)')(cons(2,i,j),j=1,NY)
		    write(23,'(1600F14.7)')(cons(3,i,j),j=1,NY)
		
		enddo
 	
	close(21)
	close(22)
	close(23)
	close(24)

	open(unit=25, file="soln.txt",action="write",status="replace")
	do i=1,NX
	 do j=1,NY
	 
	  ! cell center

	  write(25,'(5es25.10)') x(i),y(j),cons(1,i,j),cons(2,i,j),cons(3,i,j)
	 enddo
	enddo
	close(25)



 	end subroutine output




