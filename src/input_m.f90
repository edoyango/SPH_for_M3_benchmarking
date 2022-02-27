module input_m
	
	use globvar,		only: ntotal,nvirt,ntotal_loc,nhalo_loc,nvirt_loc,parts,scale_k,maxnloc
	use globvar_para,	only: procid,numprocs,bounds_glob
	use mpi
	use param,			only: dim,irho,dxo,f,hsml,mp,np,op,pp
	use error_msg_m,	only: error_msg
	use output_m,		only: write_ini_config
	
	public:: input,virt_part
	
contains
	
	!==============================================================================================================================
	subroutine input(generate)
	! Generates initial physical particle configuration.
	! 2 cases: return only number of particles retrieved, or generating the particles
	
		implicit none
		logical,intent(in):: generate
		integer:: i,j,k,d,n,n_loc,n_loc_i,n_start,n_done
		real(f):: xi,yi
		real(f),parameter:: dx = dxo,dy = dxo, xl = 25d0, yl = 25d0
	
		select case (generate)
			
			case (.false.)
			
				ntotal = mp*np
			
			case (.true.)
				
				! how many particles to generate per process
				n_loc_i = ceiling(dble(ntotal)/numprocs)
				if (procid.eq.numprocs-1) then
					n_loc = ntotal - (numprocs-1)*n_loc_i
				else
					n_loc = n_loc_i
				end if
				n_start = procid*n_loc_i + 1
				n_done = n_start + n_loc_i - 1
				
				! stopping program if array bounds are exceeded
				if ( (procid.eq.0) .and. (n_loc.gt.maxnloc) ) call error_msg(1,1)
		
				! intitial setup
				n = 0
				ntotal_loc = 0
				do i = 1,mp
					do j = 1,np
						n = n + 1 ! tracking total number of particles generated
						! Only generating particles assigned to process
						if ( (n.ge.n_start) .and. (n.le.n_done) ) then
							ntotal_loc = ntotal_loc + 1
							parts(ntotal_loc)%indglob = n
							parts(ntotal_loc)%indloc = ntotal_loc
							parts(ntotal_loc)%x(1) = (i-0.5d0)*dx
							parts(ntotal_loc)%x(2) = (j-0.5d0)*dy
							parts(ntotal_loc)%vx(:) = 0d0
							parts(ntotal_loc)%itype = 1
							parts(ntotal_loc)%rho = irho
							parts(ntotal_loc)%p = 0d0
						end if
					end do
				end do
				
				call write_ini_config
			
		end select
		
	end subroutine input
	
	!==============================================================================================================================
	subroutine virt_part(generate)
	! Generates the virtual particle configuration. Can change over time or remain static
	! 2 cases: return only number of particles retrieved, or generating the particles
		
		implicit none
		integer:: i,k,d,n
		real(f):: xi,yi,xmin_loc,ymin_loc,xmax_loc,ymax_loc
		logical,intent(in):: generate
		real(f),parameter:: dx = dxo, dy = dxo, xmin = 0d0, ymin = 0d0, xmax = 75d0, ymax = 40d0
		
		select case (generate)
			
			case (.false.)
				
				nvirt = 2*op + 2*pp
				
			case (.true.)
				
				xmin_loc = bounds_glob(1,procid+1) - scale_k*hsml
				ymin_loc = bounds_glob(2,procid+1) - scale_k*hsml
				xmax_loc = bounds_glob(3,procid+1) + scale_k*hsml
				ymax_loc = bounds_glob(4,procid+1) + scale_k*hsml
				
				nvirt_loc = 0
				n = ntotal ! counter used to track particle indices
				
				! lower boundary
				do i = 1,op
					n = n + 1
					xi = xmin + (i-0.5d0)*dx
					yi = ymin - 0.5d0*dy
					if ( (xi.ge.xmin_loc) .and. (yi.ge.ymin_loc) .and. (xi.le.xmax_loc) .and. (yi.le.ymax_loc) ) then
						nvirt_loc = nvirt_loc + 1
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%indglob = n
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%indloc = ntotal_loc+nhalo_loc+nvirt_loc
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%itype = -1
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(1) = xi
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(2) = yi
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%vx(:) = 0d0
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%rho = irho
					end if
				end do
				
				! upper boundary
				do i = 1,op
					n = n + 1
					xi = xmin + (i-0.5d0)*dx
					yi = ymax - 1.5d0*dy
					if ( (xi.ge.xmin_loc) .and. (yi.ge.ymin_loc) .and. (xi.le.xmax_loc) .and. (yi.le.ymax_loc) ) then
						nvirt_loc = nvirt_loc + 1
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%indglob = n
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%indloc = ntotal_loc+nhalo_loc+nvirt_loc
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%itype = -1
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(1) = xi
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(2) = yi
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%vx(:) = 0d0
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%rho = irho
					end if
				end do
				
				! left boundary
				do i = 1,pp
					n = n + 1
					xi = xmin - 0.5d0*dx
					yi = ymin + (i-1.5d0)*dy
					if ( (xi.ge.xmin_loc) .and. (yi.ge.ymin_loc) .and. (xi.le.xmax_loc) .and. (yi.le.ymax_loc) ) then
						nvirt_loc = nvirt_loc + 1
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%indglob = n
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%indloc = ntotal_loc+nhalo_loc+nvirt_loc
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%itype = -1
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(1) = xi
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(2) = yi
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%vx(:) = 0d0
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%rho = irho
					end if
				end do
				
				! right boundary
				do i = 1,pp
					n = n + 1
					xi = xmax + 0.5d0*dx
					yi = ymin + (i-1.5d0)*dy
					if ( (xi.ge.xmin_loc) .and. (yi.ge.ymin_loc) .and. (xi.le.xmax_loc) .and. (yi.le.ymax_loc) ) then
						nvirt_loc = nvirt_loc + 1
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%indglob = n
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%indloc = ntotal_loc+nhalo_loc+nvirt_loc
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%itype = -1
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(1) = xi
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(2) = yi
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%vx(:) = 0d0
						parts(ntotal_loc+nhalo_loc+nvirt_loc)%rho = irho
					end if
				end do
				
		end select
						
	end subroutine virt_part

end module input_m