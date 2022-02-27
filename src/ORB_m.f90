module ORB_m

	use globvar,		only: scale_k
	use globvar_para,	only: procid,numprocs,pincell_ORB,ierr,MPI_ftype,repartition_mode,node_cax,node_cut
	use mpi
	use param,			only: f,dim,hsml
	
	public:: ORB
	private:: particle_grid,ORB_bounds,P_trim,potential_neighbour_process_search,subdomain_neighbour
	
contains	
	!==============================================================================================================================
	subroutine ORB
	! Container subroutine for the bulk of the ORB algorithm, including the initial exchange of physical and halo particles
		use globvar,		only: parts,itimestep,ntotal,ntotal_loc,nhalo_loc,t_graph,t_dist
		use globvar_para,	only: n_reorients,n_parts,mintstep_bn_part,maxtstep_bn_part,mintstep_bn_reorient,maxtstep_bn_reorient,&
							box_ratio_previous,maxnode,bounds_glob,node_cax,node_segment,repartition_mode,prev_part_tstep,&
							prev_reorient_tstep,nphys_send,nphys_recv,nhalo_send,nhalo_recv,halotype_indexed,&
							haloupdatetype_indexed,prev_load,n_process_neighbour
		use param_para,		only: dcell_ORB,ORBcheck1,ORBcheck2,box_ratio_threshold
		
		use ORB_sr_m,		only: ORB_sendrecv_diffuse,ORB_sendrecv_halo
		use input_m,		only: virt_part
		
		implicit none
		integer:: i,j,k,ngridx(dim),current_load,nphys_recv_all,request_phys(2*numprocs),request_halo(2*numprocs),&
		status(MPI_STATUS_SIZE,4*numprocs),nphys_send_all,procrange_ini(2),tree_layers,gridind_ini(dim,2),diffusedepth,&
		searchrange_ini(2),n_request
		real(f):: bounds_out(2*dim),dcell,mingridx_ini(dim),maxgridx_ini(dim),current_to_previous(dim),box_ratio_current(3)
		
		!allocating partitioning arrays and initialising diagnostic variables -------------------------------------------------------------
		t_graph = t_graph - MPI_WTIME() ! commence timing of ORB algoirthm
		if (itimestep.eq.1) then
			n_reorients = 0
			n_parts = 0
			mintstep_bn_part = HUGE(1)
			maxtstep_bn_part = 0
			mintstep_bn_reorient = HUGE(1)
			maxtstep_bn_reorient = 0
			box_ratio_previous(:) = HUGE(1.d0)
			tree_layers = CEILING(LOG(DBLE(numprocs))/LOG(2.d0))
			maxnode = 2*2**tree_layers-1
			allocate( bounds_glob(2*dim,numprocs),&
					node_cax(maxnode),&
					node_segment(maxnode) )
		endif
		
		! Boundary Determiniation Algorithm ---------------------------------------------------------------------------------------
		repartition_mode = 1 !initially assumes no partition
		dcell = hsml*dcell_ORB
		
		! only checks if boundary needs updating every 50-100 time-steps
		if ( (itimestep.eq.1) .or. ((itimestep-prev_part_tstep.ge.ORBcheck1) .and. (mod(itimestep-prev_part_tstep,ORBcheck2).eq.0)) ) then
		
			! checking if change in partilces on current process > 5%
			current_load = ntotal_loc
			if (current_load.gt.prev_load+0.05d0*DBLE(ntotal)/DBLE(numprocs)) then
				repartition_mode = 2
			endif
		
			call MPI_ALLREDUCE(MPI_IN_PLACE,repartition_mode,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
		
			if ( (repartition_mode.gt.1) .or. (itimestep.eq.1) ) then
			
				n_parts = n_parts + 1
				if (itimestep.ne.1) then
					if (itimestep-prev_part_tstep.lt.mintstep_bn_part) mintstep_bn_part = itimestep-prev_part_tstep
					if (itimestep-prev_part_tstep.gt.maxtstep_bn_part) maxtstep_bn_part = itimestep-prev_part_tstep
				endif
				prev_part_tstep = itimestep
			
				call particle_grid( ngridx,dcell,mingridx_ini,maxgridx_ini )
			
				box_ratio_current(1) = DBLE(ngridx(2))/DBLE(ngridx(1))
				
				current_to_previous(1) = box_ratio_current(1)/box_ratio_previous(1)
				if ( (current_to_previous(1).gt.1.d0+box_ratio_threshold) .or. (current_to_previous(1).lt.1.d0/(1.d0+box_ratio_threshold)) ) &
					repartition_mode = 3
				
				!partition summary info
				if (repartition_mode.eq.3) then 
					box_ratio_previous(:) = box_ratio_current(:)
					if (itimestep-prev_reorient_tstep.gt.maxtstep_bn_reorient) maxtstep_bn_reorient = itimestep-prev_reorient_tstep
					if (itimestep-prev_reorient_tstep.lt.mintstep_bn_reorient) mintstep_bn_reorient = itimestep-prev_reorient_tstep
					prev_reorient_tstep = itimestep
					n_reorients = n_reorients + 1
				endif
			
				! determine subdomain boundaries using particle distribution
				gridind_ini(:,1) = 1
				gridind_ini(:,2) = ngridx(:)
				procrange_ini(1) = 0
				procrange_ini(2) = numprocs-1
				bounds_out = ORB_bounds( gridind_ini,numprocs,1,procrange_ini,ntotal,dcell,mingridx_ini,maxgridx_ini )
	
				call subdomain_neighbour
				
				! Updating sizes of select arrays to account for potential changes in neighbour list size
				if (itimestep.ne.1) deallocate(nphys_send,nphys_recv,nhalo_send,nhalo_recv,halotype_indexed,haloupdatetype_indexed)
				allocate( nphys_send(n_process_neighbour),&
						nphys_recv(n_process_neighbour),&
						nhalo_send(n_process_neighbour),&
						nhalo_recv(n_process_neighbour),&
						halotype_indexed(n_process_neighbour),&
						haloupdatetype_indexed(n_process_neighbour) )
			
			endif
		
		endif
		t_graph = t_graph + MPI_WTIME( ) ! conclude timing of ORB algorithm
	
		
		! Particle distribution (physical, halo) ----------------------------------------------------------------------------------
		t_dist = t_dist - MPI_WTIME() ! commence timing of particle distribution
		! physical particle distribution
		
		diffusedepth = 0
		searchrange_ini(:) = (/1,ntotal_loc/)
		i = ORB_sendrecv_diffuse( diffusedepth,searchrange_ini,n_request,request_phys,nphys_recv_all )
		
		! halo particle distribution
		call ORB_sendrecv_halo( request_phys,request_halo,nphys_recv_all,n_request )
		
		! update virtual particles
		call virt_part(.true.)
		
		if (repartition_mode.gt.1) prev_load = ntotal_loc
		
		do i = ntotal_loc+1,ntotal_loc+nhalo_loc
			parts(i)%indloc = i
			parts(i)%itype = 2
		end do
		
		! wait for halo particle distribution to complete
		call MPI_WAITALL(n_request,request_halo(1:n_request),status(:,1:n_request),ierr)
		
		t_dist = t_dist + MPI_WTIME( ) ! conclude timing of particle distribution
		
	end subroutine ORB
	
	!==============================================================================================================================
	subroutine particle_grid( ngridx,dcell,mingridx,maxgridx )
	! Subroutine to create a uniform rectangular grid with square cells, and counting the number of particles contained within each
	! cell. Each MPI process holds a global copy of the entire grid, in preperation for ORB
		
		use globvar,	only: parts,ntotal_loc
		
		implicit none
		integer:: i,j,k,d,n,icell,jcell,n_nonzerocells,n_nonzerocells_perprocess(numprocs),n_nonzerocells_total,pid,ngridx(dim),&
				displ(numprocs),icellmin,icellmax,jcellmin,jcellmax,icellrange,jcellrange,request(2),status(MPI_STATUS_SIZE,2)
		integer:: sendcount,recvcount(numprocs)
		real(f):: maxx(dim),minx(dim),mingridx(dim),maxgridx(dim),dcell
		integer,allocatable:: Plist_loc(:,:),Plist_all(:,:)
		
		! Local max, min, in each direction ---------------------------------------------------------------------------------------
		minx(:) = parts(1)%x(:)
		maxx(:) = parts(1)%x(:)
		
		do i = 2,ntotal_loc
			do d = 1,dim
				if (parts(i)%x(d).gt.maxx(d)) maxx(d) = parts(i)%x(d)
				if (parts(i)%x(d).lt.minx(d)) minx(d) = parts(i)%x(d)
			enddo
		enddo
		
		! Global max, min, in each direction --------------------------------------------------------------------------------------
		call MPI_IALLREDUCE(maxx(:),maxgridx(:),dim,MPI_ftype,MPI_MAX,MPI_COMM_WORLD,request(1),ierr) !finding max over all processes
		call MPI_IALLREDUCE(minx(:),mingridx(:),dim,MPI_ftype,MPI_MIN,MPI_COMM_WORLD,request(2),ierr) !finding min over all processes
		
		call MPI_WAITALL(2,request,status,ierr)
		
		mingridx(:) = mingridx(:) - dcell
		maxgridx(:) = maxgridx(:) + dcell
		
		! Number of grid cells and adjusting max extent in each direction ---------------------------------------------------------
		ngridx(:) = int((maxgridx(:)-mingridx(:))/dcell) + 1
		maxgridx(:) = mingridx(:) + dcell*ngridx(:)
		
		! Reducing search area by locating indices in which local particles are contained within
		icellmin =  int((minx(1)-mingridx(1))/dcell) + 1
		icellmax =  int((maxx(1)-mingridx(1))/dcell) + 1
		jcellmin =  int((minx(2)-mingridx(2))/dcell) + 1
		jcellmax =  int((maxx(2)-mingridx(2))/dcell) + 1
		icellrange = icellmax - icellmin + 1
		jcellrange = jcellmax - jcellmin + 1
		n = MIN(icellrange*jcellrange,ntotal_loc)
		
		! Creating list of non-zero grid cells ------------------------------------------------------------------------------------
		allocate( pincell_ORB(ngridx(1),ngridx(2)),&
				Plist_loc(dim+1,n) )
				
		pincell_ORB(icellmin:icellmax,jcellmin:jcellmax) = 0
		
		! Counting particles in each cell
		n_nonzerocells = 0
		do i = 1,ntotal_loc
			icell = int((parts(i)%x(1)-mingridx(1))/dcell) + 1
			jcell = int((parts(i)%x(2)-mingridx(2))/dcell) + 1
			if ( pincell_ORB(icell,jcell).eq.0) then
				n_nonzerocells = n_nonzerocells + 1
				Plist_loc(1,n_nonzerocells) = icell
				Plist_loc(2,n_nonzerocells) = jcell
			endif
			pincell_ORB(icell,jcell) = pincell_ORB(icell,jcell) + 1
		enddo
		
		! Exchanging info on how many non-zero cells each process has
		call MPI_IALLGATHER(n_nonzerocells,1,MPI_INTEGER,n_nonzerocells_perprocess,1,MPI_INTEGER,MPI_COMM_WORLD,request(1),ierr)
		
		! Populating Plist_loc with non-zero cell entries
		do i = 1,n_nonzerocells
			Plist_loc(3,i) = Pincell_ORB(Plist_loc(1,i),Plist_loc(2,i))
		enddo
		
		call MPI_WAIT(request(1),status(:,1),ierr)
		
		n_nonzerocells_total = 0
		do pid = 1,numprocs
			displ(pid) = n_nonzerocells_total
			n_nonzerocells_total = n_nonzerocells_total + n_nonzerocells_perprocess(pid)
		enddo
		
		allocate( Plist_all(dim+1,n_nonzerocells_total) )
		
		! Collecting all process' Plist_loc
		displ = (dim+1)*displ
		sendcount = (dim+1)*n_nonzerocells
		recvcount = (dim+1)*n_nonzerocells_perprocess(:)
		call MPI_IALLGATHERV(Plist_loc,sendcount,MPI_INTEGER,Plist_all,recvcount,displ,MPI_INTEGER,MPI_COMM_WORLD,request(1),ierr)
		
		pincell_ORB(:,:) = 0
		
		call MPI_WAIT(request(1),status(:,1),ierr)
		
		! populating the number of particles per cell
		do i = 1,n_nonzerocells_total
			pincell_ORB(Plist_all(1,i),Plist_all(2,i)) = pincell_ORB(Plist_all(1,i),Plist_all(2,i)) + Plist_all(3,i)
		enddo
		
		deallocate( Plist_loc,Plist_all )
	
	end subroutine particle_grid
	
	!==============================================================================================================================
	recursive function ORB_bounds( gridind_in,nprocs_in,node_in,procrange_in,ntotal_in,dcell,mingridx_in,maxgridx_in ) &
		result(bounds_out)
	! Recursive function that performs the 'bisection' part of the ORB algorithm
	
		use globvar_para,	only: leaf_node,bounds_glob
		use param_para,		only: bound_extend
	
		implicit none
		integer,intent(in):: gridind_in(dim,2),node_in,nprocs_in,procrange_in(2),ntotal_in
		real(f),intent(in):: mingridx_in(dim),maxgridx_in(dim),dcell
		integer:: i,j,k,d,node_out,gridind_out(dim,2),nprocs_out,ntotal_out,procrange_out(2),n_p,cax,np_per_node,pincol,&
		procrange_lo(2),procrange_hi(2),ngridx_trim(dim),A(2)
		real(f):: t1,t2,bounds_out(2*dim)

		!determining cut axis. 1 = x, 2 = y ---------------------------------------------------------------------------------------
		if (repartition_mode.eq.3) then
			call P_trim(gridind_in,ngridx_trim)
			cax = 1
			if (ngridx_trim(2) .gt. ngridx_trim(1)) cax = 2
			node_cax(node_in) = cax
		else
			cax = node_cax(node_in)
		endif

		!cut location -------------------------------------------------------------------------------------------------------------
		i = gridind_in(cax,1) - 1
		n_p = 0
		np_per_node = ceiling(real(nprocs_in)/2)/real(nprocs_in)*ntotal_in
		do while (n_p .lt. np_per_node)
			i = i + 1
			pincol = 0
			if (cax .eq. 1) then
				do j = gridind_in(2,1),gridind_in(2,2)
					pincol = pincol + pincell_ORB(i,j)
				enddo
			elseif (cax.eq.2) then
				do j = gridind_in(1,1),gridind_in(1,2)
					pincol = pincol + pincell_ORB(j,i)
				enddo
			endif
			n_p = n_p + pincol
		enddo
	
		if ( (np_per_node-(n_p-pincol).lt.n_p-np_per_node) ) then
			i = i - 1
			n_p = n_p - pincol
		endif

		!saving output information ------------------------------------------------------------------------------------------------
		procrange_lo(1) = procrange_in(1)
		procrange_lo(2) = procrange_in(1) + ceiling(real(nprocs_in)/2)-1
		procrange_hi(1) = procrange_lo(2) + 1
		procrange_hi(2) = procrange_lo(1) + nprocs_in - 1
		gridind_out(:,1) = gridind_in(:,1)
		gridind_out(:,2) = gridind_in(:,2)
		
		if (procid.le.procrange_lo(2)) then
			node_out = 2*node_in
			procrange_out = procrange_lo
			ntotal_out = n_p
			gridind_out(cax,1) = gridind_in(cax,1)
			gridind_out(cax,2) = i
		else
			node_out = 2*node_in + 1
			procrange_out = procrange_hi
			ntotal_out = ntotal_in - n_p
			gridind_out(cax,1) = i + 1
			gridind_out(cax,2) = gridind_in(cax,2)
		endif
		
		nprocs_out = procrange_out(2) - procrange_out(1) + 1
		
		!travelling to child node/saving boundary information ---------------------------------------------------------------------
		if (nprocs_out .ne. 1) then
			bounds_out = ORB_bounds( gridind_out,nprocs_out,node_out,procrange_out,ntotal_out,dcell,mingridx_in,maxgridx_in )
		else
			leaf_node = node_out
			
			bounds_out(1:dim) = mingridx_in(:) + (gridind_out(:,1)-1)*dcell
			bounds_out(dim+1:2*dim) = mingridx_in(:) + gridind_out(:,2)*dcell
			
			if (repartition_mode.eq.3) call potential_neighbour_process_search(leaf_node)
			
			if (node_cut(1).eq.0) bounds_out(3) = bounds_out(3) + bound_extend*scale_k*hsml
			if (node_cut(2).eq.0) bounds_out(1) = bounds_out(1) - bound_extend*scale_k*hsml
			if (node_cut(3).eq.0) bounds_out(4) = bounds_out(4) + bound_extend*scale_k*hsml
			if (node_cut(4).eq.0) bounds_out(2) = bounds_out(2) - bound_extend*scale_k*hsml
			
			call MPI_ALLGATHER(bounds_out,2*dim,MPI_ftype,bounds_glob,2*dim,MPI_ftype,MPI_COMM_WORLD,ierr)
			
			deallocate( pincell_ORB )
		
		endif
  
	end function ORB_bounds

	!==============================================================================================================================
	subroutine P_trim(gridind_in,ngridx_trim)
	! Trims particle-in-cell grid so as to obtain minimal bounding boxes to obtain accurate cut axis orientations
	
		implicit none
		integer,intent(in):: gridind_in(dim,2)
		integer,intent(out):: ngridx_trim(dim)
		integer:: i,j,k,d,new_start_index(dim),new_end_index(dim),i0,i1,j0,j1
		
		i0 = gridind_in(1,1)
		i1 = gridind_in(1,2)
		j0 = gridind_in(2,1)
		j1 = gridind_in(2,2)
	
		!Trimming input grid ------------------------------------------------------------------------------------------------------
		!reducing x-length of grid
		!finding new start x-index
		do i = i0,i1
			do j = j0,j1
				if (pincell_ORB(i,j).ne.0) then
					new_start_index(1) = i
					goto 1
				endif
			enddo
		enddo

		!finding new end x-index
1 		do i = i1,new_start_index(1),-1
			do j = j0,j1
				if (pincell_ORB(i,j).ne.0) then
					new_end_index(1) = i
					goto 2
				endif
			enddo
		enddo

		!reducing y-length of grid
		!finding new start y-index
2 		do j = j0,j1
			do i = new_start_index(1),new_end_index(1)
				if (pincell_ORB(i,j).ne.0) then
					new_start_index(2) = j
					goto 3
				endif
			enddo
		enddo

		!finding new end y-index
3 		do j = j1,new_start_index(2),-1
			do i = new_start_index(1),new_end_index(1)
				if (pincell_ORB(i,j).ne.0) then
					new_end_index(2) = j
					goto 4
				endif
			enddo
		enddo

4 		ngridx_trim(:) = new_end_index(:) - new_start_index(:) + 1

	end subroutine P_trim
		
	!==============================================================================================================================
	subroutine potential_neighbour_process_search(ID_node)
	! Used to determine which nodes are located at the edge of the domain in each direction.
	! IGNORE BELOW
	! Finds any potential neighbours by exploiting the binary tree data structure
	! First finds nodes in the tree that define the edge of the local process. Then finds processes that also use that node as an 
	! edge but on the opposite side. E.g. if local process has node 2 as an east edge, potential east neighbours of the local 
	! process are the processes that use node 2 as a west edge.
		
		implicit none
		integer:: i,j,k,d,pid,ID_node,start_node,blah,direction_skip
		
		!Initialization
		if (.not.allocated(node_cut)) allocate( node_cut(2*dim) )
		
		!finding nodes that define edge for current process		  
		node_cut(:) = 0
		do while (ID_node.ne.1)
		
			if (mod(ID_node,2).eq.0) then
				ID_node = ID_node/2
				if ( (node_cax(ID_node).eq.1) .and. (node_cut(1).eq.0) ) node_cut(1) = ID_node !right face of process defined by 'parent'
				if ( (node_cax(ID_node).eq.2) .and. (node_cut(3).eq.0) ) node_cut(3) = ID_node !top face of process defined by 'parent'
			else
				ID_node = (ID_node-1)/2
				if ( (node_cax(ID_node).eq.1) .and. (node_cut(2).eq.0) ) node_cut(2) = ID_node !left of process defined by 'parent'
				if ( (node_cax(ID_node).eq.2) .and. (node_cut(4).eq.0) ) node_cut(4) = ID_node !bottom of process defined by 'parent'
			endif
		
		enddo
	
	end subroutine potential_neighbour_process_search
	
	!==============================================================================================================================
	subroutine subdomain_neighbour
	!creates list of adjacent subdomains for the local subdomain by searching potential neighbours and seeing if they overlap
		
		use globvar_para,	only: proc_neighbour_list,n_process_neighbour,bounds_glob
		
		implicit none
		integer:: i,j,k,d,pid
		real(f):: xmin_loc,xmax_loc,xmin_rem,xmax_rem
		
		!Initialization
		if (repartition_mode.eq.3) then
			if(allocated(proc_neighbour_list)) deallocate(proc_neighbour_list)
			allocate( proc_neighbour_list(numprocs) )
		endif
		
		!Checking if potential neighbours overlap with local process. Overlap = adjacent.
		n_process_neighbour = 0
		do pid = 1,numprocs
			if (pid.ne.procid+1) then
				do d = 1,dim
					xmin_loc = bounds_glob(d,procid+1) - 0.5d0*scale_k*hsml
					xmax_loc = bounds_glob(dim+d,procid+1) + 0.5d0*scale_k*hsml
					xmin_rem = bounds_glob(d,pid) - 0.5d0*scale_k*hsml
					xmax_rem = bounds_glob(dim+d,pid) + 0.5d0*scale_k*hsml
					if ( (xmax_rem.lt.xmin_loc) .or. (xmin_rem.gt.xmax_loc) ) goto 1
				enddo
				n_process_neighbour = n_process_neighbour + 1
				proc_neighbour_list(n_process_neighbour) = pid - 1
			endif
	1 	enddo
	
	end subroutine subdomain_neighbour
	
end module ORB_m
