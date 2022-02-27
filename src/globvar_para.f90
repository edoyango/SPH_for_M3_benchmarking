module globvar_para
! A module that contains global variables that are needed for the parallel scheme
	
	use datatypes,	only: particles
	use param,		only: f
	
	implicit none
	
	! basic MPI variables
	integer:: procid,numprocs,ierr
	
	! derived MPI tpes
	integer:: parttype,halotype,haloupdatetype,MPI_ftype
	integer,allocatable:: halotype_indexed(:),haloupdatetype_indexed(:)

	! particle send/recv arrays -------------------------------------------------------------------------------------------------------
	integer,allocatable:: nphys_send(:),nphys_recv(:)
	integer,allocatable:: nhalo_send(:),nhalo_recv(:)
	type(particles),allocatable:: PhysPackSend(:,:)
	
	!ORB variables --------------------------------------------------------------------------------------------------------------------
	integer:: maxnode,n_process_neighbour,leaf_node,repartition_mode
	integer,allocatable:: node_cax(:),pincell_ORB(:,:),proc_neighbour_list(:),node_cut(:),node_segment(:),halo_pindex(:,:)
	real(f),allocatable:: bounds_glob(:,:)
	
	!Partition frequency variables
	integer:: prev_part_tstep,mintstep_bn_part,maxtstep_bn_part,n_parts
	integer:: prev_reorient_tstep,mintstep_bn_reorient,maxtstep_bn_reorient,n_reorients
	integer:: prev_load
	real(f):: box_ratio_previous(3)

end module globvar_para