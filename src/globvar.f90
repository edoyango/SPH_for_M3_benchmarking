module globvar
	
	use param, only: f
	use datatypes, only: particles,interactions
	
	implicit none
	
	! particle array
	type(particles),allocatable,target,public:: parts(:)
	
	! interaction array
	type(interactions),allocatable,public:: pairs(:)
	
	!global 1D variables
	integer,public:: ntotal,nvirt,ntotal_loc,nhalo_loc,nvirt_loc
	integer,public:: maxn,maxinter,maxnloc
	integer,public:: niac
	integer,public:: itimestep,maxtimestep,save_step,print_step
	real(f),public:: time
	
	real(f),public:: scale_k
	
	!timing
	real(f),public:: t_graph,t_dist,cputime,output_time,test_time
	
	public:: allocateGlobalArrays,deallocateGlobalArrays
	
! subroutines to allocate and deallocate global arrays
contains

	!==============================================================================================================================
	subroutine allocateGlobalArrays

		implicit none
		
		maxnloc = ntotal+nvirt
		maxinter = 11*maxnloc
		
		allocate( parts(maxnloc) )
		allocate( pairs(maxinter) )
				
	end subroutine allocateGlobalArrays
	
	!==============================================================================================================================
	subroutine deallocateGlobalArrays
	
		implicit none
		
		deallocate( parts )
		deallocate( pairs )
	
	end subroutine deallocateGlobalArrays
	
end module globvar

