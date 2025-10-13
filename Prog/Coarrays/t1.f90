program hello_parallel
implicit none
integer :: my_rank, num_procs

! Define a coarray variable
integer, codimension[*] :: counter

! Get the rank (ID) and number of processors
my_rank = this_image()
num_procs = num_images()

! Write out a message from each rank
write(*,*) "Hello from processor", my_rank,"out of",num_procs

! Synchronize to ensure all 'write' executions are done
sync all

! Increment the counter on the first rank
if (my_rank == 1) counter[1] = counter[1] + 1

! Again, ensure all processes are synchronized
sync all

! Print the incremented value from rank 1
if (my_rank == 1) then
write(*,*) "The counter value is now", counter[1]
endif

end program hello_parallel
