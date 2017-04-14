program nautilus_outputs

use numerical_types
use iso_fortran_env
use utilities
use nautilus_main

implicit none

! Locals
character(len=80) :: filename_output
integer :: species, output, idx_1D ! index for loops
logical :: isDefined

character(len=80) :: output_format !< format used to output data

real(double_precision), dimension(:,:,:), allocatable :: abundances_out !< abundances over time for each species. (nb_outputs, nb_species)

real(double_precision), dimension(:), allocatable :: time !< Simulation time [s]
real(double_precision), dimension(:,:), allocatable :: gas_temperature_out !< [K]
real(double_precision), dimension(:,:), allocatable :: dust_temperature_out !< [K]
real(double_precision), dimension(:,:), allocatable :: density !< [part/cm^3] 
real(double_precision), dimension(:,:), allocatable :: visual_extinction_out !< visual extinction [mag]
real(double_precision), dimension(:), allocatable :: x_rate !< X ionisation rate [s-1]

! Initialise all variables from the global_variables module. Only some of them are used here.
call initialisation()

! We calculate the total number of outputs by checking for each file if it exist or not.
nb_outputs = 0
isDefined = .true.
do while(isDefined)
  nb_outputs = nb_outputs + 1
  write(filename_output, '(a,i0.6,a)') 'abundances.',nb_outputs,'.out'
  inquire(file=filename_output, exist=isDefined)

enddo
nb_outputs = nb_outputs - 1

write(*,'(a,i0)') 'Spatial resolution: ', spatial_resolution
write(*,'(a,i0)') 'Number of time outputs: ', nb_outputs
write(*,'(a,i0)') 'Number of species: ', nb_species

! We allocate the output arrays
allocate(time(nb_outputs))
allocate(gas_temperature_out(spatial_resolution, nb_outputs))
allocate(dust_temperature_out(spatial_resolution, nb_outputs))
allocate(density(spatial_resolution, nb_outputs))
allocate(visual_extinction_out(spatial_resolution, nb_outputs))
allocate(x_rate(nb_outputs))

allocate(abundances_out(nb_outputs, nb_species, spatial_resolution))


! The next write will be written in the same line
write(*,'(a)', advance='no') 'Reading unformatted outputs...'
! We read output files
do output=1,nb_outputs
  write(filename_output, '(a,i0.6,a)') 'abundances.',output,'.out'

  open(10, file=filename_output, status='old', form='unformatted')
  read(10) time(output)
  read(10) gas_temperature_out(1:spatial_resolution, output), dust_temperature_out(1:spatial_resolution, output), &
           density(1:spatial_resolution, output), &
           visual_extinction_out(1:spatial_resolution, output), x_rate(output)
  read(10) abundances_out(output,1:nb_species, 1:spatial_resolution)
  close(10)
enddo
! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Reading unformatted outputs... Done'

! Test if the folder exists
inquire(file='ab', exist=isDefined)

! We create the folder 'ab' if he doesn't exists.
if (.not.isDefined) then
  call system("mkdir ab")
end if

! Remove all existing *.ab if needed. Will return a warning in standard output if nothing exists
call system("rm ab/*.ab")

!####################################################@@
! This part is to write one file per species, each line being one output time
!####################################################@@

if (spatial_resolution.gt.1) then
  write(filename_output, '(a,a,a)') 'ab/space.ab'
  open(10, file=filename_output)
  write(10,'(a)') '! Spatial points [AU]'
  
  do idx_1D=1, spatial_resolution
    write(10,'(es13.6e2)') grid_sample(idx_1D) / AU
  enddo
  
  close(10)
endif

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Writing *.ab ASCII files in ab/...'
! We write ASCII output file, one file per species
do species=1, nb_species
  write(filename_output, '(a,a,a)') 'ab/', trim(species_name(species)), '.ab'
  open(10, file=filename_output)
  write(10,'(a)') '! time [year]; Each column is the abundance (relative to H) [number ratio] for several spatial positions'
  
  write(output_format, *) '(es10.3e2,',spatial_resolution,'(es13.6e2," "))'
  do output=1, nb_outputs
    write(10,output_format) time(output)/YEAR, abundances_out(output, species, 1:spatial_resolution)
  enddo
  close(10)
enddo
! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Writing output files in ab/... Done'

!####################################################@@
! This part is to write one file per output time, each line being one species
!####################################################@@

!~ ! The next write will be written in the same line
!~ write(*,'(a)', advance='no') 'Writing output ASCII files...'
!~ ! We write ASCII output file, one file per output
!~ do output=1, nb_outputs
!~   write(filename_output, '(a,i0.5,a)') 'ab/abundances.', output, '.ab'
!~   open(10, file=filename_output)
!~   write(10,'(a,es10.2e2, a)') '! time =', time(output) / YEAR, ' years'
!~   write(10,'(a)') '! species name ; Abundance'
!~   do species=1, nb_species
!~     write(10,*) species_name(species), abundances_out(output, species, 1:spatial_resolution)
!~   enddo
!~   close(10)
!~ enddo
!~ ! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
!~ write(*,'(a,a)') achar(13), 'Writing output files... Done'

! Test if the folder exists
inquire(file='struct', exist=isDefined)

! We create the folder 'ab' if he doesn't exists.
if (.not.isDefined) then
  call system("mkdir struct")
end if

! Remove all existing *.ab if needed. Will return a warning in standard output if nothing exists
call system("rm struct/*.struct")

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Writing *.struct ASCII files in struct/...'
! We write ASCII output file, one file per species
do idx_1D=1, spatial_resolution
  write(filename_output, '(a,i0.5,a)') 'struct/output.', idx_1D, '.struct'
  open(10, file=filename_output)
  write(10,'(a)') '! time     ; gas temperature ; dust temperature&
                   & ; log10(H2 density)  ; visual extinction ; x ionization rate'
  write(10,'(a)') '!  [year]  ;       [K]       ;         [K]     &
                   & ; [log10(part/cm^3)] ;           [mag]   ;       [s-1] '
  
  do output=1, nb_outputs
    write(10,'(6(es10.3e2," "))') time(output)/YEAR, gas_temperature_out(idx_1D, output), dust_temperature_out(idx_1D, output), &
           density(idx_1D, output), visual_extinction_out(idx_1D, output), x_rate(output)
  enddo
  close(10)
enddo
! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Writing structure output files in struct/... Done'

end program nautilus_outputs