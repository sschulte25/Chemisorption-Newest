!******************************************************************************
! MODULE: input_output
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that contains all the routines linked to reading input files
!! or writing output files. \n\n
!!
!! Input files tends to be named *.in
!! Output files tends to be named *.out
!! Temporary files that are overwritten at each timestep are named *.tmp
!
!******************************************************************************

module input_output

use iso_fortran_env


implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Read all reactions both for gas phase and grain surface and order them
! by ITYPE
! TODO: explain the reordering process
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_reactions()

use global_variables

implicit none

! Locals
integer :: i,k,j, jk

character(len=80) :: filename !< name of the file to be read
character(len=200) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction

logical :: isDefined

! Variables for the unordered reaction file
character(len=11), dimension(MAX_COMPOUNDS,nb_gas_phase_reactions) :: SYMBOLUO1 
real(double_precision), dimension(nb_gas_phase_reactions) :: AUO1,BUO1,CUO1 
integer, dimension(nb_gas_phase_reactions) :: itypeUO1,Tmin1,Tmax1,FORMULA1,NUM1 

character (len=11), dimension(MAX_COMPOUNDS,nb_surface_reactions) :: SYMBOLUO2 
real(double_precision), dimension(nb_surface_reactions) :: AUO2,BUO2,CUO2 
integer, dimension(nb_surface_reactions) :: itypeUO2,Tmin2,Tmax2,FORMULA2,NUM2 

character (len=11), dimension(MAX_COMPOUNDS,nb_reactions) :: SYMBOLUO
real(double_precision), dimension(nb_reactions) :: AUO,BUO,CUO
integer, dimension(nb_reactions) :: itypeUO,TminUO,TmaxUO,FORMULAUO,NUMUO

character(len=80) :: line_format !< format of one line of gas_reaction.in or grain_reaction.in

! Definition of the line format, common to gas_reaction.in and grain_reaction.in
write(line_format, '(a,i0,a,i0,a)') '(', MAX_REACTANTS, 'A11,1x,', MAX_PRODUCTS, 'A11,3D11.3,23x,I3,2i7,i3,i6)'

! Reading list of reaction for gas phase
filename = 'gas_reactions.in'
inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  
  j = 0
  do
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
      comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    if (line.ne.'') then
      j = j + 1
      read(line, line_format)  (SYMBOLUO1(I,J),I=1,MAX_COMPOUNDS),AUO1(J),BUO1(J),CUO1(J), &
ITYPEUO1(J),Tmin1(j),Tmax1(j),FORMULA1(J),NUM1(J)
    
    end if
  end do
  close(10)
  
else
  write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
  call exit(1)
end if

! Reading list of reaction for grains
filename = 'grain_reactions.in'
inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  
  j = 0
  do
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
      comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    if (line.ne.'') then
      j = j + 1
      PRINT *, line
      read(line, line_format)  (SYMBOLUO2(I,J),I=1,MAX_COMPOUNDS),AUO2(J),BUO2(J),CUO2(J), &
ITYPEUO2(J),Tmin2(j),Tmax2(j),FORMULA2(J),NUM2(J)
    
    end if
  end do
  close(10)
  
else
  write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
  call exit(1)
end if

! putting everything back into the big tables
do I=1,nb_gas_phase_reactions 
  do k=1,MAX_COMPOUNDS
    SYMBOLUO(k,I)=SYMBOLUO1(k,I)
  enddo
  AUO(I)=AUO1(I)
  BUO(I)=BUO1(I)
  CUO(I)=CUO1(I)
  ITYPEUO(I)=ITYPEUO1(I)
  TminUO(I) = Tmin1(I)
  TmaxUO(I) = Tmax1(I)
  FORMULAUO(I) = FORMULA1(I)
  NUMUO(I) = NUM1(I)
enddo

do I=1,nb_surface_reactions 
  do k=1,MAX_COMPOUNDS
    SYMBOLUO(k,nb_gas_phase_reactions+I)=SYMBOLUO2(k,I)
  enddo
  AUO(nb_gas_phase_reactions+I)=AUO2(I)
  BUO(nb_gas_phase_reactions+I)=BUO2(I)
  CUO(nb_gas_phase_reactions+I)=CUO2(I)
  ITYPEUO(nb_gas_phase_reactions+I)=ITYPEUO2(I)
  TminUO(nb_gas_phase_reactions+I) = Tmin2(I)
  TmaxUO(nb_gas_phase_reactions+I) = Tmax2(I)
  FORMULAUO(nb_gas_phase_reactions+I) = FORMULA2(I)
  NUMUO(nb_gas_phase_reactions+I) = NUM2(I)
enddo

! Reorder reaction file entries with ITYPE
jk=1
do i=0,MAX_NUMBER_REACTION_TYPE
  do j=1,nb_reactions
    if (itypeuo(j).eq.i) then
      REACTION_COMPOUNDS_NAMES(:,jk)=SYMBOLUO(:,j)     
      RATE_A(jk)=AUO(j)
      RATE_B(jk)=BUO(j)
      RATE_C(jk)=CUO(j)
      REACTION_TYPE(jk)=itypeuo(j)
      REACTION_TMIN(jk) = dble(TminUO(j))
      REACTION_TMAX(jk) = dble(TmaxUO(j))
      RATE_FORMULA(jk) = FORMULAUO(J)
      REACTION_ID(jk) = NUMUO(j)
      jk=jk+1
    endif
  enddo
enddo


if (jk.ne.nb_reactions+1) then
  write(Error_unit,*) 'Some reaction was not found by the reorder process'
  write(Error_unit,*) jk,'=/',nb_reactions+1 
  call exit(4)
endif

!       replace the species names by blanks for non chemical species                                                                        
do j=1,nb_reactions-1
  do i=1,MAX_COMPOUNDS
    select case(REACTION_COMPOUNDS_NAMES(i,j))
      case ('CR', 'CRP', 'Photon')
        REACTION_COMPOUNDS_NAMES(i,j) = '           '
    end select
  enddo

enddo 


return
end subroutine read_reactions

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read all species information from gas_species.in and grain_species.in
!! files.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_species()

use global_variables
use iso_fortran_env

implicit none

! Locals
integer :: i,k

character(len=80) :: filename !< name of the file to be read
character(len=80) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction

logical :: isDefined

! Variables for the unordered reaction file
character(len=11), dimension(nb_species_for_gas) :: gas_species_label 
integer, dimension(nb_species_for_gas) :: ICG1 
integer, dimension(NB_PRIME_ELEMENTS, nb_species_for_gas) :: gas_species_composition 

character(len=11), dimension(nb_species_for_grain) :: surface_species_label 
integer, dimension(nb_species_for_grain) :: ICG2 
integer, dimension(NB_PRIME_ELEMENTS, nb_species_for_grain) :: grain_species_composition 


! Reading list of species for gas phase
filename = 'gas_species.in'
inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  
  i = 0
  do
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
      comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    if (line.ne.'') then
      i = i + 1
      read(line, '(A11,i3,13(I3))')  gas_species_label(I),ICG1(I),(gas_species_composition(K,I),K=1,NB_PRIME_ELEMENTS) 
    
    end if
  end do
  close(10)
  
else
  write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
  call exit(1)
end if


! Reading list of species for grain surface
filename = 'grain_species.in'
inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  
  i = 0
  do
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
      comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    if (line.ne.'') then
      i = i + 1
      read(line, '(A11,i3,13(I3))')  surface_species_label(I),ICG2(I),(grain_species_composition(K,I),K=1,NB_PRIME_ELEMENTS) 
    
    end if
  end do
  close(10)
  
else
  write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
  call exit(1)
end if


! putting everything back into the big tables
do I=1,nb_species_for_gas 
  species_name(I)=gas_species_label(I)
  SPECIES_CHARGE(I)=ICG1(I)
  do k=1,NB_PRIME_ELEMENTS
    species_composition(K,I)=gas_species_composition(K,I)
  enddo
enddo
PRINT *, "nb_species_for_grain=",nb_species_for_grain
do I=1,nb_species_for_grain 
  IF ( TRIM(gas_species_label(I)) .EQ. "CH3O")  THEN
    PRINT *, gas_species_label(I)
  END IF
  species_name(nb_species_for_gas+I)=surface_species_label(I)
  SPECIES_CHARGE(nb_species_for_gas+I)=ICG2(I)
  do k=1,NB_PRIME_ELEMENTS
    species_composition(K,nb_species_for_gas+I)=grain_species_composition(K,I)
  enddo
  IF ( TRIM(species_name(nb_species_for_gas+I)) .EQ. "JC2H6")  THEN
    PRINT *, nb_species_for_gas+I
  END IF
enddo

return
end subroutine read_species

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read simulation parameters from the file parameters.in
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_parameters_in()
use global_variables

implicit none

character(len=80) :: filename = 'parameters.in' !< name of the file in which parameters are stored
character(len=80) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction

logical :: isParameter, isDefined
character(len=80) :: identificator, value
!------------------------------------------------------------------------------

inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  
  do
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
      comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    call get_parameter_value(line, isParameter, identificator, value)
      
    if (isParameter) then
      select case(identificator)
      ! Solver
      case('relative_tolerance', 'RTOL') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') RELATIVE_TOLERANCE
      
      !Switches      
      case('is_3_phase')
        read(value, '(i2)') is_3_phase

      case('is_chem')
        read(value, '(i2)') is_chem

      case('is_dust_1D')
        read(value, '(i2)') is_dust_1D

      case('is_grain_reactions', 'IDUST') ! The old name is kept for compatibility reasons
        read(value, '(i2)') IS_GRAIN_REACTIONS

      case('is_h2_adhoc_form')
        read(value, '(i2)') IS_H2_ADHOC_FORM
      
      case('preliminary_test')
        read(value, '(i2)') IS_TEST
      
      case('is_absorption_h2')
        read(value, '(i2)') is_absorption_h2

      case('is_absorption_co')
        read(value, '(i2)') is_absorption_co

      case('is_absorption_n2')
        read(value, '(i2)') is_absorption_n2

      case('is_photodesorb')
      read(value, '(i2)') is_photodesorb

      case('is_crid')
      read(value, '(i2)') is_crid

      case('is_er_cir')
      read(value, '(i2)') is_er_cir

      case('grain_tunneling_diffusion', 'IGRQM') ! The old name is kept for compatibility reasons
        read(value, '(i2)') GRAIN_TUNNELING_DIFFUSION
      
      case('modify_rate_flag', 'IMODH') ! The old name is kept for compatibility reasons
        read(value, '(i2)') MODIFY_RATE_FLAG
      
      case('conservation_type', 'ICONS') ! The old name is kept for compatibility reasons
        read(value, '(i2)') CONSERVATION_TYPE
        
      case('structure_type')
        read(value, *) STRUCTURE_TYPE
      
      ! 1D definitions
      case('spatial_resolution')
        read(value, '(i4)') spatial_resolution
              
      ! Gas phase
      case('initial_gas_density', 'XNT0') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') initial_gas_density
      
      case('initial_gas_temperature', 'TEMP0') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') initial_gas_temperature
      
      case('initial_visual_extinction', 'TAU0') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') INITIAL_VISUAL_EXTINCTION
      
      case('cr_ionisation_rate', 'ZETA0') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') CR_IONISATION_RATE
      
      case('x_ionisation_rate', 'ZETAX') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') X_IONISATION_RATE
      
      case('uv_flux', 'UVGAS') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') UV_FLUX
      
      ! Grain
      case('initial_dust_temperature', 'DTEMP0') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') initial_dust_temperature
      
      case('initial_dtg_mass_ratio', 'DTOGM') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') initial_dtg_mass_ratio
      
      case('sticking_coeff_neutral', 'STICK0') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') sticking_coeff_neutral
      
      case('sticking_coeff_positive', 'STICKP') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') sticking_coeff_positive
      
      case('sticking_coeff_negative', 'STICKN') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') sticking_coeff_negative
      
      case('grain_density', 'RHOD') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') GRAIN_DENSITY
      
      case('grain_radius', 'RD') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') GRAIN_RADIUS
        
      case('diffusion_barrier_thickness', 'ACM') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') DIFFUSION_BARRIER_THICKNESS
      
      case('surface_site_density', 'SNS') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') SURFACE_SITE_DENSITY
        
      case('diff_binding_ratio_surf')
        read(value, '(e12.6)') DIFF_BINDING_RATIO_SURF

      case('diff_binding_ratio_chem')
        read(value, '(e12.6)') DIFF_BINDING_RATIO_CHEM
      
      case('diff_binding_ratio_mant')
        read(value, '(e12.6)') DIFF_BINDING_RATIO_MANT
      
      case('chemical_barrier_thickness', 'ACT') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') CHEMICAL_BARRIER_THICKNESS
      
      case('cr_peak_grain_temp', 'TSMAX') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') CR_PEAK_GRAIN_TEMP
      
      case('cr_peak_duration', 'CRT') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') CR_PEAK_DURATION
      
      case('Fe_ionisation_rate', 'CRFE') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') FE_IONISATION_RATE
      
      case('vib_to_dissip_freq_ratio', 'ARRK') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') VIB_TO_DISSIP_FREQ_RATIO

     case('ED_H2')
     read(value, '(e12.6)') ED_H2

      ! Outputs
      case('nb_outputs')
        read(value, '(i4)') NB_OUTPUTS
      
      case('start_time')
        read(value, '(e12.6)') START_TIME
      
      case('stop_time')
        read(value, '(e12.6)') STOP_TIME
      
      case('output_type')
        read(value, *) OUTPUT_TYPE
      
      ! Initial abundances
      case('minimum_initial_abundance')
        read(value, '(e12.6)') MINIMUM_INITIAL_ABUNDANCE
      
      case('is_structure_evolution')
        read(value, '(i2)') IS_STRUCTURE_EVOLUTION
      
      case('grain_temperature_type')
        read(value, *) GRAIN_TEMPERATURE_TYPE
         
      case default
        write(*,*) 'Warning: An unknown parameter has been found'
        write(*,*) "identificator='", trim(identificator), "' ; value(s)='", trim(value),"'"
      end select
    end if
  end do
  close(10)
  
else
  write (*,*) 'Warning: The file "parameters.in" does not exist. Default values have been used'
end if

START_TIME = START_TIME * YEAR
STOP_TIME = STOP_TIME * YEAR

if ((STRUCTURE_TYPE.eq.'0D').and.(spatial_resolution.ne.1)) then
  write(Error_unit,'(a,i0,a)') 'Error: In 0D, we must have one point (spatial_resolution=',spatial_resolution,')'
  call exit(22)
endif

if ((STRUCTURE_TYPE.ne.'0D').and.(spatial_resolution.eq.1)) then
  write(Error_unit,'(a,i0,a)') 'Error: In 1D, we must have more than one point (spatial_resolution=',spatial_resolution, ')'
  call exit(22)
endif

if ((STRUCTURE_TYPE.ne.'0D').and.(IS_STRUCTURE_EVOLUTION.ne.0)) then
  write(Error_unit,*) 'Error: In 1D, structure evolution is currently not supported. If is_structure_evolution=1, you must have'
  write(Error_unit,*) '       structure_type="0D"'
  
  call exit(23)
endif

return
end subroutine read_parameters_in

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief subroutine that write the simulation parameters into the file 'parameters.out'
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_parameters()
use global_variables

  implicit none
  
  character(len=80) :: filename = 'parameters.in'
  character(len=80) :: cp_command
  
  ! Copy the old file into a *.bak version just in case
  write(cp_command, '(5a)') 'cp ', trim(filename), ' ', trim(filename), '.bak'
  call system(cp_command)
  
  open(10, file=filename)
  write(10,'(a)') "!# ------------------------------------------------"
  write(10,'(a)') "!# Parameter file for various properties of the disk."
  write(10,'(a)') "!# ------------------------------------------------"
  write(10,'(a)') "!# blanck line or with spaces will be skipped."
  write(10,'(a)') "!# In fact, the only lines that matter are non commented lines with a"
  write(10,'(a)') "!# '=' character to distinguish the identificator and the value(s)"
  write(10,'(a)') "!# (each value must be separated with at least one space."
  write(10,'(a)') "!# Line must not be longer than 80 character, but comments can be far"
  write(10,'(a)') "!# bigger than that, even on line with a parameter to read."
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*   Switch 2/3 phase model  *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,i0,a)') 'is_3_phase = ', is_3_phase, ' ! 0: 2 phase, 1: 3 phase'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*          Switches         *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,i0,a)') 'preliminary_test = ', IS_TEST, ' ! Will or will not do comprehensive tests &
  &before the simulation. Switch it off when lauching thousands or simulations'
  write(10,'(a,i0,a)') 'is_structure_evolution = ', IS_STRUCTURE_EVOLUTION, ' ! If 1, physical structure properties evolve with &
                        &time, values come from structure_evolution.dat file that must exists'
  write(10,'(a,a,a)') 'grain_temperature_type = ', trim(GRAIN_TEMPERATURE_TYPE), ' ! fixed, gas, table_evolv, table_1D or computed'
  write(10,'(a)') '! fixed: Tgrain = initial_dust_temperature. ;'
  write(10,'(a)') '! gas: Tgrain = Tgas ; '
  write(10,'(a)') '! table_evolv: Tgrain is interpolated from structure_evolution.dat data file (5th optional column) ; '
  write(10,'(a)') '! table_1D: Tgrain is read in the 1D_static. dat file (5th column) ; '
  write(10,'(a)') '! computed: calculated from uv_flux and visual extinction by radiative equilibrium'
  write(10,'(a,i0,a)') 'is_dust_1D = ', is_dust_1D, ' ! Reading the grain abundance and the NH/AV factor in the 1D_static.dat file & 
                      &(mostly for disks)'
  write(10,'(a,i0,a)') 'is_chem = ', is_chem, ' ! Switches Chemisorption on and off'
  write(10,'(a,i0,a)') 'is_grain_reactions = ', IS_GRAIN_REACTIONS, ' ! Accretion, grain surface reactions'
  write(10,'(a,i0,a)') 'is_h2_adhoc_form = ', IS_H2_ADHOC_FORM, ' ! Ad hoc formation of H2 on grain surfaces (1=activated)'
  write(10,'(a,i0,a)') 'is_absorption_h2 = ', is_absorption_h2, ' ! H2 self-shielding from Lee & Herbst (1996) (1=activated)'
  write(10,'(a,i0,a)') 'is_absorption_co = ', is_absorption_co, ' ! CO self-shielding. (1: Lee & Herbst (1996), &
                                                                  2: Visser et al. (2009)'
  write(10,'(a,i0,a)') 'is_absorption_n2 = ', is_absorption_n2, ' ! N2 self-shielding from Li et al. (2013) (1=activated)'
  write(10,'(a,i0,a)') 'is_photodesorb = ', is_photodesorb, &
' ! Switch to turn on the photodesorption of ices (default yield is 1e-4)'
  write(10,'(a,i0,a)') 'is_crid = ', is_crid, &
' ! Switch to turn on the CRID (cosmic rays induced diffusion) mechanism'
  write(10,'(a,i0,a)') 'is_er_cir = ', is_er_cir, &
' ! Switch to turn on Eley-Rideal and Complex Induced Reaction mechanisms (default=0: desactivated)'
  write(10,'(a,i0,a)') 'grain_tunneling_diffusion = ', GRAIN_TUNNELING_DIFFUSION, &
  ' ! 0=thermal; For H,H2: 1=QM1; 2=QM2; 3=choose fastest'
  write(10,'(a,i0,a)') 'modify_rate_flag = ', MODIFY_RATE_FLAG, ' ! 1=modify H; 2=modify H,H2, 3=modify all, -1=H+H only'
  write(10,'(a,i0,a)') 'conservation_type = ', CONSERVATION_TYPE, ' ! 0=only e- conserved; 1=elem #1 conserved, 2=elem #1 & #2, etc'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*      1D and diffusion     *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!(diffusion is for species, not the structure)"
  write(10,'(a)') ""
  write(10,'(a,a,a)') 'structure_type = ', trim(STRUCTURE_TYPE), ' ! 0D, 1D_diff, 1D_no_diff'
  write(10,'(a,i0,a)') 'spatial_resolution = ', spatial_resolution, &
    ' ! If 1, we are in 0D, else, we are in 1D, with diffusion between gas boxes. Number of lines in 1D.'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*    Gas phase parameters   *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'initial_gas_density = ', initial_gas_density, ' ! initial gas density [part/cm-3]'
  write(10,'(a,es10.3e2,a)') 'initial_gas_temperature = ', initial_gas_temperature, ' ! initial gas temperature [K]'
  write(10,'(a,es10.3e2,a)') 'initial_visual_extinction = ', INITIAL_VISUAL_EXTINCTION, ' ! initial visual extinction'
  write(10,'(a,es10.3e2,a)') 'cr_ionisation_rate = ', CR_IONISATION_RATE, ' ! cosmic ray ionisation rate [s-1] (standard=1.3e-17)'
  write(10,'(a,es10.3e2,a)') 'x_ionisation_rate = ', X_IONISATION_RATE, ' ! Ionisation rate due to X-rays [s-1]'
  write(10,'(a,es10.3e2,a)') 'uv_flux = ', UV_FLUX, ' ! Scale factor for the UV flux, in unit of the reference flux (1.=nominal)'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*      Grain parameters     *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'initial_dust_temperature = ', initial_dust_temperature, ' ! initial dust temperature [K]&
                              & when grain_temperature_type=fixed'
  write(10,'(a,es10.3e2,a)') 'initial_dtg_mass_ratio = ', initial_dtg_mass_ratio, ' ! dust-to-gas ratio by mass'
  write(10,'(a,es10.3e2,a)') 'sticking_coeff_neutral = ', sticking_coeff_neutral, ' ! sticking coeff for neutral species'
  write(10,'(a,es10.3e2,a)') 'sticking_coeff_positive = ', sticking_coeff_positive, ' ! sticking coeff for positive species'
  write(10,'(a,es10.3e2,a)') 'sticking_coeff_negative = ', sticking_coeff_negative, ' ! sticking coeff for negative species'
  write(10,'(a,es10.3e2,a)') 'grain_density = ', GRAIN_DENSITY, ' ! mass density of grain material'
  write(10,'(a,es10.3e2,a)') 'grain_radius = ', grain_radius, ' ! grain radius [cm]'
  write(10,'(a,es10.3e2,a)') 'diffusion_barrier_thickness = ', DIFFUSION_BARRIER_THICKNESS, ' ! Barrier thickness [cm]'
  write(10,'(a,es10.3e2,a)') 'surface_site_density = ', SURFACE_SITE_DENSITY, ' ! site density on one grain [cm-2]'
  write(10,'(a,es10.3e2,a)') 'diff_binding_ratio_surf = ', DIFF_BINDING_RATIO_SURF, &
                             ' ! Ratio used to compute the DIFFUSION_BARRIER from the BINDING_ENERGY if not known (surface species)'
  write(10,'(a,es10.3e2,a)') 'diff_binding_ratio_chem = ', DIFF_BINDING_RATIO_CHEM, &
                             ' ! Ratio used to compute DIFFUSION_BARRIER from the BINDING_ENERGY if not known (chemisorbed species)'
  write(10,'(a,es10.3e2,a)') 'diff_binding_ratio_mant = ', DIFF_BINDING_RATIO_MANT, &
                             ' ! Ratio used to compute the DIFFUSION_BARRIER from the BINDING_ENERGY if not known  (mantle species)'
  write(10,'(a,es10.3e2,a)') 'chemical_barrier_thickness = ', CHEMICAL_BARRIER_THICKNESS, &
                             ' ! grain reaction activation energy barrier width. [cm]'
  write(10,'(a,es10.3e2,a)') 'cr_peak_grain_temp = ', CR_PEAK_GRAIN_TEMP, ' ! peak grain temperature [K] (CR heating)'
  write(10,'(a,es10.3e2,a)') 'cr_peak_duration = ', CR_PEAK_DURATION, ' ! duration [s] of peak grain temperature'
  write(10,'(a,es10.3e2,a)') 'Fe_ionisation_rate = ', FE_IONISATION_RATE, ' ! (cosmic) Fe-ion--grain encounter [s-1 grain-1] '
  write(10,'(a,es10.3e2,a)') '!! (for 0.1 micron grain) For cosmic photo desorptions, only Fe-ions are efficient to heat grains. '
  write(10,'(a,es10.3e2,a)') 'vib_to_dissip_freq_ratio = ', VIB_TO_DISSIP_FREQ_RATIO, &
                             ' ! [no unit] The ratio of the surface-molecule bond frequency to the frequency at'
  write(10,'(a)') '!! which energy is lost to the grain surface. Used for the RRK (Rice Ramsperger-Kessel) desorption mechanism'
  write(10,'(a)') '!! (see Garrod el al. 2007 for more). Assumed to be 1% by default.'
  write(10,'(a,es10.3e2,a)') 'ED_H2 = ', ED_H2, &
                            ' ! H2 binding energy over itself. Used for the desorption encounter mechanism. in K. '
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*  Integration and Outputs  *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'start_time = ', START_TIME/YEAR, ' ! [yrs] first output time'
  write(10,'(a,es10.3e2,a)') 'stop_time = ', STOP_TIME/YEAR, ' ! [yrs] last output time'
  write(10,'(a,i0,a)') 'nb_outputs = ', NB_OUTPUTS, ' ! Total number of outputs (used for linear or log spaced outputs)'
  write(10,'(a,a,a)') 'output_type = ', trim(OUTPUT_TYPE), ' ! linear, log'
  write(10, '(a)') '! linear: Output times are linearly spaced'
  write(10, '(a)') '! log   : Outputs times are log-spaced'
  write(10,'(a,es10.3e2, a)') 'relative_tolerance = ',RELATIVE_TOLERANCE, ' ! Relative tolerance of the solver'
  write(10,'(a,es10.3e2,a)') 'minimum_initial_abundance = ', MINIMUM_INITIAL_ABUNDANCE, ' ! default minimum initial &
                             &fraction abundance'
  close(10)
  
end subroutine write_parameters

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Write information in an output file, notably the commit and branch of the compiled binary
!! Also show if the current version had uncommitted modification that can't
!! be traced.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_general_infos()
use global_variables
use git_infos

  implicit none
  
  character(len=80) :: filename = 'info.out'
  
  open(10, file=filename)
  write(10,'(a)') '!----------------------------'
  write(10,'(a)') '!     Nautilus Version      |'
  write(10,'(a)') '!----------------------------'
  write(10,'(a,a)') 'branch = ', trim(branch)
  write(10,'(a,a)') 'commit = ', trim(commit)
  write(10,'(a,a)') '!', trim(modifs)
  write(10,'(a)') ""
  write(10,'(a)') '!----------------------------'
  write(10,'(a)') '!      General infos        |'
  write(10,'(a)') '!----------------------------'
  write(10,'(a,i0)') 'Maximum number of non-zeros values in jacobian = ', nb_nonzeros_values
  close(10)
  
end subroutine write_general_infos

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write the list of species and the corresponding index in an
!! output file 'species.out'.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_species()
use global_variables

implicit none

! Locals
integer :: i

open(10, file='species.out')
! Write 'ggo_spec.d': 5 columns of numbered species=====================
write(10,'(5(I4,")",1X,A11,1X))') (I,species_name(I),I=1,nb_species)
close(10)

return
end subroutine write_species

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Write the list of elemental species with their abundances and mass.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_elemental_abundances(filename, el_abundances)

use global_variables

implicit none

! Input
character(len=*), intent(in) :: filename !< [in] the name of the output file
real(double_precision), intent(in), dimension(NB_PRIME_ELEMENTS) :: el_abundances !< [in] Elemental abundances, either initial or current ones

! Locals
integer :: i

open(10, file=filename)
write(10, '(a)') '! Species name ; abundance (relative to H) ; mass (AMU)'
do i=1, NB_PRIME_ELEMENTS
  
  write(10, '(a,es10.4e2, f8.3)') species_name(PRIME_ELEMENT_IDX(i)), el_abundances(i), elemental_mass(i)
  
enddo
close(10)

return
end subroutine write_elemental_abundances


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief subroutine that try to split the line in two part, given a 
!! separator value (set in parameter of the subroutine)
!
!> @warning The first character of the parameter value MUST NOT be a string. All spaces will be truncated
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_parameter_value(line, isParameter, id, value)

  implicit none
  
  ! Input
  character(len=80), intent(in) :: line !< [in] Input line in which we want to retrieve a parameter and its value
  
  ! Output
  logical, intent(out) :: isParameter !< [out] a boolean to say whether or not there is a parameter on this line. 
!! i.e if there is an occurence of the separator in the input line
  character(len=80), intent(out) :: id !< [out] the name of the parameter
  character(len=80), intent(out) :: value !< [out] a string that contains the value(s) associated with the parameter name. 
!!         Note that a special attention is given to the fact that the first character of 'value' must NOT be a 'space'
  
  ! Local
  character(len=1), parameter :: SEP = '=' ! the separator of a parameter line
  character(len=1) :: first_character
  integer :: id_first_char
  integer :: sep_position ! an integer to get the position of the separator

  !------------------------------------------------------------------------------

  sep_position = index(line, SEP)
  
  if (sep_position.ne.0) then
    isParameter = .true.
    id = line(1:sep_position-1)
    
    id_first_char = sep_position +1
    first_character = line(id_first_char:id_first_char)
    do while (first_character.eq.' ')
      id_first_char = id_first_char +1
      first_character = line(id_first_char:id_first_char)
    end do
    value = line(id_first_char:)
  else
    isParameter = .false.
  end if

end subroutine get_parameter_value

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read abundances from abundances.in file. All abundances not defined
!! here will have the default value MINIMUM_INITIAL_ABUNDANCE
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_abundances()
! Writes 1D outputs
use global_variables

implicit none

! Locals
character(len=80), parameter :: filename='abundances.in'
integer :: i, j

real(double_precision), allocatable, dimension(:) :: temp_abundances
character(len=11), allocatable, dimension(:) :: temp_names
integer :: nb_lines

character(len=80) :: line
character(len=1), parameter :: comment_character = '!' ! character that will indicate that the rest of the line is a comment
integer :: comment_position ! the index of the comment character on the line. if zero, there is none on the current string
integer :: error ! to store the state of a read instruction

logical :: isParameter, isDefined
character(len=80) :: identificator, value

call get_linenumber(filename, nb_lines)

allocate(temp_abundances(nb_lines))
allocate(temp_names(nb_lines))

temp_abundances(1:nb_lines) = 0.d0
temp_names(1:nb_lines) = ''

  !------------------------------------------------------------------------------
  
inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  i = 1
  do 
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
      comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    call get_parameter_value(line, isParameter, identificator, value)
      
    if (isParameter) then
      read(value, '(e12.6)') temp_abundances(I)
      read(identificator, *) temp_names(I)
      i = i + 1
    end if
  enddo
  close(10)
endif

! We check if all species in abundances.in exists in the simulation
do j=1,nb_lines
  error = 1
  do i=1,nb_species
    if (temp_names(j).eq.species_name(i)) then
      error = 0 ! The species exist
    endif
  enddo
  
  if (error.eq.1) then
    write(*,*) j
    write(*,*) temp_names
    write(Error_Unit,*) 'Input species "', trim(temp_names(j)), '" in "', trim(filename), '" do not match those in reaction file'
    call exit(2)
  endif
enddo

! Set initial abundances================================================
abundances(1:nb_species, 1:spatial_resolution) = MINIMUM_INITIAL_ABUNDANCE

! Initial abundance for one species is assumed to be the same throughout the 1D structure initially
do i=1,nb_species
  do j=1,nb_lines
    if ((species_name(i).EQ.temp_names(j)).AND.(temp_abundances(j).NE.0.D0)) then
      abundances(i,1:spatial_resolution)=temp_abundances(j)
    endif
  enddo
enddo

deallocate(temp_abundances)
deallocate(temp_names)

return
end subroutine read_abundances

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read name and mass of all 'basic' elements (basic brick for molecules
!! such as H, He and so on).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_element_in()
! Writes 1D outputs
use global_variables

implicit none

! Locals
character(len=80), parameter :: filename='element.in'
integer :: i

character(len=80) :: line
character(len=1), parameter :: comment_character = '!' ! character that will indicate that the rest of the line is a comment
integer :: comment_position ! the index of the comment character on the line. if zero, there is none on the current string
integer :: error ! to store the state of a read instruction

logical :: isDefined

character(len=80) :: filename_gas = 'gas_species.in'
character(len=80) :: filename_grain = 'grain_species.in'

integer :: nb_columns_gas
integer :: nb_columns_grain

inquire(file=filename, exist=isDefined)

call get_linenumber(filename, NB_PRIME_ELEMENTS)

! We get the number of prime elements
nb_columns_gas = get_nb_columns(filename_gas)

if ((nb_columns_gas - 2).ne.NB_PRIME_ELEMENTS) then
  write (Error_unit,'(a,i0,a,a,a,i0,a)') 'The number of prime elements is different in "element.in" (', NB_PRIME_ELEMENTS, &
  ') and "', trim(filename_gas), '" (', nb_columns_gas-2, ') .'
  call exit(6)
endif

nb_columns_grain = get_nb_columns(filename_grain)

if ((nb_columns_grain - 2).ne.NB_PRIME_ELEMENTS) then
  write (Error_unit,'(a,i0,a,a,a,i0,a)') 'The number of prime elements is different in "element.in" (', NB_PRIME_ELEMENTS, &
  ') and "', trim(filename_grain), '" (', nb_columns_grain-2, ') .'
  call exit(6)
endif

! We allocate global variables
allocate(element_name(NB_PRIME_ELEMENTS))
allocate(elemental_mass(NB_PRIME_ELEMENTS))

element_name(1:NB_PRIME_ELEMENTS) = ''
elemental_mass(1:NB_PRIME_ELEMENTS) = 0.d0

if (isDefined) then

  open(10, file=filename, status='old')
  i = 1
  do 
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
      comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    if (line.ne.'') then
      read(line, '(a, f8.3)') element_name(i), elemental_mass(i)
      i = i + 1
    endif
  enddo
  close(10)
endif

! Other operations are done once species from gas and grain are read, because we need the indexes of all prime element in those arrays

return
end subroutine read_element_in

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant & Christophe Cossou
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write all abundances for all species in an output file at each
!! output time. The total number of output files related to abundances 
!! will be equal to the number of timestep, not equally spaced in time.\n\n
!! Output filename is of the form : abundances.000001.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_current_output(index)
! Writes 1D outputs
use global_variables

implicit none

! Input
integer, intent(in) :: index !<[in] The reference index of the current output

! Locals
character(len=80) :: filename_output

write(filename_output, '(a,i0.6,a)') 'abundances.',index,'.out'


open(UNIT=35, file=filename_output, form='unformatted')

write(35) current_time
write(35) gas_temperature(1:spatial_resolution), dust_temperature(1:spatial_resolution), &
          H_number_density(1:spatial_resolution), visual_extinction(1:spatial_resolution), X_IONISATION_RATE
write(35) abundances(1:nb_species, 1:spatial_resolution)
close(35)

return
end subroutine write_current_output

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write rates of all chemical reactions for the current timestep.
!! The total number of files will be equal to the total number of timesteps, the
!! routine being called at the end of each timestep.\n\n
!! Output filename is of the form : rates.000001.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_current_rates(index)

use global_variables

implicit none

! Input
integer, intent(in) :: index !<[in] The reference index of the current output

! Locals
character(len=80) :: filename_output

write(filename_output, '(a,i0.6,a)') 'rates.',index,'.out'

open(45, file=filename_output, form='unformatted')

write(45) reaction_rates_1D 

close(45)

return 
end subroutine write_current_rates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author
!> Valentine Wakelam
!
!> @date 2014
!
! DESCRIPTION:
!> @brief Write the H2 and CO column density computed by the model - used 
!   for the self-shielding of H2 and CO from the UV photons
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine write_H2_CO_col_dens(index)

use global_variables

implicit none

! Input
integer, intent(in) :: index !<[in] The reference index of the current output

! Locals
character(len=80) :: filename_output
integer :: i

write(filename_output, '(a,i0.6,a)') 'col_dens.',index,'.out'

open(55, file=filename_output)

! Header
write(55,'(50a)') 'H column density (cm-2) H2 column density (cm-2)  CO column density (cm-2)  N2 column density (cm-2)'

do i=1,spatial_resolution
    write(55,*) NH_z(i),NH2_z(i),NCO_z(i), NN2_z(i)
enddo

close(55)

return
end subroutine write_H2_CO_col_dens


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write the total chemical composition of the actual timestep in 
!! a file whose name is given as an input parameter. This allow to use the
!! same routine to write input, temporary and output files
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_abundances(filename)

use global_variables

implicit none

! Input
character(len=*), intent(in) :: filename !<[in] the name of the output file

! Locals
integer :: i
character(len=80) :: line_format

open(13, file=filename)

! Header
write(13,'(5(a,es10.3e2),a)') '!time =', current_time, ' s ; density = ', H_number_density, &
' part/cm^3 ; temperature=', gas_temperature,' K ; visual extinction = ', visual_extinction, &
' [mag] ; CR ionisation rate = ',CR_IONISATION_RATE,' s-1'

! To adapt the format in function of the 1D number of points
write(line_format,'(a,i0,a)') '(a," = ",', spatial_resolution, '(es12.5e2))'

do i=1,nb_species
  write(13,line_format) trim(species_name(i)), abundances(i,1:spatial_resolution)
enddo

close(13)

return
end subroutine write_abundances

end module input_output
