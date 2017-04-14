module ode_solver

use shielding

implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2003
!
! DESCRIPTION: 
!> @brief Initialize reactants and products of all reactions
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine set_chemical_reactants()
use global_variables

implicit none

! Locals
integer :: I,J,L

integer :: no_species

no_species = nb_species + 1

! By default, non existing reactants (dummy species) will be assigned (nb_species+1)
REACTION_COMPOUNDS_ID(1:MAX_COMPOUNDS, 1:nb_reactions) = no_species

do I=1,nb_reactions
  do J=1,nb_species

    do L=1,MAX_COMPOUNDS
      if (REACTION_COMPOUNDS_NAMES(L,I).EQ.species_name(J)) then
        REACTION_COMPOUNDS_ID(L,I) = J
      endif
    enddo

  enddo
enddo   

return
end subroutine set_chemical_reactants

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Count the number of non-zeros elements in each line of the jacobian
!! to dimension the arrays used in ODEPACK.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine count_nonzeros()
use global_variables
implicit none

! Locals
integer :: i

! Dummy parameters for restricted call of get_jacobian
real(double_precision), dimension(nb_species) :: DUMMYPDJ, DUMMYY
integer IDUMMY
integer, parameter :: dummy_n = 3
real(double_precision), parameter :: dummy_t = 0.d0
real(double_precision), dimension(dummy_n) :: dummy_ian, dummy_jan

integer :: max_nonzeros, NUMBERJAC

! TODO comment not needed anymore maybe include it in the routine doxygen documentation
! Dimension of the work arrays for the solver 
! The number of non zero values is checked with the testjac flag
! nb_nonzeros_values should be around the largest printed value

! Forced initialisation of global variables that will be needed, especially for the 'set_constant_rates' part. We donc care about specific values, 
!! all that counts is that we can retrieve the number of non-zeros elements.

max_nonzeros = 0

dummyy(1:nb_species) = 1.d-5

call set_constant_rates()
call set_dependant_rates(dummyy)
if(is_3_phase.eq.1) call set_dependant_rates_3phase(dummyy)


do IDUMMY=1,nb_species
  call get_jacobian(n=dummy_n, t=dummy_t, y=dummyy,j=idummy,ian=dummy_ian, jan=dummy_jan, pdj=dummypdj)
    
  NUMBERJAC=0
  do i=1,nb_species
    if (dummypdj(i).ne.0.d0) then
      NUMBERJAC = NUMBERJAC + 1
    endif
  enddo

  if (NUMBERJAC.gt.max_nonzeros) then
    max_nonzeros = NUMBERJAC
  endif
  
enddo

nb_nonzeros_values = max_nonzeros

return
end subroutine count_nonzeros

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou & Franck Hersant
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Computes a column (for the J-th species) of the chemical jacobian. 
!!\n Documentation from ODEPACK:\nuser-supplied subroutine defining the
!!\n ODE system.  The system must be put in the first-order
!!\n form dy/dt = f(t,y), where f is a vector-valued function
!!\n of the scalar t and the vector y.  Subroutine F is to
!!\n compute the function f.  It is to have the form
!!\n      SUBROUTINE F (NEQ, T, Y, YDOT)
!!\n      DOUBLE PRECISION T, Y(*), YDOT(*)
!!\n where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!!\n is output.  Y and YDOT are arrays of length NEQ.
!!\n Subroutine F should not alter y(1),...,y(NEQ).
!!\n F must be declared External in the calling program.
!!\n 
!!\n Subroutine F may access user-defined quantities in
!!\n NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
!!\n (dimensioned in F) and/or Y has length exceeding NEQ(1).
!!\n See the descriptions of NEQ and Y below.
!!\n 
!!\n If quantities computed in the F routine are needed
!!\n externally to DLSODES, an extra call to F should be made
!!\n for this purpose, for consistent and accurate results.
!!\n If only the derivative dy/dt is needed, use DINTDY instead.
!
!> @warning Even if N, T, IAN, and JAN are not actually used, they are needed
!! because ODEPACK need a routine with a specific format, and specific inputs
!! and outputs
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_jacobian(N, T, Y, J, IAN, JAN, PDJ)
use global_variables
implicit none

! Inputs
integer, intent(in) :: N !<[in] number of first order ODEs.
integer, intent(in) :: J !<[in] Index representing the J-th column of the jacobian
real(double_precision), intent(in) :: T !<[in] the initial value of the independent variable t.
real(double_precision), intent(in), dimension(N) :: IAN !<[in] structure descriptor array of size N + 1.
real(double_precision), intent(in), dimension(N) :: JAN !<[in] structure descriptor array of size NNZ.
!!\n IAN and JAN together describe the sparsity
!!\n structure of the Jacobian matrix, as used by
!!\n DLSODES when MITER = 1 or 2.
!!\n JAN contains the row indices of the nonzero
!!\n locations, reading in columnwise order, and
!!\n IAN contains the starting locations in JAN of
!!\n the descriptions of columns 1,...,NEQ, in
!!\n that order, with IAN(1) = 1.  Thus for each
!!\n j = 1,...,NEQ, the row indices i of the
!!\n nonzero locations in column j are
!!\n i = JAN(k),  IAN(j) .le. k .lt. IAN(j+1).
!!\n Note that IAN(NEQ+1) = NNZ + 1.
!!\n (If MOSS = 0, IAN/JAN may differ from the
!!\n input IA/JA because of a different ordering
!!\n in each column, and added diagonal entries.)
real(double_precision), intent(in), dimension(nb_species) :: Y !< [in] abundances

! Outputs
real(double_precision), intent(out), dimension(nb_species) :: PDJ !<[out] J-th column of df/dy

! Locals
integer :: no_species

real(double_precision), dimension(nb_species+1) :: PDJ2
integer :: i
integer :: reactant1_idx, reactant2_idx, reactant3_idx, product1_idx, product2_idx, product3_idx, product4_idx, product5_idx
integer :: reaction_idx ! The index of a given reaction

! Temp values to increase speed
real(double_precision) :: H_number_density_squared ! H_number_density*H_number_density, to gain speed
real(double_precision) :: tmp_value ! To optimize speed, temporary variable is created to avoid multiple calculation of the same thing

no_species=nb_species+1 ! Index corresponding to no species (meaning that there is no 3rd reactant for instance

H_number_density_squared = actual_gas_density * actual_gas_density

PDJ2(1:nb_species+1) = 0.d0

do i=1,nb_reactions_using_species(j)
  reaction_idx = relevant_reactions(i, j) ! j being the species index, given as a parameter

  reactant1_idx = REACTION_COMPOUNDS_ID(1, reaction_idx)
  reactant2_idx = REACTION_COMPOUNDS_ID(2, reaction_idx)
  reactant3_idx = REACTION_COMPOUNDS_ID(3, reaction_idx)
  
  product1_idx = REACTION_COMPOUNDS_ID(4, reaction_idx)
  product2_idx = REACTION_COMPOUNDS_ID(5, reaction_idx)
  product3_idx = REACTION_COMPOUNDS_ID(6, reaction_idx)
  product4_idx = REACTION_COMPOUNDS_ID(7, reaction_idx)
  product5_idx = REACTION_COMPOUNDS_ID(8, reaction_idx)
  
  ! if statements are written in a specific order to increase speed. The goal is to test first the most probable event, and 
  !! then always go to 'else' statement, not to test if we have already found our case. One, then two bodies reactions are the most 
  !! abundants reactions.

  ! One reactant only
  if (reactant2_idx.eq.no_species) then
    if (reactant1_idx.eq.J) then 
      tmp_value = reaction_rates(reaction_idx)
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(product5_idx) = PDJ2(product5_idx) + tmp_value
      PDJ2(reactant1_idx) = PDJ2(reactant1_idx) - tmp_value
    endif
  
  ! Two bodies reaction
  else if (reactant3_idx.eq.no_species) then
    if (reactant1_idx.eq.J) then 
      tmp_value = reaction_rates(reaction_idx) * Y(reactant2_idx) * actual_gas_density
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(product5_idx) = PDJ2(product5_idx) + tmp_value
      PDJ2(reactant1_idx) = PDJ2(reactant1_idx) - tmp_value
      PDJ2(reactant2_idx) = PDJ2(reactant2_idx) - tmp_value
    endif

    if (reactant2_idx.eq.J) then 
      tmp_value = reaction_rates(reaction_idx) * Y(reactant1_idx) * actual_gas_density
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(product5_idx) = PDJ2(product5_idx) + tmp_value
      PDJ2(reactant1_idx) = PDJ2(reactant1_idx) - tmp_value
      PDJ2(reactant2_idx) = PDJ2(reactant2_idx) - tmp_value
    endif
  
  ! Three bodies reaction
  else
    if (reactant1_idx.eq.J) then 
      tmp_value = reaction_rates(reaction_idx) * Y(reactant2_idx) * Y(reactant3_idx) * H_number_density_squared
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(product5_idx) = PDJ2(product5_idx) + tmp_value
      PDJ2(reactant1_idx) = PDJ2(reactant1_idx) - tmp_value
      PDJ2(reactant2_idx) = PDJ2(reactant2_idx) - tmp_value
      PDJ2(reactant3_idx) = PDJ2(reactant3_idx) - tmp_value
    endif

    if (reactant2_idx.eq.J) then 
      tmp_value = reaction_rates(reaction_idx) * Y(reactant1_idx) * Y(reactant3_idx) * H_number_density_squared
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(product5_idx) = PDJ2(product5_idx) + tmp_value
      PDJ2(reactant1_idx) = PDJ2(reactant1_idx) - tmp_value
      PDJ2(reactant2_idx) = PDJ2(reactant2_idx) - tmp_value
      PDJ2(reactant3_idx) = PDJ2(reactant3_idx) - tmp_value
    endif

    if (reactant3_idx.eq.J) then 
      tmp_value = reaction_rates(reaction_idx) * Y(reactant1_idx) * Y(reactant2_idx) * H_number_density_squared
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(product5_idx) = PDJ2(product5_idx) + tmp_value
      PDJ2(reactant1_idx) = PDJ2(reactant1_idx) - tmp_value
      PDJ2(reactant2_idx) = PDJ2(reactant2_idx) - tmp_value
      PDJ2(reactant3_idx) = PDJ2(reactant3_idx) - tmp_value
    endif

  endif

enddo

PDJ(1:nb_species)=PDJ2(1:nb_species)

return
end subroutine get_jacobian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read the full list of reactions and retrieve, for each species
!! the list of reactions involving it. This will be used in get_jacobian 
!! to increase speed. 
!!\n\n
!! max_reactions_same_species is set here. nb_reactions_using_species and relevant_reactions array are set here.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine init_relevant_reactions()
use global_variables
implicit none

! Locals
integer, dimension(nb_reactions, nb_species+1) :: is_species_used ! For each species, tell which reactions use it or not 
!! (0 if not used, 1 if used at least once)

integer :: reaction, species, idx
integer :: reactant1_idx, reactant2_idx, reactant3_idx

is_species_used(1:nb_reactions, 1:nb_species+1) = 0

do reaction=1,nb_reactions

  reactant1_idx = REACTION_COMPOUNDS_ID(1, reaction)
  reactant2_idx = REACTION_COMPOUNDS_ID(2, reaction)
  reactant3_idx = REACTION_COMPOUNDS_ID(3, reaction)
  
  is_species_used(reaction, reactant1_idx) = 1
  is_species_used(reaction, reactant2_idx) = 1
  is_species_used(reaction, reactant3_idx) = 1

enddo

! We get the total number of reactions in which each species can be involved
! We skip the 'nb_species+1' species that is only a fake species for "no reactant"
nb_reactions_using_species(1:nb_species) = sum(is_species_used(1:nb_reactions, 1:nb_species), 1)

! What is the maximum number of reactions involving one particular species?
max_reactions_same_species = maxval(nb_reactions_using_species(1:nb_species))

allocate(relevant_reactions(max_reactions_same_species, nb_species))

relevant_reactions(1:max_reactions_same_species, 1:nb_species) = 0 ! For the extra elements (because not all species 
!! will have 'max_reactions' reactions involving it).

! For each species, we get the references of reactions that have it as a reactant. The number of reactions is different for each species 
!! Thus, at least one species will have a full line of meaningfull indexes. The other will have the rest of their line completed by zeros.
do species=1,nb_species
  idx = 1
  do reaction=1, nb_reactions
    if (is_species_used(reaction, species).eq.1) then
      relevant_reactions(idx, species) = reaction
      idx = idx + 1
    endif
  enddo
enddo


return
end subroutine init_relevant_reactions

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Initialize working arrays used for ODEPACK. We must know 
!! in advance the maximum number of non-zeros elements. 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine initialize_work_arrays()
use global_variables
implicit none

lrw = 20 + 3 * nb_nonzeros_values*nb_species + 21 * nb_species
liw = 31 + 3 * nb_nonzeros_values*nb_species + 21 * nb_species

allocate(iwork(liw))
allocate(rwork(lrw))

iwork(1:liw) = 0
rwork(1:lrw) = 0.d0

return
end subroutine initialize_work_arrays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2008
!
! DESCRIPTION: 
!> @brief Calculates the position of non-zero values in the jacobian. 
!! Set the global variables iwork and rwork
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine set_work_arrays(Y)
use global_variables
implicit none

! Inputs
real(double_precision), intent(in), dimension(nb_species) :: Y !< [in] abundances

! Locals
integer :: i,j,k
real(double_precision), dimension(nb_species) :: PDJ
integer :: NNZ

! Dummy parameters for restricted call of get_jacobian
integer, parameter :: dummy_n = 3
real(double_precision), parameter :: dummy_t = 0.d0
real(double_precision), dimension(dummy_n) :: dummy_ian, dummy_jan


! For IA and JA
integer, dimension(nb_species+1) :: IA !< For each species, the number (minus 1) of non-zeros values for the previous species
integer, dimension(liw) :: JA !< List the non-zeros values. For each non-zeros values, tell to what species it correspond 

call set_constant_rates()
call set_dependant_rates(Y)
if(is_3_phase.eq.1) call set_dependant_rates_3phase(Y)


! Initialize work arrays

iwork(1:liw) = 0
rwork(1:lrw) = 0.d0
IWORK(5) = 5
RWORK(6) = 3.154D14
IWORK(6) = 10000
IWORK(7) = 2

if (.not.(first_step_done)) then
  IWORK(6)=2000
endif

k=1

do j=1,nb_species
  call get_jacobian(n=dummy_n, t=dummy_t, y=Y,j=J,ian=dummy_ian, jan=dummy_jan, pdj=PDJ)

  IA(j)=k

  do i=1,nb_species
    if (abs(PDJ(i)).gt.1.d-99) then
      JA(k)=i
      k=k+1
    endif
  enddo

enddo

IA(nb_species+1)=k

NNZ=IA(nb_species+1)-1
iwork(30+1:30+nb_species+1)=IA(1:nb_species+1)
iwork(31+nb_species+1:31+nb_species+NNZ)=JA(1:NNZ)

return
end subroutine set_work_arrays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2008
!
! DESCRIPTION: 
!> @brief Computes the chemical evolution. In particular, calculate all the
!! derivatives for abundances.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_temporal_derivatives(N,T,Y,YDOT)
use global_variables

implicit none

! Inputs
integer, intent(in) :: N !<[in] number of first order ODEs.
real(double_precision), intent(in), dimension(nb_species) :: Y !< [in] abundances
real(double_precision), intent(in) :: T !<[in] Not used by the code, but needed for the ODEPACK call format expected for FCHEM in dlsodes

! Outputs
real(double_precision), intent(out), dimension(nb_species) :: YDOT !<[out] derivative of the abundances

! Locals
integer :: no_species
real(double_precision), dimension(nb_species+1) :: YD2
real(double_precision), dimension(nb_species+1) :: YDTMP1,YDTMP2
integer :: i
integer :: reactant1_idx, reactant2_idx, reactant3_idx, product1_idx, product2_idx, product3_idx, product4_idx, product5_idx
real(double_precision) :: rate

call set_dependant_rates(y)


no_species = nb_species + 1

ydot(1:nb_species) = 0.d0
yd2(1:nb_species) = 0.d0
YDTMP1(1:nb_species) = 0.d0
YDTMP2(1:nb_species) = 0.d0

! The differential equations are calculated in a loop here
do I=1,nb_reactions

  reactant1_idx = REACTION_COMPOUNDS_ID(1, i)
  reactant2_idx = REACTION_COMPOUNDS_ID(2, i)
  reactant3_idx = REACTION_COMPOUNDS_ID(3, i)

  product1_idx = REACTION_COMPOUNDS_ID(4, i)
  product2_idx = REACTION_COMPOUNDS_ID(5, i)
  product3_idx = REACTION_COMPOUNDS_ID(6, i)
  product4_idx = REACTION_COMPOUNDS_ID(7, i)
  product5_idx = REACTION_COMPOUNDS_ID(8, i)

  ! One reactant only
  if (reactant2_idx.eq.no_species) then
    RATE = reaction_rates(I) * Y(reactant1_idx)  
  else
    if (reactant3_idx.eq.no_species) then
      ! Two bodies reactions
      RATE = reaction_rates(I) * Y(reactant1_idx) * Y(reactant2_idx) * actual_gas_density
    else 
      ! Three bodies reactions
      RATE = reaction_rates(I)*Y(reactant1_idx) * Y(reactant2_idx) * Y(reactant3_idx) * actual_gas_density * actual_gas_density
    endif
  endif

  ! Used to compute the net rate of change in the total surface material
  IF ((REACTION_TYPE(i).ne.40).and.(REACTION_TYPE(i).ne.41)) THEN
     YD2(product1_idx) = YD2(product1_idx) + RATE
     YD2(product2_idx) = YD2(product2_idx) + RATE
     YD2(product3_idx) = YD2(product3_idx) + RATE
     YD2(product4_idx) = YD2(product4_idx) + RATE
     YD2(product5_idx) = YD2(product5_idx) + RATE

     YD2(reactant1_idx) = YD2(reactant1_idx) - RATE
     YD2(reactant2_idx) = YD2(reactant2_idx) - RATE
     YD2(reactant3_idx) = YD2(reactant3_idx) - RATE

     YDTMP1(product1_idx) = YDTMP1(product1_idx) + RATE
     YDTMP1(product2_idx) = YDTMP1(product2_idx) + RATE
     YDTMP1(product3_idx) = YDTMP1(product3_idx) + RATE
     YDTMP1(product4_idx) = YDTMP1(product4_idx) + RATE
     YDTMP1(product5_idx) = YDTMP1(product5_idx) + RATE

     YDTMP2(reactant1_idx) = YDTMP2(reactant1_idx) - RATE
     YDTMP2(reactant2_idx) = YDTMP2(reactant2_idx) - RATE
     YDTMP2(reactant3_idx) = YDTMP2(reactant3_idx) - RATE

  ENDIF

enddo

IF(is_3_phase.eq.1) then
! Compute the net rate of change in the total surface material
rate_tot_acc = 0.0D+00
rate_tot_des = 0.0D+00
DO i = nb_gaseous_species+1,nb_species
   IF((species_name(i)(1:1).ne.'K').and.(species_name(i)(1:1).ne.'J').and.&
      (species_name(i)(1:1).ne.'B')) THEN
      PRINT*, "Problem in the computation of rate_tot in get_temporal_derivatives..."
      stop
   ENDIF
   IF((species_name(i)(1:1).eq.'J').or.(species_name(i)(1:1).eq.'B')) THEN
       rate_tot_acc = rate_tot_acc + YDTMP1(i)
       rate_tot_des = rate_tot_des + YDTMP2(i)
   ENDIF
ENDDO

rate_tot = rate_tot_acc + rate_tot_des

call set_dependant_rates_3phase(Y)

do i=type_id_start(40),type_id_stop(40)
  reactant1_idx = REACTION_COMPOUNDS_ID(1, i)
  reactant2_idx = REACTION_COMPOUNDS_ID(2, i)
  reactant3_idx = REACTION_COMPOUNDS_ID(3, i)

  product1_idx = REACTION_COMPOUNDS_ID(4, i)
  product2_idx = REACTION_COMPOUNDS_ID(5, i)
  product3_idx = REACTION_COMPOUNDS_ID(6, i)
  product4_idx = REACTION_COMPOUNDS_ID(7, i)
  product5_idx = REACTION_COMPOUNDS_ID(8, i)

  RATE = reaction_rates(I) * Y(reactant1_idx)

  YD2(product1_idx) = YD2(product1_idx) + RATE
  YD2(product2_idx) = YD2(product2_idx) + RATE
  YD2(product3_idx) = YD2(product3_idx) + RATE
  YD2(product4_idx) = YD2(product4_idx) + RATE
  YD2(product5_idx) = YD2(product5_idx) + RATE

  YD2(reactant1_idx) = YD2(reactant1_idx) - RATE
  YD2(reactant2_idx) = YD2(reactant2_idx) - RATE
  YD2(reactant3_idx) = YD2(reactant3_idx) - RATE

enddo

do i=type_id_start(41),type_id_stop(41)
  reactant1_idx = REACTION_COMPOUNDS_ID(1, i)
  reactant2_idx = REACTION_COMPOUNDS_ID(2, i)
  reactant3_idx = REACTION_COMPOUNDS_ID(3, i)

  product1_idx = REACTION_COMPOUNDS_ID(4, i)
  product2_idx = REACTION_COMPOUNDS_ID(5, i)
  product3_idx = REACTION_COMPOUNDS_ID(6, i)
  product4_idx = REACTION_COMPOUNDS_ID(7, i)
  product5_idx = REACTION_COMPOUNDS_ID(8, i)

  RATE = reaction_rates(I) * Y(reactant1_idx)

  YD2(product1_idx) = YD2(product1_idx) + RATE
  YD2(product2_idx) = YD2(product2_idx) + RATE
  YD2(product3_idx) = YD2(product3_idx) + RATE
  YD2(product4_idx) = YD2(product4_idx) + RATE
  YD2(product5_idx) = YD2(product5_idx) + RATE

  YD2(reactant1_idx) = YD2(reactant1_idx) - RATE
  YD2(reactant2_idx) = YD2(reactant2_idx) - RATE
  YD2(reactant3_idx) = YD2(reactant3_idx) - RATE

enddo
endif

YDOT(1:nb_species) = YD2(1:nb_species)

return
end subroutine get_temporal_derivatives

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Valentine Wakelam
!
!> @date 2012
!
! DESCRIPTION: 
!> @brief Constant reactions rates with respect to abundances. 
!! Reactions coefficient formally dependent on the abundances Y are 
!! computed in a companion subroutine: set_dependant_rates
!! Grain surface reactions, self shielding, etc...
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine set_constant_rates()

  use global_variables 
  
  implicit none

  ! Locals
  real(double_precision) :: T300, TI, TSQ
  integer :: nsta, nfin
  integer :: k, j, w, m, n
  integer, dimension(10) :: indice
  real(double_precision), dimension(10) :: distmin, distmax
  real(double_precision) :: activ, activ2, p1, mu
  real(double_precision) :: barr, barr2
  real(double_precision) :: tunneling_rate
  real(double_precision) :: kappa

  T300 = actual_gas_temp / 300.d0
  TI = 1.0d00 / actual_gas_temp
  TSQ = SQRT(actual_gas_temp)
  EVAPORATION_RATES_H2=0.D0

  !------------------------------------------------------------
  !
  !      --- GAS PHASE REACTIONS
  !
  !------------------------------------------------------------


  !-------------------------------------------------------------
  !  (0) Gas phase reactions with grains (s-1)
  !-------------------------------------------------------------

  do J=type_id_start(0),type_id_stop(0)
    reaction_rates(J) = RATE_A(J) * (T300**RATE_B(J))
  enddo


  !-------------------------------------------------------------
  !  (10-11) Ad-hoc formation of h2 on grains when
  !          is_grain_reactions = 0 is_h2_adhoc_form = 1
  !
  !  Comments:
  !     - In the case is_grain_reactions = 0 we still need the formation
  !     of H on the grains this is done with the XH species
  !     through those reactions type
  !     - When the adhoc h2 formation is activated 1/2 of
  !     the adsorbed H are availables for grain reactions
  !     (other than h2 formation) and 1/2 for the formation of h2
  !-------------------------------------------------------------

  if ((IS_GRAIN_REACTIONS.eq.0).OR.(IS_H2_ADHOC_FORM.eq.1)) then
  
    do J=type_id_start(10),type_id_stop(10)
      reaction_rates(J) = RATE_A(J) * 1.186D7 * exp(225.D0 / actual_gas_temp)**(-1) * GTODN / actual_gas_density
      ! when the adhoc h2 formation is activated 1/2 of the adsorbed H are availables for grain reactions (other than h2 formation) and
      ! 1/2 for the formation of h2
      IF(IS_H2_ADHOC_FORM.eq.1) reaction_rates(J) = 0.5D+00 * reaction_rates(J)
    enddo
    
    do J=type_id_start(11),type_id_stop(11)       
      reaction_rates(J) = RATE_A(J) * (T300**RATE_B(J)) * actual_gas_density / GTODN
    enddo       
    
  endif


  !-------------------------------------------------------------
  !  (1) Photodissoc/ionisation with cosmic and X rays (s-1)
  !-------------------------------------------------------------

  do J=type_id_start(1),type_id_stop(1)
    reaction_rates(J) = RATE_A(J) * (CR_IONISATION_RATE + X_IONISATION_RATE)
  enddo


  !-------------------------------------------------------------
  !  (4-8) Bimolecular gas phase reactions
  !        Several possible formula:
  !         - Kooij (modified Arhenius) - formula = 3
  !         - Ionpol1 - formula = 4
  !         - Ionpol2 - formula = 5
  !-------------------------------------------------------------

  W=1
  NSTA=0
  NFIN=0
  distmin(:)=9999.d0
  distmax(:)=9999.d0
  do J=4,8
    if ((type_id_start(J).NE.0).AND.(NSTA.EQ.0)) NSTA=J
    if (type_id_stop(J).NE.0) NFIN=J
  enddo
  do J=type_id_start(NSTA),type_id_stop(NFIN)

    !---------------  Kooij (modified Arhenius) - formula = 3
    if (RATE_FORMULA(J).eq.3) then

      reaction_rates(J) = RATE_A(J) * (T300**RATE_B(J)) * EXP(-RATE_C(J) / actual_gas_temp)

      ! Check for temperature bounderies
      if (actual_gas_temp.LT.REACTION_TMIN(J)) then
        reaction_rates(J) = RATE_A(J) * ((REACTION_TMIN(J) / 300.D0)**RATE_B(J)) * EXP(-RATE_C(J) / REACTION_TMIN(J))
      endif
      if (actual_gas_temp.GT.REACTION_TMAX(J)) then
        reaction_rates(J) = RATE_A(J) * ((REACTION_TMAX(J) / 300.D0)**RATE_B(J)) * EXP(-RATE_C(J) / REACTION_TMAX(J))
      endif

      ! Check for the presence of several rate coefficients present in the network for the
      ! the same reaction
      if (REACTION_ID(J+1).EQ.REACTION_ID(J)) then
        INDICE(W) = J
        distmin(w) = REACTION_TMIN(j) - actual_gas_temp
        distmax(w) = actual_gas_temp - REACTION_TMAX(j)
        W = W + 1
      endif

      if ((REACTION_ID(J+1).NE.REACTION_ID(J)).AND.(W.NE.1)) then

        INDICE(W) = J
        distmin(w) = REACTION_TMIN(j) - actual_gas_temp
        distmax(w) = actual_gas_temp - REACTION_TMAX(j)

        do M=1,W
          N = INDICE(M)
          !if(IT==1) write(*,*) N,M, REACTION_COMPOUNDS_NAMES(:,N), REACTION_TMIN(N), REACTION_TMAX(N),distmin(M),distmax(M)
          if (actual_gas_temp.LT.REACTION_TMIN(N)) reaction_rates(N) = 0.d0
          if (actual_gas_temp.GT.REACTION_TMAX(N)) reaction_rates(N) = 0.d0
        enddo

        if (maxval(reaction_rates(indice(1:w))).lt.1.d-99) then

          if (minval(abs(distmin)).lt.minval(abs(distmax))) then
            N=indice(minloc(abs(distmin),dim=1))
            reaction_rates(N)=RATE_A(N)*((REACTION_TMIN(N)/300.D0)**RATE_B(N))*EXP(-RATE_C(N)/REACTION_TMIN(N))
          else
            N=indice(minloc(abs(distmax),dim=1))
            reaction_rates(N)=RATE_A(N)*((REACTION_TMAX(N)/300.D0)**RATE_B(N))*EXP(-RATE_C(N)/REACTION_TMAX(N))
          endif
        endif

        W=1
        INDICE(:)=0
        distmin(:)=9999.d0
        distmax(:)=9999.d0
      endif
    endif

    !---------------  Ionpol1 - formula = 4
    if (RATE_FORMULA(J).EQ.4) then
      reaction_rates(J)=RATE_A(J)*RATE_B(J)*(0.62d0+0.4767d0*RATE_C(J)*((300.D0/actual_gas_temp)**0.5))

      ! Check for temperature bounderies
      if (actual_gas_temp.LT.REACTION_TMIN(J)) then
        reaction_rates(J)=RATE_A(J)*RATE_B(J)*(0.62d0+0.4767d0*RATE_C(J)*((300.D0/REACTION_TMIN(J))**0.5))
      endif
      if (actual_gas_temp.GT.REACTION_TMAX(J))  then
        reaction_rates(J)=RATE_A(J)*RATE_B(J)*(0.62d0+0.4767d0*RATE_C(J)*((300.D0/REACTION_TMAX(J))**0.5))
      endif

      ! Check for the presence of several rate coefficients present in the network for the
      ! the same reaction
      if (REACTION_ID(J+1).EQ.REACTION_ID(J)) then
        INDICE(W)=J
        distmin(w) = REACTION_TMIN(j) - actual_gas_temp
        distmax(w) = actual_gas_temp - REACTION_TMAX(j)
        W = W + 1
      endif

      if ((REACTION_ID(J+1).NE.REACTION_ID(J)).AND.(W.NE.1)) then

        INDICE(W)=J
        distmin(w) = REACTION_TMIN(j) - actual_gas_temp
        distmax(w) = actual_gas_temp - REACTION_TMAX(j)

        do M=1,W
          N=INDICE(M)
          if (actual_gas_temp.LT.REACTION_TMIN(N)) reaction_rates(N) = 0.d0
          if (actual_gas_temp.GT.REACTION_TMAX(N)) reaction_rates(N) = 0.d0
        enddo

        if (maxval(reaction_rates(indice(1:w))).lt.1.d-99) then
          if (minval(abs(distmin)).lt.minval(abs(distmax))) then
            N=indice(minloc(abs(distmin),dim=1))
            reaction_rates(N)=RATE_A(N)*RATE_B(N)*(0.62d0+0.4767d0*RATE_C(N)*((300.D0/REACTION_TMIN(N))**0.5))
          else
            N=indice(minloc(abs(distmax),dim=1))
            reaction_rates(N)=RATE_A(N)*RATE_B(N)*(0.62d0+0.4767d0*RATE_C(N)*((300.D0/REACTION_TMAX(N))**0.5))
          endif
        endif

        W=1
        INDICE(:)=0
        distmin(:)=9999.d0
        distmax(:)=9999.d0
      endif
    endif

    !---------------  Ionpol2 - formula = 5
    if (RATE_FORMULA(J).EQ.5) then
      
      ! Check for temperature boundaries and apply the formula in consequence
      if (actual_gas_temp.LT.REACTION_TMIN(J)) then
        reaction_rates(J) = RATE_A(J)*RATE_B(J)*((1.d0+0.0967*RATE_C(J)*(300.D0/REACTION_TMIN(J))**0.5)+&
                           (RATE_C(J)**2*300.d0/(10.526*REACTION_TMIN(J))))
      else if (actual_gas_temp.GT.REACTION_TMAX(J)) then
        reaction_rates(J) = RATE_A(J)*RATE_B(J)*((1.d0+0.0967*RATE_C(J)*(300.D0/REACTION_TMAX(J))**0.5)+&
                           (RATE_C(J)**2*300.d0/(10.526*REACTION_TMAX(J))))
      else
        reaction_rates(J) = RATE_A(J)*RATE_B(J)*((1.d0+0.0967*RATE_C(J)*(300.D0/actual_gas_temp)**0.5)+&
                           (RATE_C(J)**2*300.d0/(10.526*actual_gas_temp)))
      endif

      ! Check for the presence of several rate coefficients present in the network for the
      !! the same reaction
      if (REACTION_ID(J+1).EQ.REACTION_ID(J)) then
        INDICE(W)=J
        distmin(w)=REACTION_TMIN(j)-actual_gas_temp
        distmax(w)=actual_gas_temp-REACTION_TMAX(j)
        W = W + 1
      endif

      if ((REACTION_ID(J+1).NE.REACTION_ID(J)).AND.(W.NE.1)) then

        INDICE(W)=J
        distmin(w)=REACTION_TMIN(j)-actual_gas_temp
        distmax(w)=actual_gas_temp-REACTION_TMAX(j)

        do M=1,W
          N=INDICE(M)
          if (actual_gas_temp.LT.REACTION_TMIN(N)) reaction_rates(N)= 0.d0
          if (actual_gas_temp.GT.REACTION_TMAX(N)) reaction_rates(N)= 0.d0
        enddo

        if (maxval(reaction_rates(indice(1:w))).lt.1.d-99) then
          if (minval(abs(distmin)).lt.minval(abs(distmax))) then
            N=indice(minloc(abs(distmin),dim=1))
            reaction_rates(N)=RATE_A(N)*RATE_B(N)*((1.d0+0.0967*RATE_C(N)*(300.D0/REACTION_TMIN(N))**0.5) + &
                              (RATE_C(N)**2*300.d0/(10.526*REACTION_TMIN(N))))
          else
            N=indice(minloc(abs(distmax),dim=1))
            reaction_rates(N)=RATE_A(N)*RATE_B(N)*((1.d0+0.0967*RATE_C(N)*(300.D0/REACTION_TMAX(N))**0.5) + &
                              (RATE_C(N)**2*300.d0/(10.526*REACTION_TMAX(N))))
          endif
        endif

        W=1
        INDICE(:)=0
        distmin(:)=9999.d0
        distmax(:)=9999.d0
      endif
    endif
  enddo



  !------------------------------------------------------------
  !
  !      --- GRAIN REACTIONS
  !
  !------------------------------------------------------------


  if (IS_GRAIN_REACTIONS.NE.0) then


  !-------------------------------------------------------------
  !  (15) Thermal evaporation (s-1)
  !-------------------------------------------------------------

  do J=type_id_start(15),type_id_stop(15)
    reaction_rates(J)=RATE_A(J)*branching_ratio(J)*VIBRATION_FREQUENCY(reactant_1_idx(J))*EXP(-BINDING_ENERGY(reactant_1_idx(J))/&
                      actual_dust_temp)
    EVAPORATION_RATES(reactant_1_idx(J))=EVAPORATION_RATES(reactant_1_idx(J))+reaction_rates(J)
    if (species_name(reactant_1_idx(J)).EQ.YJH2) EVAPORATION_RATES_H2=RATE_A(J)*branching_ratio(J)*&
          VIBRATION_FREQUENCY(reactant_1_idx(J))*EXP(-ED_H2/actual_dust_temp)
  !If ((REACTION_COMPOUNDS_NAMES(1,J)).EQ.'G          ') then
      !PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(2,J)
 
      !PRINT *, "The reaction rate (ITYPE 15) for reaction:",J,'is=',reaction_rates(J)
      !END IF
      ! call exit()
  enddo




  !-------------------------------------------------------------
  !  (16) Cosmic-ray evaporation (s-1)
  !-------------------------------------------------------------

  do J=type_id_start(16),type_id_stop(16)
    reaction_rates(J) = RATE_A(J) * branching_ratio(J) * ((CR_IONISATION_RATE + X_IONISATION_RATE) / 1.3D-17) &
    * VIBRATION_FREQUENCY(reactant_1_idx(J)) * FE_IONISATION_RATE * CR_PEAK_DURATION &
    * EXP(-BINDING_ENERGY(reactant_1_idx(J)) / CR_PEAK_GRAIN_TEMP)

    if (grain_radius.gt.1.D-4) reaction_rates(J)=0.D0

    EVAPORATION_RATES(reactant_1_idx(J))=EVAPORATION_RATES(reactant_1_idx(J))+reaction_rates(J)
    if (species_name(reactant_1_idx(J)).EQ.YJH2) EVAPORATION_RATES_H2=EVAPORATION_RATES_H2+RATE_A(J) * branching_ratio(J) * &
        ((CR_IONISATION_RATE + X_IONISATION_RATE) / 1.3D-17) * VIBRATION_FREQUENCY(reactant_1_idx(J)) * &
        FE_IONISATION_RATE * CR_PEAK_DURATION * EXP(-ED_H2 / CR_PEAK_GRAIN_TEMP)

  enddo


  !-------------------------------------------------------------
  !  (98) Test the storage of H2S under a refractory form
  !-------------------------------------------------------------

!  if (type_id_start(98).ne.0) then
!    do j=type_id_start(98),type_id_stop(98)
!      reaction_rates(J)=RATE_A(J)*(T300**RATE_B(J))*EXP(-RATE_C(J)*TI)
!    enddo
!  endif


  !-------------------------------------------------------------
  !  Set diffusion and evaporation rates (s-1)
  !-------------------------------------------------------------

  do K=1,nb_species
     THERMAL_HOPING_RATE(K)=VIBRATION_FREQUENCY(K)*EXP(-DIFFUSION_BARRIER(K)/actual_dust_temp)/nb_sites_per_grain
     CR_HOPING_RATE(K)=VIBRATION_FREQUENCY(K)*EXP(-DIFFUSION_BARRIER(K)/CR_PEAK_GRAIN_TEMP)/nb_sites_per_grain * &
     (CR_IONISATION_RATE / 1.3D-17) * FE_IONISATION_RATE * CR_PEAK_DURATION
     if (is_crid.NE.0) THERMAL_HOPING_RATE(K)=VIBRATION_FREQUENCY(K)*EXP(-DIFFUSION_BARRIER(K)/actual_dust_temp) &
                       /nb_sites_per_grain + CR_HOPING_RATE(K)
  enddo


  !-------------------------------------------------------------
  !  (31) Complex reactions
  !       - Reactions like JC-X->JCX and JCH-X->JCHX
  !       - Choose fastest between thermal hooping and tunneling (0.4 < a < 0.5 Angstrom)
  !         a is set to 1 Angtrom for the moment (ACM in nls_control.d)
  !-------------------------------------------------------------

  IF(type_id_start(31).ne.0) THEN
     DO J=type_id_start(31),type_id_stop(31)
        tunneling_rate = TUNNELING_RATE_TYPE_2(reactant_1_idx(J))/nb_sites_per_grain
        BARR=1.0d0
        ! --------- Calculate activation energy barrier multiplier
        if(ACTIVATION_ENERGY(J).GE.1.0D-40) then
           ACTIV=ACTIVATION_ENERGY(J)/actual_dust_temp
           ! ------------ Choose fastest of classical or tunnelling
           if (ACTIV.GT.SURF_REACT_PROBA(J)) ACTIV=SURF_REACT_PROBA(J)
           BARR=EXP(-ACTIV)
        endif
        IF(tunneling_rate.GT.THERMAL_HOPING_RATE(reactant_1_idx(J))) THEN
           reaction_rates(J)=RATE_A(J)*branching_ratio(J)*tunneling_rate * BARR
        ELSE
           reaction_rates(J)=RATE_A(J)*branching_ratio(J)*THERMAL_HOPING_RATE(reactant_1_idx(J))*BARR
        ENDIF
        IF(is_er_cir.eq.0) reaction_rates(J) = 0.0d+00
     ENDDO
  ENDIF

  endif


  !-------------------------------------------------------------
  !       When dust is turned off, zero all dust rates
  !-------------------------------------------------------------

  if ((IS_GRAIN_REACTIONS.EQ.0).AND.(.not.(first_step_done))) then
    do J=type_id_start(14),type_id_stop(99)
      reaction_rates(J)=0.d0
      branching_ratio(J)=0.d0
    enddo
  endif

  return
  end subroutine set_constant_rates


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Valentine Wakelam
!
!> @date 2012
!
! DESCRIPTION: 
!> @brief Set Reactions coefficient formally dependent on the abundances Y. 
!! Grain surface reactions, self shielding, etc...
!
! EVAPORATION_RATES_TEMPO was created because this subroutine is called several times 
! during the run of one model
!
! Chemisorption added by Sean Schulte
! 2017
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine set_dependant_rates(Y)

  use global_variables
  
  implicit none
  
  ! Inputs
  real(double_precision), intent(in), dimension(nb_species) :: Y !< [in] abundances (relative to H) [number ratio]

  ! Locals
  real(double_precision) :: ACTIV !< [no unit] ratio between the energy barrier (in K) and the temperature of the grain (ACTIV)
  real(double_precision) :: BARR !< Probability to cross the barrier (interrim calculation variable) TODO unit?
  real(double_precision) :: DIFF !< TODO description/units ???
  real(double_precision) :: ACTIVCR !< [no unit] ratio between the energy barrier (in K) and the peak temperature of the grain when hit by a cosmic ray
  real(double_precision) :: BARRCR !< TODO description/units ???
  real(double_precision) :: DIFFCR !< TODO description/units ???
  real(double_precision) :: PROBH2H2 ! probability for the encounter desorption process.
  real(double_precision) :: TETABIS,TETABIS1,TETABIS2,TETABIS3
  real(double_precision) :: T300, TI, TSQ
  real(double_precision) :: YMOD1, YMOD2
  integer :: IMOD1 !< modify rate flag for reactant 1
  integer :: IMOD2 !< modify rate flag for reactant 2 
  integer :: j, l, m
  real(double_precision) :: tunneling_rate
  real(double_precision) :: LH_RT !< langmuir hinchelwood rejection term

  
  REAL(double_precision) :: ab_tot  !< Sum of all abundances on grain (surface+mantle)
  REAL(double_precision) :: ab_surf !< Sum of all abundances on grain surface
  REAL(double_precision) :: ab_chem !< Sum of all abundances on chemisorption layer
  REAL(double_precision) :: ab_mant !< Sum of all abundances on grain mantle
  real(double_precision) :: ab_lay  !< Abundance in 1 layer
  real(double_precision) :: sumlaysurf
  real(double_precision) :: sumlaymant
  
  REAL(double_precision) :: MLAY   !< Number of layers from which species can desorb
  REAL(double_precision) :: MLAY_CHEM   !< Number of layers from which species can desorb (Chemisorption)
  REAL(double_precision) :: SUMLAY !< Total number of layers on the grain surface
  REAL(double_precision) :: UVCR   !< Scaling factor for CR generated UV
                                   !! The reference used is 1.3x10^-17 s-1
  REAL(double_precision) :: cond,stick

  real(double_precision) :: abH2O   !< H2O abundance on grains surface
  real(double_precision) :: abNH3   !< NH3 abundance on grains surface
  real(double_precision) :: abCO2   !< CO2 abundance on grains surface
  real(double_precision) :: abCH4   !< CH4 abundance on grains surface
  real(double_precision) :: abCH3OH !< CH3OH abundance on grains surface
  real(double_precision) :: abCO    !< CO abundance on grains surface
  real(double_precision) :: abBCO   !< BCO abundance on grains surface
  real(double_precision) :: abBH    !< H abundance on grains surface (Chemisorption)
  real(double_precision) :: abH    !< H abundance on grains surface 
  real(double_precision) :: abBH2   !< H2 abundance on grains surface (Chemisorption)
  real(double_precision) :: abH2   !< H2 abundance on grains surface
  real(double_precision) :: abBC    !< C abundance on grains surface (Chemisorption)
  real(double_precision) :: abC    !< C abundance on grains surface
  real(double_precision) :: abBCH   !< C abundance on grains surface (Chemisorption)
  real(double_precision) :: abCH   !< C abundance on grains surface
  real(double_precision) :: abBCH2  !< C abundance on grains surface (Chemisorption)
  real(double_precision) :: abCH2  !< C abundance on grains surface 
  real(double_precision) :: abBCH3  !< C abundance on grains surface (Chemisorption)
  real(double_precision) :: abCH3  !< C abundance on grains surface
  real(double_precision) :: abBCH4  !< CH4 abundance on grains surface (Chemisorption)
  real(double_precision) :: abBC2H4 !< C2H4 abundance on grains surface (Chemisorption)
  real(double_precision) :: abBOH   !< OH abundance on grains surface (Chemisorption)
  real(double_precision) :: abBO    !< O abundance on grains surface (Chemisorption)
  real(double_precision) :: abBH2O  !< H2O abundance on grains surface (Chemisorption)
  real(double_precision) :: prob_reac, prob_deso, prob_diff

  ! For the bilinear interpolation of the shielding function of CO and N2
  real(double_precision) :: int
  real(double_precision) :: denom
  real(double_precision) :: x1, x2
  real(double_precision) :: xp, yp
  real(double_precision) :: y1, y2
  real(double_precision) :: r1, r2
  real(double_precision) :: f11, f21, f12, f22




  T300=actual_gas_temp/300.d0
  TI=1.0d00/actual_gas_temp
  TSQ=SQRT(actual_gas_temp)

  EVAPORATION_RATES_TEMPO(:) = 0.D0
  EVAPORATION_RATES_TEMPO_H2 = 0.D0

  cond=PI*grain_radius*grain_radius*SQRT(8.0d0*K_B/PI/AMU)

  do J=1,nb_species
      ACC_RATES_PREFACTOR(J)=COND*STICK_SPEC(J)/SQRT(SPECIES_MASS(J))
  ENDDO


  ab_tot  = 0.0D+00
  ab_surf = 0.0D+00
  ab_chem = 0.0D+00
  ab_mant = 0.0D+00
  abCO    = 0.0D+00
  abBCO   = 0.0D+00
  abH2O   = 0.0D+00
  abNH3   = 0.0D+00
  abCH4   = 0.0D+00
  abCH3OH  = 0.0D+00
  abBH    = 0.0D+00
  abH    = 0.0D+00
  abBH2   = 0.0D+00
  abH2   = 0.0D+00
  abBCO   = 0.0D+00
  abBH2O  = 0.0D+00
  abBC    = 0.0D+00
  abC    = 0.0D+00
  abBCH4  = 0.0D+00
  abBCH3  = 0.0D+00
  abCH3  = 0.0D+00
  abBCH2  = 0.0D+00
  abCH2  = 0.0D+00
  abBCH   = 0.0D+00
  abCH   = 0.0D+00
  abBO    = 0.0D+00
  abBOH   = 0.0D+00
  abBC2H4 = 0.0D+00

  DO J = nb_gaseous_species+1,nb_species
       ab_tot = ab_tot + Y(j)
       IF (species_name(j)(1:1).eq."J") ab_surf = ab_surf + Y(j)
       IF (species_name(j)(1:1).eq."B") ab_chem = ab_chem +  Y(j)
       IF (species_name(j)(1:1).eq."K") ab_mant = ab_mant + Y(j)
       IF(species_name(J).EQ.'JCO        ') abCO  = Y(J)
       IF(species_name(J).EQ.'BCO        ') abBCO  = Y(J)
       IF(species_name(J).EQ.'JH2O       ') abH2O = Y(J)
       IF(species_name(J).EQ.'JNH3       ') abNH3 = Y(J)
       IF(species_name(J).EQ.'JCO2       ') abCO2 = Y(J)
       IF(species_name(J).EQ.'JCH4       ') abCH4 = Y(J)
       IF(species_name(J).EQ.'JCH3OH     ') abCH3OH = Y(J)
       IF(species_name(J).EQ.'BC         ') abBC  = Y(J)
       IF(species_name(J).EQ.'JC         ') abC  = Y(J)
       IF(species_name(J).EQ.'BO         ') abBO = Y(J)
       IF(species_name(J).EQ.'BH         ') abBH  = Y(J)
       IF(species_name(J).EQ.'JH         ') abH  = Y(J)
       IF(species_name(J).EQ.'BH2        ') abBH2 = Y(J)
       IF(species_name(J).EQ.'JH2        ') abH2 = Y(J)
       IF(species_name(J).EQ.'BCH        ') abBCH = Y(J)
       IF(species_name(J).EQ.'JCH        ') abCH = Y(J)
       IF(species_name(J).EQ.'BCH2       ') abBCH2 = Y(J)
       IF(species_name(J).EQ.'JCH2       ') abCH2 = Y(J)
       IF(species_name(J).EQ.'BCH3       ') abBCH3 = Y(J)
       IF(species_name(J).EQ.'JCH3       ') abCH3 = Y(J)
       IF(species_name(J).EQ.'BCH4       ') abBCH4  = Y(J)
       IF(species_name(J).EQ.'BC2H4      ') abBC2H4 = Y(J)
       IF(species_name(J).EQ.'BOH        ') abBOH = Y(J)
       IF(species_name(J).EQ.'BH2O       ') abBH2O = Y(J)
    
  ENDDO

  ab_lay = nb_sites_per_grain/GTODN

  MLAY = 2.d0
  MLAY_CHEM = 1.d0
  SUMLAY = ab_tot*GTODN/nb_sites_per_grain
  sumlaysurf = ab_surf*GTODN/nb_sites_per_grain
  sumlaymant = ab_mant*GTODN/nb_sites_per_grain
  
  UVCR = 1.300d-17 / CR_IONISATION_RATE
 



  !------------------------------------------------------------
  !
  !      --- GAS PHASE REACTIONS
  !
  !------------------------------------------------------------


  !-------------------------------------------------------------
  !  (2) Gas phase photodissociations/ionisations by UV induced 
  !      by cosmic rays.
  !      Now it takes into account the abundance of H2 but it 
  !      still assumes that the albedo of the grains is 0.5, 
  !      hence the 0.5 in the equation. Check Wakelam et al. (2012)
  !      to understand the formula. This is now time dependent.
  !-------------------------------------------------------------

  do J=type_id_start(2),type_id_stop(2)
      reaction_rates(J) = RATE_A(J) * (CR_IONISATION_RATE + X_IONISATION_RATE)*Y(INDH2)/0.5D0
  enddo


  !-------------------------------------------------------------
  !  (3) Gas phase photodissociations/ionisations by UV of H2, CO
  !      and N2
  !      The treatment for the H2 and CO photodissociation was
  !      corrected for larger and smaller values of H2, CO column
  !      densities and Av not in the data tables.
  !-------------------------------------------------------------
  
  do J=type_id_start(3),type_id_stop(3)
    reaction_rates(J)=RATE_A(J)*EXP(-RATE_C(J)*actual_av)*UV_FLUX

    ! MODIFY THE H2, CO and N2 PHOTODISSOCIATION if flags are activated

    ! ====== H2 self-shielding
    if ((REACTION_COMPOUNDS_NAMES(1,J).EQ.YH2).and.(is_absorption_h2.EQ.1)) then
       TETABIS=1.D0


       ! ======= Linear extrapolation of the shielding factors
       do L=1,NL1-1
         if ((N1H2(L).LE.NH2).AND.(N1H2(L+1).GE.NH2)) then
           TETABIS=T1H2(L)+(NH2-N1H2(L))*(T1H2(L+1)-T1H2(L))/(N1H2(L+1)-N1H2(L))
         endif
       enddo
       if (NH2.GT.N1H2(NL1)) TETABIS = T1H2(NL1)

       reaction_rates(J)=2.54D-11*TETABIS

       reaction_rates(J)=reaction_rates(J)*UV_FLUX
    endif

    ! ====== the CO self-shielding (2 prescription)
    if (REACTION_COMPOUNDS_NAMES(1,J).EQ.YCO) then

      ! --- Self shielding of CO from Lee et al. (1996) ---
      if(is_absorption_co.eq.1) then

        TETABIS1=1.D0
        TETABIS2=1.D0
        TETABIS3=1.D0

        ! Linear extrapolation of the three shileding factors
        do L=1,NL2-1
          if ((N2CO(L).LE.NCO).AND.(N2CO(L+1).GE.NCO))  &
          TETABIS2=T2CO(L)+(NCO-N2CO(L))*(T2CO(L+1)-T2CO(L))&
          /(N2CO(L+1)-N2CO(L))
        enddo

        do L=1,NL3-1
          if ((N2H2(L).LE.NH2).AND.(N2H2(L+1).GE.NH2)) then
             TETABIS1 = T2H2(L) + (NH2 - N2H2(L)) * (T2H2(L+1) - T2H2(L)) / (N2H2(L+1) - N2H2(L))
          endif
          if ((AV2(L).LE.actual_av).AND.(AV2(L+1).GE.actual_av)) then
             TETABIS3 = T2AV(L) + (actual_av - AV2(L)) * (T2AV(L+1) - T2AV(L)) / (AV2(L+1) - AV2(L))
          endif
        enddo

        ! Saturate the rate coefficient if necessary (when density or Av are out of
        ! the shielding array, from the photodiss files)

        if (NCO.GT.N2CO(NL2)) TETABIS2 = T2CO(NL2)
        if (NH2.GT.N2H2(NL3)) TETABIS1 = T2H2(NL3)
        if (actual_av.GT.AV2(NL3))  TETABIS3 = T2AV(NL3)

        reaction_rates(J)=1.03D-10*TETABIS1*TETABIS2*TETABIS3

        reaction_rates(J)=reaction_rates(J)*UV_FLUX

      ! --- Self shielding of CO from Visser et al. (2009) ---
      elseif(is_absorption_co.eq.2) then
         int = 1.0D+00
         if (nco .gt. nco_co(47)) nco = nco_co(47) ! nco = max(nco_co) if greater than max(nco_co)
         if (nh2 .gt. nh2_co(42)) nh2 = nh2_co(42) ! nh2 = max(nh2_co) if greater than max(nh2_co)
         DO l=1,46
            IF((nco_co(l).le.nco).AND.(nco_co(l+1).ge.nco)) THEN
             DO m=1,41
                IF((NH2_co(m).le.nh2).AND.(NH2_co(m+1).ge.nh2)) THEN

                  ! Bilinear interpolation of the CO shielding function
                  xp = nh2
                  yp = nco

                  x1 = nh2_co(m)
                  x2 = nh2_co(m+1)

                  y1 = nco_co(l)
                  y2 = nco_co(l+1)

                  f11 = shielding_function_co(m,l)
                  f21 = shielding_function_co(m+1,l)
                  f12 = shielding_function_co(m,l+1)
                  f22 = shielding_function_co(m+1,l+1)

                  denom = x2 - x1
                  r1 = (x2 - xp) * f11 + (xp - x1) * f21
                  r1 = r1 / denom
                  r2 = (x2 - xp) * f12 + (xp - x1) * f22
                  r2 = r2 / denom

                  denom = y2 - y1
                  int = (y2 - yp) * r1 + (yp - y1) * r2
                  int = int / denom

                ENDIF
             ENDDO
            ENDIF
         ENDDO
         reaction_rates(J) = reaction_rates(J) * int
      endif
    endif

    ! ====== N2 self-shielding
    if ((REACTION_COMPOUNDS_NAMES(1,J).EQ.YN2).and.(is_absorption_n2.eq.1)) then
       int = 1.0D+00
       if (nn2 .gt. nn2_n2(46)) nn2 = nn2_n2(46) ! nn2 = max(nn2_n2) if greater than max(nn2_n2)
       if (nh2 .gt. nh2_n2(46)) nh2 = nn2_n2(46) ! nh2 = max(nh2_n2) if greater than max(nh2_n2)
       DO l=1,45
          IF((nn2_n2(l).le.nn2).AND.(nn2_n2(l+1).ge.nn2)) THEN
           DO m=1,45
              IF((nh2_n2(m).le.nh2).AND.(nh2_n2(m+1).ge.nh2)) THEN

                ! Bilinear interpolation of the N2 shielding function
                xp = nh2
                yp = nn2

                x1 = nh2_n2(m)
                x2 = nh2_n2(m+1)

                y1 = nn2_n2(l)
                y2 = nn2_n2(l+1)

                f11 = shielding_function_n2_nh1e20(m,l)
                f21 = shielding_function_n2_nh1e20(m+1,l)
                f12 = shielding_function_n2_nh1e20(m,l+1)
                f22 = shielding_function_n2_nh1e20(m+1,l+1)

                denom = x2 - x1
                r1 = (x2 - xp) * f11 + (xp - x1) * f21
                r1 = r1 / denom
                r2 = (x2 - xp) * f12 + (xp - x1) * f22
                r2 = r2 / denom

                denom = y2 - y1
                int = (y2 - yp) * r1 + (yp - y1) * r2
                int = int / denom

              ENDIF
           ENDDO
          ENDIF
       ENDDO
       reaction_rates(J) = reaction_rates(J) * int
    endif
  enddo



  !------------------------------------------------------------
  !
  !      --- GRAIN REACTIONS
  !
  !------------------------------------------------------------


  if (IS_GRAIN_REACTIONS.NE.0) then


  !-------------------------------------------------------------
  !  (99) Accretion rates on grain surfaces
  !-------------------------------------------------------------

  do J=type_id_start(99),type_id_stop(99)
    ! ========= Set accretion rates
    if((species_name(reactant_1_idx(J)).eq.YH).or.(species_name(reactant_1_idx(J)).eq.YH2).or.&
       (species_name(reactant_1_idx(J)).eq.'N2         ').or.&
       (species_name(reactant_1_idx(J)).eq.'CO         ').or.&
       (species_name(reactant_1_idx(J)).eq.'O2         ').or.&
       (species_name(reactant_1_idx(J)).eq.'CH4        ').or.&
       (species_name(reactant_1_idx(J)).eq.'CO2        ')) &
        call sticking_special_cases(j,sumlay)
    ACCRETION_RATES(reactant_1_idx(J)) = ACC_RATES_PREFACTOR(reactant_1_idx(J)) * TSQ * Y(reactant_1_idx(J)) * actual_gas_density
    ! when the adhoc h2 formation is activated 1/2 of the adsorbed H are availables for grain reactions (other than h2 formation) and
    ! 1/2 for the formation of h2
    IF((IS_H2_ADHOC_FORM.eq.1).AND.(species_name(reactant_1_idx(J)).eq.YH)) THEN
       ACCRETION_RATES(reactant_1_idx(J)) = 0.5D+00 * ACCRETION_RATES(reactant_1_idx(J))
    ENDIF
    ! When Eley-Rideal and complex induced reaction are activated we must be carreful on how accretions rate are computed
    IF(is_er_cir.ne.0) THEN

       IF(REACTION_COMPOUNDS_NAMES(1,J) == "C          " .AND.REACTION_COMPOUNDS_NAMES(4,J) == "JC         ") THEN
          IF(ab_surf.le.ab_lay) THEN
             ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J)) * &
                                                 (1.0D+00-(abH2O+abCO2+abNH3+abCH3OH+abCH4)/ab_lay)
          ELSE
             ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J)) * &
                                                 (1.0D+00-(abH2O+abCO2+abNH3+abCH3OH+abCH4)/ab_surf)
          ENDIF
       ENDIF
      ! PRINT *, " (1.0D+00-(abH2O+abCO2+abNH3+abCH3OH+abCH4)/ab_lay) for reaction:"&
       !       ,J,'is=', (1.0D+00-(abH2O+abCO2+abNH3+abCH3OH+abCH4)/ab_lay)
     !PRINT *, '(abBCO+abCO+abBC+abC+abBCH+abCH+abBCH2+abCH2+abBCH3+abCH3)='&
      !      , (abBCO+abCO+abBC+abC+abBCH+abCH+abBCH2+abCH2+abBCH3+abCH3)
     !PRINT *, 'ab_lay', ab_lay
     !PRINT *, 'ab_chem', ab_chem

       

       IF(REACTION_COMPOUNDS_NAMES(1,J) == "CH         " .AND.REACTION_COMPOUNDS_NAMES(4,J) == "JCH        ") THEN
          IF(ab_surf.le.ab_lay) THEN
             ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J)) * &
                                                  (1.0D+00-(abH2O+abNH3+abCH3OH)/ab_lay)
          ELSE
             ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J)) * &
                                                  (1.0D+00-(abH2O+abNH3+abCH3OH)/ab_surf)
          ENDIF
       ENDIF

       IF(REACTION_COMPOUNDS_NAMES(1,J) == "O          " .AND.REACTION_COMPOUNDS_NAMES(4,J) == "JO         ") THEN
          IF(ab_surf.le.ab_lay) THEN
             ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J))*(1.0D+00-abCO/ab_lay)
          ELSE
             ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J))*(1.0D+00-abCO/ab_surf)
          ENDIF
       ENDIF

    ENDIF

    ! in the case of ITYPE 99, reactant_2_idx(J) is not the second reactant of the reaction but the result of the 
    ! adsorption
    ACCRETION_RATES(reactant_2_idx(J)) = ACCRETION_RATES(reactant_1_idx(J))
    reaction_rates(J) = RATE_A(J) * branching_ratio(J) * ACCRETION_RATES(reactant_1_idx(J)) / Y(reactant_1_idx(J)) / GTODN
      !PRINT *, REACTION_COMPOUNDS_NAMES(1,J)
      !PRINT *, "The reaction rate (ITYPE 99) for reaction:",J,'is=',reaction_rates(J)
     ! PRINT *, "The branching ratio (ITYPE 88) for reaction:",J,'is=',branching_ratio(J)
     !PRINT *, "The ACCRETION_RATES (ITYPE 88) for reaction:",J,"is=", ACCRETION_RATES(reactant_1_idx(J))
     !PRINT *, "The ACC_RATES_PREFACTOR2(ITYPE 88) for reaction:",J,"is=", ACC_RATES_PREFACTOR2(reactant_1_idx(J))
  enddo
  

  !-------------------------------------------------------------
  !  (66) Photodesorption by external UV
  !       1.d8 is I_ISRF-FUV from
  !       Oberg et al. 2007, ApJ, 662, 23
  !-------------------------------------------------------------

  if (type_id_start(66).NE.0) then
      do J = type_id_start(66),type_id_stop(66)
          !---- Used for all species
          reaction_rates(J) = RATE_A(J) / SURFACE_SITE_DENSITY * UV_FLUX * 1.d8 * EXP(-2. * actual_av)
          !---- Specific cases
          !call photodesorption_special_cases(J,SUMLAY)
          !---- If there is more than MLAY on the grain surface, then we take into account that only
          !     the upper layers can photodesorb: this is done by assigning a reducing factor to the rate coefficient
          if(SUMLAY.GE.MLAY) reaction_rates(J) = reaction_rates(J) * MLAY / SUMLAY
          if (is_photodesorb.Eq.0) reaction_rates(J) = 0.D0
          EVAPORATION_RATES_TEMPO(reactant_1_idx(J))=EVAPORATION_RATES(reactant_1_idx(J))+reaction_rates(J)
          if (species_name(reactant_1_idx(J)).EQ.YJH2) EVAPORATION_RATES_TEMPO_H2=EVAPORATION_RATES_H2+reaction_rates(J)
      enddo
  endif


  !-------------------------------------------------------------
  !  (67) Photodesorption by CR generated UV
  !-------------------------------------------------------------

  if (type_id_start(67).NE.0) then
      do J = type_id_start(67),type_id_stop(67)
          !---- Used for all species
          reaction_rates(J) = RATE_A(J) / SURFACE_SITE_DENSITY * 1.d4 * UVCR
          !PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(4,J)  
          !PRINT *, "reaction_rates(J) (ITYPE 67) for reaction:",J,'is=',reaction_rates(J)
          !call photodesorption_special_cases(J,SUMLAY)
          !---- If there is more than MLAY on the grain surface, then we take into account that only
          !     the upper layers can photodesorb: this is done by assigning a reducing factor to the rate coefficient
          if(SUMLAY.GE.MLAY) reaction_rates(J) = reaction_rates(J) * MLAY / SUMLAY
          if (is_photodesorb.Eq.0) reaction_rates(J) = 0.D0
          EVAPORATION_RATES_TEMPO(reactant_1_idx(J))=EVAPORATION_RATES_TEMPO(reactant_1_idx(J))+reaction_rates(J)
          if (species_name(reactant_1_idx(J)).EQ.YJH2) then
          EVAPORATION_RATES_TEMPO_H2=EVAPORATION_RATES_H2+1.D30*reaction_rates(J)
          
          !PRINT *, "reaction_rates(J) (ITYPE 67) for reaction:",J,'is=',reaction_rates(J)
          endif
      enddo
  endif

!-------------------------------------------------------------
 
  !-------------------------------------------------------------
  !  (17) Photodiss by Cosmic rays on grain surfaces (s-1)
  !-------------------------------------------------------------

   do J=type_id_start(17),type_id_stop(17)
      reaction_rates(J)=RATE_A(J)*(CR_IONISATION_RATE + X_IONISATION_RATE)*Y(INDH2)/0.5D0
   enddo


  !-------------------------------------------------------------
  !  (18) Photodiss by Cosmic rays on grain surfaces (s-1)
  !-------------------------------------------------------------

   do J=type_id_start(18),type_id_stop(18)
      reaction_rates(J)=RATE_A(J)*(CR_IONISATION_RATE+X_IONISATION_RATE)*Y(INDH2)/0.5D0
! If ((REACTION_COMPOUNDS_NAMES(1,J).EQ.'J          ')) then
  !PRINT *, REACTION_COMPOUNDS_NAMES(1,J)
  !PRINT *, "The reaction rate(J) (ITYPE 18) for reaction:",J,'is=',reaction_rates(J)
  ! END IF
   enddo



  !-------------------------------------------------------------
  !  (14) Grain surface reactions
  !-------------------------------------------------------------

  do J=type_id_start(14),type_id_stop(14)
    IMOD1=0
    IMOD2=0
    BARR=1.0d0
    !PRINT *, 'actual_dust_temp=', actual_dust_temp

    ! --------- Thermal hopping diffusion method
    DIFFUSION_RATE_1(J)=THERMAL_HOPING_RATE(reactant_1_idx(J))
    DIFFUSION_RATE_2(J)=THERMAL_HOPING_RATE(reactant_2_idx(J))

    ! --------- Check for JH,JH2, and JO
    if (REACTION_COMPOUNDS_NAMES(1,J).EQ.YJH)  IMOD1=1
    if (REACTION_COMPOUNDS_NAMES(1,J).EQ.YJH2) IMOD1=2
    if (REACTION_COMPOUNDS_NAMES(1,J).EQ.YJO)  IMOD1=3
    if (REACTION_COMPOUNDS_NAMES(2,J).EQ.YJH)  IMOD2=1
    if (REACTION_COMPOUNDS_NAMES(2,J).EQ.YJH2) IMOD2=2
    if (REACTION_COMPOUNDS_NAMES(2,J).EQ.YJO)  IMOD2=3
    if (REACTION_COMPOUNDS_NAMES(1,J).EQ.YBH)  IMOD1=1
    if (REACTION_COMPOUNDS_NAMES(1,J).EQ.YBH2) IMOD1=2
    if (REACTION_COMPOUNDS_NAMES(1,J).EQ.YBO)  IMOD1=3
    if (REACTION_COMPOUNDS_NAMES(2,J).EQ.YBH)  IMOD2=1
    if (REACTION_COMPOUNDS_NAMES(2,J).EQ.YBH2) IMOD2=2
    if (REACTION_COMPOUNDS_NAMES(2,J).EQ.YBO)  IMOD2=3

    ! --------- QM for JH,JH2 only - others are too heavy
    if (IMOD1+IMOD2.NE.0) then
      ! ------------ QM1 - Tunnelling (if it's faster than thermal)
      if (GRAIN_TUNNELING_DIFFUSION.EQ.1) then
        tunneling_rate = TUNNELING_RATE_TYPE_1(reactant_1_idx(J))/nb_sites_per_grain
        if ((IMOD1.NE.0).AND.(tunneling_rate.GT.DIFFUSION_RATE_1(J))) then
          DIFFUSION_RATE_1(J)=tunneling_rate
        endif
        tunneling_rate = TUNNELING_RATE_TYPE_1(reactant_2_idx(J))/nb_sites_per_grain
        if ((IMOD2.NE.0).AND.(tunneling_rate.GT.DIFFUSION_RATE_2(J))) then
          DIFFUSION_RATE_2(J)=tunneling_rate
        endif
      endif
      ! ------------ QM2 - Tunnelling: use estimated width of lowest energy band (if it's faster than thermal)
      if (GRAIN_TUNNELING_DIFFUSION.EQ.2) then
        tunneling_rate = TUNNELING_RATE_TYPE_2(reactant_1_idx(J))/nb_sites_per_grain
        if ((IMOD1.NE.0).AND.(tunneling_rate.GT.DIFFUSION_RATE_1(J))) then
          DIFFUSION_RATE_1(J)=tunneling_rate
        endif
        tunneling_rate = TUNNELING_RATE_TYPE_2(reactant_2_idx(J))/nb_sites_per_grain
        if ((IMOD2.NE.0).AND.(tunneling_rate.GT.DIFFUSION_RATE_2(J))) then
          DIFFUSION_RATE_2(J)=tunneling_rate
        endif
      endif
      ! ------------ QM3 - Fastest out of thermal, QM1, QM2 rates
      if (GRAIN_TUNNELING_DIFFUSION.EQ.3) then
        if (IMOD1.NE.0) then
          tunneling_rate = TUNNELING_RATE_TYPE_1(reactant_1_idx(J))/nb_sites_per_grain
          if (tunneling_rate.GT.DIFFUSION_RATE_1(J)) then
            DIFFUSION_RATE_1(J) = tunneling_rate
          endif
          tunneling_rate = TUNNELING_RATE_TYPE_2(reactant_1_idx(J))/nb_sites_per_grain
          if (tunneling_rate.GT.DIFFUSION_RATE_1(J)) then
            DIFFUSION_RATE_1(J) = tunneling_rate
          endif
        endif
        if (IMOD2.NE.0) then
          tunneling_rate = TUNNELING_RATE_TYPE_1(reactant_2_idx(J))/nb_sites_per_grain
          if (tunneling_rate.GT.DIFFUSION_RATE_2(J))  then
            DIFFUSION_RATE_2(J) = tunneling_rate
          endif
          tunneling_rate = TUNNELING_RATE_TYPE_2(reactant_2_idx(J))/nb_sites_per_grain
          if (tunneling_rate.GT.DIFFUSION_RATE_2(J))  then
            DIFFUSION_RATE_2(J) = tunneling_rate
          endif
        endif
      endif
     ! ------------ QM2 - Tunnelling - case of atomic oxygen O : use estimated width of lowest energy band (if it's faster than thermal)
      if (GRAIN_TUNNELING_DIFFUSION.EQ.4) then
        tunneling_rate = TUNNELING_RATE_TYPE_2(reactant_1_idx(J))/nb_sites_per_grain
        if ((IMOD1.eq.3).AND.(tunneling_rate.GT.DIFFUSION_RATE_1(J))) then
          DIFFUSION_RATE_1(J)=tunneling_rate
        endif
        tunneling_rate = TUNNELING_RATE_TYPE_2(reactant_2_idx(J))/nb_sites_per_grain
        if ((IMOD1.eq.3).AND.(tunneling_rate.GT.DIFFUSION_RATE_2(J))) then
          DIFFUSION_RATE_2(J)=tunneling_rate
        endif
      endif
     endif

    BARR=1.0d0
    ! --------- Calculate activation energy barrier multiplier
    if (ACTIVATION_ENERGY(J).GE.1.0D-40) then
        ACTIV = ACTIVATION_ENERGY(J) / actual_dust_temp
        ! ------------ Choose fastest of classical or tunnelling
        if (ACTIV.GT.SURF_REACT_PROBA(J)) then
            ACTIV = SURF_REACT_PROBA(J)
        endif
        BARR=EXP(-ACTIV)
        if(is_reac_diff==1) then
           prob_reac = max(VIBRATION_FREQUENCY(reactant_1_idx(J)),VIBRATION_FREQUENCY(reactant_2_idx(J))) * &
                       EXP(-ACTIV)
           prob_deso = EVAPORATION_RATES_TEMPO(reactant_1_idx(J)) + &
                       EVAPORATION_RATES_TEMPO(reactant_2_idx(J))
           IF(ANY(REACTION_COMPOUNDS_NAMES(:,j)(1:1).eq.'K')) prob_deso = 0.0D+00
           prob_diff = (DIFFUSION_RATE_1(J) + DIFFUSION_RATE_2(J)) * nb_sites_per_grain

           barr = prob_reac+prob_deso+prob_diff
           barr = prob_reac / barr

        endif
    endif

    ! Modified rate have sense only for surface species
    IF(any(REACTION_COMPOUNDS_NAMES(:,j)(1:1).eq."J")) THEN

    ! --------- Modify according to MODIFY_RATE_FLAG switch:
    if (MODIFY_RATE_FLAG.NE.0) then
      ! ------------ if H+H->H2 is only modified rxn:
      if ((MODIFY_RATE_FLAG.EQ.-1).AND.(IMOD1.NE.1.OR.IMOD2.NE.1)) then
        IMOD1=0
        IMOD2=0
      endif

      ! ------------ if only H is modified:
      if ((MODIFY_RATE_FLAG.EQ.1).AND.(IMOD1.NE.1)) IMOD1=0
      if ((MODIFY_RATE_FLAG.EQ.1).AND.(IMOD2.NE.1)) IMOD2=0

      ! ------------ Set to modify all rates, if selected (just atoms)
      if (MODIFY_RATE_FLAG.EQ.3) then
        if ((REACTION_COMPOUNDS_NAMES(1,J).EQ.YJH).OR.&
        (REACTION_COMPOUNDS_NAMES(1,J).EQ.'JHe        ').OR.&
        (REACTION_COMPOUNDS_NAMES(1,J).EQ.'JC         ').OR.&
        (REACTION_COMPOUNDS_NAMES(1,J).EQ.'JN         ').OR.&
        (REACTION_COMPOUNDS_NAMES(1,J).EQ.'JO         ').OR.&
        (REACTION_COMPOUNDS_NAMES(1,J).EQ.'JS         ').OR.&
        (REACTION_COMPOUNDS_NAMES(1,J).EQ.'JSi        ').OR.&
        (REACTION_COMPOUNDS_NAMES(1,J).EQ.'JFe        ').OR.&
        (REACTION_COMPOUNDS_NAMES(1,J).EQ.'JNa        ').OR.&
        (REACTION_COMPOUNDS_NAMES(1,J).EQ.'JMg        ').OR.&
        (REACTION_COMPOUNDS_NAMES(1,J).EQ.'JP         ').OR.&
        (REACTION_COMPOUNDS_NAMES(1,J).EQ.'JF         ').OR.&
        (REACTION_COMPOUNDS_NAMES(1,J).EQ.'JCl        ')) IMOD1=3
        if ((REACTION_COMPOUNDS_NAMES(2,J).EQ.YJH).OR.&
        (REACTION_COMPOUNDS_NAMES(2,J).EQ.'JHe        ').OR.&
        (REACTION_COMPOUNDS_NAMES(2,J).EQ.'JC         ').OR.&
        (REACTION_COMPOUNDS_NAMES(2,J).EQ.'JN         ').OR.&
        (REACTION_COMPOUNDS_NAMES(2,J).EQ.'JO         ').OR.&
        (REACTION_COMPOUNDS_NAMES(2,J).EQ.'JS         ').OR.&
        (REACTION_COMPOUNDS_NAMES(2,J).EQ.'JSi        ').OR.&
        (REACTION_COMPOUNDS_NAMES(2,J).EQ.'JFe        ').OR.&
        (REACTION_COMPOUNDS_NAMES(2,J).EQ.'JNa        ').OR.&
        (REACTION_COMPOUNDS_NAMES(2,J).EQ.'JMg        ').OR.&
        (REACTION_COMPOUNDS_NAMES(2,J).EQ.'JP         ').OR.&
        (REACTION_COMPOUNDS_NAMES(2,J).EQ.'JF         ').OR.&
        (REACTION_COMPOUNDS_NAMES(2,J).EQ.'JCl        ')) IMOD2=3
      endif

      ! ------------ Modify rates (DIFFUSION_RATE_1 & DIFFUSION_RATE_2) according to their own evap/acc rates
      YMOD1=Y(reactant_1_idx(J))
      YMOD2=Y(reactant_2_idx(J))

      call modify_specific_rates(J,IMOD1,IMOD2,BARR,YMOD1,YMOD2)
    endif
    endif

    DIFF = DIFFUSION_RATE_1(J) + DIFFUSION_RATE_2(J)

    reaction_rates(J) = RATE_A(J) * branching_ratio(J) * BARR * DIFF * GTODN / actual_gas_density
      !IF (branching_ratio(J).eq.0) THEN
      !PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(2,J),&
       !        REACTION_COMPOUNDS_NAMES(4,J), REACTION_COMPOUNDS_NAMES(5,J)  
      !PRINT *, "RATE_A(J)",RATE_A(J)
      !PRINT *, "branching_ratio(J)",branching_ratio(J)
      !PRINT *, "BARR",BARR
      !PRINT *, "DIFF",DIFF
      !PRINT *, "GTODN / actual_gas_density",GTODN/actual_gas_density
      !PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(2,J) 
     ! PRINT *, "reaction_rates(J) (ITYPE 14) for reaction:",J,'is=',reaction_rates(J)
     ! ENDIF
    ! If the number of mantle layer is > 1, we consider that t(diff) (the time required by a species to 
    ! scan the entire grain sites) is given by the Number of sites on a layer time the number of layer time t(hop)
    IF(any(REACTION_COMPOUNDS_NAMES(:,j)(1:1).eq."K").and.sumlaymant.gt.1.0d0) THEN
      reaction_rates(J) = reaction_rates(J) / sumlaymant
    ENDIF

    ! H2 formation by LH mechanism is turned off when the ad hoc formation of H2 is activated
    IF ((reaction_compounds_names(1,J).EQ.YJH).AND.(reaction_compounds_names(2,J).EQ.YJH)) then
       IF(IS_H2_ADHOC_FORM.eq.1) reaction_rates(J) = 0.0D+00
    ENDIF

    ! "Encounter desorption" process for JH2 (Hincelin et al. 2014,A&A)
    ! The reaction JH2+JH2->JH2+H2 must be in the grain_reactions.in file to be accounted
    ! in practice the dominant processes are really the thermal hoping and thermal desorption of H2.
    if ((reaction_compounds_names(1,J).EQ.YJH2).AND.(reaction_compounds_names(2,J).EQ.YJH2)) then
       PROBH2H2=EVAPORATION_RATES_TEMPO_H2/(VIBRATION_FREQUENCY(reactant_1_idx(J))*EXP(-diff_binding_ratio_surf*ED_H2/&
                actual_dust_temp)/nb_sites_per_grain+ EVAPORATION_RATES_TEMPO_H2)
       reaction_rates(J)=reaction_rates(J)*PROBH2H2
    endif

    ! reaction_rates(J)=0.D0
  enddo

 

 
  !-------------------------------------------------------------
  !  (19-20) Photodissociations by UV photons on grain surfaces
  !   - 19 : Photodissociations by UV photons on grain surfaces
  !   - 20 : Photodissociations by UV photons on grain surfaces
  !          (when the gas-phase equivalent of the product is an
  !          ion)
  !-------------------------------------------------------------

  do J=type_id_start(19),type_id_stop(20)
    reaction_rates(J) = RATE_A(J) * EXP(-RATE_C(J) * actual_av) * UV_FLUX
   !PRINT *, REACTION_COMPOUNDS_NAMES(1,J)
   !PRINT *, "The reaction rate(J) (ITYPE 19 or 20) for reaction:",J,'is=',reaction_rates(J)
  enddo

  
  !-------------------------------------------------------------
  !  (30) Direct formation process with the incoming atom/molecule:
  !       Eley-Rideal process
  !-------------------------------------------------------------
  
  IF(type_id_start(30).ne.0) THEN

     DO J=type_id_start(30),type_id_stop(30)
        IF(ab_surf.le.ab_lay) THEN
           reaction_rates(J)=RATE_A(J)*branching_ratio(J)*ACC_RATES_PREFACTOR(reactant_1_idx(J))*TSQ/GTODN/ab_lay
        ELSE
           reaction_rates(J)=RATE_A(J)*branching_ratio(J)*ACC_RATES_PREFACTOR(reactant_1_idx(J))*TSQ/GTODN/ab_surf
        ENDIF
        IF(is_er_cir.eq.0) reaction_rates(J) = 0.0d+00
     ENDDO

  ENDIF

  
  
   !------------------------------------------------------------
  !
  !      --- GRAIN REACTIONS (CHEMISORPTION)
  !
  !------------------------------------------------------------

    if (is_chem.NE.0) then
  !-------------------------------------------------------------
  !  (88) Accretion rates on grain surfaces (Chemisorption)
  !-------------------------------------------------------------


  do J=type_id_start(88),type_id_stop(88)
    ! ========= Set accretion rates
    call sticking_special_cases_Chemisorption(J)
   
    ACCRETION_RATES(reactant_1_idx(J)) = ACC_RATES_PREFACTOR2(reactant_1_idx(J)) * TSQ * Y(reactant_1_idx(J)) * actual_gas_density
    ! when the adhoc h2 formation is activated 1/2 of the adsorbed H are availables for grain reactions (other than h2 formation) and
    ! 1/2 for the formation of h2
    IF((IS_H2_ADHOC_FORM.eq.1).AND.(species_name(reactant_1_idx(J)).eq.YH)) THEN
       ACCRETION_RATES(reactant_1_idx(J)) = 0.5D+00 * ACCRETION_RATES(reactant_1_idx(J))
    ENDIF

    !Calculating langmuir hinchelwood rejection terms

    IF(REACTION_COMPOUNDS_NAMES(1,j).eq."CO         ") THEN
        ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J)) * &
         (1.0d0 - (abBCO+abBC+abBCH+abBCH2+abBCH3)/ab_lay)
    !PRINT *, " (1.0d0 -(abBCO+abBC+abBCH+abBCH2+abBCH3)/ab_lay) for reaction:"&
     !       ,J,'is=', (1.0d0 -(abBCO+abBC+abBCH+abBCH2+abBCH3)/ab_lay)
    !PRINT *, " (0.5d0 -0.5d0*(abBCO+abCO+abBC+abC+abBCH+abCH+abBCH2+abCH2+abBCH3+abCH3)*1/ab_chem) for reaction:"&
           ! ,J,'is=', (0.5d0 - 0.5d0*(abBCO+abCO+abBC+abC+abBCH+abCH+abBCH2+abCH2+abBCH3+abCH3)*1/ab_chem)
     !PRINT *, '(abBCO+abCO+abBC+abC+abBCH+abCH+abBCH2+abCH2+abBCH3+abCH3)='&
            !, (abBCO+abCO+abBC+abC+abBCH+abCH+abBCH2+abCH2+abBCH3+abCH3)
     !PRINT *, 'ab_lay', ab_lay
    ! PRINT *, 'ab_chem', ab_chem
    ELSEIF (REACTION_COMPOUNDS_NAMES(1,j).eq."H2         ") THEN
        ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J)) * &
        ((1.0d0 -(abBH)/ab_lay))**2
    !PRINT *, "(1.0d0 -(abBH)/ab_lay)**2 for reaction:",J,'is=',((1.0d0 -(abBH)/ab_lay))**2
    !PRINT *, "(1.0d0 -(abBH+abH)/ab_lay)**2 for reaction:",J,'is=',(1.0d0 -(abBH+abH)*1/ab_lay)**2
    !PRINT *, "(1.0d0 -(abBH+abH)/ab_chem)**2 for reaction:",J,'is=',(1.0d0 -(abBH+abH)*1/ab_chem)**2
    
    ELSEIF (REACTION_COMPOUNDS_NAMES(1,j).eq."H          ") THEN
       ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J)) * &
        (1.0d0 -(abBH+abBH2)/ab_lay)
    !PRINT *, "(1.0d0 -(abBH+abBH2)/ab_lay) for reaction:",J,'is=',(1.0d0 -(abBH+abBH2)/ab_lay)
    ELSE
        ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J)) 
    ENDIF

     ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J)) * (1.0d0 - ab_chem/ab_lay)

    !PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(4,J)
    !PRINT *, 'ACCRETION_RATES(reactant_1_idx(J)) for reaction:',J,'is=', ACCRETION_RATES(reactant_1_idx(J))
    
    !IF(REACTION_COMPOUNDS_NAMES(1,j).eq."CO         ") THEN
     !    ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J)) * (0.5d0 - 0.5d0*ab_chem/ab_lay)
    !PRINT *, " (0.5d0 -0.5d0*(abBCO+abBC+abBCH+abBCH2+abBCH3)/ab_lay) for reaction:"&
           ! ,J,'is=', (0.5d0 -0.5d0*(abBCO+abBC+abBCH+abBCH2+abBCH3)/ab_lay)
    !PRINT *, " (0.5d0 -0.5d0*(abBCO+abCO+abBC+abC+abBCH+abCH+abBCH2+abCH2+abBCH3+abCH3)*1/ab_chem) for reaction:"&
           ! ,J,'is=', (0.5d0 - 0.5d0*(abBCO+abCO+abBC+abC+abBCH+abCH+abBCH2+abCH2+abBCH3+abCH3)*1/ab_chem)
     !PRINT *, '(abBCO+abCO+abBC+abC+abBCH+abCH+abBCH2+abCH2+abBCH3+abCH3)='&
            !, (abBCO+abCO+abBC+abC+abBCH+abCH+abBCH2+abCH2+abBCH3+abCH3)
     !PRINT *, 'ab_lay', ab_lay
    ! PRINT *, 'ab_chem', ab_chem
    !ELSE
     !   ACCRETION_RATES(reactant_1_idx(J)) = ACCRETION_RATES(reactant_1_idx(J)) * (1.0d0 - ab_chem/ab_lay)
    !PRINT *, "(1.0d0 -(abBH)/ab_lay)**2 for reaction:",J,'is=',(1.0d0 -(abBH)/ab_lay)**2
    !PRINT *, "(1.0d0 -(abBH+abH)/ab_lay)**2 for reaction:",J,'is=',(1.0d0 -(abBH+abH)*1/ab_lay)**2
    !PRINT *, "(1.0d0 -(abBH+abH)/ab_chem)**2 for reaction:",J,'is=',(1.0d0 -(abBH+abH)*1/ab_chem)**2
   ! ENDIF
    
    !IF (ab_chem.GE.ab_lay) ACCRETION_RATES(reactant_1_idx(J))= 0.0d0
   
    
 
    ! in the case of ITYPE 88, reactant_2_idx(J) is not the second reactant of the reaction but the result of the 
    ! adsorption
    ACCRETION_RATES(reactant_2_idx(J)) = ACCRETION_RATES(reactant_1_idx(J))
    reaction_rates(J) = RATE_A(J) *branching_ratio(J) * ACCRETION_RATES(reactant_1_idx(J)) / Y(reactant_1_idx(J)) / GTODN
    
    !PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(4,J)  
    !PRINT *, "The ACC_RATES_PREFACTOR2(ITYPE 88) for reaction:",J,"is=", ACC_RATES_PREFACTOR2(reactant_1_idx(J))
    !PRINT *, "The branching ratio (ITYPE 88) for reaction:",J,'is=',branching_ratio(J)
    !PRINT *, 'Y(reactant_1_idx(J))=', Y(reactant_1_idx(J))
    !PRINT *, 'RATE_A(J)=', RATE_A(J)
   
   ! PRINT *, "(1.0d0 -(abBH+abH+abH2+abBH2)*1/ab_lay) for reaction:",J,'is=',(1.0d0 -(abBH+abH+abH2+abBH2)*1/ab_lay)
    !PRINT *, "SUMLAY/MLAY for reaction:",J,'is=',SUMLAY/MLAY
    !PRINT *, "ab2CO+abCO+ab2C+abC+ab2CH+abCH+ab2CH2+abCH2+ab2CH3+abCH3 for reaction:"&
    !,J,'is=',ab2CO+abCO+ab2C+abC+ab2CH+abCH+ab2CH2+abCH2+ab2CH3+abCH3
     !PRINT *, REACTION_COMPOUNDS_NAMES(1,J)
   ! PRINT *, "The reaction rate (ITYPE 88) for reaction:",J,'is=',reaction_rates(J)
    !PRINT *, "The branching ratio (ITYPE 88) for reaction:",J,'is=',branching_ratio(J)
    ! PRINT *, "The ACCRETION_RATES (ITYPE 88) for reaction:",J,"is=", ACCRETION_RATES(reactant_1_idx(J))
    !PRINT *, "The ACC_RATES_PREFACTOR2(ITYPE 88) for reaction:",J,"is=", ACC_RATES_PREFACTOR2(reactant_1_idx(J))
  enddo

!-------------------------------------------------------------
  !  (25) Photodiss by Cosmic rays on grain surfaces (s-1) (Chemisorption)
  !-------------------------------------------------------------

   do J=type_id_start(25),type_id_stop(25)
      reaction_rates(J)=RATE_A(J)*(CR_IONISATION_RATE+X_IONISATION_RATE)*Y(INDH2)/0.5D0
! If ((REACTION_COMPOUNDS_NAMES(1,J).EQ.'J          ')) then
  !PRINT *, REACTION_COMPOUNDS_NAMES(1,J)
  !PRINT *, "The reaction rate(J) (ITYPE 25) for reaction:",J,'is=',reaction_rates(J)
  ! END IF
   enddo
   ! --------------------------------------------------------------------------
   !   - 24 :Photodissociations by UV photons on grain surfaces (Chemisorption)
   !---------------------------------------------------------------------------

  do J=type_id_start(24),type_id_stop(24)
    reaction_rates(J) = RATE_A(J) * EXP(-RATE_C(J) * actual_av) * UV_FLUX
   !PRINT *, REACTION_COMPOUNDS_NAMES(1,J)
   !PRINT *, "The reaction rate(J) (ITYPE 24) for reaction:",J,'is=',reaction_rates(J)
  enddo

  !--------------------------------------------------------------------------
  !   -29: Dissociation (Chemisorption)
  !--------------------------------------------------------------------------

  do J=type_id_start(29),type_id_stop(29)
   reaction_rates(J)=RATE_A(J)*VIBRATION_FREQUENCY(j)*EXP(-DISSOCIATION_ENERGY(j)/&
                     actual_dust_temp)
      !PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(4,J),REACTION_COMPOUNDS_NAMES(5,J) 
      !PRINT *, "The reaction rate (ITYPE 29) for reaction:",J,'is=',reaction_rates(J)
      !PRINT *, "The dissociation energy is=",DISSOCIATION_ENERGY(j)
  enddo


  !  (27) Photodesorption by external UV (Chemisorption)
  !       1.d8 is I_ISRF-FUV from
  !       Oberg et al. 2007, ApJ, 662, 23
  !-------------------------------------------------------------

  if (type_id_start(27).NE.0) then
      do J = type_id_start(27),type_id_stop(27)
          !---- Used for all species
          reaction_rates(J) = RATE_A(J) / SURFACE_SITE_DENSITY * UV_FLUX * 1.d8 * EXP(-2. * actual_av)
          !PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(4,J)  
          !PRINT *, "reaction_rates(J) (ITYPE 27) for reaction:",J,'is=',reaction_rates(J)
          !---- Specific cases
          !call photodesorption_special_cases(J,SUMLAY)
          !---- If there is more than MLAY on the grain surface, then we take into account that only
          !     the upper layers can photodesorb: this is done by assigning a reducing factor to the rate coefficient
          if(SUMLAY.GE.MLAY_CHEM) reaction_rates(J) = reaction_rates(J) * MLAY_CHEM / SUMLAY
          ! from parameters in
          if (is_photodesorb.Eq.0) reaction_rates(J) = 0.D0
          EVAPORATION_RATES_TEMPO(reactant_1_idx(J))=EVAPORATION_RATES(reactant_1_idx(J))+reaction_rates(J)
          !if (species_name(reactant_1_idx(J)).EQ.YJH2) EVAPORATION_RATES_TEMPO_H2=EVAPORATION_RATES_H2+reaction_rates(J)
         ! PRINT *, "RATE_A(J)=",RATE_A(J)
         ! PRINT *, "UV_FLUX=",UV_FLUX
         ! PRINT *, "SURFACE_SITE_DENSITY=",SURFACE_SITE_DENSITY
         ! PRINT *, "actual_av=",actual_av
         ! PRINT *, "EXP(-2. * actual_av)=",EXP(-2. * actual_av)
         ! PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(4,J)  
         !PRINT *, "reaction_rates(J) (ITYPE 27) for reaction:",J,'is=',reaction_rates(J)
         ! PRINT *, "EVAPORATION_RATES_TEMPO(ITYPE 27) for reaction:",J,'is=',EVAPORATION_RATES_TEMPO(reactant_1_idx(J))
      enddo
  endif


  !-------------------------------------------------------------
  !  (28) Photodesorption by CR generated UV (Chemisorption)
  !-------------------------------------------------------------

  if (type_id_start(28).NE.0) then
      do J = type_id_start(28),type_id_stop(28)
      !PRINT *, 'IN ITYPE 28' 
          !---- Used for all species
          reaction_rates(J) = RATE_A(J) / SURFACE_SITE_DENSITY * 1.d4 * UVCR
          ! PRINT *, "reaction_rates(J) (ITYPE 28) for reaction:",J,'is=',reaction_rates(J)
          ! PRINT *, 'UVCR=', UVCR
          ! PRINT *, 'SURFACE_SITE_DENSITY',SURFACE_SITE_DENSITY
          !call photodesorption_special_cases(J,SUMLAY)
          !---- If there is more than MLAY on the grain surface, then we take into account that only
          !     the upper layers can photodesorb: this is done by assigning a reducing factor to the rate coefficient
          if(SUMLAY.GE.MLAY_CHEM) reaction_rates(J) = reaction_rates(J) * MLAY_CHEM / SUMLAY
          ! from parameters in
          if (is_photodesorb.Eq.0) reaction_rates(J) = 0.D0
          EVAPORATION_RATES_TEMPO(reactant_1_idx(J))=EVAPORATION_RATES_TEMPO(reactant_1_idx(J))+reaction_rates(J)
         ! if (species_name(reactant_1_idx(J)).EQ.YJH2) then
          !EVAPORATION_RATES_TEMPO_H2=EVAPORATION_RATES_H2+1.D30*reaction_rates(J)
         ! PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(4,J)  
          !PRINT *, "reaction_rates(J) (ITYPE 28) for reaction:",J,'is=',reaction_rates(J)
          !PRINT *, "EVAPORATION_RATES_TEMPO(ITYPE 28) for reaction:",J,'is=',EVAPORATION_RATES_TEMPO(reactant_1_idx(J))
         ! PRINT *, 'UVCR=', UVCR
         ! PRINT *, 'SURFACE_SITE_DENSITY',SURFACE_SITE_DENSITY
         ! endif
      enddo
  endif


   !--------------------------------------------------------------------
   !(23): High temperature grain surface reaction (set grain_tunneling_diffusion to 2, QM2 or 0 for no tunneling through diffusion barrier)
   !--------------------------------------------------------------------

   do J=type_id_start(23),type_id_stop(23)
    IMOD1=0
    IMOD2=0
    BARR=1.0d0
    !PRINT *, 'actual_dust_temp=', actual_dust_temp

    ! --------- Thermal hopping diffusion method
    DIFFUSION_RATE_1(J)=THERMAL_HOPING_RATE(reactant_1_idx(J))
    DIFFUSION_RATE_2(J)=THERMAL_HOPING_RATE(reactant_2_idx(J))

    ! --------- Check for JH,JH2,BH,BH2 and JO
    if (REACTION_COMPOUNDS_NAMES(1,J).EQ.YJH)  IMOD1=1
    if (REACTION_COMPOUNDS_NAMES(1,J).EQ.YJH2) IMOD1=2
    if (REACTION_COMPOUNDS_NAMES(1,J).EQ.YJO)  IMOD1=3
    if (REACTION_COMPOUNDS_NAMES(2,J).EQ.YJH)  IMOD2=1
    if (REACTION_COMPOUNDS_NAMES(2,J).EQ.YJH2) IMOD2=2
    if (REACTION_COMPOUNDS_NAMES(2,J).EQ.YJO)  IMOD2=3
    if (REACTION_COMPOUNDS_NAMES(1,J).EQ.YBH)  IMOD1=1
    if (REACTION_COMPOUNDS_NAMES(1,J).EQ.YBH2) IMOD1=2
    if (REACTION_COMPOUNDS_NAMES(2,J).EQ.YBH)  IMOD2=1
    if (REACTION_COMPOUNDS_NAMES(2,J).EQ.YBH2) IMOD2=2

    ! --------- QM for JH,JH2,BH,BH2 only - others are too heavy
    if (IMOD1+IMOD2.NE.0) then

      ! ------------ QM2 - Tunnelling: use estimated width of lowest energy band (if it's faster than thermal)
      if (GRAIN_TUNNELING_DIFFUSION.EQ.2) then
        tunneling_rate = TUNNELING_RATE_TYPE_2(reactant_1_idx(J))/nb_sites_per_grain
        if ((IMOD1.NE.0).AND.(tunneling_rate.GT.DIFFUSION_RATE_1(J))) then
          DIFFUSION_RATE_1(J)=tunneling_rate
        endif
        tunneling_rate = TUNNELING_RATE_TYPE_2(reactant_2_idx(J))/nb_sites_per_grain
        if ((IMOD2.NE.0).AND.(tunneling_rate.GT.DIFFUSION_RATE_2(J))) then
          DIFFUSION_RATE_2(J)=tunneling_rate
        endif
      endif
     endif

    BARR=1.0d0
    ! --------- Calculate activation energy barrier multiplier
    if (ACTIVATION_ENERGY(J).GE.1.0D-40) then
        ACTIV = ACTIVATION_ENERGY(J) / actual_dust_temp
        ! ------------ Choose fastest of classical or tunnelling
        if (ACTIV.GT.SURF_REACT_PROBA(J)) then
            ACTIV = SURF_REACT_PROBA(J)
        endif
        BARR=EXP(-ACTIV)
        if(is_reac_diff==1) then
           prob_reac = max(VIBRATION_FREQUENCY(reactant_1_idx(J)),VIBRATION_FREQUENCY(reactant_2_idx(J))) * &
                       EXP(-ACTIV)
           IF(ANY(REACTION_COMPOUNDS_NAMES(:,j)(1:1).eq.'K')) prob_deso = 0.0D+00
           prob_diff = (DIFFUSION_RATE_1(J) + DIFFUSION_RATE_2(J)) * nb_sites_per_grain

           BARR = prob_reac+prob_diff
           BARR = prob_reac / BARR

        endif
    endif

  
    DIFF = DIFFUSION_RATE_1(J) + DIFFUSION_RATE_2(J)


    reaction_rates(J) = RATE_A(J) * branching_ratio(J) * BARR * DIFF * GTODN / actual_gas_density
      !PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(2,J),&
       !        REACTION_COMPOUNDS_NAMES(4,J), REACTION_COMPOUNDS_NAMES(5,J)
      !PRINT *, "RATE_A(J)",RATE_A(J)
      !PRINT *, "branching_ratio(J)=",branching_ratio(J)
      !PRINT *, "BARR",BARR
      !PRINT *, "DIFF",DIFF
      !PRINT *, "GTODN / actual_gas_density",GTODN/actual_gas_density
      !PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(2,J), REACTION_COMPOUNDS_NAMES(4,J)  
      !PRINT *, "reaction_rates(J) (ITYPE 23) for reaction:",J,'is=',reaction_rates(J)
    
    ! If the number of mantle layer is > 1, we consider that t(diff) (the time required by a species to 
    ! scan the entire grain sites) is given by the Number of sites on a layer time the number of layer time t(hop)
    !IF(any(REACTION_COMPOUNDS_NAMES(:,j)(1:1).eq."K").and.sumlaymant.gt.1.0d0) THEN
     ! reaction_rates(J) = reaction_rates(J) / sumlaymant
    !ENDIF

    ! H2 formation by LH mechanism is turned off when the ad hoc formation of H2 is activated
    IF ((reaction_compounds_names(1,J).EQ.YJH).AND.(reaction_compounds_names(2,J).EQ.YJH)) then
       IF(IS_H2_ADHOC_FORM.eq.1) reaction_rates(J) = 0.0D+00
    ENDIF
    IF ((reaction_compounds_names(1,J).EQ.YBH).AND.(reaction_compounds_names(2,J).EQ.YBH)) then
       IF(IS_H2_ADHOC_FORM.eq.1) reaction_rates(J) = 0.0D+00
    ENDIF

    ! "Encounter desorption" process for JH2 (Hincelin et al. 2014,A&A)
    ! The reaction JH2+JH2->JH2+H2 must be in the grain_reactions.in file to be accounted
    ! in practice the dominant processes are really the thermal hoping and thermal desorption of H2.
    if ((reaction_compounds_names(1,J).EQ.YJH2).AND.(reaction_compounds_names(2,J).EQ.YJH2)) then
       PROBH2H2=EVAPORATION_RATES_TEMPO_H2/(VIBRATION_FREQUENCY(reactant_1_idx(J))*EXP(-diff_binding_ratio_surf*ED_H2/&
                actual_dust_temp)/nb_sites_per_grain+ EVAPORATION_RATES_TEMPO_H2)
       reaction_rates(J)=reaction_rates(J)*PROBH2H2
    endif

    ! reaction_rates(J)=0.D0
  enddo
  
  !-------------------------------------------------------------
  !  (26) Thermal evaporation (s-1) (Chemisorption)
  !-------------------------------------------------------------

  do J=type_id_start(26),type_id_stop(26)
    reaction_rates(J)=RATE_A(J)*branching_ratio(J)*VIBRATION_FREQUENCY(reactant_1_idx(J))*EXP(-BINDING_ENERGY(reactant_1_idx(J))/&
                      actual_dust_temp)
    EVAPORATION_RATES(reactant_1_idx(J))=EVAPORATION_RATES(reactant_1_idx(J))+reaction_rates(J)
    !if (species_name(reactant_1_idx(J)).EQ.YJH2) EVAPORATION_RATES_H2=RATE_A(J)*branching_ratio(J)*&
     !     VIBRATION_FREQUENCY(reactant_1_idx(J))*EXP(-ED_H2/actual_dust_temp)
  !If ((REACTION_COMPOUNDS_NAMES(1,J)).EQ.'G          ') then
      !PRINT *, REACTION_COMPOUNDS_NAMES(1,J), REACTION_COMPOUNDS_NAMES(2,J) 
      !PRINT *, "The reaction rate (ITYPE 26) for reaction:",J,'is=',reaction_rates(J)
      !PRINT *, "The binding energy is=",BINDING_ENERGY(reactant_1_idx(J))
      !END IF
      ! call exit()
  enddo

 endif
endif
  ! Continually time-dependent gas phase rates============================
  ! H2 formation
  ! branching_ratio(1) and branching_ratio(2) are zero if IS_GRAIN_REACTIONS=1
  ! cf GRAINRATE
  ! VW Fev 2012 - this process has been removed
  !      do j=type_id_start(0),type_id_stop(0)
  !      if ((REACTION_COMPOUNDS_NAMES(1,J).eq.YH).and.(REACTION_COMPOUNDS_NAMES(2,j).eq.YH)) then
  !      reaction_rates(j)=branching_ratio(j)*A(j)*(T300**B(j))*GTODN/actual_gas_density/Y(reactant_1_idx(j))
  !      endif
  !      enddo

  ! if rate acoefficients are too small, put them to 0
  where (reaction_rates.lt.MINIMUM_RATE_COEFFICIENT) reaction_rates = 0.d0

!do j=1,nb_reactions
!write(*,*) REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS,j), reaction_rates(j)
!enddo
!stop

    return
    end subroutine set_dependant_rates
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author
!> Valentine Wakelam
!
!> @date 2012
!
! DESCRIPTION:
!> @brief Set Reactions coefficient formally dependent on the abundances Y.
!! Grain surface reactions, self shielding, etc...
!
! EVAPORATION_RATES_TEMPO was created because this subroutine is called several times
! during the run of one model
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine set_dependant_rates_3phase(Y)

  use global_variables

  implicit none

  real(double_precision), intent(in), dimension(nb_species) :: Y !< [in] abundances (relative to H) [number ratio]

  REAL(double_precision) :: ab_surf !< Sum of all abundances on grain surface
  REAL(double_precision) :: ab_mant !< Sum of all abundances on grain mantle
  real(double_precision) :: ab_lay  !< Abundance in 1 layer

  real(double_precision) :: sumlaysurf
  real(double_precision) :: sumlaymant

  real(double_precision) :: rate_tot_tmp1,rate_tot_tmp2

  REAL(double_precision) :: alpha_3_phase
  real(double_precision) :: rate_tot_mant_surf, rate_tot_surf_mant,rate_mant_to_surf, rate_surf_to_mant
  integer :: j

  ab_surf = 0.0d0
  ab_mant = 0.0d0
  DO J = nb_gaseous_species+1,nb_species
     IF (species_name(j)(1:1).eq."J") ab_surf = ab_surf + Y(j)
     IF (species_name(j)(1:1).eq."B") ab_surf = ab_surf + Y(j)
     IF (species_name(j)(1:1).eq."K") ab_mant = ab_mant + Y(j)
  ENDDO

  sumlaysurf = ab_surf*GTODN/nb_sites_per_grain
  sumlaymant = ab_mant*GTODN/nb_sites_per_grain

!-------------------------------------------------------------
!  (40-41) 3 phase model
!          - 40 : Surface to Mantle
!          - 41 : Mantle to Surface
!  The mantle to surface reactions need to be calculated first
!  for the swapping rates
!-------------------------------------------------------------

  if ((type_id_start(40).NE.0).AND.(type_id_start(41).NE.0)) then

     rate_tot_mant_surf = 0.d0
     rate_tot_surf_mant = 0.d0

     rate_tot_tmp1 = 0.0d0
     rate_tot_tmp2 = 0.0d0

     ! Mantle to Surface
     do J=type_id_start(41),type_id_stop(41)

        ! From Garrod & Pauly 2011
        alpha_3_phase = ab_mant / ab_surf
        IF(alpha_3_phase.ge.1.0d0) alpha_3_phase = 1.0d0

        reaction_rates(j) = -  alpha_3_phase * rate_tot_des / ab_mant

        ! The swapping rates are added to the rate here
        if (y(reactant_1_idx(J)).gt.1.0d-99) THEN
           rate_mant_to_surf = y(reactant_1_idx(J))*THERMAL_HOPING_RATE(reactant_1_idx(J))
           ! --- Based on Garrod 2013 (I don't really understande the multiplicative factor sumlaysurf here...)
           !if(sumlaysurf.le.sumlaymant) rate_mant_to_surf = rate_mant_to_surf * sumlaysurf / sumlaymant
           if(sumlaymant.ge.1.0d0) rate_mant_to_surf = rate_mant_to_surf / sumlaymant

           rate_tot_mant_surf = rate_tot_mant_surf + rate_mant_to_surf
           rate_mant_to_surf = rate_mant_to_surf / y(reactant_1_idx(J))
        ELSE
           rate_mant_to_surf = 0.0d0
        ENDIF

        reaction_rates(j) = reaction_rates(j) + rate_mant_to_surf

        if (is_3_phase.eq.0) reaction_rates(j) = 0.0d0

        rate_tot_tmp2 = rate_tot_tmp2 + reaction_rates(j)*y(reactant_1_idx(j))

  enddo


  ! Surface to Mantle
  do J=type_id_start(40),type_id_stop(40)

     alpha_3_phase = ab_surf*GTODN/nb_sites_per_grain
     alpha_3_phase = alpha_3_phase / nb_active_lay

     reaction_rates(j) = alpha_3_phase * rate_tot_acc / ab_surf

     ! The swapping rates are added to the rate here
     if (y(reactant_1_idx(J)).gt.1.0d-99) THEN
        rate_surf_to_mant = y(reactant_1_idx(J)) * rate_tot_mant_surf / ab_surf
        rate_tot_surf_mant = rate_tot_surf_mant + rate_surf_to_mant
        rate_surf_to_mant = rate_surf_to_mant / y(reactant_1_idx(J))
     ELSE
        rate_surf_to_mant = 0.0d0
     ENDIF

     reaction_rates(j) = reaction_rates(j) + rate_surf_to_mant

     if (is_3_phase.eq.0) reaction_rates(j) = 0.0d0

     rate_tot_tmp1 = rate_tot_tmp1 + reaction_rates(j)*y(reactant_1_idx(j))

  enddo
  !PRINT*, "Acc/Des  ", rate_tot_tmp1/rate_tot_tmp2, rate_tot_acc/rate_tot_des, sumlaysurf
  endif
  return
end subroutine set_dependant_rates_3phase

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author
!> Franck Hersant
!
!> @date 2003
!
! DESCRIPTION:
!> @brief Modify some reaction rates on the grain surface (itype=14) in
!! various conditions. Test to estimates the fastest process and replace them.
!
! EVAPORATION_RATES contains all the processes of evaporation. If you add a
! new process do not forget to add it to this table.
! DIFFUSION_RATE_1 and DIFFUSION_RATE_2 contains the CRID if the flag for CRID is 
! one in the parameter file.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine modify_specific_rates(J,IMOD1,IMOD2,BARR,YMOD1,YMOD2)
  use global_variables
  
  implicit none

  ! Inputs
  integer, intent(in) :: J !<[in] index of a given reaction
  integer, intent(in) :: IMOD1 !<[in] modify rate flag for reactant 1
  integer, intent(in) :: IMOD2 !<[in] modify rate flag for reactant 2
  real(double_precision), intent(in) :: BARR !<[in] TODO Description/units ???
  real(double_precision), intent(in) :: YMOD1 !<[in] Abundance (relative to H) [number ratio] for reactant 1
  real(double_precision), intent(in) :: YMOD2 !<[in] Abundance (relative to H) [number ratio] for reactant 2
  
  ! Locals
  real(double_precision) :: TESTREF1, TESTREF2, TESTNUM
  integer :: PICK

  EVAP_OVER_ACC_RATIO_1(J) = 0.d0
  EVAP_OVER_ACC_RATIO_2(J) = 0.d0


!if ((species_name(reactant_1_idx(J)).EQ.'JH         ').AND.(species_name(reactant_2_idx(J)).EQ.'JC         ')) then
!write(*,*) 'ACCRETION_RATES of H (in modified rates)',ACCRETION_RATES(reactant_1_idx(J))
!write(*,*) 'ACCRETION_RATES of C (in modified rates)',ACCRETION_RATES(reactant_2_idx(J))
!write(*,*) 'Activation barrier',BARR
!write(*,*) 'EVAPORATION_RATES of H (in modified rates)',EVAPORATION_RATES_TEMPO(reactant_1_idx(J))
!write(*,*) 'EVAPORATION_RATES of C (in modified rates)',EVAPORATION_RATES_TEMPO(reactant_2_idx(J))
!write(*,*) 'DIFFUSION RATES of H (before modified rates)',DIFFUSION_RATE_1(J)
!write(*,*) 'DIFFUSION RATES of C (before modified rates)',DIFFUSION_RATE_2(J)
!endif

 ! --- Check value of x = t_acc/t_evap
  ! EVAPORATION_RATES = 1/t_evap
  ! ACCRETION_RATES = 1/t_acc
  if (ACCRETION_RATES(reactant_1_idx(J)).GT.0.d0) then
    EVAP_OVER_ACC_RATIO_1(J) = EVAPORATION_RATES_TEMPO(reactant_1_idx(J)) / ACCRETION_RATES(reactant_1_idx(J))
  endif
  if (ACCRETION_RATES(reactant_2_idx(J)).GT.0.d0) then
    EVAP_OVER_ACC_RATIO_2(J)=EVAPORATION_RATES_TEMPO(reactant_2_idx(J))/ACCRETION_RATES(reactant_2_idx(J))
  endif
  ! Hence x = 0 if t_evap or t_acc = 0

  ! --- Assign max rates

  if (BARR.EQ.1.0d0) then
    if (IMOD1.NE.0) then
      if (EVAP_OVER_ACC_RATIO_1(J).LT.1.0d0) then ! accretion dominates
        if (DIFFUSION_RATE_1(J).GT.ACCRETION_RATES(reactant_1_idx(J))) then
          DIFFUSION_RATE_1(J) = ACCRETION_RATES(reactant_1_idx(J))
        endif
      else ! evaporation dominates
        if (DIFFUSION_RATE_1(J).GT.EVAPORATION_RATES_TEMPO(reactant_1_idx(J))) then
          DIFFUSION_RATE_1(J) = EVAPORATION_RATES_TEMPO(reactant_1_idx(J))
        endif
      endif
    endif

    if (IMOD2.NE.0) then
      if (EVAP_OVER_ACC_RATIO_2(J).LT.1.0d0) then ! accretion dominates
        if (DIFFUSION_RATE_2(J).GT.ACCRETION_RATES(reactant_2_idx(J))) then
          DIFFUSION_RATE_2(J) = ACCRETION_RATES(reactant_2_idx(J))
        endif
      else ! evaporation dominates
        if (DIFFUSION_RATE_2(J).GT.EVAPORATION_RATES_TEMPO(reactant_2_idx(J))) then
          DIFFUSION_RATE_2(J) = EVAPORATION_RATES_TEMPO(reactant_2_idx(J))
        endif
      endif
    endif
  endif

  ! --- Species rate to compare chosen by fastest diffusion rate
  if (BARR.NE.1.0d0) then
    PICK=0

    TESTREF1=ACCRETION_RATES(reactant_1_idx(J))
    if (EVAP_OVER_ACC_RATIO_1(J).GE.1.0d0) TESTREF1=EVAPORATION_RATES_TEMPO(reactant_1_idx(J))
    TESTREF2=ACCRETION_RATES(reactant_2_idx(J))
    if (EVAP_OVER_ACC_RATIO_2(J).GE.1.0d0) TESTREF2=EVAPORATION_RATES_TEMPO(reactant_2_idx(J))

    if (DIFFUSION_RATE_1(J).GE.DIFFUSION_RATE_2(J)) then
      TESTNUM=(DIFFUSION_RATE_1(J)+DIFFUSION_RATE_2(J))*BARR*YMOD2*GTODN
      if (YMOD2*GTODN.LT.1.0d0) TESTNUM=(DIFFUSION_RATE_1(J)+DIFFUSION_RATE_2(J))*BARR
      if (TESTNUM.GT.TESTREF1) PICK=1
    endif
    if (DIFFUSION_RATE_2(J).GT.DIFFUSION_RATE_1(J)) then
      TESTNUM=(DIFFUSION_RATE_1(J)+DIFFUSION_RATE_2(J))*BARR*YMOD1*GTODN
      if (YMOD1*GTODN.LT.1.0d0) TESTNUM=(DIFFUSION_RATE_1(J)+DIFFUSION_RATE_2(J))*BARR
      if (TESTNUM.GT.TESTREF2) PICK=2
    endif

    if (PICK.EQ.1) then
      DIFFUSION_RATE_1(J)=TESTREF1/BARR/YMOD2/GTODN
      if (YMOD2*GTODN.LT.1.0d0) DIFFUSION_RATE_1(J)=TESTREF1/BARR
      DIFFUSION_RATE_2(J)=0.d0
    endif

    if (PICK.EQ.2) then
      DIFFUSION_RATE_2(J)=TESTREF2/BARR/YMOD1/GTODN
      if (YMOD1*GTODN.LT.1.0d0) DIFFUSION_RATE_2(J)=TESTREF2/BARR
      DIFFUSION_RATE_1(J)=0.d0
    endif

  endif

return
  end subroutine modify_specific_rates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Treat special cases of photodesorption such as CO2, CO, H2O, CH3OH and N2
!!\n Data from Oberg et al. 2009:
!!\n      a - A&A, 496, 281-293
!!\n      b - ApJ, 693, 1209-1218
!!\n      c - A&A, 504, 891-913
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine photodesorption_special_cases(J,SUMLAY)

use global_variables

implicit none

! Inputs 
integer, intent(in) :: J !<[in] index of a given reaction
real(double_precision), intent(in) :: SUMLAY !<[in] TODO description/units ???

! Locals
real(double_precision) :: LTD
real(double_precision) :: fH2O
character(len=11) :: reactant1, product1, product2 !< species names for first reactant, and first and second products of a given reaction

reactant1 = REACTION_COMPOUNDS_NAMES(1,J)
product1 = REACTION_COMPOUNDS_NAMES(MAX_REACTANTS+1,J)
product2 = REACTION_COMPOUNDS_NAMES(MAX_REACTANTS+2,J)

! TODO are extra spaces in names really necessary? Test that.
!------ Photodesorption of CO2: photodesorbs as either CO2 or CO
if (actual_dust_temp.LE.3.5E+01) then
    if (product1 == 'CO2        ') then
        reaction_rates(J) = reaction_rates(J) * 1.2E-03 * (1.d0 - EXP(-SUMLAY/2.9E+00) ) / RATE_A(J)
    endif
    if((product1 == 'CO         ' .AND. product2 == 'O          ') .OR. &
       (product1 == 'O          ' .AND. product2 == 'CO         ')) then
        reaction_rates(J) = reaction_rates(J) * 1.1E-03 * (1.d0 - EXP(-SUMLAY/4.6E+00) ) / RATE_A(J)
    endif
elseif(actual_dust_temp.GT.3.5E+01) then
    if (product1 == 'CO2        ') then
        reaction_rates(J) = reaction_rates(J) * 2.2E-03 * (1.d0 - EXP(-SUMLAY/5.8E+00) ) / RATE_A(J)
    endif
    if((product1 == 'CO         ' .AND. product2 == 'O          ') .OR. &
       (product1 == 'O          ' .AND. product2 == 'CO         ')) then
       reaction_rates(J) = reaction_rates(J) * 2.2E-04 * SUMLAY / RATE_A(J)
    endif
endif
!------ Photodesorption of CO
if(reactant1 == 'JCO        ' .AND. product1 == 'CO         ') then
   reaction_rates(J) = reaction_rates(J) * ( 2.7E-03 - 1.7E-04 * (actual_dust_temp - 15E+00)) / RATE_A(J)
endif
!------ Photodesorption of N2
if(reactant1 == 'JN2        ' .AND. product1 == 'N2         ') then
   reaction_rates(J) = reaction_rates(J) * 4.0E-04 / RATE_A(J)
endif
!------ Photodesorption of CH3OH
if(reactant1 == 'JCH3OH     ' .AND. product1 == 'CH3OH      ') then
   reaction_rates(J) = reaction_rates(J) * 2.1E-03 / RATE_A(J)
endif
!------ Photodesorption of H2O
if(reactant1 == 'JH2O       ') then
   LTD = 6.0E-01 + 2.4E-02 * actual_dust_temp
   fH2O = 4.2E-01 + 2.0E-03 * actual_dust_temp
   reaction_rates(J) = reaction_rates(J) * 1.0E-03 * (1.3E+00 + 3.2E-02 * actual_dust_temp) &
           & * (1.d0 - EXP(-SUMLAY/LTD) ) * fH2O / RATE_A(J)
   if((product1 == 'OH         ' .AND. product2 == 'H          ') .OR. &
      (product1 == 'H          ' .AND. product2 == 'OH         ')) then
      reaction_rates(J) = reaction_rates(J) * (1.0E+00 - fH2O)/fH2O
   endif
endif

end subroutine photodesorption_special_cases

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author
!> Maxime Ruaud
!
!> @date 2014
!
! DESCRIPTION:
!> @brief Treat H and H2 sticking coefficient from Chaabouni et al. 2012
!!\n url: http://cdsads.u-strasbg.fr/abs/2012A%26A...538A.128C
!!\n A smooth transition between bare grains and ices is computed
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine sticking_special_cases(J,SUMLAY)

use global_variables

implicit none


! Inputs
integer, intent(in) :: J !<[in] index of a given reaction
real(double_precision), intent(in) :: SUMLAY

! Local
real(double_precision) :: cond
real(double_precision) :: stick, stick_ice, stick_bare


!PRINT *, 'in sticking_special_cases'
!call exit ()

cond=PI*grain_radius*grain_radius*SQRT(8.0d0*K_B/PI/AMU)
!print *, 'cond=', cond
if (species_name(reactant_1_idx(J)).eq.YH)  then
   stick_bare = (1.0d+00 + 2.5d+00*actual_gas_temp/25.0d+00)/&
                (1.0d+00 + actual_gas_temp/25.0d+00)**2.5d+00
   stick_ice  = (1.0d+00 + 2.5d+00*actual_gas_temp/52.0d+00)/&
                (1.0d+00 + actual_gas_temp/52.0d+00)**2.5d+00
endif

if (species_name(reactant_1_idx(J)).eq.YH2) then
   stick_bare = 0.95d+00 * (1.0d+00 + 2.5d+00*actual_gas_temp/56.0d+00)/&
                           (1.0d+00 + actual_gas_temp/56.0d+00)**2.5d+00
   stick_ice  = 0.76d+00 * (1.0d+00 + 2.5d+00*actual_gas_temp/87.0d+00)/&
                           (1.0d+00 + actual_gas_temp/87.0d+00)**2.5d+00
endif

! When the grain coverage is less than 1 ML we use the "silicates" expression and slowely tends to the "ASW ice" expression
if(sumlay.le.1.0d+00) then
  stick = (1.0d+00-sumlay)*stick_bare + sumlay*stick_ice
else
! When the grain coverage is more than 1 ML we use the "ASW ice" expression
  stick = stick_ice
endif


    ! --- Evaluate sticking coeff and accretion rate factor for each species using Acharyya 2016 paper
    !STICK=0.d0
   
 
    ! PRINT *, 'initial_gas_temperature=',initial_gas_temperature
    ! PRINT *, 'initial_dust_temperature=',initial_dust_temperature
    ! PRINT *, 'actual_gas_temp=',actual_gas_temp
    ! PRINT *, 'actual_dust_temp=',actual_dust_temp

    !gas_temperature=actual_gas_temp 
    !dust_temperature=actual_dust_temp 
    !actual_gas_temp = initial_dust_temperature
    !actual_dust_temp = initial_dust_temperature
 if (species_name(reactant_1_idx(J)).eq."N2         ")  then
  STICK=.5*(1-TANH(1.2d-01*(actual_dust_temp-4.3d-02*BINDING_ENERGY(reactant_2_idx(J)))))
  !PRINT *, 'species_name=', species_name(reactant_1_idx(J))
  !PRINT *, "STICK=",STICK
 endif
! PRINT *, "STICK=",STICK
 if (species_name(reactant_1_idx(J)).eq."CO         ")  then
  !PRINT *, 'species_name=', species_name(reactant_1_idx(J))
  !call exit()
  STICK=.5*(1-TANH(8.0d-02*(actual_dust_temp-4.0d-02*BINDING_ENERGY(reactant_2_idx(J)))))
  !PRINT *, 'species_name=', species_name(reactant_1_idx(J))
  !PRINT *, "STICK=",STICK
  !print *, '8.0d-02=',8.0d-02
  !PRINT *, 'BINDING_ENERGY(reactant_2_idx(J))', BINDING_ENERGY(reactant_2_idx(J))
  !call exit()
 endif
! PRINT *, "STICK=",STICK
 if (species_name(reactant_1_idx(J)).eq."O2         ")  then
  STICK=.5*(1-TANH(1.7d-01*(actual_dust_temp-4.2d-02*BINDING_ENERGY(reactant_2_idx(J)))))
 endif
! PRINT *, "STICK=",STICK
 if (species_name(reactant_1_idx(J)).eq."CH4         ")  then
  STICK=.5*(1-TANH(1.8d-01*(actual_dust_temp-4.5d-02*BINDING_ENERGY(reactant_2_idx(J)))))
 endif
! PRINT *, "STICK=",STICK
 if (species_name(reactant_1_idx(J)).eq."CO2         ")  then
  STICK=.5*(1-TANH(8.2d-02*(actual_dust_temp-4.4d-02*BINDING_ENERGY(reactant_2_idx(J)))))
 endif
! PRINT *, "STICK=",STICK
        
  
       
ACC_RATES_PREFACTOR(reactant_1_idx(J)) = COND*STICK/SQRT(SPECIES_MASS(reactant_1_idx(J)))
!PRINT *, 'ACC_RATES_PREFACTOR(reactant_1_idx(J)) for reaction:',(reactant_1_idx(J)),'is=', ACC_RATES_PREFACTOR(reactant_1_idx(J))
end subroutine sticking_special_cases

! ======================================================================

subroutine sticking_special_cases_Chemisorption(J)

use global_variables

implicit none

! Inputs
integer, intent(in) :: J !<[in] index of a given reaction

! Local
real(double_precision) :: cond
real(double_precision) :: STICK2,FAC1,FAC2,FAC3


STICK2=0.d0
FAC1=0.d0
FAC2=0.d0
FAC3=0.d0

cond=PI*grain_radius*grain_radius*SQRT(8.0d0*K_B/PI/AMU)
!print *, 'cond=', cond

   IF (actual_gas_temp .GE. 10) THEN 

       !PRINT *, 'actual_gas_temp=',actual_gas_temp
       !PRINT *, 'actual_dust_temp=',actual_dust_temp
       !PRINT *, 'gas_temperature=',gas_temperature
       !PRINT *, 'dust_temperature=',dust_temperature
         FAC1=SQRT(.4*((actual_gas_temp+actual_dust_temp)/100.0))
        ! PRINT *, "FAC1=", FAC1
       
         FAC2= .2*(actual_gas_temp)/100.0 
         FAC3=.8*((actual_gas_temp)/ 100.0)*.8*((actual_gas_temp)/ 100.0)
         STICK2=  1/(1 + FAC1 + FAC2 + FAC3)
        ! PRINT *, "STICK2=", STICK2
       
        
        ELSE ! ( actual_gas_temp .LT. 100 ))
        STICK2=0.d0
         
   
   END IF

      ! print *, "stick global=", STICK_GLOBAL, "stick=", STICK, "stick2 global=", STICK2_GLOBAL, "stick2=", STICK2
      ! call exit()
    ! PRINT *, "FAC1=",FAC1,"FAC2=",FAC2, "FAC3=",FAC3, "STICK2=",STICK2
     !ACC_RATES_PREFACTOR2(I)=COND*STICK2/SQRT(SPECIES_MASS(I))
   !  ACC_RATES_PREFACTOR(I)=COND*STICK/SQRT(SPECIES_MASS(I))
   !  PRINT *, 'ACC_RATES_PREFACTOR(I) for reaction:',I,'is=', ACC_RATES_PREFACTOR(I)

!ACC_RATES_PREFACTOR(reactant_1_idx(J)) = COND*STICK/SQRT(SPECIES_MASS(reactant_1_idx(J)))
!PRINT *, 'ACC_RATES_PREFACTOR(reactant_1_idx(J)) for reaction:',(reactant_1_idx(J)),'is=', ACC_RATES_PREFACTOR(reactant_1_idx(J))
ACC_RATES_PREFACTOR2(reactant_1_idx(J))= COND*STICK2/SQRT(SPECIES_MASS(reactant_1_idx(J)))
!PRINT *, 'ACC_RATES_PREFACTOR2(reactant_1_idx(J)) for reaction:',(reactant_1_idx(J)),'is=', ACC_RATES_PREFACTOR2(reactant_1_idx(J))
end subroutine sticking_special_cases_Chemisorption

! ==================================================================

end module ode_solver
