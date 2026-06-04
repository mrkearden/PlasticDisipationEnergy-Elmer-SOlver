!------------------------------------------------------------------------------
  SUBROUTINE plastic(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
!------------------------------------------------------------------------------
    USE Types
  USE DefUtils
    IMPLICIT NONE

    REAL(KIND=dp), INTENT(INOUT) :: STRESS(NTENS)
    ! Requirement for Elmer: At the time of calling the Cauchy stress T_n before
    ! the time/load increment is given
    ! Requirement for umat:  The stress T_{n+1}^{(k)} corresponding to the 
    ! current approximation of the strain increment (DSTRAN) must be returned. 
    ! If the strain increment is defined to be zero in the beginning of the
    ! nonlinear iteration, Elmer will generate a candidate for the strain increment
    ! by assuming purely elastic increment characterized by DDSDDE.

    REAL(KIND=dp), INTENT(INOUT) :: STATEV(NSTATEV)
    ! Requirement for Elmer: The state variables Q_n as specified at the 
    ! previous time/load level for converged solution are given.
    ! Requirement for umat:  The state variables Q_{n+1}^{(k)} corresponding to 
    ! the current approximation of the strain increment must be returned. If 
    ! convergence is attained, these values will be saved and associated with the 
    ! converged solution (cf. the input values)

    REAL(KIND=dp), INTENT(OUT) :: DDSDDE(NTENS,NTENS)
    ! The derivative of (Cauchy) stress response function with respect to the 
    ! strain evaluated for the current approximation must be returned

    REAL(KIND=dp), INTENT(INOUT) :: SSE, SPD, SCD
    ! Requirement for Elmer: Provide specific strain energy (sse), plastic 
    ! dissipation (spd) and creep dissipation (scd) at the previous time/load 
    ! level (these are supposed to be declared to be state variables)
    ! Requirement for umat:  The values of the energy variables corresponding to 
    ! the current approximation may be returned

    REAL(KIND=dp), INTENT(OUT) :: rpl
    ! The mechanical heating power (volumetric)

    REAL(KIND=dp), INTENT(OUT) :: ddsddt(NTENS), drplde(NTENS), drpldt

    REAL(KIND=dp), INTENT(IN) :: STRAN(NTENS)
    ! This gives the strains before the time/load increment.
    ! The strain can be computed from the deformation gradient, so this
    ! argument can be considered to be redundant. Elmer provides
    ! this information anyway. Abaqus assumes that the logarithmic strain 
    ! is used, but Elmer may also use other strain measures.

    REAL(KIND=dp), INTENT(IN) :: DSTRAN(NTENS)
    ! The current candidate for the strain increment to obtain the current 
    ! candidate for the stress. In principle this could be computed from the 
    ! deformation gradient; cf. the variable stran.

    REAL(KIND=dp), INTENT(IN) :: TIME(2)
    ! Both entries give time before the time/load increment (the time for the last
    ! converged solution

    REAL(KIND=dp), INTENT(IN) :: DTIME
    ! The time increment

    REAL(KIND=dp), INTENT(IN) :: TEMP
    ! Temperature before the time/load increment

    REAL(KIND=dp), INTENT(IN) :: dtemp
    ! Temperature increment associated wíth the time/load increment. Currently
    ! Elmer assumes isothermal conditions during the load increment.

    REAL(KIND=dp), INTENT(IN) :: predef(1), dpred(1)
    ! These are just dummy variables for Elmer

    CHARACTER(len=80), INTENT(IN) :: CMNAME
    ! The material model name

    INTEGER, INTENT(IN) :: NDI
    ! The number of direct stress components

    INTEGER, INTENT(IN) :: NSHR
    ! The number of the engineering shear strain components

    INTEGER, INTENT(IN) :: NTENS 
    ! The size of the array containing the stress or strain components

    INTEGER, INTENT(IN) :: NSTATEV
    ! The number of state variables associated with the material model

    REAL(KIND=dp), INTENT(IN) :: PROPS(NPROPS)
    ! An array of material constants

    INTEGER, INTENT(IN) :: NPROPS
    ! The number of the material constants

    REAL(KIND=dp), INTENT(IN) :: coords(3)
    ! The coordinates of the current point could be specified

    REAL(KIND=dp), INTENT(IN) :: drot(3,3)
    ! No support for keeping track of rigid body rotations 
    ! (the variable is initialized to the identity)

    REAL(KIND=dp), INTENT(INOUT) :: pnewdt
    ! Currently, suggesting a new size of time increment does not make any impact

    REAL(KIND=dp), INTENT(IN) :: celent
    ! The element size is not yet provided by Elmer

    REAL(KIND=dp), INTENT(IN) :: DFRGRD0(3,3)
    ! The deformation gradient before the time/load increment (at the previous 
    ! time/load level for converged solution)

    REAL(KIND=dp), INTENT(IN) :: DFRGRD1(3,3)
    ! The deformation gradient corresponding to the current approximation
    ! (cf. the return value of STRESS variable) 

    INTEGER, INTENT(IN) :: NOEL
    ! The element number

    INTEGER, INTENT(IN) :: NPT
    ! The integration point number

    INTEGER, INTENT(IN) :: layer, kspt, kstep, kinc
    ! kstep and kinc could be provided to give information on the incrementation
    ! procedure
!------------------------------------------------------------------------------
    ! Local variables:
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt,beta
  LOGICAL :: Found
  LOGICAL :: Transient
 TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp):: density, Cp,ctime, stime
  TYPE(Variable_t),POINTER :: VarLoadWkg, VarLoadWm3
  INTEGER :: Visit = 0, nF, i, n, dim, NSIZE,totnpt, stat, j, k
 REAL(KIND=dp) :: nu, E,E2, LambdaLame, MuLame,hslope,yield
real(kind=dp), allocatable, save :: kstress(:,:)
real(kind=dp), allocatable, save :: elstrain(:,:)
 real(kind=dp),allocatable,save :: plstrain(:,:)
 real(kind=dp),allocatable,save :: dlstrain(:,:)
 real(kind=dp),allocatable,save ::  pwi(:,:)
real(kind=dp),allocatable, save:: cumwork(:,:)
 real(kind=dp),allocatable,save :: heatinc(:,:)
 real(kind=dp),allocatable,save ::  htot(:)
 real(kind=dp),allocatable,save ::  cumheat(:)
 real (kind=dp)::pstress(ntens),pstran(ntens)
 real(kind=dp),allocatable,save :: Trise(:)
 real(kind=dp),allocatable,save :: cumT(:)
 INTEGER::kstp=0,totel,eout
  CHARACTER(LEN=1) :: c=',',header
  TYPE(ValueList_t), POINTER :: Material, SolverParams
  TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
SAVE Visit, stat,stime

    ! Get Young's modulus and the Poisson ratio:

!+++++++++++++++++++++++++++++++++++++
IF (VISIT .EQ. 0) THEN

 ! Mesh => Solver % Mesh
 ! dim = CoordinateSystemDimension()
!  SolverParams => GetSolverParams()
     open(unit=91,file='PDEOutput.csv',status='unknown')
     write(91,*) 'time',c,'noel',c,'npt',c,'ntens',c,'stress',c,'strain'&
,c,'elstrain',c,'plstrain',c,'dlstrain',c,'pwi','c,cumwork',c,'heatinc'&
,c,'cumheat',c,'Trise',c,'cumT'
close(91)
visit = 1
stat = 0
stime = 0
totel = 0
totnpt = 0
    allocate(kstress(10000,ntens))
    allocate(plstrain(10000,ntens))
    allocate(elstrain(10000,ntens))
    allocate(dlstrain(10000,ntens))
    allocate(pwi(10000,ntens))
    allocate(cumwork(10000,ntens))
    allocate(heatinc(10000,ntens))
    allocate(htot(10000))
    allocate(cumheat(10000))
    allocate(Trise(10000))
    allocate(cumT(10000))

end if
!++++++++++++++++++++++++++++++++++++

    E = Props(1)
    nu = Props(2)
    hslope = Props(3)
    yield = Props(4)
    density = Props(5)
    Cp = Props(6)
    beta = Props(7)
!    Totnpt = Props(8)
    eout = Props(8)

ctime =time(1) - stime
stime = time(1)
if (ctime .eq. 0.0) then
 if (noel.gt.totel) totel=totel+1
if (npt.gt.totnpt) totnpt=totnpt+1
end if

if (ctime .gt. 0.0) then
stat = 1
kstp = kstp + 1
end if

E2 = E

! Check stress for yield

do i=1,ntens

     if (abs(stress(i)) .gt. yield) then
       E = hslope

    end if
end do


    LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    MuLame = E / (2.0d0 * (1.0d0 + nu))

    ddsdde = 0.0d0
    ddsdde(1:ndi,1:ndi) = LambdaLame
    DO i=1,ntens
      ddsdde(i,i) = ddsdde(i,i) + MuLame
    END DO
    DO i=1,ndi
      ddsdde(i,i) = ddsdde(i,i) + MuLame
    END DO
    !
    pstress = stress
    pstran = stran
    stress = stress + MATMUL(ddsdde,dstran)
   kstress(kstep,:) = stress
if (stat.eq.1) then

   htot(kstp) = 0.0
   do i =1,ntens 
   elstrain(kstp,i) = pstress(i)/E2
if (elstrain(kstp,i).gt.pstran(i)) elstrain(kstp,i) = pstran(i)
   plstrain(kstp,i) = pstran(i)-elstrain(kstp,i)
    if (plstrain(kstp,i).lt.0.0) plstrain(kstp,i)=0.0
    dlstrain(kstp,i) = plstrain(kstp,i)-plstrain(kstp-1,i)
    if (dlstrain(kstp,i).lt.0.0) dlstrain(kstp,i)=0.0
  pwi(kstp,i) = dlstrain(kstp,i)* (kstress(kstep-1,i) + kstress(kstep,i)) 
  cumwork(kstp,i) = cumwork(kstp-1,i) + pwi(kstp,i)
  heatinc(kstp,i) = pwi(kstp,i) * 1000000.0 * beta
  end do

   do i =1,ntens 
      htot(kstp) = htot(kstp) + heatinc(kstp,i)
   end do

    cumheat(kstp) = cumheat(kstp-1) + htot(kstp)
    Trise(kstp) = htot(kstp)/(density*Cp)
    cumT(kstp) = cumT(kstp-1) + Trise(kstp)

if (noel.eq.eout) then
   do i=1,ntens
   open(unit=91,file='PDEOutput.csv',status='unknown',position='append')
    write(91,fmt=100) time(1),c,noel,c,npt,c,i,c,pstress(i),c,pstran(i),c&
,elstrain(kstp,i),c,plstrain(kstp,i),c,dlstrain(kstp,i),c,pwi(kstp,i),c,cumwork(kstp,i),c&
,htot(kstp),c,cumheat(kstp),c,Trise(kstp),c,cumT(kstp)
   close(91)
   end do
end if
100 FORMAT (F8.6,A1,3(i8,A1),10(E13.5,A1),E13.5)
  if ((noel .eq.totel).and.(npt.eq.totnpt)) then 
    stat=0
    ctime = 0.0
  end if

endif

!------------------------------------------------------------------------------
  END SUBROUTINE plastic
!------------------------------------------------------------------------------
