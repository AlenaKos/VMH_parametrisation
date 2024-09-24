program VMN_parcel_ensemble
    !parcel model for ice nucleation forced by artificial GW term sampled from the ICON/MS-GWaM results
    !double moment scheme integration and parameterisation are available 

    implicit none

    ! define constants and parameters for the run 
    INTEGER :: unit_nr, i_write, Nt
    INTEGER :: t1, t2, clock_rate, clock_max ! timers 

    REAL*8 :: Tf = 60000, t = 0.0 !,& !final time and start tie of the calculation 
    REAL*8, dimension(:), allocatable :: yn, tn, vn!, ksi_2_y, ksi_2_v
    REAL*8, dimension(:,:), allocatable :: y_ice
    LOGICAL :: exception_flag !method_switch, three_vars_switch, two_vars_switch, test_vm, exception_flag

    REAL*8 :: m0, m_to_m0 
    REAL*8 ::  n_max, n_aprox_when0
    INTEGER :: Approach_case
    REAL*8 :: t_start, t_final, dt
    REAL*8, dimension(3) :: In_conditions, res
    REAL*8, dimension(:,:), allocatable :: rr, r_writing
    REAL*8 ::  n_rand, m_rand, w00_msgwam
    CHARACTER(len = 150) :: filename_out,filename_main, filename_Sq
    REAL*8, dimension(:), allocatable :: w_hat_rand, omega_rand, phi_rand, a_rand, &
    mom_flux_out, w00_out

    ! define parameters for ensemble calcultaion. 
    ! number of tests - number_of_tests, starting number of test - nteststart,
    !  number of test for detailed output - n_out, number of nucleation events -n_nucleations, 
    ! number of gravity waves - n_GWs

    ! cases for constructing the fit: 1GW, 2GW, 5Gw, 6GW, 10GW and all this cases with the present background updraft
    ! cases for validation: 3GW, 8GW, 14GW
    ! cases for validation with different time steps: 6GW with w0 background updraft

    ! number_of_tests is usually =< n_datasize, 
    INTEGER :: nteststart = 2, number_of_tests = 500000, n_tests, n_out =58451, n_nucleations=0, n_GWs = 6 ! used for verefication 
    
    !parameters for GW and Ice nucleation constants
    REAL*8, parameter :: &
    gg = 9.8,&
    Li = 2.8e6,&
    Cp = 1005,&
    Rv = 461,&
    PI=4.D0*DATAN(1.D0) ,&
    J = 4.9 * 10**4     ,& ! 1/kg /s
    B = 337.0           ,&
    Sc = 1.5            ,&
    phi =0              ,& 
    Temp = 210.0        ,&
    epsilon0 = 0.62     ,&
    m_mean = 1.0e-12    ,&
    m_min = 1e-16       ,& 
    p_00 = 30000        ,& ! Pa
    p_si = 1            ,&!Pa for T_00 = 210, generally p_si(T)
    eps_small = 0.01          !small parameter for case 6, if eps_small = 1 -> case 3

    ! define bounds for initial conditions 
    REAL*8, parameter :: &
    n_up_bound = log(1.0e7), &
    n_low_bound = log(1.0e-4), &
    m_up_bound = log(1.0e-12),&
    m_low_bound = log(1.0e-16),&
    w_up_bound = 0.08,&
    w_low_bound = 0.0,&
    ! om_up_bound = 0.001,&
    ! om_low_bond = 0.004,&
    om_up_bound = 0.02,&
    om_low_bond = 1.0e-4,&
    phi_up_bound = PI,&
    phi_low_bound = 0
    INTEGER :: nnn, n_datasize
    INTEGER,allocatable :: seed(:)
    REAL*8 :: rhouw, maxmomflux = 1.0*1e-3, & ! set the maximum momentum flux allowed for superposition of GWs
     rhouw_spgw
    REAL*8, dimension(100) :: random_number_IC
    REAL*8,dimension(:),allocatable :: w_sample, omega_sample , vertvel_ls_sample
    character(len = 200), dimension(3) :: filename_Ftdata


    ! Data for forcing construction in the initial data set (used for the fit)
    !define path to the data from ICON
    filename_Ftdata(1)='../icon_visualisation/w_hat_18sample.dat'
    filename_Ftdata(2) = '../icon_visualisation/omega_18sample.dat'
    filename_Ftdata(3) = '../icon_visualisation/vertvel_18sample.dat'
    n_datasize=  588574!593600

    !Data for validation Data sets 
    ! filename_Ftdata(1)='../icon_visualisation/w_hat_21sample.dat'
    ! filename_Ftdata(2) = '../icon_visualisation/omega_21sample.dat'
    ! filename_Ftdata(3) = '../icon_visualisation/vertvel_21sample.dat'
    ! n_datasize=  584513 ! fixed length of the data from the ICON, R2B5 reolution and slice from h=8 to h=14 km
    

    allocate(w_sample(n_datasize))
    allocate(omega_sample(n_datasize))
    allocate(vertvel_ls_sample(n_datasize))

!    read data from msgwam for w_hat and omega 
    CALL outputdata(w_sample, omega_sample, vertvel_ls_sample,n_datasize,filename_Ftdata)
    ! print*, w_sample(30523+1), omega_sample(30523+1) ! used previously to check reading procedures
   

    call random_seed(size=nnn)
    allocate(seed(nnn))
    ! seed = 1  !for DS construction
    ! seed = 500  !
    seed = 123456789    ! for test of the parameterisation
    call random_seed(put=seed)

    allocate(w_hat_rand(n_GWs))
    allocate(omega_rand(n_GWs))
    allocate(phi_rand(n_GWs))
    allocate(a_rand(n_GWs))
    allocate(rr(14+3*n_GWs-3,number_of_tests))
    allocate(mom_flux_out(number_of_tests))
    allocate(w00_out(number_of_tests))
!Randomly choose initial conditions for n, q(or m) and forcing term: w, omega, phi

!loop over different tests

    !Initial conditions (n, S, q)
    In_conditions(2) = 1.4

    !Parameters for the run
    t_start = 0

    ! time of the intergration
    ! t_final = 6*500
    t_final = 9*500

    ! timestep
    ! dt = 30.0
    ! dt = 120.0
    ! dt = 60.0
    dt = 1.0

    ! Approach_case selection:
    ! 1 - constant mass reduced system
    ! 2 - constant mass parameterisation
    ! 3 - variable mass reduced system
    ! 4 - variable mass full system  
    ! 5 - parametrisation with extension for the mean mass variability

    Approach_case =4
    
    ! filename_out = "../partcel_mod_dat/VM_specific_out_test1.dat"
    ! filename_main = '../partcel_mod_dat/VM_dataset_out_dataset_test1.dat'  
    ! filename_Sq = '../partcel_mod_dat/VM_specific_Sq_test1.dat'

    filename_out = "../partcel_mod_dat/test1.dat"
    filename_main = '../partcel_mod_dat/test2.dat'  
    filename_Sq = '../partcel_mod_dat/test3.dat'


    Nt = (Tf-t_start) / dt+1

    ! allocate our solution arrays given the number of time steps
    allocate (yn(Nt))
    allocate (vn(Nt))
    allocate (tn(Nt+1))
    allocate(y_ice(3,Nt+1))

    ! set up the initial conditions.
    yn(1) = 1.0
    vn(1) = 0.0
    !initial conditions for LV model
    tn(1) = t_start


   IF (y_ice(1,1)>0 ) THEN
    m0= y_ice(3,1)/y_ice(1,1) 
ELSE
    m0 = m_mean
ENDIF

exception_flag = .TRUE.
DO n_tests=nteststart,number_of_tests
    !choose random IC for n, m(or q) and forcing term: w_hat, omega, phi
    call random_number(random_number_IC)
    ! print*, random_number_IC
    n_rand = exp(n_low_bound + random_number_IC(1) * (n_up_bound-n_low_bound))
    m_rand = exp(m_low_bound + random_number_IC(2) * (m_up_bound-m_low_bound))
    w_hat_rand = w_low_bound + random_number_IC(3:3+n_GWs-1) * (w_up_bound-w_low_bound)
    omega_rand = om_low_bond + random_number_IC(3+n_GWs:3+2*n_GWs-1) * (om_up_bound-om_low_bond)
    phi_rand = phi_low_bound + random_number_IC(3+2*n_GWs:3+3*n_GWs-1) * (phi_up_bound-phi_low_bound)
 

    rhouw_spgw = sum(p_00/287.0/Temp*&
    ((0.02**2-omega_rand(2:n_GWs)**2)/(omega_rand(2:n_GWs)**2-(1.0e-4)**2))**(0.5)*&
    w_hat_rand(2:n_GWs)**2)
    w_hat_rand(2:n_GWs) = w_hat_rand(2:n_GWs) *(maxmomflux/((rhouw_spgw)))**(0.5)


    w_hat_rand(1) = w_sample(n_tests)
    ! w_hat_rand(1) = 0 ! for cases without background updraft
    omega_rand(1) = omega_sample(n_tests)
    w00_msgwam = (abs(vertvel_ls_sample(n_tests)) )* gg * Li /Cp/Rv/Temp**2 !rescaled large scale vertical velocity
    
    ! calculate the momentum flux for created GW superposition 
    rhouw =sum(p_00/287.0/210.0*&
        ((0.02**2-omega_rand**2)/(omega_rand**2-(1.0e-4)**2))**(0.5)*w_hat_rand**2)
    ! print*, "calculated momentum flux: ", rhouw
    mom_flux_out(n_tests) = rhouw
    w00_out(n_tests) = w00_msgwam

    ! some checks and warnings for uploaded forcing term values
    IF (w_hat_rand(1)>1) THEN
        print*,"w hat > 1"
    ENDIF
    IF (rhouw>1) THEN
        print*, "Warning: momentum flux is too large", rhouw 
        ! EXIT
    ENDIF
! check for max w hat
    IF (maxval(w_hat_rand)>w_hat_rand(1) .AND. (maxval(w_hat_rand)>1)) THEN
        print*, maxval(w_hat_rand),w_hat_rand(1)
        ! EXIT
    ENDIF

    !rescale amplitude:
    a_rand = w_hat_rand * gg * Li /Cp/Rv/Temp**2

    In_conditions(1) = n_rand
    In_conditions(3) = In_conditions(1) *m_rand

    ! -------------------------------------------------------------------------------------
    ! solve the system of equations
    rr(:,n_tests) = solve_nucleation(Approach_case, &
    In_conditions,t_start,t_final,dt,a_rand,omega_rand,phi_rand,w00_msgwam,filename_out,n_tests)

    IF (n_tests>n_out+1) EXIT
    ! IF (n_nucleations>50) EXIT
ENDDO

 print*, 'new cleaned version is used'

 print*, 'number of nucleations:', n_nucleations, ' in ', number_of_tests, ' number of test'

 !write final output for ensamble calculations


 open(newunit = unit_nr, file = filename_main, status = 'replace')  
 do i_write = 1,number_of_tests  
    IF (rr(1,i_write)>0) THEN
     write(unit = unit_nr, fmt = *) rr(:,i_write), mom_flux_out(i_write), w00_out(i_write)
    ENDIF
 end do  

 close(unit = unit_nr) 
!  deallocate(seed)

!  contained functions
contains

!main function to execute the integration with prescribed approach, IC and time

FUNCTION solve_nucleation(Apr_c, In_conditions,t_start,t_final,dt,A_loc,om_loc,phi_loc,w00,filename_out,n_tests) result(res_hn)
    implicit none
    INTEGER, INTENT(IN) :: Apr_c, n_tests
    REAL*8,  INTENT(IN) :: t_start, t_final,dt
    REAL*8, dimension(3), INTENT(IN)  :: In_conditions
    CHARACTER(len = 150),INTENT(IN)  :: filename_out
    REAL*8, dimension(14+3*n_GWs-3) :: res_hn !n, S, q
    !parameters for calculations 
    REAL*8, allocatable :: q_1_fs(:), q_2_fs(:), &
    q_3_fs(:), Phi_1_fs(:), Phi_2_fs(:), y_in_local(:), &
    y_ice_2(:), y_Sq(:), y_Sq_loc(:)
    integer :: Nt, n, n_sparce
    real*8, dimension(:,:), allocatable :: y_ice, y_ice_short,YSq_out   
    real*8, dimension(:), allocatable :: tn
    real*8 :: t_time, yn_in, y_in
    real*8 :: m0, Temp = 210.0
    INTEGER :: ind_t0, ngw
    real*8 :: t0_loc, Ft0_loc, m0_loc, m0_recalc_loc, &
    n_post_loc, qt0_loc,F_mult, w00, q_post_loc, t0_loc_test
    real*8, dimension(n_GWs),INTENT(IN)  :: A_loc, om_loc, phi_loc
    real*8,dimension(2) :: Phi_S1, Phi_S2, q1S, q2S, q3S
    REAL*8, parameter :: &
    CMC = 4.3 * 10.0**(-8.)
    INTEGER :: jj
    REAL*8 :: rev_t0, m0_apr
    rev_t0 = 1/dt

    ! Initialaze outout for new test 
    t0_loc=0
    t0_loc_test = 0
    Ft0_loc=0
    qt0_loc = 0
    m0_loc=0
    m0_recalc_loc =0
    n_post_loc = 0

    t_time = t_start

    Nt = (t_final-t_start) / dt+1
    n_sparce = 10
    ! print*, n_GWs

    allocate (tn(Nt+1))
    allocate(y_ice(3,Nt+1))
    allocate(YSq_out(3,Nt+1))

    allocate(y_ice_short(4,INT((Nt+1)/n_sparce)+1))

    y_ice(1,1) = In_conditions(1)
    y_ice(2,1) = In_conditions(2)
    y_ice(3,1) = In_conditions(3)
    tn(1) = t_start
    SELECT CASE (Apr_c)
        ! selection is made for system of the equations used for ice physics
        ! Apr_c correspond to Approach_case

    CASE(1)
        !Constant mass case, full system
        print*, 'full system with constant mass (n, S, q)'
        ! working 3 vars reduced system 
        allocate(q_1_fs(3))
        allocate(q_2_fs(3))
        allocate(q_3_fs(3))
        allocate(Phi_1_fs(3))
        allocate(Phi_2_fs(3))
        allocate(y_in_local(3))
        allocate(y_ice_2(3))
        allocate(y_Sq(2))
        allocate(y_Sq_loc(2))
        y_Sq(1) =y_ice(2,1)
       y_Sq(2) = y_ice(3,1)
       yn_in = y_ice(1,1)


        IF (y_ice(1,1)>0 ) THEN
            m0= y_ice(3,1)/y_ice(1,1) 
        ELSE
            m0 = m_mean
        ENDIF
        y_in_local = y_ice(:,1)
        y_ice_2 =  y_ice(:,1)
        y_ice_short(1:3,1) = y_ice(:,1)
        y_ice_short(4,1) = t_time

        DO n = 1,Nt
            t_time = t_time+dt
        
            q_1_fs = dt * ice_nucl_full_dimentional_wT_testSPGW_constM(y_ice_2,tn(n),&
            Temp,m0,A_loc, om_loc, phi_loc)
            Phi_1_fs = y_ice_2+ q_1_fs / 3
            q_2_fs = dt * ice_nucl_full_dimentional_wT_testSPGW_constM(Phi_1_fs,tn(n)+dt/3.0,&
            Temp,m0,A_loc, om_loc, phi_loc) - 5 *q_1_fs / 9
            Phi_2_fs = Phi_1_fs + 15 *q_2_fs/16
            q_3_fs = dt * ice_nucl_full_dimentional_wT_testSPGW_constM(Phi_2_fs, tn(n)+dt*5/12,&
            Temp,&
            m0,A_loc, om_loc, phi_loc) - 153 *q_2_fs / 128
            y_in_local = Phi_2_fs + 8 * q_3_fs / 15
           
            IF (y_in_local(3)<0) THEN
                y_in_local(3) = 0
                y_in_local(1) = 0
            ENDIF
            IF (y_in_local(2)>Sc .AND. y_in_local(3) == 0 .AND. y_in_local(1)>0) THEN
                y_in_local(3) = y_in_local(1)*m_min
            ENDIF
            y_ice_2=y_in_local
            y_ice(1,n+1) = y_in_local(1)
            y_ice(2,n+1) = y_in_local(2)
            y_ice(3,n+1) = y_in_local(3)
            tn(n+1) = t_time
            F_mult=0
            DO ngw =1,n_GWs
                F_mult = F_mult + A_loc(ngw)*cos(om_loc(ngw)*t_time + phi_loc(ngw))
            ENDDO
            !detect nucleation event
            IF (t0_loc ==0 .AND. y_in_local(2)>Sc) THEN
                t0_loc = t_time
                ind_t0 = n
                
                Ft0_loc = F_mult
                qt0_loc = y_in_local(3)
                print*, 'Nucleation at time:', t0_loc, 'F(t0)=', Ft0_loc
                print*, 'case: ', n_tests
                print*, y_ice(1,n+1)
                n_nucleations=n_nucleations+1
            ELSEIF (t0_loc>0 .AND. y_in_local(2)>Sc .AND. y_ice(2,n)<Sc) THEN
                print*, 'second nucleation is detected for test case: ', n_tests, 'at the time: ', t_time
                EXIT
            ENDIF
            
            IF (MOD(n,n_sparce)==0) THEN
                y_ice_short(1:3,INT(n/n_sparce)+1) = y_ice(:,n)
                y_ice_short(4,INT(n/n_sparce)+1) = t_time
            ENDIF
        END DO
        
        IF (n_out == n_tests) THEN
        ! call write_data_3parameters_fname(Nt+1, tn(:), y_ice(1,:) , y_ice(2,:), y_ice(3,:), filename_out)
        call write_data_3parameters_fname(INT((Nt+1)/n_sparce)+1, y_ice_short(4,:), y_ice_short(1,:), &
        y_ice_short(2,:), y_ice_short(3,:),filename_out)
        print*, 'Detailed output for test number: ', n_out
        ENDIF

    
    CASE(2)
        !Constant mass case, parameterisation
        ! this part was done only for case of 1GW, to be rewritten
        print*, 'parameterisation with constant mass (n, S)'
        
!         yn_in = In_conditions(1)
!         y_in = In_conditions(2)
!         DO n = 1,Nt
!             t_time = t_time+dt
!         !calculation of parametrisation
!             q_1 = dt * S_param_dimentional_integration(y_in,tn(n),yn_in)
!             Phi_1 = y_in+ q_1 / 3
!             q_2 = dt * S_param_dimentional_integration(Phi_1,tn(n)+dt/3.0,yn_in) - 5 *q_1 / 9
!             Phi_2 = Phi_1 + 15 *q_2/16
!             q_3 = dt * S_param_dimentional_integration(Phi_2, tn(n)+dt*5.0/12.0,yn_in) - 153 *q_2 / 128
!             y_in = Phi_2 + 8 * q_3 / 15
            
!             IF ((y_in>Sc) .AND. (yn_in<n_0_parametrisation(Sc, tn(n)))) THEN
         
!                 yn_in =  2*n_0_parametrisation(Sc, tn(n))- yn_in
!             ELSE
!                 yn_in = yn_in
!             ENDIF
!             y_ice(1,n+1) = yn_in
!             y_ice(2,n+1) = y_in
!             y_ice(3,n+1) = yn_in * m_mean
!             tn(n+1) = t_time
!         END DO

! call write_data_3parameters_fname(Nt, tn(:),y_ice(1,:) , y_ice(2,:), y_ice(3,:), filename_out)

CASE(3)
!         !Variable mass, full system

        allocate(q_1_fs(3))
        allocate(q_2_fs(3))
        allocate(q_3_fs(3))
        allocate(Phi_1_fs(3))
        allocate(Phi_2_fs(3))
        allocate(y_in_local(3))
        allocate(y_ice_2(3))
        allocate(y_Sq(2))
        allocate(y_Sq_loc(2))
        y_Sq(1) =y_ice(2,1)
       y_Sq(2) = y_ice(3,1)
       yn_in = y_ice(1,1)


        IF (y_ice(1,1)>0 ) THEN
            m0= y_ice(3,1)/y_ice(1,1) 
        ELSE
            m0 = m_mean
        ENDIF
        y_in_local = y_ice(:,1)
        y_ice_2 =  y_ice(:,1)
        y_ice_short(1:3,1) = y_ice(:,1)
        y_ice_short(4,1) = t_time

        DO n = 1,Nt
            t_time = t_time+dt
        ! time integration of the reduced system

            q_1_fs = dt * ice_nucl_full_dimentional_wT_vm(y_ice_2,tn(n),Temp,m0,A_loc, om_loc, phi_loc,w00)
            Phi_1_fs = y_ice_2+ q_1_fs / 3
            q_2_fs = dt * ice_nucl_full_dimentional_wT_vm(Phi_1_fs,tn(n)+dt/3.0,Temp,m0,A_loc, om_loc, phi_loc,w00) - 5 *q_1_fs / 9
            Phi_2_fs = Phi_1_fs + 15 *q_2_fs/16
            q_3_fs = dt * ice_nucl_full_dimentional_wT_vm(Phi_2_fs, tn(n)+dt*5/12,Temp,&
            m0,A_loc, om_loc, phi_loc,w00) - 153 *q_2_fs / 128
            y_in_local = Phi_2_fs + 8 * q_3_fs / 15

            IF (y_in_local(3)<0) THEN
                y_in_local(3) = 0
        !         ! y_ice_2 = 0
                y_in_local(1) = 0
            ENDIF
            IF (y_in_local(2)>Sc .AND. y_in_local(3) == 0 .AND. y_in_local(1)>0) THEN
                y_in_local(3) = y_in_local(1)*m_min
            ENDIF
            y_ice_2=y_in_local
            y_ice(1,n+1) = y_in_local(1)
            y_ice(2,n+1) = y_in_local(2)
            y_ice(3,n+1) = y_in_local(3)

            ! additional check for n tendency and evaporation:
            IF (y_in_local(2)<0.95) THEN
                print*, 'evaporation happened, Ninit changed, therefor test is invalid'
                EXIT
            ENDIF

            tn(n+1) = t_time
            F_mult=0
            DO ngw =1,n_GWs
                F_mult = F_mult + A_loc(ngw)*cos(om_loc(ngw)*t_time + phi_loc(ngw))
            ENDDO
            F_mult = F_mult +w00
            ! F_mult = F_mult + gg * Li /Cp/Rv/Temp**2*0.1
            !detect nucleation event
            IF (t0_loc ==0 .AND. y_in_local(2)>Sc) THEN
                t0_loc = t_time
                ind_t0 = n
                
                ! DO ngw =1,n_GWs
                !     F_mult = F_mult + A_loc(ngw)*cos(om_loc(ngw)*t_time + phi_loc(ngw))
                ! ENDDO
                Ft0_loc = F_mult
                qt0_loc = y_in_local(3)
                print*, 'Nucleation at time:', t0_loc, 'F(t0)=', Ft0_loc
                print*, y_ice(1,n+1)
                print*, 'check for N grows:', y_ice(1,n+1),'>', y_ice(1,1), '?'
                print*, 'case: ', n_tests
                print*, 'IC for forcing: ', A_loc, om_loc, phi_loc
                n_nucleations=n_nucleations+1
            ELSEIF (t0_loc>0 .AND. y_in_local(2)>Sc .AND. y_ice(2,n)<Sc) THEN
                print*, 'second nucleation is detected for test case: ', n_tests, 'at the time: ', t_time
                ! EXIT
            ENDIF
            
            IF (MOD(n,n_sparce)==0) THEN
                y_ice_short(1:3,INT(n/n_sparce)+1) = y_ice(:,n)
                y_ice_short(4,INT(n/n_sparce)+1) = t_time
            ENDIF
        END DO

        IF (t0_loc>0) THEN
            n_post_loc = (y_ice(1,ind_t0+3000))
            print*, 'difference in Npost calc, as fixed N after 300 sec:', (y_ice(1,ind_t0+3000)), &
            'as local max:', maxval(y_ice(1,ind_t0:ind_t0+100*rev_t0))
            n_post_loc = maxval(y_ice(1,ind_t0-100*rev_t0:ind_t0+100*rev_t0))
            !if nucleation happened during the calculation, recalculate m0 from the system:
            m0_recalc_loc = (2*Sc*Ft0_loc /((n_post_loc)+y_ice(1,1))/(CMC  /epsilon0*Temp) /(Sc-1))**3
            print*, 'recalculated mass:', m0_recalc_loc, 'N post', (n_post_loc), 'n pre', (y_ice(1,1))
            ! n_post_loc = (y_ice(1,ind_t0+3000))

            t_time = t_start
            tn(1) = t_start
            ! YSq_out(1) = 
            YSq_out(1,1) = y_Sq(1)
            YSq_out(2,1) = y_Sq(2)
            YSq_out(3,1) = t_start
            IF (y_ice(1,1)>0 ) THEN
                    m0= y_ice(3,1)/y_ice(1,1) 
            ELSE
                    m0 = m_mean
            ENDIF
        !calculate ds, dq
            DO n = 1,Nt
                t_time = t_time+dt
                !integration of S, q independently from n
                q1S = dt * Sq_parametrisation_vm(y_Sq,tn(n),Temp,m0,yn_in,A_loc,om_loc,phi_loc,w00)
                Phi_S1 = y_Sq+ q1S / 3
                q2S = dt * Sq_parametrisation_vm(Phi_S1,tn(n)+dt/3.0,Temp,m0,yn_in,A_loc,om_loc,phi_loc,w00) - 5 *q1S / 9
                Phi_S2 = Phi_S1 + 15 *q2S/16
                q3S = dt * Sq_parametrisation_vm(Phi_S2, tn(n)+dt*5/12,Temp,m0,yn_in,A_loc,om_loc,phi_loc,w00) - 153 *q2S / 128
                y_Sq_loc = Phi_S2 + 8 * q3S / 15
                y_Sq = y_Sq_loc
                
                tn(n+1) = t_time
                YSq_out(1,n+1) = y_Sq(1)
                YSq_out(2,n+1) = y_Sq(2)
                YSq_out(3,n+1) = t_time
                IF (y_Sq(2)<0) THEN
                    ! y_Sq(2) = 0
                    y_Sq(2) = m_min*yn_in
                ENDIF
                IF (y_Sq(1)>Sc ) THEN
                    m0_loc = y_Sq(2)/yn_in
                    print*, 'nucleation at t:', t_time, ' m0 from ds/dt, dq/dt:', y_Sq(2)/yn_in
                EXIT
                ENDIF
            END DO

        ENDIF
            
 
        
        IF (n_out == n_tests) THEN
        ! call write_data_3parameters_fname(Nt+1, tn(:), y_ice(1,:) , y_ice(2,:), y_ice(3,:), filename_out)
        call write_data_3parameters_fname(INT((Nt+1)/n_sparce)+1, y_ice_short(4,:), y_ice_short(1,:), &
        y_ice_short(2,:), y_ice_short(3,:),filename_out)

        call write_data_3parameters_fname((Nt+1), YSq_out(1,:), &
        YSq_out(1,:) , YSq_out(2,:), YSq_out(3,:), &
        filename_Sq)
        print*, 'Detailed output for test number: ', n_out
        ENDIF

        
        
CASE(4)
!         !Variable mass, parameterisation original via algebraic equation
        ! print*, 'parameterisation with variable mass, corrected in the deposition term'
                    allocate(q_1_fs(3))
                    allocate(q_2_fs(3))
                    allocate(q_3_fs(3))
                    allocate(Phi_1_fs(3))
                    allocate(Phi_2_fs(3))
                    allocate(y_in_local(3))
                    allocate(y_ice_2(3))
                    allocate(y_Sq(2))
                    allocate(y_Sq_loc(2))
                    y_Sq(1) =y_ice(2,1)
                   y_Sq(2) = y_ice(3,1)
                   yn_in = y_ice(1,1)
            
            
                    IF (y_ice(1,1)>0 ) THEN
                        m0= y_ice(3,1)/y_ice(1,1) 
                    ELSE
                        m0 = m_mean
                    ENDIF
                    y_in_local = y_ice(:,1)
                    y_ice_2 =  y_ice(:,1)
                    y_ice_short(1:3,1) = y_ice(:,1)
                    y_ice_short(4,1) = t_time
        
                    
                    DO n = 1,Nt ! time loop
                        t_time = t_time+dt
                    
                        q_1_fs = dt * ice_nucl_full_dimentional_wT_vm_FS(y_ice_2,tn(n),Temp,m0,A_loc, om_loc,&
                         phi_loc,w00)
                        Phi_1_fs = y_ice_2+ q_1_fs / 3
                        q_2_fs = dt * ice_nucl_full_dimentional_wT_vm_FS(Phi_1_fs,tn(n)+dt/3.0,Temp,m0,A_loc,&
                         om_loc, phi_loc,w00) - 5 *q_1_fs / 9
                        Phi_2_fs = Phi_1_fs + 15 *q_2_fs/16
                        q_3_fs = dt * ice_nucl_full_dimentional_wT_vm_FS(Phi_2_fs, tn(n)+dt*5/12,Temp,&
                        m0,A_loc, om_loc, phi_loc,w00) - 153 *q_2_fs / 128
                        y_in_local = Phi_2_fs + 8 * q_3_fs / 15
            
                        
                        IF (y_in_local(3)<0) THEN
                            y_in_local(3) = 0
                    !         ! y_ice_2 = 0
                            y_in_local(1) = 0
                        ENDIF
                        IF (y_in_local(2)>Sc .AND. y_in_local(3) == 0 .AND. y_in_local(1)>0) THEN
                            y_in_local(3) = y_in_local(1)*m_min
                        ENDIF
                        y_ice_2=y_in_local
                        y_ice(1,n+1) = y_in_local(1)
                        y_ice(2,n+1) = y_in_local(2)
                        y_ice(3,n+1) = y_in_local(3)
                        tn(n+1) = t_time
                        F_mult=0
                        DO ngw =1,n_GWs
                            F_mult = F_mult + A_loc(ngw)*cos(om_loc(ngw)*t_time + phi_loc(ngw))
                        ENDDO
                        F_mult = F_mult + w00
                        !detect nucleation event
                        IF (t0_loc ==0 .AND. y_in_local(2)>1.5) THEN
                            t0_loc = t_time
                            ind_t0 = n
                            
                            Ft0_loc = F_mult
                            print*, 'Ft0 from full system', Ft0_loc
                            qt0_loc = y_in_local(3)
                            ! print*, 'Nucleation at time:', t0_loc, 'F(t0)=', Ft0_loc
                            print*, y_ice(1,n+1)
                            print*, 'N0 local =',  y_in_local(1)

                            print*, 'local q0 =', y_in_local(3), 'm0=', y_in_local(3)/y_in_local(1)
                            n_nucleations=n_nucleations+1
                        ! ELSEIF(y_in_local(2)>1.485) THEN
                        !     t0_loc = t_time
                        ! !     ! ind_t0 = t_time
                        !     ind_t0 = n-1
                        !     Ft0_loc = F_mult
                        ELSEIF (t0_loc>0 .AND. y_in_local(2)>1.475 .AND. y_ice(2,n)<1.475) THEN
                            print*, 'second nucleation is detected for test case: ', n_tests, 'at the time: ', t_time
                            ! EXIT
                        ENDIF
                        
                        IF (MOD(n,n_sparce)==0) THEN
                            y_ice_short(1:3,INT(n/n_sparce)+1) = y_ice(:,n)
                            y_ice_short(4,INT(n/n_sparce)+1) = t_time
                        ENDIF
                    END DO ! end time loop
            
                    IF (t0_loc>0) THEN
                        n_post_loc = (y_ice(1,ind_t0+30*dt))
                        n_post_loc = maxval(y_ice(1,ind_t0-100*rev_t0:ind_t0+20*dt*rev_t0))
                        jj = min(ind_t0 +50,int(t_final/dt) -100)
                        DO WHILE (y_ice(1,jj)<n_post_loc)
                            q_post_loc = y_ice(3,jj)
                            jj = jj+1
                        ENDDO
                        print*, "local max ",n_post_loc, q_post_loc
                        ! jj = ind_t0 +50
                        jj = min(ind_t0 +50,int(t_final/dt) -100)
                        DO WHILE ((y_ice(1,jj)-y_ice(1,jj-1)>0.1))
                            n_post_loc = (y_ice(1,jj))
                            q_post_loc = y_ice(3,jj)
                            jj = jj+1
                        ENDDO
                        jj = 0
                        print*, "dn/dt check ", n_post_loc, q_post_loc
                        print*, 'Nucleation at time:', t0_loc, 'F(t0)=', Ft0_loc
                        print*, 'case: ', n_tests
                        !if nucleation happened during the calculation, recalculate m0 from the system:
                        m0_recalc_loc = (2*Sc*Ft0_loc /((n_post_loc)+y_ice(1,1))/(CMC  /epsilon0*Temp) /(Sc-1))**3
                        print*, 'recalculated mass:', m0_recalc_loc, 'N post', (n_post_loc), 'n pre', (y_ice(1,1)), &
                        " qi_post = ", q_post_loc
                        print*, 'N0 in the middle: ', ((n_post_loc)+y_ice(1,1))/2
                       
                        t_time = t_start
                        tn(1) = t_start
                        ! YSq_out(1) = 
                        IF (y_ice(1,1)>0 ) THEN
                                m0= y_ice(3,1)/y_ice(1,1) 
                            ELSE
                                m0 = m_mean
                            ENDIF
                    !calculate ds, dq
                        DO n = 1,Nt
                            t_time = t_time+dt
                            !integration of S, q independently from n
                            q1S = dt * Sq_parametrisation_vm(y_Sq,tn(n),Temp,m0,yn_in,&
                            A_loc,om_loc,phi_loc,w00)
                            Phi_S1 = y_Sq+ q1S / 3
                            q2S = dt * Sq_parametrisation_vm(Phi_S1,tn(n)+dt/3.0,Temp,m0,yn_in, &
                            A_loc,om_loc,phi_loc,w00) - 5 *q1S / 9
                            Phi_S2 = Phi_S1 + 15 *q2S/16
                            q3S = dt * Sq_parametrisation_vm(Phi_S2, tn(n)+dt*5/12,Temp,m0,yn_in,&
                            A_loc,om_loc,phi_loc,w00) - 153 *q2S / 128
                            y_Sq_loc = Phi_S2 + 8 * q3S / 15
                            y_Sq = y_Sq_loc
                            
                            tn(n+1) = t_time
                            YSq_out(1,n+1) = y_Sq(1)
                            YSq_out(2,n+1) = y_Sq(2)
                            YSq_out(3,n+1) = t_time
                            IF (y_Sq(2)<0) THEN
                                ! y_Sq(2) = 0
                                y_Sq(2) = m_min*yn_in
                            ENDIF
                            IF (y_Sq(1)>Sc) THEN
                                m0_loc = y_Sq(2)/yn_in
                                print*, 'nucleation at t:', t_time, ' m0 from ds/dt, dq/dt:', y_Sq(2)/yn_in
                                F_mult=0
                            DO ngw =1,n_GWs
                                F_mult = F_mult + A_loc(ngw)*cos(om_loc(ngw)*t_time + phi_loc(ngw))
                            ENDDO
                            F_mult = F_mult + w00
                            IF (abs(F_mult-Ft0_loc)>0) print*, 'difference: ', F_mult, Ft0_loc
                            t0_loc_test = t_time
                            Ft0_loc = F_mult
                        
                            print*, 'Ft0 from Sq integration: ', Ft0_loc
                            print*, '---------'
                            EXIT
                            ENDIF
            
                        END DO
            
                    ENDIF

                    
                    IF (n_out == n_tests) THEN
                    ! IF (m0_recalc_loc >1e-10 .AND. t0_loc>0 .AND. t0_loc<2000) THEN
                    ! call write_data_3parameters_fname(Nt+1, tn(:), y_ice(1,:) , y_ice(2,:), y_ice(3,:), filename_out)
                    call write_data_3parameters_fname(INT((Nt+1)/n_sparce)+1, y_ice_short(4,:), y_ice_short(1,:), &
                    y_ice_short(2,:), y_ice_short(3,:),filename_out)
            
                    call write_data_3parameters_fname((Nt+1), YSq_out(1,:), &
                    YSq_out(1,:) , YSq_out(2,:), YSq_out(3,:), &
                    filename_Sq)
                    print*, 'Detailed output for test number: ', n_out
                    ENDIF
                CASE(5)
                    !variable mass parameterisation with a fix for smaller number concentrations
                    ! print*, 'parameterisation with variable mass'
                    !parametrisation 
                   allocate(q_1_fs(2))
                   allocate(q_2_fs(2))
                   allocate(q_3_fs(2))
                   allocate(Phi_1_fs(2))
                   allocate(Phi_2_fs(2))
                   allocate(y_in_local(3))
                   allocate(y_ice_2(2))
                   allocate(y_Sq(2))
                   y_ice_short(1:3,1) = y_ice(:,1)
                   
                !    y_in_local = y_ice(:,1)
                !    y_ice_2 =  y_ice(1:2,1)
                   y_Sq(1) =y_ice(2,1)
                   y_Sq(2) = y_ice(3,1)
                   yn_in = y_ice(1,1)
                   IF (y_ice(1,1)>0 ) THEN
                    m0= y_ice(3,1)/y_ice(1,1) 
                ELSE
                    m0 = m_min
                ENDIF
            
                DO n = 1,Nt
                    t_time = t_time+dt
            F_mult=0
            DO ngw =1,n_GWs
                F_mult = F_mult + A_loc(ngw)*cos(om_loc(ngw)*t_time + phi_loc(ngw))
            ENDDO
            F_mult = F_mult + w00
            !integration of S, q independently from n
            q_1_fs = dt * Sq_parametrisation_vm(y_Sq,tn(n),Temp,m0,yn_in,A_loc,om_loc,phi_loc,w00)
            Phi_1_fs = y_Sq+ q_1_fs / 3
            q_2_fs = dt * Sq_parametrisation_vm(Phi_1_fs,tn(n)+dt/3.0,Temp,m0,yn_in,A_loc,om_loc,phi_loc,w00) - 5 *q_1_fs / 9
            Phi_2_fs = Phi_1_fs + 15 *q_2_fs/16
            q_3_fs = dt * Sq_parametrisation_vm(Phi_2_fs, tn(n)+dt*5/12,Temp,m0,yn_in,A_loc,om_loc,phi_loc,w00) - 153 *q_2_fs / 128
            y_Sq_loc = Phi_2_fs + 8 * q_3_fs / 15
            y_Sq = y_Sq_loc
            
            tn(n+1) = t_time
            YSq_out(1,n+1) = y_Sq(1)
            YSq_out(2,n+1) = y_Sq(2)
            YSq_out(3,n+1) = t_time
            IF (y_Sq(2)<0) THEN
                ! y_Sq(2) = 0
                y_Sq(2) = m_min*yn_in
            ENDIF
            ! n_max = n_0_parametrisation(Sc, tn(n))
            n_max=n_0_parametrisation_simplefit(Sc,  tn(n),yn_in,Temp,A_loc,om_loc,phi_loc,w00,m0_apr)
    
            IF ((y_Sq(1)>Sc) .AND. (yn_in<n_max)) THEN
                m0_loc = y_Sq(2)/yn_in
                print*, 'nucleation at t:', t_time, ' m0 from ds/dt, dq/dt:', y_Sq(2)/yn_in
                print*, 'case: ', n_tests
                    yn_in = 2*n_0_parametrisation_simplefit(Sc,  tn(n),yn_in,Temp,A_loc,om_loc,phi_loc,w00,m0_apr)-yn_in
                    y_Sq(2) =yn_in*m0_apr
                    YSq_out(2,n+1) = yn_in*m0_apr
                    m0_recalc_loc =m0_apr
    
                    ! update output for nucleation events 
                    IF (t0_loc ==0) THEN
                        n_nucleations=n_nucleations+1
                        t0_loc = t_time
                        t0_loc_test = t_time
                        Ft0_loc = F_mult
                        ind_t0 = n
                    ENDIF
                    print*, t0_loc, Ft0_loc
                    n_post_loc = yn_in
                    
            ELSE 
                yn_in = yn_in
            ENDIF
            y_in_local(1) =yn_in
            y_ice(1,n+1) = yn_in
            y_ice(2,n+1) = y_Sq(1)
            y_ice(3,n+1) = y_Sq(2)
            tn(n+1) = t_time
        END DO !end time loop
    
        IF (t0_loc>0) THEN
            
            q_post_loc = y_ice(3,ind_t0 + 3)
            
        ENDIF
                
        IF (n_out == n_tests) THEN
            ! ! call write_data_3parameters_fname(Nt+1, tn(:), y_ice(1,:) , y_ice(2,:), y_ice(3,:), filename_out)
            ! call write_data_3parameters_fname(INT((Nt+1)/n_sparce)+1, y_ice_short(4,:), y_ice_short(1,:), &
            ! y_ice_short(2,:), y_ice_short(3,:),filename_out)
            call write_data_3parameters_fname(Nt+1, tn(:), y_ice(1,:) , y_ice(2,:), y_ice(3,:), filename_out)
        
            ! call write_data_3parameters_fname((Nt+1), YSq_out(1,:), &
            ! YSq_out(1,:) , YSq_out(2,:), YSq_out(3,:), &
            ! filename_Sq)
            print*, 'Detailed output for test number: ', n_out
            
    ENDIF  

    END SELECT

    !block for writing the time series 

    IF (t0_loc>0 .AND. t0_loc_test>0 .AND. (t0_loc_test-t0_loc)<400) THEN

    res_hn(1:3) = y_ice_short(1:3,1)
    res_hn(4:4+n_GWs-1) = A_loc
    res_hn(4+n_GWs:4+2*n_GWs-1) = om_loc
    res_hn(4+2*n_GWs:4+3*n_GWs-1) = phi_loc
    res_hn(4+3*n_GWs) = t0_loc
    res_hn(4+3*n_GWs+1) = Ft0_loc
    res_hn(4+3*n_GWs+2) = n_post_loc
    res_hn(4+3*n_GWs+3) = m0_loc
    res_hn(4+3*n_GWs+4) = m0_recalc_loc
    res_hn(4+3*n_GWs+5) = qt0_loc
    res_hn(4+3*n_GWs+6) = q_post_loc
    res_hn(4+3*n_GWs+7) = n_tests
    ELSE
        res_hn = 0
    ENDIF




END FUNCTION solve_nucleation


!Functions for right hand parts of system, S evolution equation and number concentration parameterisation


FUNCTION n_0_parametrisation_simplefit(S_res, t,n_init,Temp,A,omg,phi,w00,m0) result(n_res)
    implicit none
    REAL*8 :: S_res, t, n_init, Ft0, n_res
    REAL*8, dimension(n_GWs), INTENT(IN)  :: A, omg, phi
    REAL*8, INTENT(IN) :: Temp, w00
    REAL*8, INTENT(OUT) ::  m0
    REAL*8, dimension(2) :: nm_res
    REAL*8, parameter :: &
    CMC = 4.3 * 10.0**(-8) ,& !C0 = 4.3e-8 kg^2/3 /sK, m 0 10e-12 kg
    D_star = CMC  /epsilon0  !* n_c * T_w  
    REAL*8, parameter :: a1 =-25.325647262283194, a2=-74.5502113819399,&
    a3= 0.08935198805596141,a4= -2.9198249274795298e-05 !for corrected function with Ft0^1/3
    INTEGER :: ngw
    Ft0 =0
    DO ngw = 1,n_GWs
        Ft0 = Ft0 +A(ngw)* cos(omg(ngw) * t + phi(ngw))
    ENDDO
    Ft0 = Ft0 + w00
    m0 = exp(a1  + a2* Ft0**(1./3.)   + a3 *n_init**(1./3.) + a4*n_init*Ft0**(1./3.)) ! corrected function
    ! m0 = exp(a1  + a2* Ft0**(1./3.)   + a3 *n_init**(1./3.) + a4*n_init*Ft0**2)
    n_res =  Sc * Ft0/ (D_star * Temp *(Sc - 1))/m0**(1./3.)
    nm_res = [n_res, m0]
END FUNCTION n_0_parametrisation_simplefit


FUNCTION ice_nucl_full_dimentional_wT_vm_FS(y, t, Temp,m0,A,omg,phi,w00) result(dy)
    ! full system with superposition of GWs
    implicit none
    REAL*8, dimension(3) :: y, dy
    REAL*8, INTENT(IN) :: Temp, w00!, m0
    REAL*8 :: t, m0, mass, F_loc
    REAL*8, dimension(n_GWs), INTENT(IN)  :: A, omg, phi
    INTEGER :: ngw
    REAL*8, parameter :: &
    CMC = 4.3 * 10.0**(-8.) ,&
    qv_r = 2 * 10.0**(-5)  ,&!epsilon0 * p_si / p_00  ,&!2 * 10.0**(-5) 
    D_star = CMC  /epsilon0   

    F_loc =0
    DO  ngw = 1,n_GWs
        F_loc = F_loc +A(ngw)* cos(omg(ngw) * t + phi(ngw))
    ENDDO
    F_loc = F_loc + w00

    IF (y(1)>0 .AND. y(3)>0) THEN
        mass = y(3)/y(1)
        
        !additional boundary for small masses
        IF (mass<m_min) THEN
            mass = m_min
        ENDIF

        !RHP if mass !=0
        ! dn/dt
        dy(1) = J * exp(B * (y(2)-Sc)) ! dn/dt
        ! dS/dt
        dy(2) = - D_star* Temp * (y(2) - 1) *y(1)*mass**(1.0/3.0)+ y(2) *F_loc -(m_min)*p_00/p_si *J * exp(B * (y(2)-Sc))  !-p_00/p_si *J * exp(B * (y(2)-Sc))  !A/(Temp ** 2)  * y(2) * cos(omg * t + phi) !dS/dt
        ! dqi/dt
        dy(3) = D_star*qv_r*Temp*(y(2)-1) * y(1) *mass**(1.0/3.0) + J * exp(B * (y(2)-Sc)) *(m_min)!+ J * exp(B * (y(2)-Sc))
    ELSE
        !RHP if m=0
        mass = 0
        dy(1) = J * exp(B * (y(2)-Sc)) ! dn/dt
        dy(2) = - D_star* Temp * (y(2) - 1) *y(1)*mass**(1.0/3.0)+ y(2) *F_loc   -(m_min)*p_00/p_si *J * exp(B * (y(2)-Sc))  !-p_00/p_si*J * exp(B * (y(2)-Sc))   !A/(Temp ** 2)  * y(2) * cos(omg * t + phi) !dS/dt
        dy(3) = D_star*qv_r*Temp*(y(2)-1) * y(1) *mass**(1.0/3.0) + J * exp(B * (y(2)-Sc)) *(m_min) ! +J * exp(B * (y(2)-Sc))
    ENDIF


    END FUNCTION ice_nucl_full_dimentional_wT_vm_FS

FUNCTION ice_nucl_full_dimentional_wT_vm(y, t, Temp,m0,A,omg,phi,w00) result(dy)
    implicit none
    REAL*8, dimension(3) :: y, dy
    REAL*8, INTENT(IN) :: Temp !, m0
    REAL*8 :: t, m0, mass, F_loc, w00
    REAL*8, dimension(n_GWs), INTENT(IN)  :: A, omg, phi
    INTEGER :: ngw
    REAL*8, parameter :: &
    CMC = 4.3 * 10.0**(-8.) ,&
    qv_r = 2 * 10.0**(-5)  ,&!epsilon0 * p_si / p_00  ,&!2 * 10.0**(-5) 
    D_star = CMC  /epsilon0   

    F_loc =0
    DO ngw = 1,n_GWs
        F_loc = F_loc +A(ngw)* cos(omg(ngw) * t + phi(ngw))
    ENDDO
    F_loc = F_loc + w00
    ! F_loc = F_loc +gg * Li /Cp/Rv/Temp**2*0.1

    IF (y(1)>0 .AND. y(3)>0) THEN

        mass = y(3)/y(1)
        !additional boundary for small masses
        IF (mass<m_min) THEN
            mass = m_min
        ENDIF
        

        !RHP if mass !=0
        dy(1) = J * exp(B * (y(2)-Sc)) ! dn/dt
        dy(2) = - D_star* Temp * (y(2) - 1) *y(1)*mass**(1.0/3.0)+ y(2) *F_loc     !A/(Temp ** 2)  * y(2) * cos(omg * t + phi) !dS/dt
        dy(3) = D_star*qv_r*Temp*(y(2)-1) * y(1) *mass**(1.0/3.0)
    ELSE
        !RHP if m=0
        mass = 0
        dy(1) = J * exp(B * (y(2)-Sc)) ! dn/dt
        dy(2) = - D_star* Temp * (y(2) - 1) *y(1)*mass**(1.0/3.0)+ y(2) *F_loc      !A/(Temp ** 2)  * y(2) * cos(omg * t + phi) !dS/dt
        dy(3) = D_star*qv_r*Temp*(y(2)-1) * y(1) *mass**(1.0/3.0)   
    ENDIF
END FUNCTION ice_nucl_full_dimentional_wT_vm


FUNCTION ice_nucl_full_dimentional_wT_testSPGW_constM(y, t, Temp,m0,A,omg,phi) result(dy)
    implicit none
    REAL*8, dimension(3) :: y, dy
    REAL*8, INTENT(IN) :: Temp !, m0
    REAL*8 :: t, m0, mass, F_loc
    REAL*8, dimension(n_GWs), INTENT(IN)  :: A, omg, phi
    INTEGER :: ngw
    REAL*8, parameter :: &
    CMC = 4.3 * 10.0**(-8.) ,&
    qv_r = 2 * 10.0**(-5)  ,&!epsilon0 * p_si / p_00  ,&!2 * 10.0**(-5) 
    D_star = CMC  /epsilon0   
    
    F_loc =0
    DO ngw = 1,n_GWs
        F_loc = F_loc +A(ngw)* cos(omg(ngw) * t + phi(ngw))
    ENDDO
    
    IF (y(1)>0 .AND. y(3)>0) THEN
    
        mass = (m_mean)**(1./3.) 
    
        !RHP if mass !=0
        dy(1) = J * exp(B * (y(2)-Sc)) ! dn/dt
        dy(2) = - D_star* Temp * (y(2) - 1) *y(1)*mass**(1.0/3.0)+ y(2) *F_loc     !A/(Temp ** 2)  * y(2) * cos(omg * t + phi) !dS/dt
        dy(3) = D_star*qv_r*Temp*(y(2)-1) * y(1) *mass**(1.0/3.0)
    ELSE
           !RHP if m=0
        mass = 0
        dy(1) = J * exp(B * (y(2)-Sc)) ! dn/dt
        dy(2) = - D_star* Temp * (y(2) - 1) *y(1)*mass**(1.0/3.0)+ y(2) *F_loc      !A/(Temp ** 2)  * y(2) * cos(omg * t + phi) !dS/dt
        dy(3) = D_star*qv_r*Temp*(y(2)-1) * y(1) *mass**(1.0/3.0)   
    ENDIF
    END FUNCTION ice_nucl_full_dimentional_wT_testSPGW_constM

FUNCTION Sq_parametrisation_vm(y, t, Temp,m0,n,A,omg,phi,w00) result(dy)
    implicit none
    REAL*8, dimension(2) :: y, dy
    REAL*8, INTENT(IN) :: Temp !, m0
    REAL*8 :: t, m0, n, mass, F_loc, w00
    INTEGER :: ngw
    REAL*8, dimension(n_GWs), intent(IN) :: A, omg, phi
    REAL*8, parameter :: &
    CMC = 4.3 * 10.0**(-8.) ,&
    qv_r = 2 * 10.0**(-5.) ,&
    D = CMC  /epsilon0  
    F_loc = 0
    DO ngw = 1,n_GWs
        F_loc = F_loc +A(ngw)* cos(omg(ngw) * t + phi(ngw))
    ENDDO
    F_loc = F_loc +w00
     
    
    IF (y(2)>0 .AND. n>0) THEN    
        mass = y(2)/n
  
        !additional boundary for small masses
        IF (mass<m_min) THEN
            mass = m_min
        ENDIF
    
        dy(1) = - D* Temp * (y(1) - 1) *n*mass**(1.0/3.0)+ y(1) *F_loc     !A/(Temp ** 2)  * y(2) * cos(omg * t + phi) !dS/dt
        dy(2) = D*qv_r*Temp*(y(1)-1) * n *mass**(1.0/3.0)
    ELSE
        mass = 0
    
        dy(1) = - D* Temp * (y(1) - 1) *n*mass**(1.0/3.0)+ y(1) * F_loc     !A/(Temp ** 2)  * y(2) * cos(omg * t + phi) !dS/dt
        dy(2) = D*qv_r*Temp*(y(1)-1) * n *mass**(1.0/3.0)
    ENDIF
    
END FUNCTION Sq_parametrisation_vm


! subroutine to wrie data into files
subroutine write_data_3parameters_fname(Nt, tn,nn, ss, qq, filename)
    implicit none

    real*8, dimension(:), intent(in) :: nn, ss, qq, tn
    integer, intent(in) :: Nt
    integer :: i, unit_nr
    character(len = 150) :: filename
    ! print *, qq

    ! open(newunit = unit_nr, file = 'data_Ice_new.dat', status = 'replace')  
    open(newunit = unit_nr, file = filename, status = 'replace')  
    do i = 1,Nt  
        write(unit = unit_nr, fmt = *) tn(i), nn(i), ss(i), qq(i)
    end do  

    close(unit = unit_nr) 

end subroutine write_data_3parameters_fname

! subroutine to read data from ICON outputs 
subroutine outputdata(w, om, w_ls,n,filename_Ftdata)
    implicit none
    integer :: n != 521956 ! Number of floats to read, read all the data available anyways 
    real*8 :: data(n), omega(n), hatw(n)
    real*8 :: w(n), om(n), w_ls(n) 
    integer :: status
    integer :: unit = 10 ! Unit number for the file
    integer :: i
    character(len = 200),dimension(3) :: filename_Ftdata
    
    open(unit, file=filename_Ftdata(1), access='stream', &
    form='unformatted', status='old', action='read', iostat=status)
    if (status /= 0) then
        print *, 'Error opening file'
        stop
    end if
    
    read(unit) hatw
    
    close(unit)


    open(unit, file=filename_Ftdata(2), access='stream', &
    form='unformatted', status='old', action='read', iostat=status)
    if (status /= 0) then
        print *, 'Error opening file'
        stop
    end if
    
    read(unit) omega
    
    close(unit)

    open(unit, file=filename_Ftdata(3), access='stream', &
    form='unformatted', status='old', action='read', iostat=status)
    if (status /= 0) then
        print *, 'Error opening file'
        stop
    end if
    
    read(unit) w_ls
    
    close(unit)
    

   !  test the data correctness for bounds :
   DO i=1,n
      IF (hatw(i)>4) print*,'hat_w more than 4'
      IF (hatw(i)<0) print*, 'hat_w >0'
      IF (omega(i)>0.02 ) print*, 'omega>0.02'
      IF (omega(i)<1.0e-4) print*, 'omega<1e-4'
   ENDDO

   print*, 'check for correctness of the readed data: ',hatw(1000), omega(1050)
   print*, hatw(30523+1)
   
   w = hatw
   om = omega

   print*, "check after the read", w(30523+1)

end subroutine outputdata

end program VMN_parcel_ensemble