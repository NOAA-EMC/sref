PROGRAM lbc_perturb_wrf

!  USE module_mpiio

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  INTEGER :: io_form
  INTEGER, PARAMETER :: io_netcdf = 2

  CHARACTER (LEN=256) :: ctl_fname, p_fname, n_fname
  INTEGER :: ctl_id, p_id, n_id

  INTEGER :: rcode, istat
  INTEGER :: id_var
  INTEGER :: dyn_opt, id_dynopt, attnum
  INTEGER :: im, jm, lm, lm_ctl, im_u, jm_v
  INTEGER :: im_id, jm_id, lm_id
  INTEGER :: i,k,j,l
  REAL :: rmse_t, en_limit, en_scale
  LOGICAL :: diff_pres

  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: PRES,H,T,Q,U,V
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: PRES_p,H_p,T_p,Q_p,U_p,V_p
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: PRES_n,H_n,T_n,Q_n,U_n,V_n

  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dH, dH_ctl
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dT, dT_ctl
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dQ, dQ_ctl
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dU, dU_ctl
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dV, dV_ctl

  INTEGER :: iret

!-------------------------------------------------------------------------------

  io_form = 2

  ctl_fname="fort.31"
  p_fname="fort.32"
  n_fname="fort.33"

!
! CTL
!
  SELECT CASE ( io_form  )
  CASE ( IO_NETCDF )

      CALL check( nf_open(ctl_fname, NF_NOWRITE, ctl_id) )

      CALL check( nf_inq_dimid(ctl_id, 'west_east', im_id) )
      CALL check( nf_inq_dimlen(ctl_id, im_id, im) )
      CALL check( nf_inq_dimid(ctl_id, 'south_north', jm_id) )
      CALL check( nf_inq_dimlen(ctl_id, jm_id, jm) )
      CALL check( nf_inq_dimid(ctl_id, 'num_metgrid_levels', lm_id) )
      CALL check( nf_inq_dimlen(ctl_id, lm_id, lm) )

!d      CALL check( nf_inq_attid(ctl_id, id_dynopt, "DYN_OPT", attnum) )
      CALL check( nf_get_att_int(ctl_id, NF_GLOBAL, "DYN_OPT", dyn_opt) )

      im_u = im
      if (dyn_opt ==2) im_u = im + 1
      jm_v = jm
      if (dyn_opt ==2) jm_v = jm + 1

      lm_ctl = lm
      print *, "CTL im,jm,lm ", im,jm,lm,im_u,jm_v
      allocate (PRES(im,jm,lm))
      allocate (H(im,jm,lm))
      allocate (T(im,jm,lm))
      allocate (Q(im,jm,lm))
      allocate (U(im_u,jm,lm))
      allocate (V(im,jm_v,lm))

      allocate (dH_ctl(im,jm,lm))
      allocate (dT_ctl(im,jm,lm))
      allocate (dQ_ctl(im,jm,lm))
      allocate (dU_ctl(im_u,jm,lm))
      allocate (dV_ctl(im,jm_v,lm))

      CALL check( nf_inq_varid (ctl_id, "PRES", id_var) )
      CALL check( nf_get_var_real(ctl_id, id_var, PRES) )
      print *, "CTL pressure levels ", PRES(1,1,:)

      CALL check( nf_inq_varid (ctl_id, "GHT", id_var) )
      CALL check( nf_get_var_real(ctl_id, id_var, H) )

      CALL check( nf_inq_varid (ctl_id, "TT", id_var) )
      CALL check( nf_get_var_real(ctl_id, id_var, T) )

      CALL check( nf_inq_varid (ctl_id, "RH", id_var) )
      CALL check( nf_get_var_real(ctl_id, id_var, Q) )

      CALL check( nf_inq_varid (ctl_id, "UU", id_var) )
      CALL check( nf_get_var_real(ctl_id, id_var, U) )

      CALL check( nf_inq_varid (ctl_id, "VV", id_var) )
      CALL check( nf_get_var_real(ctl_id, id_var, V) )

      CALL check( nf_close(ctl_id) )
    !
    ! p
    !
      CALL check( nf_open(p_fname, NF_NOWRITE, p_id) )

      CALL check( nf_inq_dimid(p_id, 'west_east', im_id) )
      CALL check( nf_inq_dimlen(p_id, im_id, im) )
      CALL check( nf_inq_dimid(p_id, 'south_north', jm_id) )
      CALL check( nf_inq_dimlen(p_id, jm_id, jm) )
      CALL check( nf_inq_dimid(p_id, 'num_metgrid_levels', lm_id) )
      CALL check( nf_inq_dimlen(p_id, lm_id, lm) )

      im_u = im
      if (dyn_opt ==2) im_u = im + 1
      jm_v = jm
      if (dyn_opt ==2) jm_v = jm + 1

      print *, "MEMBER im,jm,lm ", im,jm,lm
      allocate (PRES_p(im,jm,lm))
      allocate (H_p(im,jm,lm))
      allocate (T_p(im,jm,lm))
      allocate (Q_p(im,jm,lm))
      allocate (U_p(im_u,jm,lm))
      allocate (V_p(im,jm_v,lm))

      CALL check( nf_inq_varid (p_id, "PRES", id_var) )
      CALL check( nf_get_var_real(p_id, id_var, PRES_p) )
      print *, "MEMBER pressure levels ", PRES_p(1,1,:)

      CALL check( nf_inq_varid (p_id, "GHT", id_var) )
      CALL check( nf_get_var_real(p_id, id_var, H_p) )

      CALL check( nf_inq_varid (p_id, "TT", id_var) )
      CALL check( nf_get_var_real(p_id, id_var, T_p) )

      CALL check( nf_inq_varid (p_id, "RH", id_var) )
      CALL check( nf_get_var_real(p_id, id_var, Q_p) )

      CALL check( nf_inq_varid (p_id, "UU", id_var) )
      CALL check( nf_get_var_real(p_id, id_var, U_p) )

      CALL check( nf_inq_varid (p_id, "VV", id_var) )
      CALL check( nf_get_var_real(p_id, id_var, V_p) )

      CALL check( nf_close(p_id) )

    !
    ! n
    !
      CALL check( nf_open(n_fname, NF_NOWRITE, n_id) )

      CALL check( nf_inq_dimid(p_id, 'west_east', im_id) )
      CALL check( nf_inq_dimlen(p_id, im_id, im) )
      CALL check( nf_inq_dimid(p_id, 'south_north', jm_id) )
      CALL check( nf_inq_dimlen(p_id, jm_id, jm) )
      CALL check( nf_inq_dimid(p_id, 'num_metgrid_levels', lm_id) )
      CALL check( nf_inq_dimlen(p_id, lm_id, lm) )

      im_u = im
      if (dyn_opt ==2) im_u = im + 1
      jm_v = jm
      if (dyn_opt ==2) jm_v = jm + 1

      print *, "c00 im,jm,lm ", im,jm,lm
      allocate (PRES_n(im,jm,lm))
      allocate (H_n(im,jm,lm))
      allocate (T_n(im,jm,lm))
      allocate (Q_n(im,jm,lm))
      allocate (U_n(im_u,jm,lm))
      allocate (V_n(im,jm_v,lm))

      CALL check( nf_inq_varid (n_id, "PRES", id_var) )
      CALL check( nf_get_var_real(n_id, id_var, PRES_n) )
      print *, "c00 pressure levels ", PRES_p(1,1,:)

      CALL check( nf_inq_varid (n_id, "GHT", id_var) )
      CALL check( nf_get_var_real(n_id, id_var, H_n) )

      CALL check( nf_inq_varid (n_id, "TT", id_var) )
      CALL check( nf_get_var_real(n_id, id_var, T_n) )

      CALL check( nf_inq_varid (n_id, "RH", id_var) )
      CALL check( nf_get_var_real(n_id, id_var, Q_n) )

      CALL check( nf_inq_varid (n_id, "UU", id_var) )
      CALL check( nf_get_var_real(n_id, id_var, U_n) )

      CALL check( nf_inq_varid (n_id, "VV", id_var) )
      CALL check( nf_get_var_real(n_id, id_var, V_n) )

      CALL check( nf_close(n_id) )

  CASE DEFAULT
    write(0,*)'unknown io_form ', io_form
    stop
  END SELECT


  allocate(dH(im,jm,lm))
  allocate(dT(im,jm,lm))
  allocate(dQ(im,jm,lm))
  allocate(dU(im_u,jm,lm))
  allocate(dV(im,jm_v,lm))

  dH = H_p - H_n
  dT = T_p - T_n
  dQ = Q_p - Q_n
  dU = U_p - U_n
  dV = V_p - V_n

  do l=1,lm
     write(*,*)l,'min/max H_p',minval(H_p(:,:,l)), maxval(H_p(:,:,l))
     write(*,*)l,'min/max H_n',minval(H_n(:,:,l)), maxval(H_n(:,:,l))
     write(*,*)l,'min/max T_p',minval(T_p(:,:,l)), maxval(T_p(:,:,l))
     write(*,*)l,'min/max T_n',minval(T_n(:,:,l)), maxval(T_n(:,:,l))
     write(*,*)l,'min/max Q_p',minval(Q_p(:,:,l)), maxval(Q_p(:,:,l))
     write(*,*)l,'min/max Q_n',minval(Q_n(:,:,l)), maxval(Q_n(:,:,l))
     write(*,*)l,'min/max U_p',minval(U_p(:,:,l)), maxval(U_p(:,:,l))
     write(*,*)l,'min/max U_n',minval(U_n(:,:,l)), maxval(U_n(:,:,l))
     write(*,*)l,'min/max V_p',minval(V_p(:,:,l)), maxval(V_p(:,:,l))
     write(*,*)l,'min/max V_n',minval(V_n(:,:,l)), maxval(V_n(:,:,l))

     write(*,*)l,'min/max dH',minval(dH(:,:,l)), maxval(dH(:,:,l))
     write(*,*)l,'min/max dT',minval(dT(:,:,l)), maxval(dT(:,:,l))
     write(*,*)l,'min/max dQ',minval(dQ(:,:,l)), maxval(dQ(:,:,l))
     write(*,*)l,'min/max dU',minval(dU(:,:,l)), maxval(dU(:,:,l))
     write(*,*)l,'min/max dV',minval(dV(:,:,l)), maxval(dV(:,:,l))
  end do

  diff_pres = .false.
  if (lm == lm_ctl) then
     do l=1,lm
        write(*,*) ' PRES_p(1,1,l) /= PRES(1,1,l) ',  PRES_p(1,1,l),  PRES(1,1,l)
        if ( PRES_p(1,1,l) /= PRES(1,1,l) ) then
          diff_pres = .true.
          exit
        end if
     end do
  end if

  write(*,*)' lm, lm_ctl, diff_pres ',  lm, lm_ctl, diff_pres 

  if (lm /= lm_ctl .or. diff_pres ) then
  CALL interp_press2press_lin(PRES_p, PRES,          &
                              dH,dH_ctl,lm,          &
                              .FALSE.,.TRUE.,.FALSE., & ! extrap, ignore_lowest, t_field
                              1,im,1,jm,1,lm_ctl,    &
                              1,im,1,jm,1,lm_ctl )

  CALL interp_press2press_lin(PRES_p, PRES,          &
                              dT,dT_ctl,lm,          &
                              .FALSE.,.TRUE.,.FALSE., & ! extrap, ignore_lowest, t_field
                              1,im,1,jm,1,lm_ctl,    &
                              1,im,1,jm,1,lm_ctl )

  CALL interp_press2press_lin(PRES_p, PRES,          &
                              dQ,dQ_ctl,lm,          &
                              .FALSE.,.TRUE.,.FALSE., & ! extrap, ignore_lowest, t_field
                              1,im,1,jm,1,lm_ctl,    &
                              1,im,1,jm,1,lm_ctl )

  CALL interp_press2press_lin(PRES_p, PRES,          &
                              dU,dU_ctl,lm,          &
                              .FALSE.,.TRUE.,.FALSE., & ! extrap, ignore_lowest, t_field
                              1,im,1,jm,1,lm_ctl,    &
                              1,im_u,1,jm,1,lm_ctl )

  CALL interp_press2press_lin(PRES_p, PRES,          &
                              dV,dV_ctl,lm,          &
                              .FALSE.,.TRUE.,.FALSE., & ! extrap, ignore_lowest, t_field
                              1,im,1,jm,1,lm_ctl,    &
                              1,im,1,jm_v,1,lm_ctl )
  else
     dH_ctl = dH
     dT_ctl = dT
     dQ_ctl = dQ
     dU_ctl = dU
     dV_ctl = dV
  end if

  do l=1,lm_ctl
     write(*,*)l,'min/max dH_ctl',minval(dH_ctl(:,:,l)), maxval(dH_ctl(:,:,l))
     write(*,*)l,'min/max dT_ctl',minval(dT_ctl(:,:,l)), maxval(dT_ctl(:,:,l))
     write(*,*)l,'min/max dQ_ctl',minval(dQ_ctl(:,:,l)), maxval(dQ_ctl(:,:,l))
     write(*,*)l,'min/max dU_ctl',minval(dU_ctl(:,:,l)), maxval(dU_ctl(:,:,l))
     write(*,*)l,'min/max dV_ctl',minval(dV_ctl(:,:,l)), maxval(dV_ctl(:,:,l))
  end do

  H = H + dH_ctl
  T = T + dT_ctl
  Q = Q + dQ_ctl
  U = U + dU_ctl
  V = V + dV_ctl

  do i=1,im
   do j=1,jm
    do l=1,lm
! do l=1,lm_ctl
     if(Q(i,j,l).lt.0.0) Q(i,j,l)=0.0000000001
     if(Q(i,j,l).gt.0.03) Q(i,j,l)=0.03
    enddo
   enddo
  enddo

  SELECT CASE ( io_form  )
  CASE ( IO_NETCDF )

      CALL check( nf_open(ctl_fname, NF_WRITE, ctl_id) )

      CALL check( nf_inq_varid (ctl_id, "GHT", id_var) )
      CALL check( nf_put_var_real(ctl_id, id_var, H) )
      CALL check( nf_inq_varid (ctl_id, "TT", id_var) )
      CALL check( nf_put_var_real(ctl_id, id_var, T) )
      CALL check( nf_inq_varid (ctl_id, "RH", id_var) )
      CALL check( nf_put_var_real(ctl_id, id_var, Q) )
      CALL check( nf_inq_varid (ctl_id, "UU", id_var) )
      CALL check( nf_put_var_real(ctl_id, id_var, U) )
      CALL check( nf_inq_varid (ctl_id, "VV", id_var) )
      CALL check( nf_put_var_real(ctl_id, id_var, V) )

      CALL check( nf_close(ctl_id) )

  CASE DEFAULT
    write(0,*)'unknown io_form ', io_form
    stop
  END SELECT

  print *, "End of lbc_perturb_wrf"

END PROGRAM lbc_perturb_wrf

  SUBROUTINE check(status)
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    INTEGER, INTENT (IN) :: status
    
    IF (status /= nf_noerr) THEN 
      PRINT *, trim(nf_strerror(status))
      STOP "Stopped"
    END IF
  END SUBROUTINE check  

  SUBROUTINE interp_press2press_lin(press_in,press_out, &
                                    data_in, data_out,generic          &
     &,                             extrapolate,ignore_lowest,TFIELD    &
     &,                             IDS,IDE,JDS,JDE,KDS,KDE             &
     &,                             IMS,IME,JMS,JME,KMS,KME)

    ! Interpolates data from one set of pressure surfaces to
    ! another set of pressures

    INTEGER                            :: IDS,IDE,JDS,JDE,KDS,KDE
    INTEGER                            :: IMS,IME,JMS,JME,KMS,KME
    INTEGER                            :: generic

    REAL, INTENT(IN)                   :: press_in(IDS:IDE,JDS:JDE,generic)
    REAL, INTENT(IN)                   :: press_out(IDS:IDE,JDS:JDE,KDS:KDE)
    REAL, INTENT(IN)                   :: data_in(IMS:IME,JMS:JME,generic)
    REAL, INTENT(OUT)                  :: data_out(IMS:IME,JMS:JME,KMS:KME)
    LOGICAL, INTENT(IN)                :: extrapolate, ignore_lowest, TFIELD
    LOGICAL                            :: col_smooth

    INTEGER                            :: i,j
    INTEGER                            :: k,kk
    REAL                               :: desired_press
    REAL                               :: dvaldlnp,dlnp,tadiabat,tiso

    REAL, PARAMETER                    :: ADIAFAC=9.81/1004.
    REAL, PARAMETER                    :: TSTEXTRAPFAC=.0065


!d    data_out(:,:,:) = -99999.9

    IF (ignore_lowest) then
       LMIN=2
    ELSE
       LMIN=1
    ENDIF

    print *, ' JDS, JDE ', JDS, JDE
    print *, ' IDS, IDE ', IDS, IDE
    DO j = JDS, JDE
      DO i = IDS, IDE

       col_smooth=.false.

        output_loop: DO k = KDS,KDE

          desired_press = press_out(i,j,k)

        if (K .gt. KDS) then
        if (TFIELD .and. col_smooth .and. desired_press .lt. press_in(i,j,LMIN) &
                                    .and. press_out(i,j,k-1) .gt. press_in(i,j,LMIN)) then
          MAX_SMOOTH=K
        endif
        endif

! keep track of where the extrapolation begins

          IF (desired_press .GT. press_in(i,j,LMIN)) THEN
           IF (TFIELD .and. K .eq. 1  .and. (desired_press - press_in(i,j,LMIN)) .gt. 3000.) then
            col_smooth=.TRUE.   ! due to large extrapolation distance
           ENDIF


            IF ((desired_press - press_in(i,j,LMIN)).LT. 50.) THEN ! 0.5 mb
               data_out(i,j,k) = data_in(i,j,LMIN)
            ELSE
              IF (extrapolate) THEN
                ! Extrapolate downward because desired P level is below
                ! the lowest level in our input data.  Extrapolate using simple
                ! 1st derivative of value with respect to ln P for the bottom 2
                ! input layers.

                ! Add a check to make sure we are not using the gradient of
                ! a very thin layer

                if (TFIELD) then
                  tiso=0.5*(data_in(i,j,1)+data_in(i,j,2))
                endif


                IF ( (press_in(i,j,LMIN)-press_in(i,j,LMIN+1)) .GT. 500.) THEN ! likely isobaric data
                  dlnp     = log(press_in(i,j,LMIN))-log(press_in(i,j,LMIN+1))
                  dvaldlnp = (data_in(i,j,LMIN) - data_in(i,j,LMIN+1)) / dlnp
                ELSE                                                           ! assume terrain following
                  dlnp     = log(press_in(i,j,LMIN))-log(press_in(i,j,LMIN+5))
                  dvaldlnp = (data_in(i,j,LMIN) - data_in(i,j,LMIN+5)) / dlnp
                ENDIF
                data_out(i,j,k) = data_in(i,j,LMIN) + dvaldlnp * &
                               ( log(desired_press)-log(press_in(i,j,LMIN)) )

        if (TFIELD .and. data_out(i,j,k) .lt. tiso-0.2) then

! restrict slope to -1K/10 hPa
          dvaldlnp=max(dvaldlnp, -1.0/ &
                                log( press_in(i,j,LMIN) / &
                                   ( press_in(i,j,LMIN)-1000.)  ))

          data_out(I,J,K)= data_in(i,j,LMIN) + dvaldlnp * &
                               ( log(desired_press)-log(press_in(i,j,LMIN)) )

        elseif (TFIELD .and. data_out(i,j,k) .gt. tiso+0.2) then

! restrict slope to +0.8K/10 hPa
          dvaldlnp=min(dvaldlnp, 0.8/ &
                                log( press_in(i,j,LMIN) / &
                                   ( press_in(i,j,LMIN)-1000.)  ))

          data_out(I,J,K)= data_in(i,j,LMIN) + dvaldlnp * &
                               ( log(desired_press)-log(press_in(i,j,LMIN)) )

         endif

              ELSE
                data_out(i,j,k) = data_in(i,j,LMIN)
              ENDIF
            ENDIF
          ELSE IF (desired_press .LT. press_in(i,j,generic)) THEN
            IF ( (press_in(i,j,generic) - desired_press) .LT. 10.) THEN
               data_out(i,j,k) = data_in(i,j,generic)
            ELSE
              IF (extrapolate) THEN
                ! Extrapolate upward
                IF ((press_in(i,j,generic-1)-press_in(i,j,generic)).GT.50.) THEN
                  dlnp    =log(press_in(i,j,generic))-log(press_in(i,j,generic-1))
                  dvaldlnp=(data_in(i,j,generic)-data_in(i,j,generic-1))/dlnp
                ELSE
                  dlnp    =log(press_in(i,j,generic))-log(press_in(i,j,generic-2))
                  dvaldlnp=(data_in(i,j,generic)-data_in(i,j,generic-2))/dlnp
                ENDIF
                data_out(i,j,k) =  data_in(i,j,generic) + &
                  dvaldlnp * (log(desired_press)-log(press_in(i,j,generic)))
              ELSE
                data_out(i,j,k) = data_in(i,j,generic)
              ENDIF
            ENDIF
          ELSE
            ! We can trap between two levels and linearly interpolate

            input_loop:  DO kk = LMIN, generic-1
              IF (desired_press .EQ. press_in(i,j,kk) )THEN
                data_out(i,j,k) = data_in(i,j,kk)
                EXIT input_loop
              ELSE IF ( (desired_press .LT. press_in(i,j,kk)) .AND. &
                        (desired_press .GT. press_in(i,j,kk+1)) ) THEN

!       do trapped in lnp

         dlnp = log(press_in(i,j,kk)) - log(press_in(i,j,kk+1))
         dvaldlnp = (data_in(i,j,kk)-data_in(i,j,kk+1))/dlnp
         data_out(i,j,k) = data_in(i,j,kk+1)+ &
                           dvaldlnp*(log(desired_press)-log(press_in(i,j,kk+1)))

                EXIT input_loop
              ELSE IF (desired_press .EQ. press_in(i,j,kk+1) )THEN
                data_out(i,j,k) = data_in(i,j,kk+1)
                EXIT input_loop
              ENDIF

            ENDDO input_loop
          ENDIF
        ENDDO output_loop

        if (col_smooth) then
       do K=max(KDS,MAX_SMOOTH-4),MAX_SMOOTH+4
       data_out(I,J,K)=0.5*(data_out(I,J,K)+data_out(I,J,K+1))
       enddo
        endif

      ENDDO
    ENDDO
  END SUBROUTINE interp_press2press_lin

SUBROUTINE wrf_abort
  STOP
END SUBROUTINE wrf_abort
