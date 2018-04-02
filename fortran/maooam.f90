
!  maooam.f90
!
!> Fortran 90 implementation of the modular arbitrary-order ocean-atmosphere 
!> model MAOOAM.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM maooam 
  USE params, only: ndim, dt, tw, t_trans, t_run, writeout
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE tl_ad_tensor, only: init_tltensor
  USE tl_ad_integrator, only: init_tl_ad_integrator,tl_step
  USE stat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew       !< State variable in the model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dy     !<  Initial perturbation
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tlm    !< Updated perturbation
  REAL(KIND=8) :: t=0.D0,t_tlm                             !< Time variable
  REAL(KIND=8) :: t_up

  INTEGER :: i,j

  PRINT*, 'Model MAOOAM v1.2'
  PRINT*, 'Loading information...'

  CALL init_aotensor    ! Compute the tensor

  CALL init_tltensor    ! Compute the tensor

  CALL load_IC          ! Load the initial condition

  CALL init_integrator  ! Initialize the integrator

  CALL init_tl_ad_integrator  ! Initialize the TL & AD integrator

  t_up=dt/t_trans*100.D0

  IF (writeout) OPEN(10,file='evol_field_tlm_takuma_ic_2.dat')

  ALLOCATE(X(0:ndim),Xnew(0:ndim))

  ALLOCATE(dy(0:ndim,0:ndim),tlm(0:ndim,0:ndim))

  X=IC

  PRINT*, 'Starting the transient time evolution...'

  DO WHILE (t<t_trans)
     CALL step(X,t,dt,Xnew)
     X=Xnew
     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
  END DO

  PRINT*, 'Starting the time evolution...'

  CALL init_stat
  
  t=0.D0

!  IF (writeout) WRITE(10,*) t,X(1:ndim)

  dy = 0.d00
  DO i=0,ndim
     dy(i,i) = 1.0d0
  END DO


  DO WHILE (t<t_run)
     !print*, t 
     CALL step(X,t,dt,Xnew)
     DO j=0,ndim
        t=t-dt  
        CALL tl_step(dy(:,j),X,t,dt,tlm(:,j))
     END DO
     dy = tlm
     X=Xnew
     
     IF (mod(t,tw)<dt) THEN
        IF (writeout) WRITE(10,*) t,X(1:ndim),tlm(1:ndim,1:ndim)
        CALL acc(X)
         dy = 0.d00
         DO i=0,ndim
            dy(i,i) = 1.0d0
         END DO
         print*,"write",t
     END IF
!     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
 print*, t
 END DO

  PRINT*, 'Evolution finished.'
  print*, Xnew
  IF (writeout) CLOSE(10)

DEALLOCATE (X,Xnew,tlm)
END PROGRAM maooam 
