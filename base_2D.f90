!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                KMC Simulation                !!
!!   1D morphology with different materials     !!
!!  Probability rates given by FRET mechanism   !!
!! Random exciton distribution and annihilation !!
!!					        !!
!! 09/01/2022                                   !!
!! Last update: 11/02/2022                      !!
!!         Laura Simonassi                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module precision

integer,parameter                 ::  digits = 15
integer,parameter                 ::  rp = selected_real_kind(digits)
integer,parameter                 ::  ip = selected_int_kind(10)

end module precision

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module var
use precision
integer(ip)                          :: matC,matR,matL,matU,matD,matUR,matDR,matUL,matDL
integer(ip), parameter               :: NR=1000, NE=1, NM=2             !round, excitons, materials
integer(ip), dimension(NE)           :: iev                              !Index of live excitons. 1-live 0-flourished
real(rp), dimension(NE)              :: pos_x,pos_y,pos_z,pos_hx,pos_hy,pos_hz,pos_xi,pos_yi,pos_zi,dt
real(rp), parameter                  :: xi=-300.0_rp          
real(rp), parameter                  :: xf= 300.0_rp
real(rp), parameter                  :: yi=-300.0_rp           
real(rp), parameter                  :: yf= 300.0_rp
real(rp), parameter                  :: zi= 0.0_rp           
real(rp), parameter                  :: zf= 0.0_rp
real(rp)                             :: dt_min 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(rp)                             :: pi    = dacos(-1.0d0)
real(rp), parameter                  :: kb    = 1.38064852_rp*10E-24        !SI
real(rp), parameter                  :: emass = 9.10938356_rp*10E-32        !kg
real(rp), parameter                  :: temp  = 300.0_rp                    !K
real(rp), parameter                  :: alpha = 1.15_rp*0.53_rp             !cte
real(rp), parameter                  :: r     = 10.0_rp                     !angstrom
real(rp), dimension(NM)              :: moment                              !atomic unit
real(rp), dimension(NM)              :: temi                                !seconds
real(rp), dimension(NM)              :: k_emi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(rp)                             :: kf_R,kf_L,kf_U,kf_D,kf_UR,kf_UL,kf_DR,kf_DL, kf_tot
real(rp)                             :: prob_R,prob_L,prob_U,prob_D,prob_UR,prob_UL,prob_DR,prob_DL  
real(rp), dimension(2,2)             :: rf 
end module var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program kmc
use var
implicit none
integer(ip)                       :: j, k
integer                           :: seed(33)
real(rp)                          :: probm, t!, x 
real(rp)                          :: num, poshx,poshy           
!call system_clock(seed(1))
seed(1)=1
call random_seed(put=seed)

!PTh centro (pos 2)
moment(1)= 0!3.79_rp                   
moment(2)= 0!2.02_rp                 
temi(1) = 0!5000_rp*1E-12    
temi(2) = 2260.95185812_rp*1E-12 
rf(1,1)=0!42.0_rp 
rf(2,2)=28.6_rp   
rf(1,2)=0!25.7_rp 
rf(2,1)=0!43.5_rp  
k_emi(1)=1.0_rp/temi(1)
k_emi(2)=1.0_rp/temi(2)

!PPV centro (pos 2)
!moment(1)= 2.02_rp               
!moment(2)= 3.79_rp                  
!temi(1) = 5000_rp*1E-12    
!temi(2) = 5000_rp*1E-12 
!rf(1,1)=28.6_rp 
!rf(2,2)=42.0_rp   
!rf(1,2)=43.5_rp 
!rf(2,1)=25.7_rp
!k_emi(1)=1.0_rp/temi(1)
!k_emi(2)=1.0_rp/temi(2)

open(10, file = "Simulation.txt", status='unknown')
do j=1,NR
    call exciton_distribution
    t = 0.0_rp
   ! x = 0.0_rp
    iev=1
    loopt: do
 !       call annihilation(t, j)
        pos_hx=pos_x
        pos_hy=pos_y
        pos_hz=pos_z
        dt_min=minval(dt)

        do k=1,NE
            if (all(iev .eq. 0)) exit loopt
            if (iev(k) .ne. 0) then
                call morphology(pos_x(k)  ,pos_y(k)  ,pos_z(k),matC)
                call morphology(pos_x(k)-r,pos_y(k)  ,pos_z(k),matL)
                call morphology(pos_x(k)+r,pos_y(k)  ,pos_z(k),matR)
                call morphology(pos_x(k)  ,pos_y(k)+r,pos_z(k),matU)
                call morphology(pos_x(k)  ,pos_y(k)-r,pos_z(k),matD)
                call morphology(pos_x(k)+r,pos_y(k)+r,pos_z(k),matUR)
                call morphology(pos_x(k)+r,pos_y(k)-r,pos_z(k),matDR)           
                call morphology(pos_x(k)-r,pos_y(k)+r,pos_z(k),matUL)
                call morphology(pos_x(k)-r,pos_y(k)-r,pos_z(k),matDL)

!                call gaussian(x)
            
                kf_R =(1.0_rp/temi(matC))*((rf(matC,matR)/(alpha*moment(matC)+r))**6)
                kf_L =(1.0_rp/temi(matC))*((rf(matC,matL)/(alpha*moment(matC)+r))**6) 
                kf_U =(1.0_rp/temi(matC))*((rf(matC,matU)/(alpha*moment(matC)+r))**6)             
                kf_D =(1.0_rp/temi(matC))*((rf(matC,matD)/(alpha*moment(matC)+r))**6)  
                kf_UR=(1.0_rp/temi(matC))*((rf(matC,matUR)/(alpha*moment(matC)+r*sqrt(2.0_rp)))**6)
                kf_UL=(1.0_rp/temi(matC))*((rf(matC,matUL)/(alpha*moment(matC)+r*sqrt(2.0_rp)))**6)
                kf_DR=(1.0_rp/temi(matC))*((rf(matC,matDR)/(alpha*moment(matC)+r*sqrt(2.0_rp)))**6)
                kf_DL=(1.0_rp/temi(matC))*((rf(matC,matDL)/(alpha*moment(matC)+r*sqrt(2.0_rp)))**6)

                kf_tot = kf_R + kf_L + kf_U + kf_D + kf_UR + kf_UL + kf_DR + kf_DL

                prob_R = kf_R/kf_tot
                prob_L = kf_L/kf_tot
                prob_U = kf_U/kf_tot
                prob_D = kf_D/kf_tot
                prob_UR=kf_UR/kf_tot
                prob_UL=kf_UL/kf_tot
                prob_DR=kf_DR/kf_tot
                prob_DL=kf_DL/kf_tot
           
                                
                call random_number(num)
            
                if (num .lt. prob_R) then
                    poshx=pos_hx(k) + r
                    poshy=pos_hy(k)
                    dt(k)=1.0_rp/kf_R
                    probm=k_emi(matC)/(k_emi(matC)+kf_R)
                end if

                if (num .ge. prob_R .and. num .lt. prob_R+prob_L ) then
                    poshx=pos_hx(k)-r
                    poshy=pos_hy(k)
                    dt(k)=1.0_rp/kf_L
                    probm=k_emi(matC)/(k_emi(matC)+kf_L)
                end if

                if (num .ge. prob_R+prob_L .and. num .lt. prob_R+prob_L+prob_U) then
                    poshx=pos_hx(k)
                    poshy=pos_hy(k)+r
                    dt(k)=1.0_rp/kf_U
                    probm=k_emi(matC)/(k_emi(matC)+kf_U)
                end if

               if (num .ge. prob_R+prob_L+prob_U .and. num .lt. prob_R+prob_L+prob_U+prob_D) then
                    poshx=pos_hx(k)
                    poshy=pos_hy(k)-r
                    dt(k)=1.0_rp/kf_D
                    probm=k_emi(matC)/(k_emi(matC)+kf_D)
                end if
               
               if (num .ge. prob_R+prob_L+prob_U+prob_D .and. num .lt. prob_R+prob_L+prob_U+prob_D+prob_UR) then
                    poshx=pos_hx(k)+r
                    poshy=pos_hy(k)+r
                    dt(k)=1.0_rp/kf_UR
                    probm=k_emi(matC)/(k_emi(matC)+kf_UR)
                end if

               if (num .ge. prob_R+prob_L+prob_U+prob_D+prob_UR .and. num .lt. prob_R+prob_L+prob_U+prob_D+prob_UR+prob_UL ) then
                    poshx=pos_hx(k)-r
                    poshy=pos_hy(k)+r
                    dt(k)=1.0_rp/kf_UL
                    probm=k_emi(matC)/(k_emi(matC)+kf_UL)
                end if

               if (num .ge. prob_R+prob_L+prob_U+prob_D+prob_UR+prob_UL .and. num .lt. prob_R+prob_L+prob_U+prob_D+prob_UR+prob_UL+prob_DR ) then
                    poshx=pos_hx(k)+r
                    poshy=pos_hy(k)-r
                    dt(k)=1.0_rp/kf_DR
                    probm=k_emi(matC)/(k_emi(matC)+kf_DR)
                end if

               if (num  .ge. prob_R+prob_L+prob_U+prob_D+prob_UR+prob_UL+prob_DR .and. num .lt. prob_R+prob_L+prob_U+prob_D+prob_UR+prob_UL+prob_DR+prob_DL) then
                    poshx=pos_hx(k)-r
                    poshy=pos_hy(k)-r
                    dt(k)=1.0_rp/kf_DL
                    probm=k_emi(matC)/(k_emi(matC)+kf_DL)
                end if

                call random_number(num)
                if(num .lt. probm) then
                    write(10,101) t, sqrt((pos_x(k)-pos_xi(k))**2+(pos_y(k)-pos_yi(k))**2+(pos_z(k)-pos_zi(k))**2), pos_x(k), pos_y(k),pos_z(k),"Fluor", matC
                    iev(k)=0
                endif

                if(num .ge. probm) then
                    pos_hx(k)=poshx
                    pos_hy(k)=poshy
                end if

               call random_number(num)
               if(num .lt. dt_min/dt(k)) then
                   pos_x(k)=pos_hx(k)
                   pos_y(k)=pos_hy(k)
               end if
           end if
       end do    
       
       t=t+dt_min           
       end do loopt
   end do
close(10)
101 format(E19.10,1x,E19.10,1x,E19.10,1x,E19.10,1x,E19.10,1x,A5,1x,I2.2)
end program kmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine exciton_distribution
use var
implicit none
integer               :: i
real(rp),dimension(2) :: m

do i=1,NE
   call random_number(m)
   pos_x(i)=m(1)*(xf-xi)+xi
   pos_xi(i)=pos_x(i)
   pos_y(i)=m(2)*(xf-xi)+xi
   pos_yi(i)=pos_y(i)
   pos_z(i)=0.0_rp
   pos_zi(i)=pos_z(i)
enddo
return
end subroutine exciton_distribution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine morphology(x,y,z,material)
use precision
use var
implicit none
integer(ip)        :: material
real(rp)           :: x,y,z

if(x .ge. xi .and. x .le. xf .and. y .ge. yi .and. y .le. yf .and. z .ge. zi .and. z .le. zf) then
    material=2
else
    material=1
end if

return
end subroutine morphology

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine annihilation(t,iround)
use var
implicit none
integer          :: i,j
integer(ip)      :: iround,material
real(rp)         :: t,d
do i=1,NE-1
   if(iev(i) .ne. 0) then
      do j=i+1,NE
         if(iev(j) .ne. 0) then                   
             call morphology(pos_x(j),pos_y(j),pos_z(j),material)
             d=sqrt((pos_x(j)-pos_x(i))**2+(pos_y(j)-pos_y(i))**2+(pos_z(j)-pos_z(i))**2)
             if (d .le. r) then 
                 iev(j)=0
                 write(10,101)t,pos_x(j),pos_y(j),pos_z(j),"Annih",material
             endif
         endif
      enddo
   endif
enddo
101 format(E19.10,1x,E19.10,1x,E19.10,1x,E19.10,1x,A5,1x,I2.2)
return
end subroutine annihilation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!subroutine gaussian(y)
!use var
!implicit none 
!real(rp)       :: rnd(2), y, sigma, kforster
!
!!kforster= (1.0_rp/temi(2))*((rf(mat2,mat1)/(alpha*moment(2) + r))**6)
!!sigma = sqrt((kb*temp)/(emass*(kforster)**2))
!sigma = 3.0E0 
!
!call random_number (rnd)
!y= sigma*sqrt(-2.0*log(rnd(1)))*cos(2.0*pi*rnd(2)) + r 
!
!return
!end subroutine gaussian
