!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                KMC Simulation                !!
!!   1D morphology with different materials     !!
!!  Probability rates given by FRET mechanism   !!
!! Random exciton distribution and annihilation !!
!!					        !!
!! 15/01/2021                                   !!
!! Last update: 11/02/2022                      !!
!!         Laura Simonassi                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module precision

integer,parameter                 ::  digits = 15
integer,parameter                 ::  rp = selected_real_kind(digits)
integer,parameter                 ::  ip = selected_int_kind(10)

end module precision

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module var
use precision
integer(ip)                          :: matC,matR,matL
integer(ip), parameter               :: NR=1000, NE=1, NM=2                                            !round, excitons, materials
integer(ip), dimension(NE)           :: iev                                                            !Index of live excitons. 1-live 0-flourished
real(rp), dimension(NE)              :: pos_x,pos_y,pos_z,pos_hx,pos_hy,pos_hz,pos_xi,pos_yi,pos_zi,dt
real(rp), parameter                  :: xi=-300.0_rp          
real(rp), parameter                  :: xf= 300.0_rp
real(rp), parameter                  :: yi= 0.0_rp           
real(rp), parameter                  :: yf= 0.0_rp
real(rp), parameter                  :: zi= 0.0_rp           
real(rp), parameter                  :: zf= 0.0_rp
real(rp)                             :: dt_min 
real(rp), parameter                  :: alpha   = 1.15_rp*0.53_rp                                       !cte
real(rp), parameter                  :: r       = 10.0_rp                                               !angstrom
real(rp), dimension(NM)              :: kf, prob, moment, temi, k_emi
real(rp), dimension(2,2)             :: rf 
end module var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program kmc
use var
implicit none
integer(ip)                       ::  j, k
integer                           :: seed(33)
real(rp)                          :: probm, t
real(rp)                          :: num, poshx,poshy,poshz           
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

!PTh centro (pos 2)
!moment(1)= 3.79_rp                   
!moment(2)= 2.02_rp                 
!temi(1) = 846.056735236_rp*1E-12    
!temi(2) = 2260.95185812_rp*1E-12 
!rf(1,1)=42.0_rp 
!rf(2,2)=28.6_rp   
!rf(1,2)=25.7_rp 
!rf(2,1)=43.5_rp  
 
!PPV centro (pos 2)
!moment(1)= 2.02_rp               
!moment(2)= 3.79_rp                  
!temi(1) = 2260.95185812_rp*1E-12    
!temi(2) = 846.056735236_rp*1E-12
!rf(1,1)=28.6_rp 
!rf(2,2)=42.0_rp   
!rf(1,2)=43.5_rp 
!rf(2,1)=25.7_rp

open(10, file = "Simulation.txt", status='unknown')
do j=1,NR
    call exciton_distribution
    t = 0.0_rp
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
                call morphology(pos_x(k)-r,pos_y(k),pos_z(k),matL)
                call morphology(pos_x(k)  ,pos_y(k),pos_z(k),matC)
                call morphology(pos_x(k)+r,pos_y(k),pos_z(k),matR)

                kf(1)= (1.0_rp/temi(matC))*((rf(matC,matL)/(alpha*moment(matC) + r))**6)
                kf(2)= (1.0_rp/temi(matC))*((rf(matC,matR)/(alpha*moment(matC) + r))**6)
            
                prob(1)=kf(1)/sum(kf)
                prob(2)=kf(2)/sum(kf)
            
                call random_number(num)
            
                if (num .lt. prob(1)) then
                    poshx=pos_hx(k) - r
                    dt(k)=1.0_rp/kf(1)
                    probm=k_emi(matC)/(k_emi(matC)+kf(1))
                end if

                if (num .ge. prob(1) .and. num .lt. prob(1)+prob(2) ) then
                    poshx=pos_hx(k)+r
                    dt(k)=1.0_rp/kf(2)
                    probm=k_emi(matC)/(k_emi(matC)+kf(2))
                end if

                call random_number(num)
                if(num .lt. probm) then
                    write(10,101) t, sqrt((pos_x(k)-pos_xi(k))**2+(pos_y(k)-pos_yi(k))**2+(pos_z(k)-pos_zi(k))**2), pos_x(k), pos_y(k),pos_z(k),"Fluor", matC
                    iev(k)=0
                endif

                if(num .ge. probm) then
                    pos_hx(k)=poshx
                end if

               call random_number(num)
               if(num .lt. dt_min/dt(k)) then
                   pos_x(k)=pos_hx(k)
               end if
           end if
       end do    
       
       t=t+dt_min           
       end do loopt
   end do
close(10)
101 format(E19.10,1x,E19.10,1x,E19.10,1x,E19.10,1x,E19.10,1x,A5,1x,I2.2)
end program kmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine exciton_distribution
use var
implicit none
integer               :: i
real(rp),dimension(2) :: m

do i=1,NE
   call random_number(m)
   pos_x(i)=m(1)*(xf-xi)+xi
   pos_xi(i)=pos_x(i)
   pos_y(i)=0.0_rp
   pos_yi(i)=pos_y(i)
   pos_z(i)=0.0_rp
   pos_zi(i)=pos_z(i)
enddo

return
end subroutine exciton_distribution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

