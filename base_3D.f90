!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                KMC Simulation                !!
!!   3D morphology with different materials     !!
!!  Probability rates given by FRET mechanism   !!
!! Random exciton distribution and annihilation !!
!!                                              !!
!! 09/01/2022                                   !!
!! Last update: 04/02/2021                      !!
!!         Laura Simonassi                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module precision
integer,parameter            ::  digits = 15
integer,parameter            ::  rp = selected_real_kind(digits)
integer,parameter            ::  ip = selected_int_kind(10)
end module precision
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module var
use precision
integer(ip)                          :: matloc
integer(ip), parameter               :: NR=1000, NE=50, NM=2                       !round, excitons, materials
integer(ip), dimension(NE)           :: iev                                        !Index of live excitons. 1-live 0-flourished
integer(rp), dimension(26)           :: mat
real(rp)                             :: dt_min,kf_tot 
real(rp), dimension(NE)              :: pos_x,pos_y,pos_z,pos_hx,pos_hy,pos_hz,dt
real(rp), dimension(26)              :: dx,dy,dz,prob,kf
real(rp), dimension(NM)              :: moment                                     !atomic unit
real(rp), dimension(NM)              :: temi                                       !seconds
real(rp), dimension(NM)              :: k_emi
real(rp), dimension(2,2)             :: rf
!!!!!!!!!!!Parameters!!!!!!!!!!!!!!!!!
real(rp), parameter                  :: alpha = 1.15_rp*0.53_rp                    !cte
real(rp), parameter                  :: r     = 10.0_rp                            !angstrom
real(rp), parameter                  :: xi=-300.0_rp, xf= 300.0_rp
real(rp), parameter                  :: yi=-300.0_rp, yf= 300.0_rp
real(rp), parameter                  :: zi=-300.0_rp, zf= 300.0_rp 
end module var
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program kmc
use var
implicit none
integer(ip)                       :: j, k, i
integer                           :: seed(33)
real(rp)                          :: probm, t!, x 
real(rp)                          :: num, poshx,poshy,poshz           
!call system_clock(seed(1))
seed(1)=1
call random_seed(put=seed)
call parameters

!PTh centro (pos 2)
moment(1)= 0.0_rp!3.79_rp                   
moment(2)= 2.02_rp                 
temi(1) = 0.0_rp!5000_rp*1E-12    
temi(2) = 5000_rp*1E-12 
rf(1,1)=0.0_rp!42.0_rp 
rf(2,2)=28.6_rp   
rf(1,2)=0.0_rp!25.7_rp 
rf(2,1)=0.0_rp!43.5_rp  
k_emi(1)=1.0_rp/temi(1)
k_emi(2)=1.0_rp/temi(2)

!!PPV centro (pos 2)
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
!
open(10, file = "Simulation.txt", status='unknown')
do j=1,NR
    call exciton_distribution
    t = 0.0_rp
   ! x = 0.0_rp
    iev=1
    loopt: do
        call annihilation(t, j)
        pos_hx=pos_x
        pos_hy=pos_y
        pos_hz=pos_z
        dt_min=minval(dt)

        do k=1,NE
            if (all(iev .eq. 0)) exit loopt
            if (iev(k) .ne. 0) then 
               call morphology(pos_x(k),pos_y(k),pos_z(k),matloc)
               do i=1,26
                  call morphology(pos_x(k)+dx(i),pos_y(k)+dy(i),pos_z(k)+dz(i),mat(i))
                  kf(i)=(1.0_rp/temi(matloc))*((rf(matloc,mat(i))/(alpha*moment(matloc)+sqrt(dx(i)**2+dy(i)**2+dz(i)**2)))**6)
               end do
!              call gaussian(x)

               kf_tot =sum(kf)
               prob=kf/kf_tot
               
               do i=2,26
                 prob(i)=prob(i)+prob(i-1) 
               end do

               call random_number(num)
            
               if (num .lt. prob(1)) then
                   poshx=pos_hx(k) + dx(1)
                   poshy=pos_hy(k) + dy(1)
                   poshz=pos_hz(k) + dz(1)
               end if
               
               do i=1,25
                   if (num .ge. prob(i) .and. num .lt. prob(i+1)) then
                       poshx=pos_hx(k)+dx(i+1)        
                       poshy=pos_hy(k)+dy(i+1)        
                       poshz=pos_hz(k)+dz(i+1)
                       dt(k)=1.0_rp/kf(i+1)
                       probm=k_emi(matloc)/(k_emi(matloc)+kf(i+1))        
                   end if 
               end do

               call random_number(num)
               if(num .lt. probm) then
                   write(10,101) t, pos_x(k), pos_y(k), pos_z(k),"Fluor", matloc
                   iev(k)=0
               endif

               if(num .ge. probm) then
                   pos_hx(k)=poshx
                   pos_hy(k)=poshy
                   pos_hz(k)=poshz
               end if

               call random_number(num)
               if(num .lt. dt_min/dt(k)) then
                   pos_x(k)=pos_hx(k)
                   pos_y(k)=pos_hy(k)
                   pos_z(k)=pos_hz(k)
               end if
            end if
        end do    
       
        t=t+dt_min           
    end do loopt
end do
close(10)
101 format(E19.10,1x,E19.10,1x,E19.10,1x,E19.10,1x,A5,1x,I2.2)
end program kmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine exciton_distribution
use var
implicit none
integer               :: i
real(rp),dimension(3) :: num

do i=1,NE
   !pos_x(i)=100.0_rp
   !pos_y(i)=0.0_rp
   !pos_z(i)=0.0_rp
   call random_number(num)
   pos_x(i)=num(1)*(xf-xi)+xi
   pos_y(i)=num(2)*(yf-yi)+yi
   pos_z(i)=num(3)*(zf-zi)+zi
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parameters
use var 
implicit none
integer           :: i,j,k,l

l=1
do i=-1,1,1
    do j=-1,1,1
        do k=-1,1,1
           if (i .ne. 0 .or. j .ne. 0 .or. k .ne. 0) then
               dx(l)=dfloat(i)*r
               dy(l)=dfloat(j)*r
               dz(l)=dfloat(k)*r
               l=l+1
           end if
        end do
    end do
end do

return
end subroutine parameters
