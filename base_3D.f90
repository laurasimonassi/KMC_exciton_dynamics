!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     Simulação em Monte Carlo para o Random Walk 2D   !!
!! Implementação de Taxas de FRET e emissão de excitons !!
!! Implementação de diferentes materiais                !!
!! Implementação de desordem posicional                 !!
!! Implementação de aniquilação                         !!
!!                                                      !!
!! 09/01/2022                                           !!
!! Última atualização: 25/01/2021                       !!
!!         Laura Simonassi                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module precision

integer,parameter                 ::  digits = 15
integer,parameter                 ::  rp = selected_real_kind(digits)
integer,parameter                 ::  ip = selected_int_kind(10)

end module precision
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module var
use precision
integer(ip)                          ::matPOO,matPPO,matPMO,matPOP,matPOM,matMOO,matMPO,matMMO,matMOP,matMOM,matOPO,matOPP,matOPM,matOMO,matOMP,matOMM,matOOO,matOOP,matOOM,matPPP,matPPM,matPMP,matPMM,matMPP,matMPM,matMMP,matMMM
integer(ip), parameter               :: NR=1000, NE=50, NM=2             !round, excitons, materials
integer(ip), dimension(NE)           :: iev                              !Index of live excitons. 1-live 0-flourished
real(rp), dimension(NE)              :: pos_x,pos_y,pos_z,pos_hx,pos_hy,pos_hz,dt
real(rp), parameter                  :: xi=-300.0_rp          
real(rp), parameter                  :: xf= 300.0_rp  
real(rp), parameter                  :: yi=-300.0_rp              
real(rp), parameter                  :: yf= 300.0_rp
real(rp), parameter                  :: zi=-300.0_rp 
real(rp), parameter                  :: zf= 300.0_rp
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
real(rp)                             ::kf_tot,kf_POO,kf_PPO,kf_PMO,kf_POP,kf_POM,kf_MOO,kf_MPO,kf_MMO,kf_MOP,kf_MOM,kf_OPO,kf_OPP,kf_OPM,kf_OMO,kf_OMP,kf_OMM,kf_OOP,kf_OOM,kf_PPP,kf_PPM,kf_PMP,kf_PMM,kf_MPP,kf_MPM,kf_MMP,kf_MMM
real(rp)                             ::prob_POO,prob_PPO,prob_PMO,prob_POP,prob_POM,prob_MOO,prob_MPO,prob_MMO,prob_MOP,prob_MOM,prob_OPO,prob_OPP,prob_OPM,prob_OMO,prob_OMP,prob_OMM,prob_OOP,prob_OOM,prob_PPP,prob_PPM,prob_PMP,prob_PMM,prob_MPP,prob_MPM,prob_MMP,prob_MMM
real(rp), dimension(2,2)             :: rf 
end module var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program rw
use var
implicit none
integer(ip)                       :: signal, j, k
integer                           :: seed(33)
real(rp)                          :: probm, t!, x 
real(rp)                          :: num, poshx,poshy,poshz           
!call system_clock(seed(1))
seed(1)=1
call random_seed(put=seed)

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
open(10, file = "Simulation_Ne_XXX.txt", status='unknown')
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
		 call morfology(pos_x(k)+r ,pos_y(k)   ,pos_z(k)   ,matPOO)
		 call morfology(pos_x(k)+r ,pos_y(k)+r ,pos_z(k)   ,matPPO)
		 call morfology(pos_x(k)+r ,pos_y(k)-r ,pos_z(k)   ,matPMO)
		 call morfology(pos_x(k)+r ,pos_y(k)   ,pos_z(k)+r ,matPOP)
		 call morfology(pos_x(k)+r ,pos_y(k)   ,pos_z(k)-r ,matPOM)
		 call morfology(pos_x(k)-r ,pos_y(k)   ,pos_z(k)   ,matMOO)
		 call morfology(pos_x(k)-r ,pos_y(k)+r ,pos_z(k)   ,matMPO)
		 call morfology(pos_x(k)-r ,pos_y(k)-r ,pos_z(k)   ,matMMO)
		 call morfology(pos_x(k)-r ,pos_y(k)   ,pos_z(k)+r ,matMOP)
		 call morfology(pos_x(k)-r ,pos_y(k)   ,pos_z(k)-r ,matMOM)
		 call morfology(pos_x(k)   ,pos_y(k)+r ,pos_z(k)   ,matOPO)
		 call morfology(pos_x(k)   ,pos_y(k)+r ,pos_z(k)+r ,matOPP)
		 call morfology(pos_x(k)   ,pos_y(k)+r ,pos_z(k)-r ,matOPM)
		 call morfology(pos_x(k)   ,pos_y(k)-r ,pos_z(k)   ,matOMO)
                 call morfology(pos_x(k)   ,pos_y(k)-r ,pos_z(k)+r ,matOMP)
		 call morfology(pos_x(k)   ,pos_y(k)-r ,pos_z(k)-r ,matOMM)
		 call morfology(pos_x(k)   ,pos_y(k)   ,pos_z(k)   ,matOOO)
		 call morfology(pos_x(k)   ,pos_y(k)   ,pos_z(k)+r ,matOOP)
		 call morfology(pos_x(k)   ,pos_y(k)   ,pos_z(k)-r ,matOOM)
		 call morfology(pos_x(k)+r ,pos_y(k)+r ,pos_z(k)+r ,matPPP)
		 call morfology(pos_x(k)+r ,pos_y(k)+r ,pos_z(k)-r ,matPPM)
		 call morfology(pos_x(k)+r ,pos_y(k)-r ,pos_z(k)+r ,matPMP)
		 call morfology(pos_x(k)+r ,pos_y(k)-r ,pos_z(k)-r ,matPMM)
		 call morfology(pos_x(k)-r ,pos_y(k)+r ,pos_z(k)+r ,matMPP)
		 call morfology(pos_x(k)-r ,pos_y(k)+r ,pos_z(k)-r ,matMPM)
		 call morfology(pos_x(k)-r ,pos_y(k)-r ,pos_z(k)+r ,matMMP)
		 call morfology(pos_x(k)-r ,pos_y(k)-r ,pos_z(k)-r ,matMMM)

!                call gaussian(x)
		
		kf_POO=(1.0_rp/temi(matOOO))*((rf(matOOO,matPOO)/(alpha*moment(matOOO)+r))**6)
		kf_PPO=(1.0_rp/temi(matOOO))*((rf(matOOO,matPPO)/(alpha*moment(matOOO)+r*sqrt(2.0_rp)))**6)
		kf_PMO=(1.0_rp/temi(matOOO))*((rf(matOOO,matPMO)/(alpha*moment(matOOO)+r*sqrt(2.0_rp)))**6)
		kf_POP=(1.0_rp/temi(matOOO))*((rf(matOOO,matPOP)/(alpha*moment(matOOO)+r*sqrt(2.0_rp)))**6)
		kf_POM=(1.0_rp/temi(matOOO))*((rf(matOOO,matPOM)/(alpha*moment(matOOO)+r*sqrt(2.0_rp)))**6)
		kf_MOO=(1.0_rp/temi(matOOO))*((rf(matOOO,matMOO)/(alpha*moment(matOOO)+r))**6)
		kf_MPO=(1.0_rp/temi(matOOO))*((rf(matOOO,matMPO)/(alpha*moment(matOOO)+r*sqrt(2.0_rp)))**6)
		kf_MMO=(1.0_rp/temi(matOOO))*((rf(matOOO,matMMO)/(alpha*moment(matOOO)+r*sqrt(2.0_rp)))**6)
		kf_MOP=(1.0_rp/temi(matOOO))*((rf(matOOO,matMOP)/(alpha*moment(matOOO)+r*sqrt(2.0_rp)))**6)
		kf_MOM=(1.0_rp/temi(matOOO))*((rf(matOOO,matMOM)/(alpha*moment(matOOO)+r*sqrt(2.0_rp)))**6)
		kf_OPO=(1.0_rp/temi(matOOO))*((rf(matOOO,matOPO)/(alpha*moment(matOOO)+r))**6)
		kf_OPP=(1.0_rp/temi(matOOO))*((rf(matOOO,matOPP)/(alpha*moment(matOOO)+r*sqrt(2.0_rp)))**6)
		kf_OPM=(1.0_rp/temi(matOOO))*((rf(matOOO,matOPM)/(alpha*moment(matOOO)+r*sqrt(2.0_rp)))**6)
		kf_OMO=(1.0_rp/temi(matOOO))*((rf(matOOO,matOMO)/(alpha*moment(matOOO)+r))**6)
		kf_OMP=(1.0_rp/temi(matOOO))*((rf(matOOO,matOMP)/(alpha*moment(matOOO)+r*sqrt(2.0_rp)))**6)
		kf_OMM=(1.0_rp/temi(matOOO))*((rf(matOOO,matOMM)/(alpha*moment(matOOO)+r*sqrt(2.0_rp)))**6)
		kf_OOP=(1.0_rp/temi(matOOO))*((rf(matOOO,matOOP)/(alpha*moment(matOOO)+r))**6)
		kf_OOM=(1.0_rp/temi(matOOO))*((rf(matOOO,matOOM)/(alpha*moment(matOOO)+r))**6)
		kf_PPP=(1.0_rp/temi(matOOO))*((rf(matOOO,matPPP)/(alpha*moment(matOOO)+r*sqrt(3.0_rp)))**6)
		kf_PPM=(1.0_rp/temi(matOOO))*((rf(matOOO,matPPM)/(alpha*moment(matOOO)+r*sqrt(3.0_rp)))**6)
		kf_PMP=(1.0_rp/temi(matOOO))*((rf(matOOO,matPMP)/(alpha*moment(matOOO)+r*sqrt(3.0_rp)))**6)
		kf_PMM=(1.0_rp/temi(matOOO))*((rf(matOOO,matPMM)/(alpha*moment(matOOO)+r*sqrt(3.0_rp)))**6)
		kf_MPP=(1.0_rp/temi(matOOO))*((rf(matOOO,matMPP)/(alpha*moment(matOOO)+r*sqrt(3.0_rp)))**6)
		kf_MPM=(1.0_rp/temi(matOOO))*((rf(matOOO,matMPM)/(alpha*moment(matOOO)+r*sqrt(3.0_rp)))**6)
		kf_MMP=(1.0_rp/temi(matOOO))*((rf(matOOO,matMMP)/(alpha*moment(matOOO)+r*sqrt(3.0_rp)))**6)
		kf_MMM=(1.0_rp/temi(matOOO))*((rf(matOOO,matMMM)/(alpha*moment(matOOO)+r*sqrt(3.0_rp)))**6)


                kf_tot = kf_POO + kf_PPO + kf_PMO + kf_POP + kf_POM + kf_MOO + kf_MPO + kf_MMO + kf_MOP + kf_MOM + kf_OPO + kf_OPP + kf_OPM + kf_OMO + kf_OMP + kf_OMM + kf_OOP + kf_OOM + kf_PPP + kf_PPM + kf_PMP + kf_PMM + kf_MPP + kf_MPM + kf_MMP + kf_MMM


		prob_POO=kf_POO/kf_tot
		prob_PPO=kf_PPO/kf_tot
		prob_PMO=kf_PMO/kf_tot
		prob_POP=kf_POP/kf_tot
		prob_POM=kf_POM/kf_tot
		prob_MOO=kf_MOO/kf_tot
		prob_MPO=kf_MPO/kf_tot
		prob_MMO=kf_MMO/kf_tot
		prob_MOP=kf_MOP/kf_tot
		prob_MOM=kf_MOM/kf_tot
		prob_OPO=kf_OPO/kf_tot
		prob_OPP=kf_OPP/kf_tot
		prob_OPM=kf_OPM/kf_tot
		prob_OMO=kf_OMO/kf_tot
		prob_OMP=kf_OMP/kf_tot
		prob_OMM=kf_OMM/kf_tot
		prob_OOP=kf_OOP/kf_tot
		prob_OOM=kf_OOM/kf_tot
		prob_PPP=kf_PPP/kf_tot
		prob_PPM=kf_PPM/kf_tot
		prob_PMP=kf_PMP/kf_tot
		prob_PMM=kf_PMM/kf_tot
		prob_MPP=kf_MPP/kf_tot
		prob_MPM=kf_MPM/kf_tot
		prob_MMP=kf_MMP/kf_tot
		prob_MMM=kf_MMM/kf_tot

                call random_number(num)
            
                if (num .lt. prob_POO) then
                    poshx=pos_hx(k) + r
                    poshy=pos_hy(k)
		    poshz=pos_hz(k)
                    dt(k)=1.0_rp/kf_POO
                    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_POO)
                end if
		
		if (num .ge. prob_POO .and. num .lt. prob_POO+prob_PPO) then
		    poshx=pos_hx(k) + r
		    poshy=pos_hy(k) + r
		    poshz=pos_hz(k) 
		    dt(k)=1.0_rp/kf_PPO
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_PPO)
		end if

		if (num .ge. prob_POO+prob_PPO .and. num .lt. prob_POO+prob_PPO+prob_PMO) then
		    poshx=pos_hx(k) + r
		    poshy=pos_hy(k) - r
		    poshz=pos_hz(k) 
		    dt(k)=1.0_rp/kf_PMO
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_PMO)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP) then
		    poshx=pos_hx(k) + r
		    poshy=pos_hy(k) 
		    poshz=pos_hz(k) + r
		    dt(k)=1.0_rp/kf_POP
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_POP)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM) then
		    poshx=pos_hx(k) + r
		    poshy=pos_hy(k) 
		    poshz=pos_hz(k) - r
		    dt(k)=1.0_rp/kf_POM
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_POM)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO) then
		    poshx=pos_hx(k) - r
		    poshy=pos_hy(k) 
		    poshz=pos_hz(k) 
		    dt(k)=1.0_rp/kf_MOO
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_MOO)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO) then
		    poshx=pos_hx(k) - r
		    poshy=pos_hy(k) + r
		    poshz=pos_hz(k) 
		    dt(k)=1.0_rp/kf_MPO
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_MPO)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO) then
		    poshx=pos_hx(k) - r
		    poshy=pos_hy(k) - r
		    poshz=pos_hz(k) 
		    dt(k)=1.0_rp/kf_MMO
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_MMO)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP) then
		    poshx=pos_hx(k) - r
		    poshy=pos_hy(k) 
		    poshz=pos_hz(k) + r
		    dt(k)=1.0_rp/kf_MOP
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_MOP)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM) then
		    poshx=pos_hx(k) - r
		    poshy=pos_hy(k) 
		    poshz=pos_hz(k) - r
		    dt(k)=1.0_rp/kf_MOM
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_MOM)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO) then
		    poshx=pos_hx(k) 
		    poshy=pos_hy(k) + r
		    poshz=pos_hz(k) 
		    dt(k)=1.0_rp/kf_OPO
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_OPO)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP) then
		    poshx=pos_hx(k) 
		    poshy=pos_hy(k) + r
		    poshz=pos_hz(k) + r 
		    dt(k)=1.0_rp/kf_OPP
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_OPP)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM) then
		    poshx=pos_hx(k) 
		    poshy=pos_hy(k) + r
		    poshz=pos_hz(k) - r 
		    dt(k)=1.0_rp/kf_OPM
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_OPM)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO) then
		    poshx=pos_hx(k) 
		    poshy=pos_hy(k) - r
		    poshz=pos_hz(k)  
		    dt(k)=1.0_rp/kf_OMO
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_OMO)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP) then
		    poshx=pos_hx(k) 
		    poshy=pos_hy(k) - r
		    poshz=pos_hz(k) + r 
		    dt(k)=1.0_rp/kf_OMP
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_OMP)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM) then
		    poshx=pos_hx(k) 
		    poshy=pos_hy(k) - r
		    poshz=pos_hz(k) - r 
		    dt(k)=1.0_rp/kf_OMM
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_OMM)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP) then
		    poshx=pos_hx(k) 
		    poshy=pos_hy(k) 
		    poshz=pos_hz(k) + r
		    dt(k)=1.0_rp/kf_OOP
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_OOP)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM) then
		    poshx=pos_hx(k) 
		    poshy=pos_hy(k) 
		    poshz=pos_hz(k) - r
		    dt(k)=1.0_rp/kf_OOM
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_OOM)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP) then
		    poshx=pos_hx(k) + r
		    poshy=pos_hy(k) + r
		    poshz=pos_hz(k) + r
		    dt(k)=1.0_rp/kf_PPP
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_PPP)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM) then
		    poshx=pos_hx(k) + r
		    poshy=pos_hy(k) + r
		    poshz=pos_hz(k) - r
		    dt(k)=1.0_rp/kf_PPM
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_PPM)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM+prob_PMP) then
		    poshx=pos_hx(k) + r
		    poshy=pos_hy(k) - r
		    poshz=pos_hz(k) + r
		    dt(k)=1.0_rp/kf_PMP
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_PMP)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM+prob_PMP .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM+prob_PMP+prob_PMM) then
		    poshx=pos_hx(k) + r
		    poshy=pos_hy(k) - r
		    poshz=pos_hz(k) - r
		    dt(k)=1.0_rp/kf_PMM
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_PMM)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM+prob_PMP+prob_PMM .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM+prob_PMP+prob_PMM+prob_MPP) then
		    poshx=pos_hx(k) - r
		    poshy=pos_hy(k) + r
		    poshz=pos_hz(k) + r
		    dt(k)=1.0_rp/kf_MPP
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_MPP)
		end if


		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM+prob_PMP+prob_PMM+prob_MPP .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM+prob_PMP+prob_PMM+prob_MPP+prob_MPM) then
		    poshx=pos_hx(k) - r
		    poshy=pos_hy(k) + r
		    poshz=pos_hz(k) - r
		    dt(k)=1.0_rp/kf_MPM
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_MPM)
		end if

		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM+prob_PMP+prob_PMM+prob_MPP+prob_MPM .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM+prob_PMP+prob_PMM+prob_MPP+prob_MPM+prob_MMP) then
		    poshx=pos_hx(k) - r
		    poshy=pos_hy(k) - r
		    poshz=pos_hz(k) + r
		    dt(k)=1.0_rp/kf_MMP
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_MMP)
		end if


		if (num .ge. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM+prob_PMP+prob_PMM+prob_MPP+prob_MPM+prob_MMP .and. num .lt. prob_POO+prob_PPO+prob_PMO+prob_POP+prob_POM+prob_MOO+prob_MPO+prob_MMO+prob_MOP+prob_MOM+prob_OPO+prob_OPP+prob_OPM+prob_OMO+prob_OMP+prob_OMM+prob_OOP+prob_OOM+prob_PPP+prob_PPM+prob_PMP+prob_PMM+prob_MPP+prob_MPM+prob_MMP+prob_MMM) then
		    poshx=pos_hx(k) - r
		    poshy=pos_hy(k) - r
		    poshz=pos_hz(k) - r
		    dt(k)=1.0_rp/kf_MMM
		    probm=k_emi(matOOO)/(k_emi(matOOO)+kf_MMM)
		end if

                call random_number(num)
                if(num .lt. probm) then
                    write(10,101) t, pos_x(k), pos_y(k), pos_z(k),"Fluor", matOOO
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
end program rw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine exciton_distribution
use var
implicit none
integer               :: i
real(rp),dimension(3) :: m

do i=1,NE
   !pos_x(i)=100.0_rp
   !pos_y(i)=0.0_rp
   !pos_z(i)=0.0_rp
   call random_number(m)
   pos_x(i)=m(1)*(xf-xi)+xi
   pos_y(i)=m(2)*(yf-yi)+yi
   pos_z(i)=m(3)*(zf-zi)+zi
enddo

return
end subroutine exciton_distribution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine morfology(x,y,z,material)
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
end subroutine morfology

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
         if(iev(j) .ne. 0) then                    !até aqui ele confere se os dois excitons estao vivos
             call morfology(pos_x(j),pos_y(j),pos_z(j),material)
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
