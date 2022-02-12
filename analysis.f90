!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     Simulação em Monte Carlo para o Random Walk 1D   !!
!! Implementação de Taxas de FRET e emissão de excitons !!
!! Implementação de diferentes materiais                !!
!! Implementação da desordem posicional com gerador de  !!
!! gaussiana                                            !!
!! Cálculo de omega prara a largura da gaussiana        !!
!! Implementação da aniquilação                         !!
!! Programa para análise de TRPL                        !!
!!                                                      !!
!! 09/06/2021                                           !!
!!         Laura Simonassi                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program analysis
implicit none
integer      :: i,NE,IOSTAT,IOS

integer,parameter                   :: n_bin=100

real(kind=8)                        :: total_sum, meanpos, variancepos,meanx, variancex,meany, variancey,meanz, variancez, diffusion, ldpos,ldx,ldy,ldz, dbin
real,dimension(:),allocatable       :: time,pos,x,y,z
real(kind=8),dimension(n_bin)       :: bound

total_sum = 0.0d0
meanpos      = 0.0d0
variancepos  = 0.0d0
diffusion = 0.0d0
ldpos        = 0.0d0
meanx      = 0.0d0
variancex  = 0.0d0
ldx        = 0.0d0
meany      = 0.0d0
variancey  = 0.0d0
ldy        = 0.0d0
meanz      = 0.0d0
variancez  = 0.0d0
ldz        = 0.0d0


open(1, file ='Simulation.txt', status= 'old')
   NE=-2
   do
      read(1,*,IOSTAT=IOS)
      NE=NE+1
      if(IOS .lt. 0)exit
   end do
close(1)

allocate(time(NE),pos(NE),x(NE),y(NE),z(NE))

open(1, file ='Simulation.txt', status= 'old')
   do i=1,NE
      read(1,*)time(i),pos(i),x(i),y(i),z(i)
   enddo
close(1)

meanpos=sum(pos)/NE
variancepos=sum((pos-meanpos)**2)/dfloat(NE)
ldpos= dsqrt(variancepos)

meanx=sum(x)/NE
variancex=sum((x-meanx)**2)/dfloat(NE)
ldx= dsqrt(variancex)

meany=sum(y)/NE
variancey=sum((y-meany)**2)/dfloat(NE)
ldy= dsqrt(variancey)

meanz=sum(z)/NE
variancez=sum((z-meanz)**2)/dfloat(NE)
ldz= dsqrt(variancez)

open (2, file='analysis.txt', status='unknown')
    write(2,*) ldpos,ldx,ldy,ldz
close(2)

!dbin=(8.0d0*1e-8)/dfloat(n_bin)                !passo do bin tornar o nbin um float dfloat()
dbin=(maxval(time))/dfloat(n_bin)            !passo do bin tornar o nbin um float dfloat()
bound(1)=minval(time)+dbin                    !bound são os intervalos
do i=1,n_bin-1
   bound(i+1)=bound(i)+dbin
enddo

open(1,file="TRPL.txt")
      write(1,10)0, count(time <= bound(1))                                  !bound(1), count(tempo <= bound(1))
   do i=1,n_bin-1
      write(1,10)bound(i), count(time <= bound(i+1))-count(time <= bound(i))
   enddo
close(1)

10 format(E19.10,i10)
end program analysis


