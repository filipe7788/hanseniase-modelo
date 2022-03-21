c  para compilar os arquivos: gfortran hansr1-modif.f -o hansr1-M.exe


c     hansr0 - programa para resolver o modelo da hanseniase sem termos de difusao

      external fct,outp
      dimension y(5)
      dimension prmt(5),aux(8,5),dery(5)
	real mu,gama,betaP,betaM,rho,sigma,omegaP,omegaM,
     *deltaP,deltaM
	character*30 saida
      common mu,gama,betaP,betaM,rho,sigma,omegaP,omegaM,
     *deltaP,deltaM

      open (unit=3,file='ehansr0.dat',status='old')
      read (3,*)tf,xint,rh,tprin,prmt(4)
	read (3,*)mu,gama,betaP,betaM,rho,sigma,omegaP,omegaM,
     *deltaP,deltaM
	read (3,*)(y(i),i=1,5)
	read (3,510)saida

      close(unit=3)

c	x = 0.

      open (unit=5,file=saida)

      write(5,*)'# resultados do programa hansr0'

      write(5,520)tf,xint,rh,prmt(4)
      write(5,530)mu,gama,betaP,betaM,rho,sigma,omegaP,omegaM,
     *deltaP,deltaM
      write(5,540)
c      write(5,*)(y(i),i=1,5)

      t = 0
      prmt(1) = t
      prmt(2) = t + rh
      prmt(3) = xint

      do i = 1,5
       dery(i) = 0.2
      enddo

100   write(5,300)prmt(1),(y(i),i=1,5),y(3)+y(4)

c     write(*,300)(prmt(i),i=1,5),(y(i),i=1,5)
c	pause

      call rkgs(prmt,y,dery,5,ihlf,fct,outp,aux)

      prmt(1)=prmt(2)
      prmt(2)=prmt(2)+rh
      if(prmt(2).gt.(tf+rh))go to 110
      go to 100


510   format (a30)
520   format(/'#  parametros '
     ./'#    tf     tint     rh       prec'
     ./ 5(e10.3,2x))
530   format(
     . '#    mu     gama    betaP     betaM         rho       sigma'
     ./ 6(e10.3,2x)
     ./'#     omegaP    omegaM       deltaP    deltaM'
     ./ 4(e10.3,2x)/)
540   format(' trajetoria')
300	format (1x,7(e14.7,2x))

110   write(5,300)prmt(1),(y(i),i=1,5)

      stop
      end


      subroutine fct(x,y,dery)
      real x,dery(5),y(5)
      	real mu,gama,betaP,betaM,rho,sigma,omegaP,omegaM,deltaP,deltaM
      common mu,gama,betaP,betaM,rho,sigma,omegaP,omegaM,deltaP,deltaM

c	x = 0.

      dery(1) = mu*(1-gama)-(betaM*y(3)+betaP*y(4))*y(1)-
     *mu*y(1)

      dery(2) =(betaM*y(3)+betaP*y(4))*y(1)-mu*y(2)-(omegaP*
     *sigma+omegaM*(1-sigma))*y(2)
      
      dery(3) = rho*sigma*omegaM*y(2)-mu*y(3)-deltaM*y(3)
      
      dery(4) = rho*(1-sigma)*omegaP*y(2)-mu*y(4)-deltaP*y(4)

      dery(5) = (1-rho)*sigma*omegaP*y(2)+(1-rho)*(1-sigma)*omegaM 
     **y(2)+deltaM*y(3)+deltaP*y(4)-mu*y(5)+mu*gama 
       
700   format(1x,3(i4,1x),4(e13.6))
      return
      end

      subroutine outp(x,ihlf)
c      subroutine outp(x,y,dery,ihlf,ndim,prmt)
c      real y(256,128),dery(256,128),prmt(5)
      if(ihlf.ne.11)go to 500
   20 write(*,400)x
  400 format(1x,'erro na solucao em t =',e15.8)
      write(*,*)'sai do outp'
  500 return
      end


      subroutine rkgs(prmt,yaux,deryau,ndim,ihlf,fct,outp,aux)
      dimension yaux(1),deryau(1),aux(8,1),a(4),b(4),c(4),prmt(5)
c	write(*,*)'entra rkgs'
      do 1 i=1,ndim
    1 aux(8,i)=.06666667*deryau(i)
      x=prmt(1)
      xend=prmt(2)
      h=prmt(3)
      prmt(5)=0.
      call fct(x,yaux,deryau)
      if(h*(xend-x))38,37,2
    2 a(1)=.5
      a(2)=.2928932
      a(3)=1.707107
      a(4)=.1666667
      b(1)=2.
      b(2)=1.
      b(3)=1.
      b(4)=2.
      c(1)=.5
      c(2)=.2928932
      c(3)=1.707107
      c(4)=.5
      do 3 i=1,ndim
      aux(1,i)=yaux(i)
      aux(2,i)=deryau(i)
      aux(3,i)=0.
    3 aux(6,i)=0.
      irec=0
      h=h+h
      ihlf=-1
      istep=0
      iend=0
    4 if((x+h-xend)*h)7,6,5
    5 h=xend-x
    6 iend=1
    7 call outp(x,ihlf)
c    7 call outp(x,yaux,deryau,irec,ndim,prmt)
c      write(*,*)'voltei do prim. outp'
      if(prmt(5))40,8,40
    8 itest=0
    9 istep=istep+1
      j=1
   10 aj=a(j)
      bj=b(j)
      cj=c(j)
      do 11 i=1,ndim
      r1=h*deryau(i)
      r2=aj*(r1-bj*aux(6,i))
      yaux(i)=yaux(i)+r2
      r2=r2+r2+r2
   11 aux(6,i)=aux(6,i)+r2-cj*r1
      if(j-4)12,15,15
   12 j=j+1
      if(j-3)13,14,13
   13 x=x+.5*h
   14 call fct(x,yaux,deryau)
      goto 10
   15 if(itest)16,16,20
   16 do 17 i=1,ndim
   17 aux(4,i)=yaux(i)
      itest=1
      istep=istep+istep-2
   18 ihlf=ihlf+1
      x=x-h
      h=.5*h
      do 19 i=1,ndim
      yaux(i)=aux(1,i)
      deryau(i)=aux(2,i)
   19 aux(6,i)=aux(3,i)
      goto 9
   20 imod=istep/2
      if(istep-imod-imod)21,23,21
   21 call fct(x,yaux,deryau)
      do 22 i=1,ndim
      aux(5,i)=yaux(i)
   22 aux(7,i)=deryau(i)
      goto 9
   23 delt=0.
      do 24 i=1,ndim
   24 delt=delt+aux(8,i)*abs(aux(4,i)-yaux(i))
      if(delt-prmt(4))28,28,25
   25 if(ihlf-10)26,36,36
   26 do 27 i=1,ndim
   27 aux(4,i)=aux(5,i)
      istep=istep+istep-4
      x=x-h
      iend=0
      goto 18
   28 call fct(x,yaux,deryau)
      do 29 i=1,ndim
      aux(1,i)=yaux(i)
      aux(2,i)=deryau(i)
      aux(3,i)=aux(6,i)
      yaux(i)=aux(5,i)
   29 deryau(i)=aux(7,i)
      call outp(x-h,ihlf)
c      call outp(x-h,yaux,deryau,ihlf,ndim,prmt)
c      write(*,*)'voltei do seg. outp'
      if(prmt(5))40,30,40
   30 do 31 i=1,ndim
      yaux(i)=aux(1,i)
   31 deryau(i)=aux(2,i)
      irec=ihlf
      if(iend)32,32,39
   32 ihlf=ihlf-1
      istep=istep/2
      h=h+h
      if(ihlf)4,33,33
   33 imod=istep/2
      if(istep-imod-imod)4,34,4
   34 if(delt-.02*prmt(4))35,35,4
   35 ihlf=ihlf-1
      istep=istep/2
      h=h+h
      goto 4
   36 ihlf=11
      call fct(x,yaux,deryau)
      goto 39
   37 ihlf=12
      goto 39
   38 ihlf=13
   39 call outp(x,ihlf)
c   39 call outp(x,yaux,deryau,ihlf,ndim,prmt)
c      write(*,*)'voltei do terc. outp'
   40 return
      end
      irec=ihlf
      if(iend)32,32,39
   32 ihlf=ihlf-1
      istep=istep/2
      h=h+h
      if(ihlf)4,33,33
   33 imod=istep/2
      if(istep-imod-imod)4,34,4
   34 if(delt-.02*prmt(4))35,35,4
   35 ihlf=ihlf-1
      istep=istep/2
      h=h+h
      goto 4
   36 ihlf=11
      call fct(x,yaux,deryau)
      goto 39
   37 ihlf=12
      goto 39
   38 ihlf=13
   39 call outp(x,ihlf)
c   39 call outp(x,yaux,deryau,ihlf,ndim,prmt)
c      write(*,*)'voltei do terc. outp'
   40