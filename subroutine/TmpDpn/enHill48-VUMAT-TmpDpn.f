c====================================================================
c          Program for 3D anisotropic hardening 
c          Anisotropic Hill48 yield function
c          Hill7 - anisotropic hardening & r-value evolution
c          by Junhe Lian
c          junhe.lian@iehk.rwth-aachen.de
c          February 2016 
c          Abaqus 6.11-1
c
c          !DO NOT DISTRIBUTE WITHOUT AUTHOR'S PERMISSION!
c====================================================================
c
       subroutine vumat(
C Read only -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     3  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     5  stressNew, stateNew, enerInternNew, enerInelasNew)
C
      include 'vaba_param.inc'
C
C All arrays dimensioned by (*) are not used in this algorithm
      dimension props(nprops), density(nblock),
     1  coordMp(nblock,*),
     2  charLength(*), strainInc(nblock,ndir+nshr),
     3  relSpinInc(*), tempOld(*),
     4  stretchOld(*), defgradOld(*),
     5  fieldOld(*), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(*),
     8  stretchNew(*), defgradNew(*), fieldNew(*),
     9  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     1  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C
C Defining numerical constants
      parameter(zero=0.,one=1.,two=2.,three=3.,four=4.,six=6.,
     1  half=0.5,threeHalves = 1.5,third = one/three,Tol=1.D-8,
     1  twoThirds=two/three,const1=sqrt(threeHalves),NTENS=6,NDI=3)
      REAL A0,e0,xn0,A45,e45,xn45,A90,e90,xn90,AEB,eEB,xnEB,r0,r45,r90,
     1  rEB,F,G,H,L,M,N,F1,G1,H1,L1,M1,N1,C1,C2,C3,C4,C5,C6,C7,C8,C9	 
C     User defined tensors
      DIMENSION dfds(6),Ydirection(6),dgds(6),DSTRESS(6),
     1 ST(6),EELAS(6),EPLAS(6),EELT(6),EPLT(6),STRESST(6),
     2 EELAST(6),EPLAST(6),DDSDDE(6,6),DSTRAN(6) 
C Define material properties constants
c ELASTIC properties
       EMOD = PROPS(1)
       ENU  = PROPS(2)
c ELASTIC STIFFNESS
       EBULK3=EMOD/(one-two*ENU)
       EG2=EMOD/(one+ENU)
       EG=EG2/two
       EG3=three*EG
       ELAM=(EBULK3-EG2)/three
	   PI= four*atan(one)
c Define dimension of flow curve and r-value input table	   
       nvalue = (nprops-17)/17	   
C----------------------------------------------------------------------C
C Computation per material point starts here
      do 9000 i = 1,nblock
C Fill in the elastic stiffness tensor
      DO 20 K1=1,NTENS
        DO 10 K2=1,NTENS
          DDSDDE(K2,K1)=zero
 10     CONTINUE
 20   CONTINUE
C
      DO 40 K1=1,NDI
        DO 30 K2=1,NDI
          DDSDDE(K2,K1)=ELAM
 30     CONTINUE
        DDSDDE(K1,K1)=EG2+ELAM
 40   CONTINUE
      DO 50 K1=NDI+1,NTENS
        DDSDDE(K1,K1)=EG
 50   CONTINUE
C  Elastic predictor and plastic corrector scheme in conjuction with the
C  Tangent Cutting Plane (TCP) algorithm 
C  Elastic predictor
      eqplastrue = stateOld(i,1)
      eqplas = eqplastrue
C Transite DSTRAN from strainInc, note diff. between UMAT and VUMAT.  
      DSTRAN(1)=strainInc(i,1)
      DSTRAN(2)=strainInc(i,2)
      DSTRAN(3)=strainInc(i,3)
      DSTRAN(4)=strainInc(i,4)*two
      DSTRAN(5)=strainInc(i,5)*two
      DSTRAN(6)=strainInc(i,6)*two
C Calculate Trial stress
      DO 70 K1=1,NTENS
          STRESST(K1)=stressOld(i,K1)+DDSDDE(K1,1)*DSTRAN(1)
     1   +DDSDDE(K1,2)*DSTRAN(2)+DDSDDE(K1,3)*DSTRAN(3)
     1   +DDSDDE(K1,4)*DSTRAN(4)+DDSDDE(K1,5)*DSTRAN(5)
     1   +DDSDDE(K1,6)*DSTRAN(6)
 70   CONTINUE
       ST1=STRESST(1)
       ST2=STRESST(2)
       ST3=STRESST(3)
       ST4=STRESST(4)
       ST5=STRESST(5)
       ST6=STRESST(6)  
C Calculate anisotropic parameters in Hill48 based on flow stresses	          
C Fetching yield stress from flow curve
       call ahard(sigmayield0,hard0,sigmayield45,sigmayield90,sigmayieldEB,
     1	  rvalue0,rvalue45,rvalue90,eqplas,PROPS(9),nvalue)	   
       sigmayield=sigmayield0
	   xhard=hard0
	   s0=sigmayield
       s45=sigmayield45
       s90=sigmayield90
       sEB=sigmayieldEB
       F=((s0**2)/(s90**2)-one+(s0**2/sEB**2))/two
       G=(one-(s0**2/s90**2)+(s0**2/sEB**2))/two
       H=(one+(s0**2/s90**2)-(s0**2/sEB**2))/two
       N=((four*s0**2/s45**2)-(s0**2/sEB**2))/two
	   L=threeHalves
	   M=threeHalves
C Calculate subvalue of the trial Hill48 equivalent stress		   
       subsigtrial= F*(ST2-ST3)**2+G*(ST3-ST1)**2
     1             +H*(ST1-ST2)**2+two*N*ST4**2
     2             +two*L*ST5**2+two*M*ST6**2
c To prevent subsigtrial =0
       if(subsigtrial.lt.Tol) subsigtrial = Tol
C Calculate the trial Hill48 equivalent stress 
       sigtrial = sqrt(subsigtrial)
c To prevent sigtrial =0	   
       if(sigtrial.lt.Tol .and. sigtrial.ge.zero) sigtrial = Tol
       if(sigtrial.gt.(zero-Tol) .and. sigtrial.le.zero) sigtrial = zero-Tol	   
C Calculate the temperature dependent function
       T = stateOld(i, 3)
       if (T.le.Tol .and. T.ge.(zero-Tol)) then
          T = T0
       endif
       fT = cT1*exp(zero-cT2*T)+cT3+cT4*exp(zero-((T-cT5)/cT6)**2)+cT7*exp(zero-((T-cT8)/cT9)**2)
C Check yielding
       radius = sigmayield*fT
	   Xhardnew= xhard*fT
       xf = sigtrial - radius
       if(xf.le.(Tol)) then            		 
          goto 1000
       endif
       dgama=zero
       iflag=zero
C Plastic corrector and Begin iterations
5000   continue
       iflag=iflag+one   
c Calculate anisotropic parameters in Hill48 based on flow stresses
C Fetching yield stress from flow curve
       call ahard(sigmayield0,hard0,sigmayield45,sigmayield90,sigmayieldEB,
     1	  rvalue0,rvalue45,rvalue90,eqplas,PROPS(9),nvalue)	   
       sigmayield=sigmayield0
	   xhard=hard0
	   s0=sigmayield
       s45=sigmayield45
       s90=sigmayield90
       sEB=sigmayieldEB
       r0=rvalue0
       r45=rvalue45
       r90=rvalue90
       F=((s0**2)/(s90**2)-one+(s0**2/sEB**2))/two
       G=(one-(s0**2/s90**2)+(s0**2/sEB**2))/two
       H=(one+(s0**2/s90**2)-(s0**2/sEB**2))/two
       N=((four*s0**2/s45**2)-(s0**2/sEB**2))/two
	   L=threeHalves
	   M=threeHalves
C Calculate subvalue of the trial Hill48 equivalent stress		   
       subsigtrial=F*(ST2-ST3)**2+G*(ST3-ST1)**2
     1             +H*(ST1-ST2)**2+two*N*ST4**2
     2             +two*L*ST5**2+two*M*ST6**2
c To prevent subsigtrial =0	 
       if(subsigtrial.lt.Tol) subsigtrial = Tol
C Calculate the trial Hill48 equivalent stress	   
       sigtrial = sqrt(subsigtrial)
c To prevent sigtrial =0	   
       if(sigtrial.lt.Tol .and. sigtrial.ge.zero) sigtrial = Tol
       if(sigtrial.gt.(zero-Tol) .and. sigtrial.le.zero) sigtrial = zero-Tol
c Calculate anisotropic parameters in Hill48 based on r-values	 
       F1=r0/(r90)/(one+r0)
       G1=one/(one+r0)
       H1=r0/(one+r0)
       L1=threeHalves
       M1=threeHalves
       N1=(r0+r90)*(one+two*r45)/two/r90/(one+r0)
C Calculate subvalue of the trial Hill48 flow potential
	   subsigmaHill=F1*(ST2-ST3)**2+G1*(ST3-ST1)**2
     1          +H1*(ST1-ST2)**2+two*N1*ST4**2
     2          +two*L1*ST5**2+two*M1*ST6**2
c To prevent subsigmaHill =0	 
       if(subsigmaHill.lt.Tol) subsigmaHill = Tol
C Calculate the trial Hill48 flow potential	   
       sigmaHill=sqrt(subsigmaHill)
c To prevent sigmaHill =0	   
       if(sigmaHill.lt.Tol .and. sigmaHill.ge.zero) sigmaHill = Tol
       if(sigmaHill.gt.(zero-Tol) .and. sigmaHill.le.zero) sigmaHill = zero-Tol		   
C Update parameters with Non-AFR and calculate sigmayield and sigmaHill at every iteration 
C Calculate the flow direction vector based on Flow potential
       dgds(1)=(H1*(ST1-ST2)-G1*(ST3-ST1))/sigmaHill
       dgds(2)=(F1*(ST2-ST3)-H1*(ST1-ST2))/sigmaHill
       dgds(3)=(G1*(ST3-ST1)-F1*(ST2-ST3))/sigmaHill    
       dgds(4)=(two*N1*ST4)/sigmaHill
       dgds(5)=(two*L1*ST5)/sigmaHill
       dgds(6)=(two*M1*ST6)/sigmaHill
C Calculate the flow direction vector based on Yield function 
       dfds(1)=(H*(ST1-ST2)-G*(ST3-ST1))/sigtrial
       dfds(2)=(F*(ST2-ST3)-H*(ST1-ST2))/sigtrial
       dfds(3)=(G*(ST3-ST1)-F*(ST2-ST3))/sigtrial 
       dfds(4)=(two*N*ST4)/sigtrial
       dfds(5)=(two*L*ST5)/sigtrial
       dfds(6)=(two*M*ST6)/sigtrial
C Calculate the flow direction in the TCP algorithm
       Ydirection(1)=dgds(1)*DDSDDE(1,1)+dgds(2)*DDSDDE(1,2)
     1              +dgds(3)*DDSDDE(1,3)+dgds(4)*DDSDDE(1,4)
     1              +dgds(5)*DDSDDE(1,5)+dgds(6)*DDSDDE(1,6)
       Ydirection(2)=dgds(1)*DDSDDE(2,1)+dgds(2)*DDSDDE(2,2)
     1              +dgds(3)*DDSDDE(2,3)+dgds(4)*DDSDDE(2,4)
     1              +dgds(5)*DDSDDE(2,5)+dgds(6)*DDSDDE(2,6)
       Ydirection(3)=dgds(1)*DDSDDE(3,1)+dgds(2)*DDSDDE(3,2)
     1              +dgds(3)*DDSDDE(3,3)+dgds(4)*DDSDDE(3,4)
     1              +dgds(5)*DDSDDE(3,5)+dgds(6)*DDSDDE(3,6)
       Ydirection(4)=dgds(1)*DDSDDE(4,1)+dgds(2)*DDSDDE(4,2)
     1              +dgds(3)*DDSDDE(4,3)+dgds(4)*DDSDDE(4,4)
     1              +dgds(5)*DDSDDE(4,5)+dgds(6)*DDSDDE(4,6)
       Ydirection(5)=dgds(1)*DDSDDE(5,1)+dgds(2)*DDSDDE(5,2)
     1              +dgds(3)*DDSDDE(5,3)+dgds(4)*DDSDDE(5,4)
     1              +dgds(5)*DDSDDE(5,5)+dgds(6)*DDSDDE(5,6)
       Ydirection(6)=dgds(1)*DDSDDE(6,1)+dgds(2)*DDSDDE(6,2)
     1              +dgds(3)*DDSDDE(6,3)+dgds(4)*DDSDDE(6,4)
     1              +dgds(5)*DDSDDE(6,5)+dgds(6)*DDSDDE(6,6)
C Calculate the first term in the TCP algorithm
c dsdgama caculation=CNN, dot product, prepare for iteration
C Calculate MCN
       dsdga=dfds(1)*Ydirection(1)+dfds(2)*Ydirection(2)
     1      +dfds(3)*Ydirection(3)+dfds(4)*Ydirection(4)
     1      +dfds(5)*Ydirection(5)+dfds(6)*Ydirection(6)
C Calculate the second term in the TCP algorithm
C Correspongding flow direction, Nf
	   radius = sigmayield*fT
	   Xhardnew= xhard*fT
	   Xdeltagamma = dsdga+Xhardnew	   
c To prevent Xdeltagamma =0	   
       if(Xdeltagamma.lt.Tol .and. Xdeltagamma.ge.zero) Xdeltagamma = Tol
       if(Xdeltagamma.gt.(zero-Tol) .and. Xdeltagamma.le.zero) Xdeltagamma=zero-Tol	   
	   xddgamma = xf/Xdeltagamma
       dgama = dgama+xddgamma
c Update trial stress with new dgama
      STRESST(1)=stressOld(i,1)+(DDSDDE(1,1)*(DSTRAN(1)-dgama*dgds(1))
     1          +DDSDDE(1,2)*(DSTRAN(2)-dgama*dgds(2))
     1          +DDSDDE(1,3)*(DSTRAN(3)-dgama*dgds(3))
     1          +DDSDDE(1,4)*(DSTRAN(4)-dgama*dgds(4))
     1          +DDSDDE(1,5)*(DSTRAN(5)-dgama*dgds(5))
     1          +DDSDDE(1,6)*(DSTRAN(6)-dgama*dgds(6)))
      STRESST(2)=stressOld(i,2)+(DDSDDE(2,1)*(DSTRAN(1)-dgama*dgds(1))
     1          +DDSDDE(2,2)*(DSTRAN(2)-dgama*dgds(2))
     1          +DDSDDE(2,3)*(DSTRAN(3)-dgama*dgds(3))
     1          +DDSDDE(2,4)*(DSTRAN(4)-dgama*dgds(4))
     1          +DDSDDE(2,5)*(DSTRAN(5)-dgama*dgds(5))
     1          +DDSDDE(2,6)*(DSTRAN(6)-dgama*dgds(6)))	  
      STRESST(3)=stressOld(i,3)+(DDSDDE(3,1)*(DSTRAN(1)-dgama*dgds(1))
     1          +DDSDDE(3,2)*(DSTRAN(2)-dgama*dgds(2))
     1          +DDSDDE(3,3)*(DSTRAN(3)-dgama*dgds(3))
     1          +DDSDDE(3,4)*(DSTRAN(4)-dgama*dgds(4))
     1          +DDSDDE(3,5)*(DSTRAN(5)-dgama*dgds(5))
     1          +DDSDDE(3,6)*(DSTRAN(6)-dgama*dgds(6)))
      STRESST(4)=stressOld(i,4)+(DDSDDE(4,1)*(DSTRAN(1)-dgama*dgds(1))
     1          +DDSDDE(4,2)*(DSTRAN(2)-dgama*dgds(2))
     1          +DDSDDE(4,3)*(DSTRAN(3)-dgama*dgds(3))
     1          +DDSDDE(4,4)*(DSTRAN(4)-dgama*dgds(4))
     1          +DDSDDE(4,5)*(DSTRAN(5)-dgama*dgds(5))
     1          +DDSDDE(4,6)*(DSTRAN(6)-dgama*dgds(6)))
      STRESST(5)=stressOld(i,5)+(DDSDDE(5,1)*(DSTRAN(1)-dgama*dgds(1))
     1          +DDSDDE(5,2)*(DSTRAN(2)-dgama*dgds(2))
     1          +DDSDDE(5,3)*(DSTRAN(3)-dgama*dgds(3))
     1          +DDSDDE(5,4)*(DSTRAN(4)-dgama*dgds(4))
     1          +DDSDDE(5,5)*(DSTRAN(5)-dgama*dgds(5))
     1          +DDSDDE(5,6)*(DSTRAN(6)-dgama*dgds(6)))
      STRESST(6)=stressOld(i,6)+(DDSDDE(6,1)*(DSTRAN(1)-dgama*dgds(1))
     1          +DDSDDE(6,2)*(DSTRAN(2)-dgama*dgds(2))
     1          +DDSDDE(6,3)*(DSTRAN(3)-dgama*dgds(3))
     1          +DDSDDE(6,4)*(DSTRAN(4)-dgama*dgds(4))
     1          +DDSDDE(6,5)*(DSTRAN(5)-dgama*dgds(5))
     1          +DDSDDE(6,6)*(DSTRAN(6)-dgama*dgds(6)))
c Update stress tensors
      ST1=STRESST(1)
      ST2=STRESST(2)
      ST3=STRESST(3)
      ST4=STRESST(4)
      ST5=STRESST(5)
      ST6=STRESST(6)		  
C Update the equivalent plastic strain	  
      eqplastrue = stateOld(i,1)
      eqplas = eqplastrue+dgama*(sigmaHill/sigtrial)
	  deqplas = dgama*(sigmaHill/sigtrial)	  
C Update anisotropic parameters in Hill48 based on flow stresses
C Fetching yield stress from flow curve
       call ahard(sigmayield0,hard0,sigmayield45,sigmayield90,sigmayieldEB,
     1	  rvalue0,rvalue45,rvalue90,eqplas,PROPS(9),nvalue)	   
       sigmayield=sigmayield0
	   xhard=hard0
	   radius=sigmayield*fT	   
	   s0=sigmayield
       s45=sigmayield45
       s90=sigmayield90
       sEB=sigmayieldEB
       F=((s0**2)/(s90**2)-one+(s0**2/sEB**2))/two
       G=(one-(s0**2/s90**2)+(s0**2/sEB**2))/two
       H=(one+(s0**2/s90**2)-(s0**2/sEB**2))/two
       N=((four*s0**2/s45**2)-(s0**2/sEB**2))/two 
	   L=threeHalves
	   M=threeHalves	   
C Update subvalue of the trial Hill48 equivalent stress	   
       subsigtrial=F*(ST2-ST3)**2+G*(ST3-ST1)**2
     1            +H*(ST1-ST2)**2+two*N*ST4**2
     2            +two*L*ST5**2+two*M*ST6**2
c To prevent subsigtrial =0	 
       if(subsigtrial.lt.Tol) subsigtrial = Tol
C Update the trial Hill48 equivalent stress 	   
       sigtrial = sqrt(subsigtrial)
c To prevent sigtrial =0	   
       if(sigtrial.lt.Tol .and. sigtrial.ge.zero) sigtrial = Tol
       if(sigtrial.gt.(zero-Tol) .and. sigtrial.le.zero) sigtrial = zero-Tol
C Update the radius	   
	   radius = sigmayield*fT
	   Xhardnew= xhard*fT  
       xf = sigtrial - radius
C Check if the stress state is located on the updated yield loci
      if(iflag.gt.(10.+Tol)) then
         goto 1000
      endif
      if(xf.le.Tol .and. xf.ge.(zero-Tol)) then
         goto 1000
      else 
         goto 5000
      endif		   

1000  continue
C     Update PEEQ		 
      stateNew(i,1) = eqplas
c      write(*,*) "stateNew(i,1)", stateNew(i,1)	 
C     Update equivalent stress
	  stateNew(i,2) = sigtrial
      Update Temperature
      stateNew(i,3) = T
      stateNew(i,4) = facT
C     Update stress
      DO 180 K1=1,NTENS          
         stressNew(i,K1)=STRESST(K1)
180   CONTINUE
C   
9000  continue
      return
      end
c
      subroutine ahard(sigmayield0,hard0,sigmayield45,sigmayield90,sigmayieldEB,
     1	  rvalue0,rvalue45,rvalue90,eqplas,table,nvalue,
     2   C1, C2, C3, C4, C5, C6, C7, C8, C9)
C
      include 'vaba_param.inc'
      dimension table(17,nvalue)
      interger :: k1
      real :: eqpl0, eqpl1, deqpl
      real :: sigmayield0a, sigmayield0b, hard45, hard90, hardEB
      real :: sigmayield45a, sigmayield45b
      real :: sigmayield90a, sigmayield90b
      real :: sigmayieldEBa, sigmayieldEBb
      real :: rvalue0a, rvalue0b, hardr0, hardr45, hardr90
      real :: rvalue45a, rvalue45b
      real :: rvalue90a, rvalue90b
      real :: C1a, C1b, hardC1
      real :: C2a, C2b, hardC2
      real :: C3a, C3b, hardC3
      real :: C4a, C4b, hardC4
      real :: C5a, C5b, hardC5
      real :: C6a, C6b, hardC6
      real :: C7a, C7b, hardC7
      real :: C8a, C8b, hardC8
      real :: C9a, C9b, hardC9
C
C     Initialize properties with default values from the last entry
      sigmayield0 = table(2, nvalue)
      sigmayield45 = table(3, nvalue)
      sigmayield90 = table(4, nvalue)
      sigmayieldEB = table(5, nvalue)
      rvalue0 = table(6, nvalue)
      rvalue45 = table(7, nvalue)
      rvalue90 = table(8, nvalue)
      C1 = table(9, nvalue)
      C2 = table(10, nvalue)
      C3 = table(11, nvalue)
      C4 = table(12, nvalue)
      C5 = table(13, nvalue)
      C6 = table(14, nvalue)
      C7 = table(15, nvalue)
      C8 = table(16, nvalue)
      C9 = table(17, nvalue)
C
      hard0 = zero
      hard45 = zero
      hard90 = zero
      hardEB = zero
      hardr0 = zero
      hardr45 = zero
      hardr90 = zero
      hardC1 = zero
      hardC2 = zero
      hardC3 = zero
      hardC4 = zero
      hardC5 = zero
      hardC6 = zero
      hardC7 = zero
      hardC8 = zero
      hardC9 = zero
C
C     If more than one entry, search table
      if (nvalue .gt. 1) then
         do k1 = 1, nvalue - 1
            eqpl0 = table(1, k1)
            eqpl1 = table(1, k1 + 1)
            if (eqplas .ge. eqpl0 .and. eqplas .le. eqpl1) then
               deqpl = eqpl1 - eqpl0
               if (deqpl .le. 1e-12) deqpl = 1e-12
C
C              Interpolate yield stresses
               sigmayield0a = table(2, k1)
               sigmayield0b = table(2, k1 + 1)
               hard0 = (sigmayield0b - sigmayield0a) / deqpl
               sigmayield0 = sigmayield0a + (eqplas - eqpl0) * hard0
C
               sigmayield45a = table(3, k1)
               sigmayield45b = table(3, k1 + 1)
               hard45 = (sigmayield45b - sigmayield45a) / deqpl
               sigmayield45 = sigmayield45a + (eqplas - eqpl0) * hard45
C
               sigmayield90a = table(4, k1)
               sigmayield90b = table(4, k1 + 1)
               hard90 = (sigmayield90b - sigmayield90a) / deqpl
               sigmayield90 = sigmayield90a + (eqplas - eqpl0) * hard90
C
               sigmayieldEBa = table(5, k1)
               sigmayieldEBb = table(5, k1 + 1)
               hardEB = (sigmayieldEBb - sigmayieldEBa) / deqpl
               sigmayieldEB = sigmayieldEBa + (eqplas - eqpl0) * hardEB
C
C              Interpolate r-values
               rvalue0a = table(6, k1)
               rvalue0b = table(6, k1 + 1)
               hardr0 = (rvalue0b - rvalue0a) / deqpl
               rvalue0 = rvalue0a + (eqplas - eqpl0) * hardr0
C
               rvalue45a = table(7, k1)
               rvalue45b = table(7, k1 + 1)
               hardr45 = (rvalue45b - rvalue45a) / deqpl
               rvalue45 = rvalue45a + (eqplas - eqpl0) * hardr45
C
               rvalue90a = table(8, k1)
               rvalue90b = table(8, k1 + 1)
               hardr90 = (rvalue90b - rvalue90a) / deqpl
               rvalue90 = rvalue90a + (eqplas - eqpl0) * hardr90
C
C              Interpolate C1 to C9
               C1a = table(9, k1)
               C1b = table(9, k1 + 1)
               hardC1 = (C1b - C1a) / deqpl
               C1 = C1a + (eqplas - eqpl0) * hardC1
C
               C2a = table(10, k1)
               C2b = table(10, k1 + 1)
               hardC2 = (C2b - C2a) / deqpl
               C2 = C2a + (eqplas - eqpl0) * hardC2
C
               C3a = table(11, k1)
               C3b = table(11, k1 + 1)
               hardC3 = (C3b - C3a) / deqpl
               C3 = C3a + (eqplas - eqpl0) * hardC3
C
               C4a = table(12, k1)
               C4b = table(12, k1 + 1)
               hardC4 = (C4b - C4a) / deqpl
               C4 = C4a + (eqplas - eqpl0) * hardC4
C
               C5a = table(13, k1)
               C5b = table(13, k1 + 1)
               hardC5 = (C5b - C5a) / deqpl
               C5 = C5a + (eqplas - eqpl0) * hardC5
C
               C6a = table(14, k1)
               C6b = table(14, k1 + 1)
               hardC6 = (C6b - C6a) / deqpl
               C6 = C6a + (eqplas - eqpl0) * hardC6
C
               C7a = table(15, k1)
               C7b = table(15, k1 + 1)
               hardC7 = (C7b - C7a) / deqpl
               C7 = C7a + (eqplas - eqpl0) * hardC7
C
               C8a = table(16, k1)
               C8b = table(16, k1 + 1)
               hardC8 = (C8b - C8a) / deqpl
               C8 = C8a + (eqplas - eqpl0) * hardC8
C
               C9a = table(17, k1)
               C9b = table(17, k1 + 1)
               hardC9 = (C9b - C9a) / deqpl
               C9 = C9a + (eqplas - eqpl0) * hardC9
C
               exit
            endif
         end do
C
C        Handle extrapolation beyond the table range
         if (eqplas .ge. table(1, nvalue)) then
            eqpl0 = table(1, nvalue - 1)
            eqpl1 = table(1, nvalue)
            deqpl = eqpl1 - eqpl0
            if (deqpl .le. 1e-12) deqpl = 1e-12
C
C           Extrapolate yield stresses
            sigmayield0a = table(2, nvalue - 1)
            sigmayield0b = table(2, nvalue)
            hard0 = (sigmayield0b - sigmayield0a) / deqpl
            sigmayield0 = sigmayield0b + (eqplas - eqpl1) * hard0
C
            sigmayield45a = table(3, nvalue - 1)
            sigmayield45b = table(3, nvalue)
            hard45 = (sigmayield45b - sigmayield45a) / deqpl
            sigmayield45 = sigmayield45b + (eqplas - eqpl1) * hard45
C
            sigmayield90a = table(4, nvalue - 1)
            sigmayield90b = table(4, nvalue)
            hard90 = (sigmayield90b - sigmayield90a) / deqpl
            sigmayield90 = sigmayield90b + (eqplas - eqpl1) * hard90
C
            sigmayieldEBa = table(5, nvalue - 1)
            sigmayieldEBb = table(5, nvalue)
            hardEB = (sigmayieldEBb - sigmayieldEBa) / deqpl
            sigmayieldEB = sigmayieldEBb + (eqplas - eqpl1) * hardEB
C
C           Extrapolate r-values
            rvalue0a = table(6, nvalue - 1)
            rvalue0b = table(6, nvalue)
            hardr0 = (rvalue0b - rvalue0a) / deqpl
            rvalue0 = rvalue0b + (eqplas - eqpl1) * hardr0
C
            rvalue45a = table(7, nvalue - 1)
            rvalue45b = table(7, nvalue)
            hardr45 = (rvalue45b - rvalue45a) / deqpl
            rvalue45 = rvalue45b + (eqplas - eqpl1) * hardr45
C
            rvalue90a = table(8, nvalue - 1)
            rvalue90b = table(8, nvalue)
            hardr90 = (rvalue90b - rvalue90a) / deqpl
            rvalue90 = rvalue90b + (eqplas - eqpl1) * hardr90
C
C           Extrapolate C1 to C9
            C1a = table(9, nvalue - 1)
            C1b = table(9, nvalue)
            hardC1 = (C1b - C1a) / deqpl
            C1 = C1b + (eqplas - eqpl1) * hardC1
C
            C2a = table(10, nvalue - 1)
            C2b = table(10, nvalue)
            hardC2 = (C2b - C2a) / deqpl
            C2 = C2b + (eqplas - eqpl1) * hardC2
C
            C3a = table(11, nvalue - 1)
            C3b = table(11, nvalue)
            hardC3 = (C3b - C3a) / deqpl
            C3 = C3b + (eqplas - eqpl1) * hardC3
C
            C4a = table(12, nvalue - 1)
            C4b = table(12, nvalue)
            hardC4 = (C4b - C4a) / deqpl
            C4 = C4b + (eqplas - eqpl1) * hardC4
C
            C5a = table(13, nvalue - 1)
            C5b = table(13, nvalue)
            hardC5 = (C5b - C5a) / deqpl
            C5 = C5b + (eqplas - eqpl1) * hardC5
C
            C6a = table(14, nvalue - 1)
            C6b = table(14, nvalue)
            hardC6 = (C6b - C6a) / deqpl
            C6 = C6b + (eqplas - eqpl1) * hardC6
C
            C7a = table(15, nvalue - 1)
            C7b = table(15, nvalue)
            hardC7 = (C7b - C7a) / deqpl
            C7 = C7b + (eqplas - eqpl1) * hardC7
C
            C8a = table(16, nvalue - 1)
            C8b = table(16, nvalue)
            hardC8 = (C8b - C8a) / deqpl
            C8 = C8b + (eqplas - eqpl1) * hardC8
C
            C9a = table(17, nvalue - 1)
            C9b = table(17, nvalue)
            hardC9 = (C9b - C9a) / deqpl
            C9 = C9b + (eqplas - eqpl1) * hardC9	  
        endif
      endif
      return
C
C Iteration ends here
      end

	  