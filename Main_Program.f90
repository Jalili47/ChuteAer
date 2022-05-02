! This program is written by Dr. Mohammadreza Jalili Ghazizadeh  (m_jalili@sbu.ac.ir) under supervision Dr. Amir reza zarrati (zarrati@aut.ac.ir)"
							



	  PROGRAM MAIN
	 !  print *, "This program is written by Dr. Mohammadreza Jalili Ghazizadeh  (m_jalili@sbu.ac.ir) under supervision Dr. Amir reza zarrati (zarrati@aut.ac.ir)"
 	 		 	   	   	
	  CALL INPUT
	  CALL INITIALC5
	  CALL JUDGE
	  STOP
	  
      END


	  
 
	  
	  SUBROUTINE INPUT	  !  Assigning the input values to the model parameters
	  PARAMETER (NY=403,NX=1)
	  DOUBLE PRECISION HB,ZOB,DX,g,tgteta,Qjudge,DZo,Qinitial,Hinitial, &
	& ZINITIAL,ZEND,XLENGTH,UstarJudge,e,Hjudge,RatioResU,Ujudge,Teta,Akapa,AKs,phiaB,phiaF,Vspilt,&
	& eshmit,xplot,dxplot,phicol,AInitialphi0,AInitialphi2,AInitialphi1 &
	& ,rowaterpair	
	  DOUBLE PRECISION Zair,Pair,relaxphi,Vjudge,Pjudge

	  CHARACTER temps*100
	  integer optionkey
	  character file1*20,file2*20,file3*20,filen*30
      COMMON/C12/NCELLI,ZINITIAL,ZEND,Qinitial
	  COMMON/C1267/XLENGTH
	  COMMON/C13/UstarJudge,e,MaxIterU,MaxIterUstar,INDEX,AKs
	  COMMON/C139/Akapa
	  COMMON/C123/tgteta,Hinitial
	  COMMON/C127/DZo
	  COMMON/C1236789/DX(NX)
	  COMMON/C1235/g,teta
	  COMMON/C1246R/Hjudge,Qjudge,Vjudge,Pjudge
      COMMON/C1246I/MaxIter
	  COMMON/C12346/RatioResU,Ujudge 
!	  COMMON/C17/XPlot,DXPlot				  
	  COMMON/C17/XPlot,DXPlot,NSecPlot
	  COMMON/C1234789R/HB,ZOB
	  COMMON/C1234789I/NCELLJ
	  common/covdif/phiaB(NY),phiaF(NY),Vspilt,phicol,eshmit
	  common/covdiB/ idboundbottm,idboundsurf,IdChanson,idVspAuto,relaxphi,IdLinear
	  common/covdiI/AInitialphi0,AInitialphi2,AInitialphi1,zjn1,zjn2,xin1,xin2
	  common/ids/Idconvdif,idfull,idnik,Idprofile
	  common/files/file1,file2,file3
	  common/inputtype2/Idtypecvin,nodata,Zair(NY),Pair(NY)
!	  common/mass/ ResQmass,Qwmass,Qfmass
	  common/optiout/idoptiout,optionkey

	  idoptiout=0
      	
!     *....INITIAL DATA SET
	  open (234,FILE='Entry\Name.txt') 
	  read(234,*) filen
	  open (200,FILE=filen) 
	  if(idoptiout.eq.1) open (10,FILE='ustar.dat')
!    data Feed
!	 *-----Initial&Geo-data'
	 read(200,*)temps
	 read(200,1900)Qinitial
	 read(200,1900)Hinitial
	 read(200,1900)teta
	 read(200,1920)optionkey
	 read(200,1920)NCELLJ
	 read(200,1900)Dx(1)
	 read(200,1900)XLENGTH
	 read(200,1900)Qjudge
	 read(200,1920)NSecPlot
	 read(200,*)temps
	 read(200,*)file1
	 read(200,*)file2
	 read(200,*)file3
	 read(200,*)temps
	 read(200,1924)idprofile
	 read(200,1905)Vjudge
	 read(200,1905)UstarJudge
	 read(200,1905)Ujudge
	 read(200,*)temps
	 read(200,1900)AKs
	 read(200,1900)Akapa
	 read(200,1924)idfull
	 read(200,*)temps
	 read(200,1925)Idconvdif
	 read(200,1900)Vspilt
	! Vspilt=Vspilt*cosd(teta)	 
     read(200,1925)idVspAuto	 
	 read(200,1925)IdChanson
	 read(200,1925)IdLinear
	 read(200,1905)relaxphi
	 read(200,1905)Pjudge
	 read(200,1904)eshmit
     read(200,1925)idboundsurf
	 read(200,1925)idboundbottm
	 read(200,1904)AInitialphi0	    ! BACKGROUND PHI
	 read(200,1904)AInitialphi1  	! WALL AERATION
     read(200,1904)AInitialphi2		! SELF AERATION
	 rowaterpair=800.0d0
	 AInitialphi0=(1.0d0/(1.0d0-AInitialphi0*(1-1.0d0/rowaterpair))-1.0d0)/(rowaterpair-1.0d0)
	 AInitialphi2=(1.0d0/(1.0d0-AInitialphi2*(1-1.0d0/rowaterpair))-1.0d0)/(rowaterpair-1.0d0)
	 AInitialphi1=(1.0d0/(1.0d0-AInitialphi1*(1-1.0d0/rowaterpair))-1.0d0)/(rowaterpair-1.0d0)
	 read(200,1900)zjn1
	 read(200,1900)zjn2
	 read(200,1900)xin1
	 read(200,1900)xin2
	 read(200,1920)nodata
	 if (nodata.ne.0.0) then
	 Idtypecvin=2
	 do ii=1,nodata
	 read(200,1900)Zair(ii)
	 read(200,1900)Pair(ii)
	 Pair(ii)=(1.0d0/(1.0d0-(Pair(ii)/100.0)*(1-1.0d0/rowaterpair))-1.0d0)/(rowaterpair-1.0d0)
	 enddo
     else
	 Idtypecvin=1
	 end if


	 if (optionkey.eq.0) write(*,*)'op=0<usual zero cell>'
	 if (optionkey.eq.1) write(*,*)'op=1<submit-13>'
	 if (optionkey.eq.3) write(*,*)'op=3<external zero cell>'
	 if (optionkey.eq.4) write(*,*)'op=4<external cell>'


!----*methode for lenght calculation
	  idnik=0	   ! idnik=0 Normal formula L=ky

	  
	  HB=Hinitial
 
	  ZINITIAL=0.0
	  ZEND=0.
	  ZOB=ZINITIAL
	  NCELLI=XLENGTH/DX(1)

	  tgteta=-TAND(Teta) 

 
	  DZo=0.0
	  INDEX=0

!---  JUDGE VALUES  	 	   
	  MaxIter=5000
	  MaxIterU=2000
	  MaxIterUstar=10000

!---- FOR SECTION PLOT SELECTION

	  DXPlot=XLENGTH/NSecPlot 
	  XPlot=0.0
!     *--- FLUID PROPERTIES
	  g=0.9806D1

!---- DATA PRINT
    	
				    OPEN (3,FILE=file1)
                    OPEN (4,FILE=file2) 	
				    open(15,FILE=file3)
if(idoptiout.eq.1)	open(17,FILE='turbulnc.dat')
open (5,FILE='zarayebup.txt') 
if(idoptiout.eq.1)	open (5,FILE='zarayebu.txt') 
if(idoptiout.eq.1)	open (6,FILE='zarayebp.txt') 
if(idoptiout.eq.1)	OPEN (8,FILE='PHI.txt') 
if(idoptiout.eq.1)	OPEN (7,FILE='TMP.txt') 

!---- eco input
!	 *-----Initial&Geo-data'
 	  write(3,901)Qinitial,Hinitial,teta
	  write(3,902)NCELLJ,Dx(1),XLENGTH
	  write(3,903)NSecPlot
	  write(3,904)g,Qjudge
	  write(3,912)Vjudge,Pjudge,idprofile
	  write(3,*)
!*-----Turbilence-data'
      write(3,905)AKs,Akapa
      write(3,906)idfull,UstarJudge,idnik
	  write(3,*)
!*-----convectiondiffusion-data'
	  write(3,907)Idconvdif
	  write(3,908)Vspilt,eshmit
	  write(3,909)idboundsurf,idboundbottm
	  write(3,910)AInitialphi0,AInitialphi2,AInitialphi1
	  write(3,911)zjn1,zjn2,xin1,xin2
	  write(3,*)
	  write(3,*)file1
	  write(3,*)file2
	  write(3,*)file3

901	  format('Qw=',D10.5,2x,'H0=',D10.5,2x,"te=",D10.5)
902	  format('NJ=',I4,2x,'DX=',D10.5,2x,"Le=",D10.5)
903	  format('NSplt=',I4)
904	  format('g=',D10.5,2x,'Qj=',D10.5)
905	  format('Ks=',D10.5,2x,'Kp=',D10.5)
906	  format('Idf=',I1,2x,'Ustarjudge=',D10.5,2x,'nikorud=',I1)
907	  format('Icondif=',I1)
908	  format('Vs=',D10.5,2x,'eshmit=',D10.5)
909	  format('Idsurface=',I1,2x,'Idbottom',I1)
910	  format('A0=',D10.5,2x,'A1=',D10.5,2x,'A2',D10.5)
911	  format('z1=',F7.2,1x,'z2=',F7.2,1x,'x1=',F7.2,1x,'x2=',F7.2)
912	  format('Vjudge=',D15.5,1x,'Phijudge=',D15.5,1x,'idprf=',I1)


1900	 format(3x,D20.5)
1904	 format(4x,D20.5)
1905	 format(10x,D20.5)
1902	 format(I2)
1920	 format(3x,I25)
1924	 format(4x,I25)
1925	 format(10x,I25)

	  RETURN
	  END



      SUBROUTINE INITIALC5	   !	It calculates the initial conditions (the first boundary condition) 
	  PARAMETER (NY=403,NX=1)
	  DOUBLE PRECISION HF,HB,ZOF,ZOB,ZF,ZB,DZF,DZB,DX,PF,PB,Y,RoF,RoB,UB,UF,g,tgteta,Qw,DZo,&
    & Qinitial,Hinitial,ZINITIAL,ZEND,XLENGTH,RatioResU,Ujudge,ViscosC,difusivity,DY,WF,WB,Yc,Re,XB,teta,&
	& Qjudge,Hjudge,UOLD1,phiaB,phiaF,Vspilt,phicol,eshmit,AInitialphi0,AInitialphi2, &
	& AInitialphi1,CAF,CAB,Zair,Pair     !,Prpr
	  DOUBLE PRECISION ResQmass,Qwmass,Qfmass,relaxphi,Vjudge,Pjudge,VraisF,VraisB

	  character file1*20,file2*20,file3*20
	  integer optionkey

	  COMMON/C12/NCELLI,ZINITIAL,ZEND,Qinitial
	  COMMON/C1267/XLENGTH
	  COMMON/C2379/UOLD1(NY),ViscosC(NY),difusivity(NY)
      COMMON/C27/XB
	  COMMON/C1235/g,teta
	  COMMON/C1246R/Hjudge,Qjudge,Vjudge,Pjudge
      COMMON/C1246I/MaxIter
	  COMMON/C1236789/DX(NX)
      COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
	  COMMON/C2356789/RoF(NY),PF(NY),WF(NY)  
	  COMMON/C2346811/QW
	  COMMON/C127/DZo
	  COMMON/C123/tgteta,Hinitial
	  COMMON/C12346/RatioResU,Ujudge
	  COMMON/C23789/UB(NY),RoB(NY),PB(NY),ZB(NY),Y(NY),DY(NY),WB(NY)
	  COMMON/C2346789/ UF(NY)
	  COMMON/C1234789R/HB,ZOB
	  COMMON/C1234789I/NCELLJ

	  common/covdif/phiaB(NY),phiaF(NY),Vspilt,phicol,eshmit
	  common/covdiB/ idboundbottm,idboundsurf,idshansen,idVspAuto,relaxphi,IdLinear
	  common/covdiI/AInitialphi0,AInitialphi2,AInitialphi1,zjn1,zjn2,xin1,xin2
	  common/condiC/CAF(NY),CAB(NY)
	  common/ncol/ ncol1
	  common/ids/Idconvdif,idfull,idnik,Idprofile
	  common/files/file1,file2,file3
	  common/inputtype2/Idtypecvin,nodata,Zair(NY),Pair(NY)
	  common/mass/ ResQmass,Qwmass,Qfmass
	  common/nokesid/idnokes
	  common/optiout/idoptiout,optionkey
	  COMMON/raisingvelosity/ VraisF(NY),VraisB(NY)

	  idnokes=0


	  ncol1=0
!     *--- GRID SET
!	  x - location
      NNODEI=NCELLI+1
      NNODEJ=NCELLJ+2

!	  z - location
	  HF=HB
	  DZB=HB/NCELLJ
	  DZF=DZB

!*--- LOCATION NODAL POINTS
      XB=0.0
	  ZOF=ZOB+DZo
	  ZF(1)=ZOF
	  ZF(2)=ZOF+DZF/2
      ZF(NNODEJ)=HF+ZOF
	  ZF(NNODEJ-1)=ZF(NNODEJ)-DZF/2
	  DO J=3,NNODEJ-2
       ZF(J)=ZF(J-1)+DZF
	  END DO
      DO J=1,NNODEJ
      ZB(J)=ZF(J)-DZo
	  END DO
  	  DO J=1,NNODEJ
      Y(J)=ZF(J)/HF
      ENDDO 
	  DY(1)=0.0
	  DY(NNODEJ)=0.0
	  DO J=2,NNODEJ-1
	  DY(J)=1./NCELLJ
	  END DO

! 	  PRINT GRID AND GEOMETRY DATA 
	  !**********---Boundary conditon for Phia*********************************
	  do j=1,NNODEJ
	  PhiaB(J)=AInitialphi0
 	  end do


!	  boundary cond
	  if (Idtypecvin.eq.1) then	  !continuos data input 
	  do j=1,NNODEJ
	  if(ZB(J).ge.zjn1.and.ZB(J).le.zjn2) then
	  phiaB(J)=AInitialphi1
	  end if
	  enddo
	  END IF 
	  
	  if (Idtypecvin.eq.2) then
	  DO J=1,NNODEJ	       
	  do i=1,nodata
	  if(ZB(J).eq.Zair(i)) PhiaB(J)=Pair(i)
	  if(ZB(J).gt.Zair(i).and.ZB(J).lt.Zair(i+1)) PhiaB(J)=Pair(i)+(Pair(i+1)-Pair(i))/(Zair(i+1)-Zair(i))*(ZB(J)-Zair(i))
	  enddo
	  enddo
	  end if
	   ! SELF AERATION Input Data 
	   if (idboundsurf.eq.1) then
	   PhiaB(NNODEJ)=0.0d0 !****************************************
	   end if		!	End of S.A.

   	  DO J=1,NNODEJ	       
	  PhiaF(J)=phiaB(J)
	  END DO


!******************************************************************
	  Phicol=0.0d0
	  DO J=1,NNODEJ	       
	  Phicol=phicol+PhiaF(J)
	  enddo

     CALL ROW_CALCULATION
     DO j=1,NNODEJ
	 ROB(J)=ROF(J)
	 END DO
     CALL C_CALCULATION
     Call PRESSURE
	 DO j=1,NNODEJ
	 CAB(J)=CAF(J)
	 PB(J)=PF(J)
	 END DO
	
	 CALL RISING_VELOCITY
     VraisF(1)=0.0d0
  !  VraisF(NNODEJ)=0.00d0
     
  	 DO J=1,NNODEJ
	 VraisB(J)=VraisF(J)
	 ENDDO


 	  DO J=1,NNODEJ	       
	  rn1p7=7.0
	   if (Idprofile.eq.1) then
	   UB(J)=(Qinitial*((1.0/rn1p7)+1.0))/Hinitial*((ZB(J)-ZOB)/Hinitial)**(1.0/rn1p7)
	   else 
	   UB(J)=Qinitial/Hinitial
	   endif
	   
	  END DO


!  *--- WATER DISCHARGE CAL. 1
	 QW=0.0	 
	 Qwmass=0.0d0
  	 DO J=2,NNODEJ-1
	 QW=QW+UB(J)*DZB*(1.0-CAB(J))/(1.0-PhiaB(J))
     Qwmass=Qwmass+UB(J)*ROB(J)*DZB
	 END DO

!**	 velocity correction to get exact Q
	 DO J=1,NNODEJ
	  UB(J)=UB(J)*Qinitial/Qw  ! end of velocity correction
     
	 END DO 

	 !  *--- WATER DISCHARGE CAL. 2
	 QW=0.0
	 Qwmass=0.0d0
  	 DO J=2,NNODEJ-1
	 QW=QW+UB(J)*DZB*(1.0-CAB(J))
!	 QW=QW+UB(J)*DZB*(1.0/(1+(PhiaB(J)*(800.0-1.0))))
	 Qwmass=Qwmass+UB(J)*ROB(J)*DZB
	 END DO
!	   next section
       DO J=1,NNODEJ
       UF(J)=UB(J)
       WB(J)=0.0
       WF(J)=WB(J)
	   ViscosC(J)=1D-3
      
	  END DO
	! Nokes example
	  if(idnokes.eq.1) then 
!	  call Umoment2
	  UF(NNODEJ)=UF(NNODEJ-1)
      UB(NNODEJ)=UB(NNODEJ-1)
	  UB(1)=0.0
	  UF(1)=0.0
	  DO J=1,NNODEJ
	  IF (UF(J)<0.0) THEN
	   WRITE (6,1076)
	   PAUSE 'WARNING 387 '
	  END IF
	  END DO
	  END IF


!** CALCUL. Re & Froud No.	for water only 
	 Yc=(Qw*Qw/g)**(1./3.)
	 Re=Qw/1e-6
	 
     WRITE(3,1000) NNODEJ,NCELLJ,DX(1),XLENGTH,ZINITIAL,ZOB,ZOF,HF,teta,tgteta,Ujudge,Qjudge,Qw
     WRITE(4,1000) NNODEJ,NCELLJ,DX(1),XLENGTH,ZINITIAL,ZOB,ZOF,HF,teta,tgteta,Ujudge,Qjudge,Qw
	! WRITE(3,1140) XB,HF,HF-HB,Yc,Re,Qw

	! WRITE(3,2019) 
	 WRITE(3,2021)
	  







1000 FORMAT(' NNODEJ='I10/,' NCELLJ='I10/,' DX(1)='F14.8/,' XLENGTH='F12.3/,' ZINITIAL='F10.3/,' ZOB='F10.3/,' ZOF=',F10.7/,' HF='F10.7/,' teta=',f10.6/,' tgteta=',f14.12/,' Ujudge',F15.9/,' Qjudge',F15.12/,' Qw=',F10.6//)
1080 FORMAT ('XB=',F8.5,'HB=',F8.5,5X,'Qw=',F8.6)	 
1090 FORMAT (' U(' ,I3,')=',F10.7,10X,' W(' ,I3,')=',F10.7)
!1140 FORMAT ('NSECTION',8X,'XB(m)',8X,'HF(m)',8X,'HF-HB',8X,'Yc',13x,'Renolds',6x,'Qw'//,5X,' 1',2X,4F15.9,F14.2,F11.7)   
1076 FORMAT (//' THE VALUE OF VELOCITY IS NEGATIVE ...')

!2019 FORMAT ('**************************************************************************************************'   /)
2021 FORMAT ('   NSECTION           X(m)           Z(m)          '   /)

	 RETURN
     END
!**************************************************************************************************************************



!*---maning
      SUBROUTINE UMOMENTUM	 ! 	  Solving the u momentum equation.
      PARAMETER(NY=403,NX=1)
	  DOUBLE PRECISION HF,HB,ZOF,ZOB,ZF,DZF,DZB,DX,PF,PB,ZB,Dh,Y,RoF,SF,RoB,UB,UF,g,tgteta,Qw,&
   &  Hinitial,UstarJudge,e,RatioResU,Ujudge,ViscosC,difusivity,DY,WF,WB,USTAR,XF,ResUtotal,Akapa,C,Crough,&
   &  AKs,USTAROLD,UOLD1,nman,cfman,Ustarman
	  DOUBLE PRECISION Fsw,Fnw,FP,FW,Dsw,Dnw,Fswt,Fnwt,Fslope,AS,AP,AN,S
	  DOUBLE PRECISION CC,UU,DD,XX,BB,ResU,teta
	  DOUBLE PRECISION ResQmass,Qwmass,Qfmass

	  integer optionkey
	  COMMON/C2379/UOLD1(NY),ViscosC(NY),difusivity(NY)
	  COMMON/C2356789/RoF(NY),PF(NY),WF(NY)
	  COMMON/C13/UstarJudge,e,MaxIterU,MaxIterUstar,INDEX,AKs
      COMMON/C139/Akapa
	  COMMON/C123/tgteta,Hinitial
	  COMMON/C1235/g,teta
	  COMMON/C2346811/QW
	  COMMON/C23789/UB(NY),RoB(NY),PB(NY),ZB(NY),Y(NY),DY(NY),WB(NY)
	  COMMON/C1236789/DX(NX)
	  COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
  	  COMMON/C2346789/ UF(NY)
      COMMON/C1234789R/HB,ZOB
	  COMMON/C1234789I/NCELLJ
      COMMON/C12346/RatioResU,Ujudge 
	  COMMON/C3467/ITER
	  COMMON/C37/SF,USTAR
      COMMON/C378/Dh(NY)
      COMMON/C3678/XF,NSECTION
	  common/mass/ ResQmass,Qwmass,Qfmass
	  common/optiout/idoptiout,optionkey
	  DIMENSION Fsw(NY),Fnw(NY),FP(NY),FW(NY),Dsw(NY),Dnw(NY),Fswt(NY),Fnwt(NY),AS(NY),AP(NY),AN(NY),S(NY),ResU(NY),Fslope(NY)
	  DIMENSION  CC(NY-2),UU(NY-3),DD(NY-3),XX(NY-2),BB(NY-2)



!	  LOOP IN HEIGHT
!      Calculation of C 
	    Crough=30.*(DZF+DZB)/2./2./Aks
		C=Crough
!	   END OF Calculation C
	   DO J=1,NNODEJ

	   Dh(j)=ZF(j)-ZB(j)	  !!

	   END DO
	  
	   DO J=2,NNODEJ-1
	   Fsw(J)=0.25*(RoF(J-1)+RoF(J)+RoB(J-1)+RoB(J))*WF(J-1)*DX(1)
	   Fnw(J)=0.25*(RoF(J)+RoF(J+1)+RoB(J)+RoB(J+1))*WF(J)*DX(1)
	   Dsw(J)=0.5*(ViscosC(J-1)+ViscosC(J))*DX(1)/(0.5*(HF+HB))/(0.5*(DY(J)+DY(J-1))) 
	   Dnw(J)=0.5*(ViscosC(J+1)+ViscosC(J))*DX(1)/(0.5*(HF+HB))/(0.5*(DY(J)+DY(J+1)))
	   FW(J)=RoB(J)*UB(J)*HB*DY(J)
	   END DO  
	   IterU=0
       II=0
	   id1=0
	   UstarOld=0D0

	   Relax=1.0d0	! Under relaxation for Ustar

!----------maning calculation
	   cfman=(6.0d0+5.75*log10(0.5*(HF+HB)/AKs))**(-2)
	   nman=(cfman*((0.5*(HF+HB))**(1.0/3.0))/g)**0.5
	   Ustarman=((g**0.5)*nman*Qw)/((0.5*(HF+HB))**(7.0/6.0))

	   DO   	   
	   II=II+1
! 	   Ustar calculation
	   Ustar=(UF(2)+UB(2))/2.*Akapa/LOG(C)
!********************************************************************88	   

5252   format('X=',F10.3,2X,'It=',I3,2X,'Usn=',D10.3,2x,'Uso=',D10.3,2x,'Uf=',D10.3)
       IF (Relax.le.1.0d0.and.II.ge.2.) Ustar=UstarOld+Relax*(Ustar-UstarOld)
!???????????????????????????????????????????????????????????????????///
	          
!      ROUGH OR SMOOTH?
	      SR=Ustar*Aks/1E-6
	      IF (SR < 3.5) THEN
	      PRINT *, "BED IS NOT ROUGH" ,"       SR=",SR

	      END IF


	     II=II+1

!***********************************************************************
	   shakhes=ABS(Ustar-UstarOld)/Ustarman




	   IF (shakhes < UstarJudge)  EXIT
	   	   	  
	   Id1=1
	   UstarOld=Ustar
	   DO J=2,NNODEJ-1
10	    FP(J)=RoF(J)*UF(J)*HF*DY(J)
	    Fswt(J)=0.25*(RoF(J-1)+RoF(J)+RoB(J-1)+RoB(J))*0.25*(UF(J-1)+UF(J)+UB(J-1)+UB(J))*0.5*(Dh(J)+Dh(J-1))
        Fnwt(J)=0.25*(RoF(J+1)+RoF(J)+RoB(J+1)+RoB(J))*0.25*(UF(J+1)+UF(J)+UB(J+1)+UB(J))*0.5*(Dh(J)+Dh(J+1))
	    Fslope(J)=0.5*(HB+HF)*0.5*(RoF(J)+RoB(J))*g*SIND(Teta)*DX(1)*DY(J)
	   END DO
  
	  fsw(2)=0.0d0
	  Fswt(2)=0.0d0	  

      Dsw(2)=ViscosC(1)*DX(1)/(0.5*(HF+HB))/(0.5*(DY(2)+DY(1)))	  
!*--- CORRECTION COEFFICIENT FOR UPER MOST CELL (1)
       Fnw(NNODEJ-1)=0.5*(RoF(NNODEJ)+RoB(NNODEJ))*WF(NNODEJ-1)*DX(1)
       Fnwt(NNODEJ-1)=0.5*(RoF(NNODEJ)+RoB(NNODEJ))*0.5*(UF(NNODEJ-1)+UB(NNODEJ))*Dh(NNODEJ)
       Dnw(NNODEJ-1)=ViscosC(NNODEJ)*DX(1)/(0.5*(HF+HB))/(0.5*(DY(NNODEJ-1)+DY(NNODEJ)))

!     *--- MAIN COEFFICIENT
       DO J=2,NNODEJ-1
	    AS(J)=Fswt(J)/4-Dsw(J)/2-Fsw(J)/4
        AN(J)=-Fnwt(J)/4-Dnw(J)/2+Fnw(J)/4
        AP(J)=FP(J)+Fswt(J)/4-Fnwt(J)/4+Dsw(J)/2+Dnw(J)/2+Fnw(J)/4-Fsw(J)/4
    
	   S(J)=(Fw(J)+Fnwt(J)/4-Fswt(J)/4+Fsw(J)/4-Fnw(J)/4-Dnw(J)/2-Dsw(J)/2)*UB(J)&
      &  +(Fnwt(J)/4+Dnw(J)/2-Fnw(J)/4)*UB(J+1)&
	  &	 +(-Fswt(J)/4+Dsw(J)/2+Fsw(J)/4)*UB(J-1)&
	  &	 +(HB*PB(J)-HF*PF(J))*DY(J)+0.5*(Dh(J)+Dh(J+1))*0.25*(PF(J+1)+PB(J+1)+PF(J)+PB(J))-0.5*(Dh(J)+Dh(J-1))*0.25*(PF(J-1)+PB(J-1)+PF(J)+PB(J))+Fslope(J)
	    END DO

!*--- CORRECTION COEFFICIENT FOR BOTTOM CELL (2)
       Dsw2=(RoB(1)+RoF(1))/2.*ABS(Ustar)*DX(1)*Akapa/2./LOG(C) 
	   AS(2)=0.0
       AP(2)=FP(2)-Fnwt(2)/4+Dsw2+Dnw(2)/2+Fnw(2)/4
  
       S(2)=(Fw(2)+Fnwt(2)/4-Fnw(2)/4-Dnw(2)/2-Dsw2)*UB(2)+&
	 & (Fnwt(2)/4+Dnw(2)/2-Fnw(2)/4)*UB(3)+(HB*PB(2)-HF*PF(2))*DY(2)+(0.5*(Dh(2)+Dh(3))/4)*(PF(3)+PB(3)+PF(2)+PB(2))-(Dh(1)/2)*(PF(1)+PB(1))+Fslope(2)
	   	 
! *--- CORRECTION COEFFICIENT FOR UPER MOST CELL (2)
	   AN(NNODEJ-1)=0.0
	   AP(NNODEJ-1)=FP(NNODEJ-1)+Fswt(NNODEJ-1)/4-Fnwt(NNODEJ-1)/2+Dsw(NNODEJ-1)/2+Fnw(NNODEJ-1)/2-Fsw(NNODEJ-1)/4

	   S(NNODEJ-1)=(Fw(NNODEJ-1)-Fswt(NNODEJ-1)/4+Fsw(NNODEJ-1)/4-Dsw(NNODEJ-1)/2)*UB(NNODEJ-1) &
	&   +(Fnwt(NNODEJ-1)/2.-Fnw(NNODEJ-1)/2)*UB(NNODEJ)&
	&   +(-Fswt(NNODEJ-1)/4+Dsw(NNODEJ-1)/2+Fsw(NNODEJ-1)/4)*UB(NNODEJ-2)  &
	&   +(HB*PB(NNODEJ-1)-HF*PF(NNODEJ-1))*DY(NNODEJ-1)+(Dh(NNODEJ)/2)*(PB(NNODEJ)+PF(NNODEJ))-(0.5*(Dh(NNODEJ-1)+Dh(NNODEJ-2))/4)*(PF(NNODEJ-2)+PB(NNODEJ-2)+PF(NNODEJ-1)+PB(NNODEJ-1))+Fslope(NNODEJ-1)

!**   CALCULATION OF ResU
	   IdResU=0
	   IF (IdResU.eq.1) THEN
	   DO J=2 ,NNODEJ-1
	    ResU(J)=AS(J)*UF(J-1)+AP(J)*UF(J)+AN(J)*UF(J+1)-S(J)
	   END DO
	   ResUtotal=0.0
	   DO J=2,NNODEJ-1
	    ResUtotal=ResUtotal+ABS(ResU(J))
	   END DO

      RatioResU=ResUtotal/1000./(Qw*Qw/Hinitial/Hinitial*Hinitial) 	 ! end of ResU

       1109 FORMAT (I7,5X,'RatioResU=',F12.9,3X,'HF=',F12.9,3X,'HB=',F12.9) 

	  IdResU=0
	  RETURN
	  END IF
	
	 
	   IF (ITER > 50000) THEN
	   if(idoptiout.eq.1) WRITE(6,1107) ITER,RatioResU
1107   FORMAT('ITR=',I6,5X,'RatioResU=',F12.10)

	   END IF

  !  ---PRINT DATA
  !	   GO FOR TDMA
	   DO J = 1,NCELLJ
	    CC(J)=AP(J+1)
	    BB(J)=S(J+1)
	    XX(J)=UF(J+1)
	   END DO
	   DO J= 1,NCELLJ-1
	    UU(J)=AN(J+1)
	    DD(J)=AS(J+2)
	   END DO


	   CALL TDMA (NCELLJ,CC,UU,DD,XX,BB) 
       DO I=1,NCELLJ
	    UF(I+1)=XX(I)
	   END DO
	   UF(NNODEJ)=UF(NNODEJ-1)
	   UF(1)=0.0
 	  ENDDO
 1008  FORMAT('HF=',F8.6,2X,'HB=',F8.6,' ZF(',I3,')=',F10.8,3X,'ZB(',I3,')=',F10.8,3X,'DZ(',I3,')=',F10.8,3X,' Dh(',I3,')=',F10.8)
    IF (INDEX > 0.) THEN


     WRITE(1,1070) (J,UF(J),J,WF(J),J,ViscosC(J),J=1,NNODEJ)
  END IF
1010  FORMAT(' ZF(',I3,')=',F10.8,10X,'Y(',I3,')=',F10.7,10X,' DY(',I3,')=',F10.8)
1020  FORMAT(' PB(',I3,')=',F10.4,10X,' PF(',I3,')=',F10.4,10X,'Dh(',I3,')=',F16.10)
1030  FORMAT ('FW(',I2,')= ',F13.7,1X,'FP(',I2,')=' ,F13.7,1X,'Fnw(' ,I2,')=',F15.7,1X,' Fsw(' ,I2,')=',F15.7) 
1040  FORMAT ('Dnw(',I2,')=',F11.8,1X,' Dsw(',I2,')=',F11.8,1X,' Fnwt(',I2,')=',F15.7,1X,' Fswt(',I2,')=',F15.7) 
1050  FORMAT (' AN(' ,I3,')=',F12.7,1X,' AP(' ,I3,')=',F12.7,1X,' AS('  ,I3,')=',F12.7,1X,' S('   ,I3,')=',F15.7) 
1060  FORMAT (/,F10.5,4X,' Ustar=' ,F10.5,4X,' Dsw2=',F10.5) 
1070  FORMAT (' UF(' ,I3,')=',F14.10,10X,' WF(' ,I3,')=',F10.6,' ViscousC(' ,I3,')=',F14.10) 
1075  FORMAT (//' THE VALUE OF Sf IS NEGATIVE ...')
1152  FORMAT ('NSEC=',I7,2X,'XF=',F18.9,2X,'HF=',F16.12,5X,'HB=',F16.8)
1077  FORMAT (/'XF=',F12.6,2X,'HF=',F8.6,2X,'HB=',F8.6,2X,'SF=',F12.6/)	 
 	  RETURN
	  END


      SUBROUTINE TDMA (NN,CC,UU,DD,XX,BB)    ! Solving the systems of equations using the Tridiagonal Matrix Algorithm (TDMA) .
	  DIMENSION  CC(NN),UU(NN-1),DD(NN-1),XX(NN),BB(NN)
      DOUBLE PRECISION CC,UU,DD,XX,BB,CO
   !   DOUBLE PRECISION NN,NCELLJ
	  DO I=2,NN
	  CO=DD(I-1)/CC(I-1)
	  CC(I)=CC(I)-CO*UU(I-1)
      BB(I)=BB(I)-CO*BB(I-1)
      END DO

!     BACKSUBSTITUDE
      XX(NN)=BB(NN)/CC(NN)
	   DO I=NN-1,1,-1
        XX(I)=(BB(I)-UU(I)*XX(I+1))/CC(I)
	   END DO
      RETURN
      END


	  SUBROUTINE DEPTH    ! Calculating flow depth and z coordinate of computational cells
 	  PARAMETER (NY=403,NX=1)
      DOUBLE PRECISION HF,HB,ZOF,ZOB,ZF,DZF,DZB,UF,HFOLD,Qw,Qf,Qjudge,Hjudge,ResH,ResQ,RatioResU,&
	& Ujudge,CAF,CAB,resdualv,RoF,PF,WF,phiaB,phiaF,Vspilt,phicol,eshmit
	  DOUBLE PRECISION ResQmass,Qwmass,Qfmass,Vjudge,Pjudge
!	  DOUBLE PRECISION NCELLI,MaxIterU,MaxIterUstar,INDEX,MaxIter,NCELLJ,NNODEJ
	  integer optionkey
	  COMMON/C2346811/QW
	  COMMON/C478/QF
	  COMMON/C2346789/UF(NY)
      COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
  	  COMMON/C46/ResH,ResQ,IdResQ
      COMMON/C1246R/Hjudge,Qjudge,Vjudge,Pjudge
      COMMON/C1246I/MaxIter
	  COMMON/C1234789R/HB,ZOB
      COMMON/C1234789I/NCELLJ
      COMMON/C3467/ITER
	  COMMON/C12346/RatioResU,Ujudge	 
	  common/covdif/phiaB(NY),phiaF(NY),Vspilt,phicol,eshmit
	  common/condiC/CAF(NY),CAB(NY)
	  common/resid/resdualv
	  COMMON/C2356789/RoF(NY),PF(NY),WF(NY)
	  common/ids/Idconvdif,idfull,idnik,Idprofile
	  common/mass/ ResQmass,Qwmass,Qfmass
	  common/optiout/idoptiout,optionkey


	  HFOLD=HF
!  *--- WATER DISCHARGE CAL.
	 QF=0.0
	 Qfmass=0.0d0
	 DO J=2,NNODEJ-1
	 QF=QF+UF(J)*DZF*ROF(J)/1000.0
	 diff=ROF(J)/1000.0-(1.0-CAF(J))/(1.0-PhiaF(J))
	 if (diff.gt.000001) then
	 write(*,*) diff
	 PAUSE " LINE 662"
	 end if

	 Qfmass=Qfmass+UF(J)*DZF*ROF(J)/(1.0-PhiaF(J))


	 END DO

     IF (IdResQ.eq.1) THEN
	 ResQ=ABS(Qw-Qf)/Qw
	 ResQmass=ABS(Qwmass-Qfmass)/Qwmass
	 IdResQ=0
	 return
	 END IF
		   

  	 HF=Qwmass/Qfmass*HFOLD

!--- CALCULATION OF Z NEW VALUES 
	  DZF=HF/NCELLJ
	  ZF(1)=ZOF
	  ZF(2)=ZOF+DZF/2.
      ZF(NNODEJ)=HF+ZOF
	  ZF(NNODEJ-1)=ZF(NNODEJ)-DZF/2
	  DO J=3,NNODEJ-2
       ZF(J)=ZF(J-1)+DZF
	  END DO


  	 
	 RETURN
	 END



	 SUBROUTINE PRESSURE	 ! Calculating nodal pressure 
   	 DOUBLE PRECISION HF,ZoF,ZF,PF,DZF,DZB,RoF,g,WF,teta
	 DOUBLE PRECISION HB,ZOB,UB,RoB,PB,ZB,Y,DY,WB
  
	 PARAMETER (NY=403,NX=1)
	 integer optionkey
	 COMMON/C2356789/RoF(NY),PF(NY),WF(NY)
	 COMMON/C1234789R/HB,ZOB
	 COMMON/C1234789I/NCELLJ
	 COMMON/C1235/g,Teta
	 COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
	 COMMON/C23789/UB(NY),RoB(NY),PB(NY),ZB(NY),Y(NY),DY(NY),WB(NY)
	 common/optiout/idoptiout,optionkey
     
 ! *---PRESSURE CALCULATION
	 PF(NNODEJ)=0.0d0
	 DO J=NNODEJ-1,1,-1
	 
	 PF(J)=PF(J+1)+(ZF(J+1)-ZF(J))*0.5*(ROF(J+1)+ROF(J))*cosd(Teta)*g
	 END DO
	 RETURN
	 END




	 SUBROUTINE JUDGE	! Comparing the convergence criteria
	 PARAMETER (NY=403,NX=1)
   	 DOUBLE PRECISION DX,Qjudge,XLENGTH,Hjudge,RatioResU,Ujudge,XF,ResH,ResQ,QW,HF,ZOF,ZF,DZF,DZB,&
   & UF,PF,ROF,WF,resphi,resdualv
     DOUBLE PRECISION ResQmass,Qwmass,Qfmass,Vjudge,Pjudge
  
	 integer optionkey
	 COMMON/C1246R/Hjudge,Qjudge,Vjudge,Pjudge
     COMMON/C1246I/MaxIter
	 COMMON/C12346/RatioResU,Ujudge 
	 COMMON/C2346811/QW
	 COMMON/C46/ResH,ResQ,IdResQ
	 COMMON/C3678/XF,NSECTION
	 COMMON/C1267/XLENGTH
	 COMMON/C1236789/DX(NX)
	 COMMON/C3467/ITER
	 COMMON/C2346789/UF(NY)
	 COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
	 COMMON/C2356789/RoF(NY),PF(NY),WF(NY)
	 common/ids/Idconvdif,idfull,idnik,Idprofile
	 common/resid/resdualv
	 common/nokesid/idnokes
	 common/visnocalc/idvis
	 common/mass/ ResQmass,Qwmass,Qfmass
	 common/optiout/idoptiout,optionkey
	 Common/C68/IdResV
	 Common/C611/IdResPhi,resphi

	 Idvis=1
	 NSECTION=1 
	 XF=DX(1)
	 numzone=0
	 nploted=0
	 call TECPLOT(numzone)
	 
	 XF=DX(1)
	 DO
	
	 WRITE(*,1104) XF/XLENGTH*100,Iter
	 ! WRITE(*,*) XF/XLENGTH*100,Iter,idoptiout

	  ResQ=Qjudge
	  ResH=Hjudge
	  RatioResU=Ujudge
	   	 Iter=0.0
if (idnokes.eq.0) then
	 DO
       	   IF (Iter > 1000)	then 
		   WRITE(*,1104) XF/XLENGTH*100,Iter,Resphi,resdualv,resq
		   end if 
		   
		   IF(Iter > MaxIter) THEN 
		   WRITE(*,1104) XF/XLENGTH*100,Iter,Resphi,resdualv,resq
			
		    CALL WARNING (MaxIter,6)  
	       ENDIF 

		    Iter=Iter+1

                    CALL MIXINGLENGTH                    
             		CALL UMOMENTUM
			   	    !   GO for calculation of ResQ	
					       IdResQ=1
					       call DEPTH	! the goal is finding ResQ no depth calculation
				    !End of calculation ResQ

if(Idconvdif.eq.1)	THEN
                  !   GO for calculation of ResPhi	
					       IdResPhi=1
					       call CONVECTION_DIFFUSION	! the goal is finding ResPhi not  calculation of Phi
				    !End of calculation ResPhi
				   ! call ResdPhi
					END IF 
					
					!   GO for calculation of ResV	
					       IdResV=1
					       call CONTINUITY	! the goal is finding ResV calculation of V
				    !End of calculation ResQ

 if (iter.eq.1399) then
alaki=0.


end if				
  IF (ResQ < Qjudge.and.resphi < Pjudge.and.resdualv.lt.Vjudge) EXIT

	                CALL DEPTH
					CALL PRESSURE					
                 	CALL CONTINUITY
if(Idconvdif.eq.1)  CALL CONVECTION_DIFFUSION 
 

 	   IF (ITER > 50000) THEN

		  PRINT *,ResQ,Qjudge,RatioResU,Ujudge
		  END IF
		    
	 END DO
                 ELSE


     end if
	  IF(XLENGTH < XF-DX(1)) EXIT 
		  NSECTION=NSECTION+1
		  CALL SUBSTITUATION 
	 
	  END DO 

	  

	 RETURN
1100 FORMAT ('SECTION NUMBER IS:',I4)
1110 FORMAT ('ITERATION NUMBER IS:',I4)
1104 FORMAT (F5.0,5x,I6,5X,F12.9,2X,F12.9,2X,F12.9)
1106 FORMAT (' ResQ MUST BE LESS THAN Qjudge, THERE IS A PROBLEM, CHECK PROGRAM',3X,'ResQ=',F12.9,3X,'Qjudge=',F15.11) 
	 END



	 SUBROUTINE SUBSTITUATION  ! Substitute "k+1" th section with  the "k" th section based on the marching technique.
     PARAMETER (NY=403,NX=1)
  	 DOUBLE PRECISION HF,HB,ZOF,ZOB,ZF,ZB,DZF,DZB,DX,PF,PB,Dh,RoF,SF,Y,RoB,UB,UF,Qf,DZo,XLENGTH,DY,&
   & WF,WB,USTAR,XF,XB,DFB,UOLD1,ViscosC,difusivity,phiaB,phiaF,Vspilt,phicol,eshmit,XPlot,DXPlot,CAF,CAB,eshmitn,DRoDY,rich
	 DOUBLE PRECISION ResQmass,Qwmass,Qfmass
	 DOUBLE PRECISION DUDY,DUDX,aLength,DVDY,DVDX,gent
	 DOUBLE PRECISION delta(NY),VraisF,VraisB
  
	 integer optionkey
	 COMMON/C2379/UOLD1(NY),ViscosC(NY),difusivity(NY)
  	 COMMON/C3678/XF,NSECTION
   	 COMMON/C2346789/UF(NY)
     COMMON/C478/QF
     COMMON/C23789/UB(NY),RoB(NY),PB(NY),ZB(NY),Y(NY),DY(NY),WB(NY)
	 COMMON/C2356789/RoF(NY),PF(NY),WF(NY)
	 COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
	 COMMON/C1236789/DX(NX)
  	 COMMON/C1234789R/HB,ZOB
	 COMMON/C1234789I/NCELLJ

	 COMMON/C1267/XLENGTH
     COMMON/C127/DZo

	 COMMON/C17/XPlot,DXPlot,NSecPlot
     COMMON/C27/XB
	 COMMON/C37/SF,USTAR
     COMMON/C378/Dh(NY)
 	 COMMON/C3467/ITER
	 common/covdif/phiaB(NY),phiaF(NY),Vspilt,phicol,eshmit
	 common/condiC/CAF(NY),CAB(NY)
	 common/turtec/DUDY(NY),DUDX(NY),aLength(NY),DVDY(NY),DVDX(NY),gent
	 common/mass/ ResQmass,Qwmass,Qfmass
	 common/optiout/idoptiout,optionkey
	 COMMON/raisingvelosity/ VraisF(NY),VraisB(NY)
	 common/eshmitno/eshmitn(NY),DRoDY(NY),rich(NY)

	        do jj2=1,NNODEJ
	if (UF(NNODEJ-1).ne.0.0) delta(jj2)=(UF(jj2))/(UF(NNODEJ-1))
	enddo

	 IF (XB >= XPLOT.or. XF >= XLENGTH) THEN

	  WRITE(4,1160)
	    IF (XF >= XLENGTH) THEN
		 DO J=1,NNODEJ
		 ! WRITE(4,1100)J,XF,ZF(J)-ZOF,UF(J),WF(J),ViscosC(J),(ViscosC(J)/1000/HF/Ustar),phiaF(j),CAF(j)*100,VraisF(J)
		  WRITE(4,1100)J,XF,ZF(J)-ZOF,UF(J),WF(J),CAF(j)*100

		  write(15,*)XF,ZF(J)-ZOF,UF(J),WF(J),PF(J),CAF(J)*100,ViscosC(J),ROF(J),delta(J),VraisF(J),eshmitn(J),DRoDY(J),rich(J)
		  if(idoptiout.eq.1) write(17,*)XF,ZF(J)-ZOF,UF(J),WF(J),aLength(J),ViscosC(J),DUDY(J),DUDX(J),DVDY(J),DVDX(J),gent		  
	     END DO
		ELSE
         DO J=1,NNODEJ
	    !  WRITE(4,1100)J,XB,ZB(J)-ZOB,UB(J),WF(J),ViscosC(J),(ViscosC(J)/1000/HF/Ustar),phiaB(j),CAB(j)*100,VraisB(J)
		WRITE(4,1100)J,XB,ZB(J)-ZOB,UB(J),WF(J),CAB(j)*100
  
		  write(15,*)XB,ZB(J)-ZOB,UB(J),WB(J),PB(J),CAB(J)*100,ViscosC(J),ROB(J),delta(J),VraisB(J),eshmitn(J),DRoDY(J),rich(J)
		  if(idoptiout.eq.1) write(17,*)XB,ZB(J)-ZOB,UB(J),WB(J),aLength(J),ViscosC(J),DUDY(J),DUDX(J),DVDY(J),DVDX(J),gent		  		 
		 END DO
   		 END IF
	    WRITE(4,*)
	 XPLOT=XPLOT+DXPLOT
	  ! WRITE(3,1150) NSECTION,XF,HF,HF-HB,ITER

	   WRITE(3,1150) NSECTION,XF,HF




	    
 END IF

!**  CALCULATION OF CRITERIA VALUE 
	 CRITICAL=ABS(Y(NNODEJ-1)/HF*(HF-HB))/DY(NNODEJ-1) !END OF CALC.
	 IF (CRITICAL > 1.0) THEN
	 	  if(idoptiout.eq.1) WRITE(6,1175)
	  PAUSE 'WARNING: The sigma-coordinate criteria has not been satisfied LINE 905 '
	 END IF

 !--- CHANGE TO NEW SECTION 	 
	 DFB=HF-HB	   ! A Linear progress assumed.
!	 DFB=0.0	   ! assumed HF=HB in the first guess
	 ZoB=ZoF
	 HB=HF
	 DzB=DzF
	 DO J=1,NNODEJ
	 ZB(J)=ZF(J)
	 UB(J)=UF(J)
	 PB(J)=PF(J)
	 WB(J)=WF(J)
!	 WF(J)=0.0
	 phiaB(J)=phiaF(J)
	 RoB(J)=RoF(J)	
	 CAB(J)=CAF(J)
	 VraisB(J)=VraisF(J)
	 END DO
	 XB=XF
!*--- CALCULATION OF NEW VALUES
	
	 HF=HF+DFB		   
	 ZoF=ZoF+DZo
     XF=XF+DX(1)
    
	  
 !*-- CALCULATION ZF
	  DZF=HF/NCELLJ
	  ZF(1)=ZOF
	  ZF(2)=ZOF+DZF/2.
      ZF(NNODEJ)=HF+ZOF
	  ZF(NNODEJ-1)=ZF(NNODEJ)-DZF/2
	  DO J=3,NNODEJ-2
      ZF(J)=ZF(J-1)+DZF
	  END DO

1080 FORMAT (' XF=',F8.5,' HF=',F8.5,5X,'Qf=',F8.6)	 
!   1100 FORMAT (I4,1X,F15.7,2X,F10.8,5X,F11.6,5x,F11.6,5X,F11.6,5X,F11.6,5X,F11.6,5X,F11.6,5X,F11.6)
1100 FORMAT (I4,1X,F15.7,2X,F10.8,5X,F11.6,5x,F11.6,8X,F11.6)
1150 FORMAT (I7,2X,2F18.9,2X,F16.12,5X,I4)
!     1160 FORMAT ('   J       X(m)        Z(m)             U(m/s)           W(m/s)          VISCOS          1/Re(U*)'   /) 
1160 FORMAT ('   J       X(m)        Z(m)             U(m/s)           W(m/s)          Air_Concentration(%)           '   /) 
1175 FORMAT (//' THE CRETERIA VALUE IS LESS THAN 1 AND THIS IS ON CONTRARY WITH CONVERGENCY CRETERIA, CHECK THE VALUES OF PARAMETERS'//)



	 RETURN 
	 END


	 SUBROUTINE WARNING (MaxIter,NSub) 	! Sending divergence warnings
	 integer optionkey
	 common/optiout/idoptiout,optionkey
     if(idoptiout.eq.1) WRITE(6,1112) MaxIter,Nsub 
	  
     PAUSE 'WARNING 1012'
	 STOP
 1112 FORMAT (//' NUMBER OF ITERATION EXCCED OF',I5,'  IN SUBROUTINE  ',I3,'  PROGRAM DOES NOT CONVERGE')

	 END




      SUBROUTINE CONTINUITY	 ! Solving the continuity equation to calculate the v parameter
      PARAMETER(NY=403 ,NX=1)
      DOUBLE PRECISION HF,HB,ZOF,ZOB,ZF,DZF,DZB,DX,PF,PB,ZB,Dh,Y,RoF,RoB,UB,UF,DY,WF,WB,XF,QF,Qw
   	  DOUBLE PRECISION FCswt,FCnwt,ASC,ANC,SC,resdualv
	  DOUBLE PRECISION ResQmass,Qwmass,Qfmass,XPlot,DXPlot
   
	  integer optionkey
	  COMMON/C2356789/RoF(NY),PF(NY),WF(NY)
	  COMMON/C23789/UB(NY),RoB(NY),PB(NY),ZB(NY),Y(NY),DY(NY),WB(NY)
      COMMON/C1236789/DX(NX)
	  COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
	  COMMON/C17/XPlot,DXPlot,NSecPlot
  	  COMMON/C2346789/ UF(NY)
      COMMON/C1234789R/HB,ZOB
 	  COMMON/C1234789I/NCELLJ

	  COMMON/C378/Dh(NY)
	  COMMON/C3678/XF,NSECTION
	  COMMON/C2346811/QW
	  COMMON/C478/QF
	  common/resid/resdualv
	  common/mass/ ResQmass,Qwmass,Qfmass
	  common/optiout/idoptiout,optionkey
	  Common/C68/IdResV

	  DIMENSION  FCswt(NY),FCnwt(NY),ASC(NY),ANC(NY),SC(NY)



      DO J=1,NNODEJ

	  Dh(j)=ZF(j)-ZB(j)	  !!

	  END DO

	  DO J=2,NNODEJ-1 
	  FCswt(J)=0.25*(RoF(J-1)+RoF(J)+RoB(J-1)+RoB(J))*0.25*(UF(J-1)+UF(J)+UB(J-1)+UB(J))*0.5*(Dh(J)+Dh(J-1))
      FCnwt(J)=0.25*(RoF(J+1)+RoF(J)+RoB(J+1)+RoB(J))*0.25*(UF(J+1)+UF(J)+UB(J+1)+UB(J))*0.5*(Dh(J)+Dh(J+1))
	  END DO

!*--- CORRECTION COEFFICIENT FOR BOUNDERIES CELLS  (1)
 	  FCnwt(NNODEJ-1)=0.5*(RoF(NNODEJ)+RoB(NNODEJ))*0.5*(UF(NNODEJ)+UB(NNODEJ))*Dh(NNODEJ)  
	  FCswt(2)=0.0		   ! end of correction
	 
	  DO J=2,NNODEJ-1
!     *--- MAIN COEFFICIENT
      ASC(J)=-0.25*(RoF(J-1)+RoF(J)+RoB(J-1)+RoB(J))*DX(1)
	  ANC(J)=0.25*(RoF(J+1)+RoF(J)+RoB(J+1)+RoB(J))*DX(1)
	  SC(J)=(HB*RoB(J)*UB(J)-HF*RoF(J)*UF(J))*DY(J)+FCnwt(J)-FCswt(J)
 
      END DO

!*--- CORRECTION COEFFICIENT FOR BOUNDERIES CELLS (2)
      ASC(2)=-0.5*(RoF(1)+RoB(1))*DX(1)
  	  ANC(NNODEJ-1)=0.5*(RoF(NNODEJ)+RoB(NNODEJ))*DX(1)

	  IF (IdResV.eq.1) then
  	  resdualv=0.0d0
	  DO J=2,NNODEJ-1
	  resdualv=resdualv+abs((SC(J)-ASC(J)*WF(J-1)-ANC(J)*WF(J)))
	  END DO
	  resdualv=resdualv/(NNODEJ-2)/Qwmass

	  IdResV=0
	  RETURN
      END IF
!*--- CALCULATION OF W 

	  DO J=2,NNODEJ-1
	  WF(J)=(SC(J)-ASC(J)*WF(J-1))/ANC(J)
	  END DO
!*--- CALCULATION OF V (W IN TRANSFORMED COORDINATE)

1011  FORMAT('ZF(',I3,')=',F14.10,2X,' ZB(',I3,')=',F14.10,2X,'Y(',I3,')=',F12.9,2X,' DY(',I3,')=',F10.8)
1150  FORMAT ('NSEC=',I5,1X,'XF=',F12.9,1X,'HF=',F14.11,1X,'HB=',F14.11,1x,'Qw=',F12.9,1X,'Qf=',F12.9)
1040  FORMAT ('FCnwt(',I3,')=',F16.10,3X,' FCswt(',I3,')=',F16.10,3X,'Dh(',I3,')=',F15.12) 
1050  FORMAT (' ANC(' ,I3,')=',F14.10,2X,' ASC('  ,I3,')=',F14.10,2X,' SC('   ,I3,')=',F16.10) 
1072  FORMAT ('WF(' ,I3,')=',F16.11,1X,' V(' ,I3,')=',F15.11,1X,' UF(' ,I3,')=',F15.11,1X,' UB(' ,I3,')=',F15.11) 
	  RETURN
	  END



      SUBROUTINE MIXINGLENGTH  ! Implement the Mixing Length turbulence model
      PARAMETER(NY=403,NX=1)
	  DOUBLE PRECISION HF,HB,ZOF,ZOB,ZF,DX,ZB,Y,RoF,RoB,UB,UF,ViscosC,difusivity,DY,WF,WB,Akapa,PF,PB,DZF,DZB,UOLD1
	  DOUBLE PRECISION DUDY,DUDX,aLength,DVDY,DVDX,gent,XF
    !  DOUBLE PRECISION NCELLI,MaxIterU,MaxIterUstar,INDEX,MaxIter,NCELLJ,NNODEJ
	  integer optionkey
	  COMMON/C2379/UOLD1(NY),ViscosC(NY),difusivity(NY)
	  COMMON/C2356789/RoF(NY),PF(NY),WF(NY)
	  COMMON/C23789/UB(NY),RoB(NY),PB(NY),ZB(NY),Y(NY),DY(NY),WB(NY)
	  COMMON/C1236789/DX(NX)
	  COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
  	  COMMON/C2346789/ UF(NY)
      COMMON/C1234789R/HB,ZOB
	  COMMON/C1234789I/NCELLJ
	  COMMON/C139/Akapa
	  common/ids/Idconvdif,idfull,idnik,Idprofile
	  common/turtec/DUDY(NY),DUDX(NY),aLength(NY),DVDY(NY),DVDX(NY),gent
	  COMMON/C3678/XF,NSECTION
	  common/optiout/idoptiout,optionkey
	  
	 
!*--- calculation of length
	  DO J=2,NNODEJ
	     if (idnik.eq.0) then 
					IF (0.5*(ZF(J)+ZB(J)).le.0.3*0.5*(HF+HB))THEN
					aLength(J)=Akapa*0.5*(ZF(J)+ZB(J))
					ELSE
					aLength(J)=0.12*0.5*(HF+HB)
     				ENDIF
 	     endif
	  
	  if (idnik.eq.1) then
	  


	    aLength(J)=0.5*(HF+HB)*(0.14-0.08*(1-((ZF(J)+ZB(J))/(HF+HB)))**2-0.06*(1-((ZF(J)+ZB(J))/(HF+HB)))**4)

    

	  end if
	 
	  if (idnik.eq.2) then
	  


	  aLength(J)=0.5*(HF+HB)*(0.412*((1.0-((ZF(J)+ZB(J))/(HF+HB)))**0.5)*(1/((ZF(J)+ZB(J))/(HF+HB))+3.14159*0.2* sin(3.14159*((ZF(J)+ZB(J))/(HF+HB))))**(-1.0)) 

	   
	  endif 


	  END DO	  
!*--- calculation of DUDX
	   DO J=2,NNODEJ-2
       DUDX(J)=((HF*UF(J)-HB*UB(J))*dY(J)-((HF-HB))*((Y(J)+Y(J+1))/2.*0.25*(UF(J)+UF(J+1)+UB(J)+UB(J+1))-(Y(J)+Y(J-1))/2.*0.25*(UF(J)+UF(J-1)+UB(J)+UB(J-1))))/(0.5*(HF+HB)*dx(1)*dY(J))
	   END DO
	   DUDX(NNODEJ-1)=((HF*UF(NNODEJ-1)-HB*UB(NNODEJ))*dY(NNODEJ-1)-(HF-HB)*(Y(NNODEJ)*(UF(NNODEJ-1)+UB(NNODEJ-1))/2.0-(Y(NNODEJ-1)+Y(NNODEJ-2))/2.*(UF(NNODEJ-1)+UF(NNODEJ-2)+UB(NNODEJ-1)+UB(NNODEJ-2))/4.0))/(0.5*(HF+HB)*dx(1)*dY(NNODEJ-1))
      !end of DUDX calculation

!*--- calculation of DUDY
	   DO J=3,NNODEJ-2
       DUDY(J)=0.5*(UF(J+1)+UB(J+1)-UF(J-1)-UB(J-1))/(HF+HB)/DY(J)
	   END DO
       DUDY(2)=0.5*(UF(3)+UB(3)+UF(2)+UB(2))/(HF+HB)/DY(2)
	   DUDY(NNODEJ-1)=(0.5*(UF(NNODEJ-1)+UB(NNODEJ))-0.25*(UF(NNODEJ-1)&
	&  +UB(NNODEJ-1)+UF(NNODEJ-2)+UB(NNODEJ-2)))/(0.5*(HF+HB))/DY(NNODEJ-1)
	 

!*--- calculation of DVDX
       DO J=2,NNODEJ-2
	   if (XF.le.DX(1)) then
	   DVDX(J)=0.0d0
	   else
	   DVDX(J)=((HF*(WF(J)+WF(J-1))*0.5-HB*(WB(J)+WB(J-1))*0.5)*dY(J)-(HF-HB)*(((Y(J)+Y(J+1))/2.*WF(J))-((Y(J)+Y(J-1))/2.*WF(J-1))))/(0.5*(HF+HB)*dx(1)*dY(J))
	   end if
	   enddo

	   DVDX(NNODEJ-1)=((HF*(WF(NNODEJ-1)+WF(NNODEJ-2))*0.5-HB*(WB(NNODEJ-1)+WB(NNODEJ-2))*0.5)*dY(NNODEJ-1)-(HF-HB)*(Y(NNODEJ)*WF(NNODEJ-1)-0.5*(Y(NNODEJ-1)+Y(NNODEJ-2))*WF(NNODEJ-2)))/(0.5*(HF+HB)*dx(1)*dY(NNODEJ-1))
	   !end of DVDX calculation

5555  format(f12.3,2x,I3,2x,d15.5,2x,d15.5,2x,d15.5)
!*--- calculation of DVDY
	  DO J=2,NNODEJ-1
	  DVDY(J)=2.0*(WF(J)-WF(J-1))/(HF+HB)/DY(J)
	  END DO

	   
 	  DO J=2,NNODEJ-1

	  if(idfull.eq.0) gent=DUDY(J)
	  if(idfull.eq.1) gent=(2*(DUDX(J)**2)+2*(DVDY(J))**2+(DUDY(J))**2)**0.5 !dv/dx
	  if(idfull.eq.2) gent=(2*(DUDX(J)**2)+(DUDY(J))**2)**0.5 !dv/dx dv/dy
	  if(idfull.eq.3) gent=(2*(DUDX(J)**2)+2*(DVDY(J))**2+(DUDY(J)+DVDX(J))**2)**0.5

	  !dv/dx ******************************************************************************

     ViscosC(J)=(RoF(J)+RoB(J))/2.*((aLength(j)**2.*gent)+1D-6)


!****************************************************88888888
	  END DO
   
IDdifusivity=0	!: Eddy diffusivity equals to eddy viscosity
!  IDdifusivity=1	!: Constant Eddy diffusivity in upper region
! PAY ATTENTION TO MESH DEPENDENCY BECAUSE OF NUMBER OF MESH !!
DO J=1,NNODEJ

  if (IDdifusivity.eq.1) then
 
  

     	IF (0.5*(ZF(J)+ZB(J)).le.0.6*0.5*(HF+HB))THEN
		   difusivity(J)=ViscosC(J)
		else
           difusivity(J)=ViscosC(0.2*NNODEJ)/5.0
		endif
 

 	   

  else
   difusivity(J)=ViscosC(J)
  endif
ENDDO




	  RETURN
	  END


!------------------------------------------------------------------------

	subroutine ROW_CALCULATION	! Calculating the air-water density
	PARAMETER (NY=403,NX=1)

	  DOUBLE PRECISION HF,ZOF,ZF,ZB,DZF,DZB,Y,DY,WB,WF,PF,PB,RoF,RoB,UB,phiaB,phiaF,Vspilt,phicol,&
	& eshmit

	  integer optionkey
	  COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
	  COMMON/C23789/UB(NY),RoB(NY),PB(NY),ZB(NY),Y(NY),DY(NY),WB(NY)
	  COMMON/C2356789/RoF(NY),PF(NY),WF(NY)  	  
	  common/covdif/phiaB(NY),phiaF(NY),Vspilt,phicol,eshmit
	  common/optiout/idoptiout,optionkey
	 	  
	  rowaterpair=800.000d0
	  	  	  
 	  DO J=1,NNODEJ
	  RoF(J)=1000.0/(1+(PhiaF(J)*(rowaterpair-1.0)))
      end do

	  return
	  end
!-----------------------------------------------------------------------

	  subroutine C_CALCULATION   !  Calculating the air volume concentration
	PARAMETER (NY=403,NX=1)

	  DOUBLE PRECISION HF,ZOF,ZF,DZF,DZB,phiaB,phiaF,Vspilt,phicol,eshmit,&
	& rowaterpair,CAF,CAB
	  integer optionkey
	  COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB	 	  
	  common/covdif/phiaB(NY),phiaF(NY),Vspilt,phicol,eshmit
	  common/condiC/CAF(NY),CAB(NY)
	  common/optiout/idoptiout,optionkey
	  
	  rowaterpair=800.000d0
	  	  	  
 	  DO J=1,NNODEJ
	  CAF(j)=PhiaF(j)*rowaterpair/(1.0d0+PhiaF(j)*(rowaterpair-1))	  
	  end do
	  return
	  end
!-----------------------------------------------------------------------

	subroutine TECPLOT(numzone)  ! Preparing data for Tecplot software
	PARAMETER (NY=403,NX=1)
	double precision ZINITIAL,ZEND,Qinitial,XF,HF,ZOF,ZF,DZF,DZB,HB,ZOB,XLENGTH,XPlot,DXPlot
	integer optionkey
 	COMMON/C12/NCELLI,ZINITIAL,ZEND,Qinitial
	COMMON/C3678/XF,NSECTION
	COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
	COMMON/C1234789R/HB,ZOB
	COMMON/C1234789I/NCELLJ
	COMMON/C1267/XLENGTH
	COMMON/C17/XPlot,DXPlot,NSecPlot
	common/optiout/idoptiout,optionkey
	
 
	numzone=numzone+1  
	write(*,*)'writing tecplot zone ...',numzone
	if(numzone.eq.1) then
        write(15,*)'TITLE = "concentration profil"'
        write(15,*)'VARIABLES=x,y,u,v,p,phia,vis,ro,delta,Vs,Sch,dr,rich'	  
	  		
	end if
    ni=NsecPlot+1
  
  	nj=NNODEJ
      
	write(15,1)numzone,ni,nj
1     format('ZONE T="',i4,'", J=',i6,', I=',i6,', F=POINT')
    if(idoptiout.eq.1) write(17,1)numzone,ni,nj
	return 
	end
!c-------------------------------------------------


    SUBROUTINE RISING_VELOCITY	! Calculating the rising velocity
    PARAMETER (NY=403,NX=1)
    Double precision HF,ZOF,ZF,DZF,DZB
	Double precision phiaB,phiaF,Vspilt,phicol,eshmit,VraisF,VraisB
	Double precision CAF,CAB,g,teta,relaxphi
	COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
	COMMON/raisingvelosity/VraisF(NY),VraisB(NY)
	COMMON/C1235/g,teta
	common/covdif/phiaB(NY),phiaF(NY),Vspilt,phicol,eshmit	
	common/condiC/CAF(NY),CAB(NY)
	common/covdiB/ idboundbottm,idboundsurf,IdChanson,idVspAuto,relaxphi,IdLinear

	
	Do J=2,NNODEJ

	if(idVspAuto.eq.1) then

!----Specifing Raising Velosity according to C air and Bobble Size

    IF ( CAF(J).LE.0.04) Then    
 	VraisF(J)=(0.06/0.04)*CAF(J)
	endif 
	IF ( CAF(J).GT.0.04.AND.CAF(J).LE.0.15) Then
	VraisF(J)=(0.25-0.06)/(0.15-0.04)*(CAF(J)-0.04)+0.06
    endif
	IF ( CAF(J).GT.0.15.AND.CAF(J).LE.0.7) Then   
 	VraisF(J)=0.25
	END IF
	IF ( CAF(J).GT.0.7.AND.CAF(J).LE.1.0) Then
	VraisF(J)=0.25
   
	endif


   	ELSE

	VraisF(J)=Vspilt

	ENDIF

!----Chanson
	if(IdChanson.eq.1) then


!	old Chanson coeficient 
	VraisF(J)=VraisF(J)*((1.0-CAF(J))**0.5)

	
	end if

	END DO
 

   	IF(IdLinear==1) then 

	
	DO J=int(NNODEJ/2.0),NNODEJ
 	VraisF(J)=VraisF(int(NNODEJ/2.0))/(HF-ZF(int(NNODEJ/2.0)))*(HF-ZF(J))

   	end do
	ENDIF
	return
	END

	! if(idVspAuto.eq.1) then

!----Specifing Raising Velosity according to C air and Bobble Size

!	VraisF(J)=((4E-05*((CAF(J)*100.00)**2)+0.4512*(CAF(J)*100.00)+ 14.27)/100.00)
	

 !  	ELSE
!	VraisF(J)=Vspilt
!	ENDIF

!*--- 
	
		   ! 11th subroutine ConvectionDiffusion

! this version
! 1-option 3 
! 2-dV/dy in turb. is omitted
! 3-PhiaB(nnj)=0.0 in initial

	 

      SUBROUTINE CONVECTION_DIFFUSION	  ! Solving the convection- diffusion equations
      PARAMETER(NY=403,NX=1)
	  DOUBLE PRECISION HF,HB,ZOF,ZOB,ZF,DZF,DZB,DX,PF,PB,ZB,Dh,Y,RoF,RoB,UB,UF,g,tgteta,&
	& Hinitial,ViscosC,difusivity,DY,WF,WB,XF,UOLD1,phiaB,phiaF,Vspilt,eshmit,eshmitn,DRoDY,rich,xplot,dxplot,phicol,sumdif
	  DOUBLE PRECISION Fsw,Fnw,FP,FW,Dsw,Dnw,Fswt,Fnwt,AS,AP,AN,S
	  DOUBLE PRECISION CC,UU,DD,XX,BB,teta,AInitialphi0,AInitialphi2,AInitialphi1,CAF,CAB
	  DOUBLE PRECISION maxdif,Phiaold,relaxphi,resphi,ReferencePhi,Qw,VraisF,VraisB
	  DOUBLE PRECISION RoSw,RoNw,VraisSw,VraisNw,RoNNODEJ,VraisNNODEJ,PhiaNNODEJ
      integer optionkey

	  COMMON/C2379/UOLD1(NY),ViscosC(NY),difusivity(NY)
	  COMMON/C2356789/RoF(NY),PF(NY),WF(NY)
      COMMON/C123/tgteta,Hinitial
	  COMMON/C1235/g,teta
	  COMMON/C23789/UB(NY),RoB(NY),PB(NY),ZB(NY),Y(NY),DY(NY),WB(NY)
	  COMMON/C1236789/DX(NX)
	  COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
  	  COMMON/C2346789/ UF(NY)
      COMMON/C1234789R/HB,ZOB
	  COMMON/C1234789I/NCELLJ
      COMMON/C2346811/QW
      COMMON/C378/Dh(NY)
      COMMON/C3678/XF,NSECTION
	  common/covdif/phiaB(NY),phiaF(NY),Vspilt,phicol,eshmit
	  common/covdiB/ idboundbottm,idboundsurf,idshansen,idVspAuto,relaxphi,IdLinear
	  common/covdiI/AInitialphi0,AInitialphi2,AInitialphi1,zjn1,zjn2,xin1,xin2
	  common/condiC/CAF(NY),CAB(NY)
	  COMMON/C17/XPlot,DXPlot,NSecPlot
	  common/optiout/idoptiout,optionkey
      COMMON/raisingvelosity/ VraisF(NY),VraisB(NY)
	  Common/C611/IdResPhi,resphi
	  common/eshmitno/eshmitn(NY),DRoDY(NY),rich(NY)

	  DIMENSION Fsw(NY),Fnw(NY),FP(NY),FW(NY),Dsw(NY),Dnw(NY),Fswt(NY),Fnwt(NY),AS(NY),AP(NY),AN(NY),S(NY)
	  DIMENSION CC(NY-2),UU(NY-3),DD(NY-3),XX(NY-2),BB(NY-2)
	  Dimension Phiaold(NY)

	  
	  Iterro=0
	  DO 	  !MAIN INTERNAL LOOP
	  Iterro=Iterro+1
	  if (iterro.gt.9999) then
	  write(*,*)'iter=',iterro

	  end if

	   Xselfaeration=11.4d0		  ! for cain 450mm fisrt xs=3.0



!      edit version for gradually self aeration
	   idgrad=1
	   if (idgrad.eq.1) then 
	       if (idboundsurf.eq.1) then
		
				Xself2=3.0
				if(XF.ge.(Xself2+Xselfaeration)) then
				PhiaF(NNODEJ)=AInitialphi2
				else
			    	if (XF.le.(Xself2+Xselfaeration).and.XF.ge.Xselfaeration) then
	     			 PhiaF(NNODEJ)=(XF-Xselfaeration)/Xself2*AInitialphi2
	    
			    	else 
				    PhiaF(NNODEJ)=0.0d0
				    PhiaB(NNODEJ)=0.0d0
				    end if
				end if
	       end if
	  
!    end of  edit version for gradually self aeration 
	   else 

	        if (idboundsurf.eq.1) then
				if(XF.ge.Xselfaeration) then
				PhiaF(NNODEJ)=AInitialphi2
				else
				PhiaF(NNODEJ)=0.0d0
				PhiaB(NNODEJ)=0.0d0
				end if
	         end if
    	endif    


	   DO J=1,NNODEJ-1
	    Dh(J)=ZF(J)-ZB(J)
	   END DO

 



  ! SCHMIDT VARIABLE CALCULATION
IdSchmit=0

 if (IdSchmit.eq.1) then 

    DO J=1,NNODEJ


    IF ( CAF(J).LE.0.8) Then    
 	eshmitn(J)=(1.0-0.6)/(0.8-0.0)*(CAF(J)-0.0)+0.6
	
    endif 
	IF ( CAF(J).GT.0.8.AND.CAF(J).LE.1.0) Then
    eshmitn(J)=(0.8-1.0)/(0.9-0.8)*(CAF(J)-0.8)+1.0
	endif


  

	ENDDO



 else 
     DO J=1,NNODEJ
		eshmitn(J)=eshmit
	 ENDDO	 
 endif


if ( eshmit.eq.0) then 
call  cal_Schmidt
end if
!end of  SHCMIDT VARIABLE CALCULATION


	   Call RISING_VELOCITY


	   VraisF(1)=0.0d0
	   VraisF(NNODEJ)=0.00d0
	   VraisB(1)=0.0d0
	   VraisB(NNODEJ)=0.00d0
	   

	   
	   DO J=2,NNODEJ-1
	 !UTILITY Coeficients
	  RoSw=0.25*(RoF(J-1)+RoF(J)+RoB(J-1)+RoB(J))
	  RoNw=0.25*(RoF(J)+RoF(J+1)+RoB(J)+RoB(J+1))
	  VraisSw=0.25*(VraisF(J-1)+VraisF(J)+VraisB(J-1)+VraisB(J))
	  VraisNw=0.25*(VraisF(J+1)+VraisF(J)+VraisB(J+1)+VraisB(J))
	  PhiaSw=0.25*(PhiaF(J)+PhiaB(J)+PhiaF(J-1)+PhiaB(J-1))
	  PhiaNw=0.25*(PhiaF(J)+PhiaB(J)+PhiaF(J+1)+PhiaB(J+1))
	 ! END OF UTILITY COEFFICIENTS
	   Fsw(J)=RoSw*(WF(J-1)+(1-PhiaSw)*VraisSw*cosd(teta))*DX(1)
	   Fnw(J)=RoNw*(WF(J)  +(1-PhiaNw)*VraisNw*cosd(teta))*DX(1)
	   Dsw(J)=0.5*(difusivity(J-1)/eshmitn(J-1)+difusivity(J)/eshmitn(J))*DX(1)/(0.5*(HF+HB))/(0.5*(DY(J)+DY(J-1)))
	   Dnw(J)=0.5*(difusivity(J+1)/eshmitn(J+1)+difusivity(J)/eshmitn(J))*DX(1)/(0.5*(HF+HB))/(0.5*(DY(J)+DY(J+1)))
	 
	  
	   FW(J)=RoB(J)*(UB(J)-(1-PhiaB(J))*VraisB(J)*sind(teta))*HB*DY(J)	   	  
       FP(J)=RoF(J)*(UF(J)-(1-PhiaF(J))*VraisF(J)*sind(teta))*HF*DY(J)	   
	   Fswt(J)=RoSw*(0.25*(UF(J-1)+UF(J)+UB(J-1)+UB(J))-(1-PhiaSw)*VraisSw*sind(teta))*0.5*(Dh(J)+Dh(J-1))  
       Fnwt(J)=RoNw*(0.25*(UF(J+1)+UF(J)+UB(J+1)+UB(J))-(1-PhiaNw)*VraisNw*sind(teta))*0.5*(Dh(J)+Dh(J+1))  


!     floating strings
 	  
	   enddo
!	  
!*--- CORRECTION COEFFICIENT FOR BOTTOM CELL (1)

	  fsw(2)=0.0d0
	  Fswt(2)=0.0d0	  
	  Dsw(2)=0.0d0 !could be assigned zero just for partial boundary at bottom (Dphi/Dy=0.0),check it 
	        	 	  
!*--- CORRECTION COEFFICIENT FOR UPER MOST CELL (1)


if (optionkey.eq.0) then	  
	   PhiaNNODEJ=0.5*(PhiaF(NNODEJ)+PhiaB(NNODEJ))	 
	   RoNNODEJ=0.5*(RoF(NNODEJ)+RoB(NNODEJ))
	   VraisNNODEJ=0.5*(VraisF(NNODEJ)+VraisB(NNODEJ))
	   Fnw(NNODEJ-1)= 0.5*(RoF(NNODEJ)+Rob(NNODEJ))*(WF(NNODEJ-1)+(1-0.5*(PhiaF(NNODEJ)+PhiaB(NNODEJ)))*(VraisNNODEJ)*cosd(teta))*DX(1)	          	   
       Dnw(NNODEJ-1)=0.5*(difusivity(NNODEJ)/eshmitn(NNODEJ)+difusivity(NNODEJ-1)/eshmitn(NNODEJ-1))*DX(1)/(0.5*(HF+HB))/(0.5*(DY(NNODEJ-1)+DY(NNODEJ)))
	   Fnwt(NNODEJ-1)=0.5*(RoF(NNODEJ)+Rob(NNODEJ))*(0.5*(UF(NNODEJ-1)+UB(NNODEJ))-(1-0.5*(PhiaF(NNODEJ)+PhiaB(NNODEJ)))*0.5*(VraisF(NNODEJ)+VraisB(NNODEJ))*sind(teta))*Dh(NNODEJ)
end if

if (optionkey.eq.3.or.optionkey.eq.4) then
	   RoNNODEJ=0.25*(RoF(NNODEJ)+RoB(NNODEJ)+RoF(NNODEJ-1)+RoB(NNODEJ-1))
	   VraisNNODEJ=0.25*(VraisF(NNODEJ)+VraisB(NNODEJ)+VraisF(NNODEJ-1)+VraisB(NNODEJ-1))
	   Fnw(NNODEJ-1)= 0.25*(RoF(NNODEJ-1)+Rob(NNODEJ-1)+RoF(NNODEJ)+Rob(NNODEJ))*(WF(NNODEJ-1)+(1-0.25*(PhiaF(NNODEJ-1)+PhiaB(NNODEJ-1)+PhiaF(NNODEJ)+PhiaB(NNODEJ)))*VraisNNODEJ*cosd(teta))*DX(1)	          	   
       Dnw(NNODEJ-1)=0.5*(difusivity(NNODEJ)/eshmitn(NNODEJ)+difusivity(NNODEJ-1)/eshmitn(NNODEJ-1))*DX(1)/(0.5*(HF+HB))/(0.5*(DY(NNODEJ-1)+DY(NNODEJ)))
       Fnwt(NNODEJ-1)=0.25*(RoF(NNODEJ-1)+Rob(NNODEJ-1)+RoF(NNODEJ)+Rob(NNODEJ))*(0.5*(UF(NNODEJ-1)+UB(NNODEJ))-(1-0.25*(PhiaF(NNODEJ-1)+PhiaB(NNODEJ-1)+PhiaF(NNODEJ)+PhiaB(NNODEJ)))*VraisNNODEJ*sind(teta))*((Dh(NNODEJ)))  
if (optionkey.eq.4) Dnw(NNODEJ-1)=Dnw(NNODEJ-1)/2.0
end if
if (optionkey.eq.1) then	  
	   PhiaNNODEJ=0.5*(PhiaF(NNODEJ)+PhiaB(NNODEJ))	 
	   RoNNODEJ=0.5*(RoF(NNODEJ)+RoB(NNODEJ))
       Fnw(NNODEJ-1)= 0.5*(RoF(NNODEJ-1)+Rob(NNODEJ-1))*(WF(NNODEJ-1)+(1-0.5*(PhiaF(NNODEJ-1)+PhiaB(NNODEJ-1)))*(0.5*(VraisF(NNODEJ-1)+VraisB(NNODEJ-1)))*cosd(teta))*DX(1)	          	   
       Dnw(NNODEJ-1)=difusivity(NNODEJ)/eshmitn(NNODEJ)*DX(1)/(0.5*(HF+HB))/(0.5*(DY(NNODEJ-1)+DY(NNODEJ)))
	   Fnwt(NNODEJ-1)=0.5*(RoF(NNODEJ)+Rob(NNODEJ))*(0.5*(UF(NNODEJ)+UB(NNODEJ))-(1-0.5*(PhiaF(NNODEJ)+PhiaB(NNODEJ)))*0.5*(VraisF(NNODEJ-1)+VraisB(NNODEJ-1))*sind(teta))*(Dh(NNODEJ))  
end if


!     *--- MAIN COEFFICIENT
      DO J=2,NNODEJ-1
	  AS(J)= Fswt(J)/4.0 -Dsw(J)/2.0 -Fsw(J)/4.0
      AN(J)=-Fnwt(J)/4.0 -Dnw(J)/2.0 +Fnw(J)/4.0
      AP(J)= FP(J)+Fswt(J)/4.0 -Fnwt(J)/4.0 +Dsw(J)/2.0 +Dnw(J)/2.0 +Fnw(J)/4.0-Fsw(J)/4.0
	  S(J)=	(Fw(J)-Fswt(J)/4.0+fnwt(J)/4.0-Fnw(J)/4.0+Fsw(J)/4.0-Dnw(J)/2.0-Dsw(J)/2.0)*phiaB(J)+(0.0-Fswt(J)/4.0+Fsw(J)/4.0+Dsw(J)/2.0)*phiaB(J-1)+(Fnwt(J)/4.0-Fnw(J)/4.0+Dnw(J)/2.0)*phiaB(J+1)	  

	  END DO


!*--- CORRECTION COEFFICIENT FOR BOTTOM CELL (2)
!     Bottom 
	  j=2
	  AS(J)= 0.0
	  
!     Derative Boundary 
      AP(J)= FP(J)+Fswt(J)/2.0 -Fnwt(J)/4.0  +Dnw(J)/2.0 +Fnw(J)/4.0-Fsw(J)/2.0
	  S(J)=	(Fw(J)+fnwt(J)/4.0-Fnw(J)/4.0-Dnw(J)/2.0)*phiaB(J)+(Fnwt(J)/4.0-Fnw(J)/4.0+Dnw(J)/2.0)*phiaB(J+1)+(-Fswt(J)/2.0+Fsw(J)/2.0)*phiaB(J-1)


! *--- CORRECTION COEFFICIENT FOR UPER MOST CELL (2)
!     Surface 
	  J=NNODEJ-1
	  AN(J)= 0.0
      if(idboundsurf.eq.0) then
!     Derative Boundary 
      AP(J)= FP(J)+Fswt(J)/4.0 -Fnwt(J)/2.0  +Dsw(J)/2.0 +Fnw(J)/2.0-Fsw(J)/4.0
	  S(J)=	(Fw(J)-Fswt(J)/4.0+Fsw(J)/4.0-Dsw(J)/2.0)*phiaB(J)+(-Fswt(J)/4.0+Fsw(J)/4.0+Dsw(J)/2.0)*phiaB(J-1)+(fnwt(J)/2.0-Fnw(J)/2.0)*phiaB(J+1)
	  else
!	  Deritchleh

if (optionkey.eq.0) then	  
      AP(J)= FP(J)+Fswt(J)/4.0 +Dsw(J)/2.0 +Dnw(J)/2.0 -Fsw(J)/4.0
	  S(J)=	(Fw(J)-fswt(J)/4.0+Fsw(J)/4.0-Dnw(J)/2.0-Dsw(J)/2.0)*phiaB(J)+(-Fswt(J)/4.0+Fsw(J)/4.0+Dsw(J)/2.0)*phiaB(J-1)+(Fnwt(J)/2.0-Fnw(J)/2.0+Dnw(J)/2.0)*phiaB(J+1)+(Fnwt(J)/2.0-Fnw(J)/2.0+Dnw(J)/2.0)*PhiaF(J+1)


end if

if (optionkey.eq.3.or.optionkey.eq.4) then	  
	  S(J)=	(Fw(J)-fswt(J)/4.0+Fnwt(J)/4.0-Fnw(J)/4.0+Fsw(J)/4.0-Dnw(J)/2.0-Dsw(J)/2.0)*phiaB(J)+(-Fswt(J)/4.0+Fsw(J)/4.0+Dsw(J)/2.0)*phiaB(J-1)+(Fnwt(J)/4.0-Fnw(J)/4.0+Dnw(J)/2.0)*phiaB(J+1)+(Fnwt(J)/4.0-Fnw(J)/4.0+Dnw(J)/2.0)*PhiaF(J+1)
end if
if (optionkey.eq.1) then	  
      AP(J)= FP(J)+Fswt(J)/4.0 +Dsw(J)/2.0 +Dnw(J)/2.0 -Fsw(J)/4.0-Fnwt(J)/2.0+Fnw(J)/2.0
	  S(J)=	(Fw(J)-fswt(J)/4.0+Fsw(J)/4.0-Dnw(J)/2.0-Dsw(J)/2.0+Fnwt(J)/2.0-Fnw(J)/2.0)*phiaB(J)+(-Fswt(J)/4.0+Fsw(J)/4.0+Dsw(J)/2.0)*phiaB(J-1)+((Dnw(J)/2.0)*phiaB(J+1))+((Dnw(J)/2.0)*phiaF(J+1))


end if



	  end if
	  
	  IF (IdResPhi.eq.1) then
	  resphi=0.0d0
      DO J=2,NNODEJ-1
      resphi=resphi+abs(AN(J)*PhiaF(J+1)+AP(J)*PhiaF(J)+AS(J)*PhiaF(J-1)-S(J))
      enddo
	  ReferencePhi=100*0.011*(Qw/Hinitial)*sind(teta)*Hinitial !Acoording to 90% air concentration at free surface
 	
	 
      resphi=resphi/ReferencePhi
	  IdResPhi=0
 	  RETURN
	  END if

  !  ---PRINT DATA
  !	  GO FOR TDMA
	  DO J = 1,NCELLJ
	  CC(J)=AP(J+1)
	  BB(J)=S(J+1)
	  XX(J)=phiaF(J+1)
	  END DO
	  DO J= 1,NCELLJ-1
	  UU(J)=AN(J+1)
	  DD(J)=AS(J+2)
	  END DO
  

	  CALL TDMA (NCELLJ,CC,UU,DD,XX,BB) 
	  DO I=1,NCELLJ	  	  
	  Phiaold(I+1)=PhiaF(I+1)
	  PhiaF(I+1)=XX(I)
	  end do
!---under relaxation
	  do j=2,NNODEJ-1
	  phiaF(j)=relaxphi*PhiaF(j)+(1.0-relaxphi)*Phiaold(j)
	  end do

	  sumdif=0.0d0
	  maxdif=0.0d0	  
      DO j=2,NNODEJ-1
	  sumdif=sumdif+abs(Phiaold(j)-phiaF(j))
	  if(abs(phiaF(j)-Phiaold(j)).gt.maxdif) then 
	  maxdif=abs(phiaF(j)-Phiaold(j))
	  idmaxdif=j
	  end if
	  END DO

	  if(PhiaF(idmaxdif).ne.0.0) then
	  maxdif=maxdif/PhiaF(idmaxdif)
	  else

	  end if
	  phiaF(1)=phiaF(2)
	  if(idboundsurf.eq.0) phiaF(NNODEJ)=phiaF(NNODEJ-1)
	  call ROW_CALCULATION	  

	  if(maxdif.lt.0.005.and.iterro.gt.1) EXIT

	  END DO
      call C_CALCULATION	  

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  	  INDEX=0
  IF (INDEX > 0) THEN
	 WRITE(1,1152) NSECTION,XF,HF,HB
     WRITE(1,1010) (J,ZF(J),J,Y(J),J,DY(J),J=1,NNODEJ) !
	 WRITE(1,1020) (J,PB(J),J,PF(J),J,Dh(J),J=1,NNODEJ)
     WRITE(1,1030) (J,FW(J),J,FP(J),J,Fnw(J),J,Fsw(J),J=2,NNODEJ-1) !
	 WRITE(1,1040) (J,Dnw(J),J,Dsw(J),J,Fnwt(J),J,Fswt(J),J=2,NNODEJ-1)
	 WRITE(1,1050) (J,AN(J),J,AP(J),J,AS(J),J,S(J),J=2,NNODEJ-1)
     WRITE(1,1070) (J,UF(J),J,WF(J),J=1,NNODEJ)
	 WRITE(1,1575) (J,phiaB(J),J,phiaF(J),J=1,NNODEJ)

  END IF
1010  FORMAT(' ZF(',I3,')=',F10.8,10X,'Y(',I3,')=',F10.7,10X,' DY(',I3,')=',F10.8)
1020  FORMAT(' PB(',I3,')=',F10.4,10X,' PF(',I3,')=',F10.4,10X,'Dh(',I3,')=',F16.10)
1030  FORMAT ('FW(',I2,')= ',F13.7,1X,'FP(',I2,')=' ,F13.7,1X,'Fnw(' ,I2,')=',F15.7,1X,' Fsw(' ,I2,')=',F15.7) 
1040  FORMAT ('Dnw(',I2,')=',F11.8,1X,' Dsw(',I2,')=',F11.8,1X,' Fnwt(',I2,')=',F15.7,1X,' Fswt(',I2,')=',F15.7) 
1050  FORMAT (' AN(' ,I3,')=',F12.7,1X,' AP(' ,I3,')=',F12.7,1X,' AS('  ,I3,')=',F12.7,1X,' S('   ,I3,')=',F15.7) 
1060  FORMAT (/,F10.5,4X,' Ustar=' ,F10.5,4X,' Dsw2=',F10.5) 
1070  FORMAT (' UF(' ,I3,')=',F14.10,10X,' WF(' ,I3,')=',F10.6) 
1575  FORMAT (' phiaB(' ,I3,')=',F14.10,10X,' phiaF(' ,I3,')=',F10.6) 
1075  FORMAT (//' THE VALUE OF Sf IS NEGATIVE ...')
1152  FORMAT ('NSEC=',I7,2X,'XF=',F18.9,2X,'HF=',F16.12,5X,'HB=',F16.8)
1077  FORMAT (/'XF=',F12.6,2X,'HF=',F8.6,2X,'HB=',F8.6,2X,'SF=',F12.6/)	 
 	  RETURN
	  END
!------------------
!---------------------------

	
    !  SUBROUTINE CONTINUITY5
    !  PARAMETER(NY=403 ,NX=1)	 
	!  DOUBLE PRECISION phiaB,phiaF,Vspilt,phicol,eshmit
    !  DOUBLE PRECISION HF,HB,ZOF,ZOB,ZF,DZF,DZB,DX,PF,PB,ZB,Dh,Y,RoF,RoB,UB,UF,DY,WF,WB,XF,QF,Qw
   	!  DOUBLE PRECISION FCswt,FCnwt,ASC,ANC,SC,resdualv
	!  DOUBLE PRECISION ResQmass,Qwmass,Qfmass,XPlot,DXPlot
	!  integer optionkey
	!  COMMON/C2356789/RoF(NY),PF(NY),WF(NY)
	!  COMMON/C23789/UB(NY),RoB(NY),PB(NY),ZB(NY),Y(NY),DY(NY),WB(NY)
    !  COMMON/C1236789/DX(NX)
	!  COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
	!  COMMON/C17/XPlot,DXPlot,NSecPlot
  	!  COMMON/C2346789/ UF(NY)
    !  COMMON/C1234789R/HB,ZOB
 	!  COMMON/C1234789I/NCELLJ
	!  COMMON/C378/Dh(NY)
	!  COMMON/C3678/XF,NSECTION
	!  COMMON/C2346811/QW
	!  COMMON/C478/QF
	!  common/resid/resdualv
	!  common/mass/ ResQmass,Qwmass,Qfmass
	!  common/optiout/idoptiout,optionkey
	!  Common/C68/IdResV
	!  common/covdif/phiaB(NY),phiaF(NY),Vspilt,phicol,eshmit
	  

	!  DIMENSION  FCswt(NY),FCnwt(NY),ASC(NY),ANC(NY),SC(NY)


    !  DO J=1,NNODEJ

	!  Dh(j)=ZF(j)-ZB(j)	  !!

	!  END DO

	!  DO J=2,NNODEJ-1 
	!  phisw=0.25*(PhiaF(J-1)+PhiaF(J)+PhiaB(J-1)+PhiaB(J))
	!  phinw=0.25*(PhiaF(J+1)+PhiaF(J)+PhiaB(J+1)+PhiaB(J))
	!  FCswt(J)=0.25*(RoF(J-1)+RoF(J)+RoB(J-1)+RoB(J))/(1.0-phisw)*0.25*(UF(J-1)+UF(J)+UB(J-1)+UB(J))*0.5*(Dh(J)+Dh(J-1))
    !  FCnwt(J)=0.25*(RoF(J+1)+RoF(J)+RoB(J+1)+RoB(J))/(1.0-phinw)*0.25*(UF(J+1)+UF(J)+UB(J+1)+UB(J))*0.5*(Dh(J)+Dh(J+1))
	!  END DO

!*--- CORRECTION COEFFICIENT FOR BOUNDERIES CELLS  (1)
	!  phinw=0.5*(PhiaB(NNODEJ)+PhiaF(NNODEJ))
 	!  FCnwt(NNODEJ-1)=0.5*(RoF(NNODEJ)+RoB(NNODEJ))/(1.0-phinw)*0.5*(UF(NNODEJ)+UB(NNODEJ))*Dh(NNODEJ)  
	!  FCswt(2)=0.0		   ! end of correction
	 
	!  DO J=2,NNODEJ-1
!     *--- MAIN COEFFICIENT
	!  phisw=0.25*(PhiaF(J-1)+PhiaF(J)+PhiaB(J-1)+PhiaB(J))
	!  phinw=0.25*(PhiaF(J+1)+PhiaF(J)+PhiaB(J+1)+PhiaB(J))
    !  ASC(J)=-0.25*(RoF(J-1)+RoF(J)+RoB(J-1)+RoB(J))*DX(1)/(1.0-phisw)
	!  ANC(J)=0.25*(RoF(J+1)+RoF(J)+RoB(J+1)+RoB(J))*DX(1)/(1.0-phinw)
	!  SC(J)=(HB*RoB(J)/(1.0-PhiaB(J))*UB(J)-HF*RoF(J)/(1.0-PhiaF(J))*UF(J))*DY(J)+FCnwt(J)-FCswt(J)
 
    !  END DO

!*--- CORRECTION COEFFICIENT FOR BOUNDERIES CELLS (2)
    !  ASC(2)=-0.5*(RoF(1)+RoB(1))*DX(1)/(1.0-0.5*(PhiaF(1)+PhiaB(1)))
  	!  ANC(NNODEJ-1)=0.5*(RoF(NNODEJ)+RoB(NNODEJ))/(1.0-0.5*(PhiaF(NNODEJ)+PhiaB(NNODEJ)))*DX(1)

	!  IF (IdResV.eq.1) then
  	!  resdualv=0.0d0
	!  DO J=2,NNODEJ-1
	!  resdualv=resdualv+abs((SC(J)-ASC(J)*WF(J-1)-ANC(J)*WF(J)))
	!  END DO
	!  resdualv=resdualv/(NNODEJ-2)/Qwmass
	!  IdResV=0
	!  RETURN
    !  END IF
!*--- CALCULATION OF W 

	!  DO J=2,NNODEJ-1
	!  WF(J)=(SC(J)-ASC(J)*WF(J-1))/ANC(J)
	!  END DO

! 1011  FORMAT('ZF(',I3,')=',F14.10,2X,' ZB(',I3,')=',F14.10,2X,'Y(',I3,')=',F12.9,2X,' DY(',I3,')=',F10.8)
! 1150  FORMAT ('NSEC=',I5,1X,'XF=',F12.9,1X,'HF=',F14.11,1X,'HB=',F14.11,1x,'Qw=',F12.9,1X,'Qf=',F12.9)
! 1040  FORMAT ('FCnwt(',I3,')=',F16.10,3X,' FCswt(',I3,')=',F16.10,3X,'Dh(',I3,')=',F15.12) 
! 1050  FORMAT (' ANC(' ,I3,')=',F14.10,2X,' ASC('  ,I3,')=',F14.10,2X,' SC('   ,I3,')=',F16.10) 
! 1072  FORMAT ('WF(' ,I3,')=',F16.11,1X,' V(' ,I3,')=',F15.11,1X,' UF(' ,I3,')=',F15.11,1X,' UB(' ,I3,')=',F15.11) 
	 ! RETURN
	 ! END
!------------------------------------



      SUBROUTINE cal_Schmidt	  ! Calculating the Schmidt number
      PARAMETER(NY=403,NX=1)
	  DOUBLE PRECISION HF,HB,ZOF,ZOB,ZF,DX,ZB,Y,RoF,RoB,UB,UF,ViscosC,difusivity,DY,WF,WB,Akapa,PF,PB,DZF,DZB,UOLD1
	  DOUBLE PRECISION DUDY,DUDX,aLength,DVDY,DVDX,gent,XF,eshmitn,DRoDY,rich,g,teta
	  integer optionkey
	  COMMON/C2379/UOLD1(NY),ViscosC(NY),difusivity(NY)
	  COMMON/C2356789/RoF(NY),PF(NY),WF(NY)
	  COMMON/C23789/UB(NY),RoB(NY),PB(NY),ZB(NY),Y(NY),DY(NY),WB(NY)
	  COMMON/C1236789/DX(NX)
	  COMMON/C23456789/NNODEJ,HF,ZOF,ZF(NY),DZF,DZB
  	  COMMON/C2346789/ UF(NY)
      COMMON/C1234789R/HB,ZOB
	  COMMON/C1234789I/NCELLJ
	  COMMON/C1235/g,teta
	  COMMON/C139/Akapa
	  common/ids/Idconvdif,idfull,idnik,Idprofile
	  common/turtec/DUDY(NY),DUDX(NY),aLength(NY),DVDY(NY),DVDX(NY),gent
	  COMMON/C3678/XF,NSECTION
	  common/optiout/idoptiout,optionkey
	  common/eshmitno/eshmitn(NY),DRoDY(NY),rich(NY)
	 
!*--- calculation of DUDX
	   DO J=2,NNODEJ-2
       DUDX(J)=((HF*UF(J)-HB*UB(J))*dY(J)-((HF-HB))*((Y(J)+Y(J+1))/2.*0.25*(UF(J)+UF(J+1)+UB(J)+UB(J+1))-(Y(J)+Y(J-1))/2.*0.25*(UF(J)+UF(J-1)+UB(J)+UB(J-1))))/(0.5*(HF+HB)*dx(1)*dY(J))
!			   12				  2		  23	 3	  2	23			 3	  3							  3	3			3	 3							 321 1	  					 1
	   END DO
	   DUDX(2)=((HF*UF(2)-HB*UB(2))*dY(2)-(HF-HB)*(Y(2)+Y(3))/2.*(UF(2)+UF(3)+UB(2)+UB(3))/4.0)/(0.5*(HF+HB)*dx(1)*dY(2))
!			   12				  2		  2		2 2 		2	 2 	 					 2    1	1    2	   2			1
	   DUDX(NNODEJ-1)=((HF*UF(NNODEJ-1)-HB*UB(NNODEJ))*dY(NNODEJ-1)-(HF-HB)*(Y(NNODEJ)*(UF(NNODEJ-1)+UB(NNODEJ-1))/2.0-(Y(NNODEJ-1)+Y(NNODEJ-2))/2.*(UF(NNODEJ-1)+UF(NNODEJ-2)+UB(NNODEJ-1)+UB(NNODEJ-2))/4.0))/(0.5*(HF+HB)*dx(1)*dY(NNODEJ-1))
!					  12							   2			  2		2 3									   3	 3						 3	  3													  3	   21
      !end of DUDX calculation

!*--- calculation of DUDY
	   DO J=3,NNODEJ-2
       DUDY(J)=0.5*(UF(J+1)+UB(J+1)-UF(J-1)-UB(J-1))/(HF+HB)/DY(J)
	   END DO
       DUDY(2)=0.5*(UF(3)+UB(3)+UF(2)+UB(2))/(HF+HB)/DY(2)
	   DUDY(NNODEJ-1)=(0.5*(UF(NNODEJ-1)+UB(NNODEJ))-0.25*(UF(NNODEJ-1)&
	&  +UB(NNODEJ-1)+UF(NNODEJ-2)+UB(NNODEJ-2)))/(0.5*(HF+HB))/DY(NNODEJ-1)
	  !end of DUDY calculation

!*--- calculation of DVDX
       DO J=2,NNODEJ-2
	   if (XF.le.DX(1)) then
	   DVDX(J)=0.0d0
	   else
	   DVDX(J)=((HF*(WF(J)+WF(J-1))*0.5-HB*(WB(J)+WB(J-1))*0.5)*dY(J)-(HF-HB)*(((Y(J)+Y(J+1))/2.*WF(J))-((Y(J)+Y(J-1))/2.*WF(J-1))))/(0.5*(HF+HB)*dx(1)*dY(J))
	   end if
	   enddo

	   DVDX(NNODEJ-1)=((HF*(WF(NNODEJ-1)+WF(NNODEJ-2))*0.5-HB*(WB(NNODEJ-1)+WB(NNODEJ-2))*0.5)*dY(NNODEJ-1)-(HF-HB)*(Y(NNODEJ)*WF(NNODEJ-1)-0.5*(Y(NNODEJ-1)+Y(NNODEJ-2))*WF(NNODEJ-2)))/(0.5*(HF+HB)*dx(1)*dY(NNODEJ-1))
!					  12   3						 3		  3							3	 2				2	  2	2						 	3                       3             21 
	   !end of DVDX calculation

5555  format(f12.3,2x,I3,2x,d15.5,2x,d15.5,2x,d15.5)
!*--- calculation of DVDY
	  DO J=2,NNODEJ-1
	  DVDY(J)=2.0*(WF(J)-WF(J-1))/(HF+HB)/DY(J)
	  END DO
	  !end of DVDX calculation

	  !*--- calculation of DRoDY
	   DO J=3,NNODEJ-1
       DRoDY(J)=0.5*(RoF(J+1)+RoB(J+1)-RoF(J-1)-RoB(J-1))/(HF+HB)/DY(J)
	   END DO
       DRoDY(2)=(0.25*(RoF(3)+RoB(3)+RoF(2)+RoB(2))-0.5*(RoF(1)+RoB(1)))/2/(HF+HB)/DY(2)
	  !end of DRoDY calculation

	   
 	  DO J=2,NNODEJ-1


	  rich(J)=-g*DRoDY(J)/(0.5*(RoF(J)+RoB(J))*((2*(DUDX(J)**2))+(DUDY(J)**2)))
!						  1    2             2 23              3 3          321  					 

	  if (rich(J).le.-0.00001.or.rich(J).ge.0.15) then 	  
	  end if

	  if (rich(J).gt.0.145) then
	  temp=rich(J) 	  
      rich(J)=0.145
	  eshmitn(J)=0.7*((1.0-rich(J))**2)/(1.0-rich(J)/0.15) 
	  rich(J)=temp 	  
	  else
	  eshmitn(J)=0.7*((1.0-rich(J))**2)/(1.0-rich(J)/0.15) 	  
	  end if

	  


	  if (eshmitn(J).lt.0.1.or.eshmitn(J).gt.2) then	  

	  end if

	  if (eshmitn(J).gt.1.0) then	  
	  eshmitn(J)=1.00d0
	  end if
	  

	  ENddo
	  eshmitn(1)=1.0
	  eshmitn(NNODEJ)=1.0
	  return
		
	END
	 




	
	


