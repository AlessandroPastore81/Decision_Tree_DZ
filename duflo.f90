program DZ
  implicit none

  integer :: Z,N
  double precision :: a(33),Y

  call  setpara0(a)


  write(*,*)' Write the neutron/proton number'
  read(*,*)N,Z

  call  EMASSDZ (a,Z,N,Y) 

  write(*,*)
  write(*,*)' Binding energy [MeV]=',Y



end program DZ






subroutine setpara0(par_inout)
  implicit none

  double precision :: par_inout(33)

  !    original DZ parameters :

  !    data par_inout/0.7043,17.7418,16.2562,37.5562,53.9017,0.4711,2.1307, &
  !       0.0210,40.5356,6.0632/


  par_inout(1)=9.0914 ! FM+
  par_inout(2)=6.3355 ! fm+
  par_inout(3)=4.5791 ! FS+
  par_inout(4)=19.8946 !fs+
  par_inout(5)= 1.7325 ! FS-
  par_inout(6)= 7.5247 ! fs-
  !

  par_inout(7)= -7.1953 !FC+
  par_inout(8)= -39.9787 ! fc+
  par_inout(9)= -0.3976  ! PM+
  par_inout(10)=  0.8131  ! pm+
  par_inout(11)= -0.7435  ! PS+
  par_inout(12)= -3.7291  ! ps+

  par_inout(13)=  -0.1305 ! PS-
  par_inout(14)=  -0.6387 ! ps-
  par_inout(15)=   0.4534 ! S3
  par_inout(16)=   2.0605 ! s3
  par_inout(17)=   0.3449 ! SQ-
  par_inout(18)=   1.4727 ! sq-

  par_inout(19)=   -1.0433 ! D3
  par_inout(20)=   0.0000  ! d3
  par_inout(21)=   5.2495  ! QQ+
  par_inout(22)=   0.0000  ! qq+
  par_inout(23)= -32.1007  ! D0
  par_inout(24)=-151.1164  ! d0

  par_inout(25)=  -4.6103 ! QQ-
  par_inout(26)= -32.4238 ! qq-
  par_inout(27)= 37.3226 ! TT segno -
  par_inout(28)= 52.1673 ! tt
  par_inout(29)=   0.9597 ! SS
  par_inout(30)=   3.0024 ! ss

  par_inout(31)=   0.6977 ! C
  par_inout(32)=   6.0390 ! P0
  par_inout(33)=  17.7960 ! P1



  return
end subroutine setpara0


!
!!
!

SUBROUTINE EMASSDZ (a,Z,X,Y)  !a parameters Z protons, X neutrons , Y Bind.Energy

  implicit none

  integer :: Z,X
  integer :: NN(2),Imax,MaxP,NDEF
  integer :: Ju,JUP(2),JUD(2),nz,nx,KK
  integer :: J,I,NOC(18,2),M,L
  integer :: NCUM,i2,ID,IP,MOC,Z2
  integer :: II,IPL,JF,N,MSS,MO,K
  !
  double precision :: A(33),DYDA(33),FYDA(33),op(2,3,2),fyd0(33)
  double precision :: onps(2),ONP(0:8,2,2),OT(0:8,2,2),OEI(2)
  double precision :: OP2(2),dei(2)
  double precision :: YM(2),OP1(2),n2(2),shell(2),sshell(2)
  double precision :: Y,vmj,vm2,T,V,r,s,Rc,Ra
  double precision :: vxz,vmr,uxz,txz,S3,alfa
  double precision :: FACN,FACZ,FAC,QN,OTX
  double precision :: rxz,QZ,qqp,qqm,opxx,qq0,qq1
  double precision :: DI1N,DI1Z,degi,degr,dnnb,DZZB

  !  c Data=1751 RMS= 0.330  (march 95)
  !  cFM+*   9.09 fm+*   6.34 FS+*   4.58 fs+*  19.89 FS-*   1.73 fs-*   7.52
  !  cFC+*  -7.20 fc+* -39.98 PM+*  -0.40 pm+*   0.81 PS+*  -0.74 ps+*  -3.73
  !  cPS-*  -0.13 ps-*  -0.64 S3 *   0.45 s3 *   2.06 SQ-*   0.34 sq-*   1.47
  !  cD3 *  -1.04 d3     0.00 QQ+*   5.25 qq+    0.00  D0* -32.10  d0*-151.12
  !  cQQ-*  -4.61 qq-* -32.42  TT* -37.32  tt* -52.17  SS*   0.96  ss*   3.00
  !  c C *   0.70 P0 *   6.04 P1 *  17.80
  !  c                    ----------------------

  IMAX=18
  MO=2
  MAXP=8
  NN(1)=X
  NN(2)=Z
  nx=x
  nz=z
  V=(X+Z)*1.d0
  T=ABS(X*1.d0-Z*1.d0)
  r=v**(1.d0/3.d0)
  s=r*r
  Rc=R*(1.d0-0.25d0*(t/v)**2)
  Ra=(rc*rc)/r

  !C******************
  DO NDEF=1,2               !  NDEF=1 spherical  NDEF=2  deformed
     ym(ndef)=0.d0
     !   NDEF= 1->SPH// 2->DEF//
     JU=0
     JUP(1)=0
     JUP(2)=0
     JUD(1)=0
     JUD(2)=0
     !C-----------
     IF(NDEF.EQ.2.and.nz.gt.50)JU=4
     DO KK=1,33
        FYDA(KK)=0.d0
        FYD0(KK)=0.d0
        dyda(KK)=0.d0
     END DO
     DO J=1,2
        DO I=1,IMAX
           NOC(I,J)=0
        ENDDO
        DO K=0,MAXP
           DO M=1,MO
              ONP(K,M,J)=0.d0
           ENDDO
        ENDDO
     ENDDO
     DO J=1,2                  !LOOP OVER FLUIDS (1=N,2=P)
        NCUM=0
        I=0
20      I=I+1
        i2=(I/2)*2
        if(i2.ne.i)then
           ID=I+1
           if(ncum.lt.nn(j))sshell(J)=1.d0     !sscouche  J
        else
           ID=I*(I-2)/4
           if(ncum.lt.nn(j))sshell(J)=2.d0     ! SSC R
        endif
        NCUM=NCUM+ID
        IF(NCUM.LE.NN(J))THEN
           NOC(I,J)=ID
           GO TO 20
        ENDIF
        shell(J)=I        !N0 de sscouche  (SSc 2-> 0 nucl)
        IP=(I-1)/2        !N0 couche HO
        MOC=NN(J)-NCUM+ID
        IF(NDEF.EQ.2) THEN
           if(i2.ne.i)then
              JUD(J)=MAX (0,JU-MOC)
              JUP(J)=0
           Else
              JUP(J)=MIN (JU,MOC)
              JUP(J)=JU
              JUD(J)=0
           ENDIF
        ENDIF
        NOC(I,J)=MOC-JUP(J)+JUD(J)
        NOC(I+1,J)=JUP(j)
        NOC(I-1,J)=NOC(I-1,J)-JUD(J)
        if(i2.ne.i)then
           OEI(J)=MOC+IP*(IP-1)-JU
           DEI(J)=IP*(IP+1)+2
        ELSE
           OEI(J)=MOC-JU
           DEI(J)=(IP+1)*(IP+2)+2
        ENDIF
        ! HERE,DEGENERACIES AND NUMBER OF ACTIVE PARTICLES  FOR  EI.
        IPL=0
        vmr=0.
        vmj=0.
        DO II=1,IMAX
           onps(j)=0.
           IP=(II-1)/2
           DEGI=(IP+1)*(IP+2)
           FAC=1.d0/SQRT(DEGI)
           IF(IP.NE.IPL)IPL=IPL+1
           if((2*IP+1).eq.II)then
              VM2=(0.5d0*IP)/(IP+1)
              degr=ip*(ip-1)
              if(ip.gt.2)then
                 vmr=(0.5d0*(ip-1))/ip
                 vmj=-1.d0/ip
                 if(noc(ii,j).le.degr)onps(j)=noc(ii,j)*vmr
                 if(noc(ii,j).gt.degr)onps(j)=degr*vmr+(noc(ii,j)-degr)*vmj
              endif
           endif
           if((2*IP+1).ne.II)VM2=-1.d0/(IP+1)      !SSc. j
           ONP(IPL,2,J)=ONP(IPL,2,J)+NOC(II,J)*VM2
           ONP(IPL,1,J)=ONP(IPL,1,J)+NOC(II,J)*FAC
           fyd0(29)=fyd0(29)+onps(j)*(onp(ipl,1,j)+onp(ipl,2,j))
        ENDDO
     ENDDO           !END OF LOOP OVER FLUIDS
     !**************

     IF(NDEF.EQ.2) THEN
        ALFA=0.d0
     ELSE
        ALFA=0.5d0
     ENDIF
     FACN=DEI(1)**ALFA
     FACZ=DEI(2)**ALFA
     DNNB=OEI(1)*(DEI(1)-OEI(1))/DEI(1)
     DZZB=OEI(2)*(DEI(2)-OEI(2))/DEI(2)
     QN=DNNB*FACN/SQRT(DEI(1))
     QZ=DZZB*FACZ/SQRT(DEI(2))
     DI1N=DNNB*(2*OEI(1)-DEI(1))*FACN*FACN/(DEI(1))
     DI1Z=DZZB*(2*OEI(2)-DEI(2))*FACZ*FACZ/(DEI(2))
     S3=DI1Z+DI1N
     QQ0=(QN+QZ)*(QN+QZ)
     QQ1=(QN-QZ)*(QN-QZ)
     qqp=qq0+qq1
     qqm=qq0-qq1
     !c------------
     DO M=1,MO
        DO I=0,MAXP
           OT(I,M,1)=ONP(I,M,1)+ONP(I,M,2)
           OT(I,M,2)=ONP(I,M,1)-ONP(I,M,2)
        ENDDO
     ENDDO
     !c------------
     DO L=1,2
        DO M=1,3
           DO J=1,2
              OP(L,M,J)=0.d0
           ENDDO
        ENDDO
     ENDDO
     DO I=0,MAXP
        DEGI=(I+1)*(I+2)
        FAC=SQRT(DEGI)
        DO N=1,2
           DO M=1,MO
              OP(1,M,N)=OP(1,M,N)+OT(I,M,N)
              OTX=OT(I,M,N)*OT(I,M,N)*FAC
              OP(2,M,N)=OP(2,M,N)+OTX
           ENDDO
        ENDDO
     ENDDO
     DO N=1,2
        OP1(N)=0.d0
        OP2(N)=0.d0
        DO M=1,MO
           OPXX=OP(1,M,N)*OP(1,M,N)
           OP(1,M,N)=OPXX
        ENDDO
     ENDDO
     DO N=1,2
        DO I=0,MAXP
           DEGI=(I+1)*(I+2)
           FAC=sqrt(degi)
           OP1(N)=OP1(N)+OT(I,1,N)/FAC
           OP2(N)=OP2(N)+OT(I,2,N)/FAC
        ENDDO
     ENDDO
     DO N=1,2
        OP(1,3,N)=OP1(N)*OP2(N)
     ENDDO
     K=-1
     DO L=1,2
        DO M=1,3
           DO N=1,2
              K=K+2
              FYD0(K)=OP(L,M,N)
           ENDDO
        ENDDO
     ENDDO
     !C-----------
     DO JF=1,17,4
        FYD0(JF)=FYD0(JF)+FYD0(JF+2)
        FYD0(JF+2)=FYD0(JF)-2.*FYD0(JF+2)
     ENDDO
     !C-----------
     fyda(1)=fyd0(1)
     fyda(3)=fyd0(5)
     fyda(5)=fyd0(7)
     fyda(7)=fyd0(9)
     fyda(9)=fyd0(13)
     fyda(11)=fyd0(17)
     fyda(13)=fyd0(19)
     FYDA(27)= -T*(T+2)/s
     fyda(29)=fyd0(29)
     !C-----------
     IF(NDEF.EQ.1)then
        FYDA(15)=S3
        FYDA(17)=QQm
     ELSE
        FYDA(19)=S3
        FYDA(21)=QQp
        FYDA(23)=16.d0-qqm
        FYDA(25)=QQm
     ENDIF
     !C-----------
     DO MSS=1,29,2
        DYDA(MSS)=FYDA(MSS)/Ra
        DYDA(MSS+1)=-DYDA(MSS)/Ra
        ym(ndef)=ym(ndef)+dyda(mss)*a(mss)+dyda(mss+1)*a(mss+1)
     ENDDO
     !c---------
     Z2=Z*(Z-1)
     DYDA(31)=(-Z2+0.76d0*Z2**(2.d0/3.d0))/Rc      !Coulomb
     !c---------
     rxz=1.d0/Ra                           !Pairing
     vxz=1.d0/v
     txz=rxz*(t/v)
     uxz=rxz-txz
     do l=1,2
        n2(l)=2*(nn(l)/2)
        dyda(32)=-rxz+dyda(32)
        if(sshell(l).eq.2.) dyda(33)=vxz+dyda(33)
        !effet de couche en 1/A
        if(n2(l).eq.nn(l))dyda(32)=uxz+dyda(32)
     enddo
     j=2                     !  Z>N
     if(nn(1).ge.nn(2))j=1   !  N>ou=Z
     k=3-j
     if(n2(j).eq.nn(j).and.n2(k).ne.nn(k))dyda(32)=-txz+dyda(32)
     !c********
     DO MSS=31,33
        ym(ndef)=ym(ndef)+dyda(mss)*a(mss)
     ENDDO
  ENDDO  !End of spherical and deformed calculations
  Y=ym(1)
  if(z.gt.50.and.ym(2).gt.ym(1)) Y=ym(2)

  write(*,*)' Various terms of the formula'
  write(*,*)' Parameter [MeV]    Term'
  do l=1,33
     write(*,*)a(l),dyda(l)
  enddo

  !epair=dyda(32)*a(32)+dyda(33)*a(33)
END SUBROUTINE EMASSDZ

