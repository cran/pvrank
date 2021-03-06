!----------------------------------------------
subroutine CRS(lam,mu2n,mu4n,Eval,n,zran,icoef)
        implicit none
        integer, parameter :: MyR = selected_real_kind(15,100)
        real (kind=MyR), intent(in) :: mu2n,mu4n,n
        integer, intent(in) :: icoef
        real (kind=MyR), intent(in out) :: zran(30,2)
        real (kind=MyR), intent(out)    :: Lam(2)
        real (kind=MyR), intent(out)    :: Eval
        real (kind=MyR) :: Half=0.5_MyR,One=1._MyR,Two=2._MyR,unif
        real (kind=MyR) :: Rg(2,3),Por(30,3),bsas(2)
        real (kind=MyR) :: ws1,ws2,ws3,ws4,ws5,ws6,nm1
        real (kind=MyR) :: tol1=0.000000000001_MyR, Crit=0.00000000000000001_MyR
        real (kind=MyR) :: fmax,fmin,fbar,fnew,test,large=1.0d99
        integer :: n3,n4,bm,jm,imas,nver,imis,att,pass,tryag,em(3)
        integer :: count,maxev
        n3=2;bsas(1)=0.;bsas(2)=0.;imas=0;nm1=n-1.0_MyR
        n4=n3+1;count=1
        nver=idint(10.0_MyR*n4)       
        maxev=idint(125000.0_Myr*sqrt(Log10(n)))
        Rg(1,1)=1.5_MyR;Rg(1,2)=3.5_MyR
        CALL seed_init()
        select case (icoef)
                case (1)                                         
                        Rg(2,1)=0.3_MyR*n;Rg(2,2)=0.6_MyR*n
                case (2)
                        Rg(2,1)=0.3_MyR*n; Rg(2,2)=0.9_MyR*n              
                case (3)                                         
                        Rg(2,1)=0.7_MyR*n;Rg(2,2)=1.2_MyR*n
                case (4)
                        Rg(2,1)=0.2_MyR*n;Rg(2,2)=0.7_MyR*n
                case (5)
                        Rg(2,1)=0.2_MyR*n;Rg(2,2)=0.6_MyR*n
                case (6)
                        Rg(2,1)=0.2_MyR*n;Rg(2,2)=0.6_MyR*n
        end select
        Rg(2,3)=Rg(2,2)-Rg(2,1);Rg(1,3)=Rg(1,2)-Rg(1,1)
        do bm=1,nver
                Zran(bm,1)=Rg(1,1)+Zran(bm,1)*Rg(1,3)
                Zran(bm,2)=Rg(2,1)+Zran(bm,2)*Rg(2,3) 
                Lam(1)=zran(bm,1);Lam(2)=zran(bm,2)
                call ObjFun(Lam,mu2n,mu4n,Ws6)
                Por(bm,1)=Lam(1);Por(bm,2)=Lam(2);Por(bm,3)=Ws6
        enddo
        fmax=-large;fmin=large
        do bm=1,nver
                if (Por(bm,n4)>fmax) then  
                        fmax=Por(bm,n4)
                        imas=bm
                endif
                if (Por(bm,n4)<fmin) then  
                        fmin=Por(bm,n4)
                        imis=bm
                endif
        enddo
91      Test=fmax-fmin
        if (Test<Crit) then
                do jm=1,n3
                        lam(jm)=Por(imis,jm)
                enddo
                Eval=Por(imis,n4)
                return
        endif 
        if (count>maxev) then
                do jm=1,n3
                        lam(jm)=Por(imis,jm)
                enddo
                !write(*,12)
                !12 format("CRS exceed maximum number of evaluations")
                return
        endif
        att=0
93      att=att+1
        if (att>4*n4) goto 95
        call sample(n4,nver,imis,em) 
        do jm=1,n3
                bsas(jm)=por(imis,jm)
        enddo
        do      bm=1,n3
                do jm=1,n3
                        bsas(jm)=bsas(jm)+Por(Em(bm),jm)
                enddo
        enddo
        do jm=1,n3
                bsas(jm)=bsas(jm)/n4
        enddo
        pass=0
        do jm=1,n3      
                lam(jm)=Two*bsas(jm)-Por(Em(n4),jm)
                if(lam(jm)<Rg(jm,1).or.lam(jm)>Rg(jm,2)) goto 93
        enddo
        call ObjFun(lam,mu2n,mu4n,fbar) 
        count=count+1
        if (fbar<fmax .and. fbar>fmin) then
                do jm=1,n3
                        Por(imas,jm)=lam(jm)
                enddo
                Por(imas,n4)=fbar
                fmax=-large
                do bm=1,nver
                        if(Por(bm,n4)>fmax) then 
                                fmax=Por(bm,n4)
                                imas=bm
                        endif
                enddo
                goto 93
        endif
        if(fbar<fmin) then 
                do jm=1,n3
                        Por(imas,jm)=por(imis,jm)
                        Por(imis,jm)=lam(jm)
                enddo
                por(imas,n4)=por(imis,n4)
                por(imis,n4)=fbar
                fmax=-large
                fmin=large
                do bm=1,nver
                        if(Por(bm,n4)>fmax) then 
                                fmax=Por(bm,n4)
                                imas=bm
                        endif                   
                        if (Por(bm,n4)<fmin) then  
                                fmin=Por(bm,n4)
                                imis=bm
                        endif
                enddo
                goto 91
        endif
95      Tryag=0
96      Tryag=Tryag+1
        if(Tryag>4*n4) goto 91
        call sample(n3,nver,imis,em)  
        do jm=1,n3
                Ws1=(Por(Em(1),jm)**2-Por(Em(2),jm)**2)*Por(imis,n4)
                Ws2=(Por(Em(2),jm)**2-Por(imis,jm)**2)*Por(Em(1),n4)
                Ws3=(Por(imis,jm)**2-Por(Em(1),jm)**2)*Por(Em(2),n4)
                Ws4=Ws1+Ws2+Ws3
                Ws1=(Por(Em(1),jm)-Por(Em(2),jm))*Por(imis,n4)
                Ws2=(Por(Em(2),jm)-Por(imis,jm))*Por(Em(1),n4)
                Ws3=(Por(imis,jm)-Por(Em(1),jm))*Por(Em(2),n4)
                Ws5=Ws1+Ws2+Ws3
                if (abs(Ws5)>Tol1) then
                        lam(jm)=Half*Ws4/Ws5 
                else
                        CALL RANDOM_NUMBER(Unif)
                        lam(jm)=unif*Por(imis,jm)+(One-Unif)*Por(imas,jm)
                endif         
                if (lam(jm)<Rg(jm,1) .or. lam(jm)>Rg(jm,2)) goto 95
        enddo
        call ObjFun(lam,mu2n,mu4n,fnew)
        count=count+1
        if (fnew<fmax .and. fnew>fmin) then 
                do jm=1,n3
                        Por(imas,jm)=lam(jm)
                enddo
                Por(imas,n4)=fnew
                fmax=-large
                do bm=1,nver
                        if(Por(bm,n4)>fmax) then 
                                fmax=Por(bm,n4)
                                imas=bm
                        endif
                enddo
                goto 96
        endif
        if(fnew<fmin) then 
                do jm=1,n3
                        Por(imas,jm)=por(imis,jm)
                        Por(imis,jm)=lam(jm)
                enddo
                por(imas,n4)=por(imis,n4)
                por(imis,n4)=fnew
                fmax=-large
                fmin=large
                do bm=1,nver
                        if(Por(bm,n4)>fmax) then 
                                fmax=Por(bm,n4)
                                imas=bm
                        endif                   
                        if (Por(bm,n4)<fmin) then  
                                fmin=Por(bm,n4)
                                imis=bm
                        endif
                enddo
        endif
        goto 96
end subroutine
subroutine seed_init 
       integer :: seed_size
       integer,allocatable :: seed(:)
       call random_seed(size=seed_size) ! find out size of seed
       allocate(seed(seed_size))
       call random_seed(get=seed) ! get system generated seed
       seed=3141569 
       call random_seed(put=seed) ! set current seed
       call random_seed(get=seed) ! get current seed
       deallocate(seed)           ! safe
end subroutine seed_init
!-.-.-.-.-.-.-.-.-.-.-.-
subroutine logamma(zx,Loz)
! Compute the logarithm of gamma function
        implicit none
        integer, parameter :: MyR = selected_real_kind(15,100)
        real (kind=MyR), intent(in) :: zx
        real (kind=MyR), intent(out) :: Loz
        real (kind=MyR) :: Xlg,X2,X1,Y,XLGE,ALR2PI,Onefiv,Half
        real (kind=MyR) :: Zero, One, Four,Twelve
        real (kind=MyR) :: R11,R12,R13,R14,R15,R16,R17,R18,R19,R21,R22,R23,R24,R25,R26,R27,R28
        real (kind=MyR) :: R29,R31,R32,R33,R34, R35,R36,R37,R38,R39,R41,R42,R43,R44,R45
        data R11,R12,R13/-2.66685511495_MyR,-24.4387534237_MyR,-21.9698958928_MyR/
        data R14,R15,R16/11.1667541262_MyR,3.13060547623_MyR,0.607771387771_MyR/
        data R17,R18,R19/11.9400905721_MyR,31.4690115749_MyR,15.2346874070_MyR/
        data R21,R22,R23/-78.3359299449_MyR,-142.046296688_MyR,137.51941641_MyR/
        data R24,R25,R26/78.6994924154_MyR,4.16438922228_MyR,47.0668766060_MyR/
        data R27,R28,R29/313.399215894_MyR,263.505074721_MyR,43.340002251_MyR/
        data R31,R32,R33/-212159.572323_MyR,230661.510616_MyR,27464.7644705_MyR/
        data R34,R35,R36/-40262.1119975_MyR,-2296.60729780_MyR, -116328.495004_MyR/
        data R37,R38,R39/-146025.937511_MyR,-24235.7409629_MyR,-570.691009324_MyR/
        data R41,R42,R43/0.279195317918525_MyR,0.4917317610505968_MyR,0.0692910599291889_MyR/
        data R44,R45/3.350343815022304_MyR,6.012459259764103_MyR/
        data XLGE,ALR2PI,Onefiv,Half/976._MyR,0.918938533204673_MyR,1.5_MyR,0.5_MyR/
        data Zero,One,Four,Twelve/0.0_MyR,1.0_MyR,4.0_MyR,12.0_MyR/
        Xlg=Zx
        Loz=Zero
        if (Xlg<Onefiv) then
                if(Xlg<Half) then
                        Loz=-log(Xlg)
                        Y=Xlg+One
                        if (Y==One) return
                else 
                        Loz=Zero
                        Y=Xlg
                        Xlg=(Xlg-Half)-Half
                endif
                Loz=Loz+Xlg*((((R15*Y+R14)*Y+R13)*Y+R12)*Y+R11)/((((Y+R19)*Y+R18)*Y+R17)*Y+R16)
                return
        endif
        if (Xlg<Four) then
                Y=(Xlg-One)-One
                Loz=Y*((((R25*Xlg+R24)*Xlg+R23)*Xlg+R22)*Xlg+R21)/((((Xlg+R29)*Xlg+R28)*Xlg+R27)*Xlg+R26)
                return
        endif
        if (Xlg<Twelve) then
                Loz=((((R35*Xlg+R34)*Xlg+R33)*Xlg+R32)*Xlg+R31)/((((Xlg+R39)*Xlg+R38)*Xlg+R37)*Xlg+R36)
                return
        endif
        Y=LOG(Xlg)
        Loz=Xlg*(Y-One)-Half*Y+ALR2PI
        if (Xlg>XLGE) return
        X1=One/Xlg
        X2=X1*X1
        Loz=Loz+X1*((R43*X2+R42)*X2+R41)/((X2+R45)*X2+R44)
        return
end subroutine
!-----------------------------------------------------------------
subroutine beta(ax,bx,betx)
! Emulator of the complete beta function
        implicit none
        integer, parameter :: MyR = selected_real_kind(15,100)
        real (kind=MyR), intent(in) :: ax,bx
        real (kind=MyR), intent(out) :: betx
        real (kind=MyR) ::  Pl1,Pl2,Pl3
        call logamma(ax,Pl1)
        call logamma(bx,Pl2)
        call logamma(ax+bx,Pl3)
        betx=Pl1+Pl2-Pl3
        betx=exp(betx)
 end subroutine
!-----------------------------------------------------------------
subroutine Sample(nsam,npop,imis,em)
! Generate a random sample without replacement 
        implicit none
        integer, parameter :: MyR = selected_real_kind(15,100)
        integer, intent(in) :: nsam,npop,imis
        integer, intent(out) :: em(3)
        real (kind=MyR) :: unif
        integer :: i,ifo,inu
101 CALL RANDOM_NUMBER(unif)
        inu=1+floor(unif*npop)
        if (inu==imis) goto 101
        ifo=1
        em(ifo)=inu
102 CALL RANDOM_NUMBER(unif)
        inu=1+floor(unif*npop)
        do i=1,ifo      
                if (inu==em(i)) goto 102
        enddo
        ifo=ifo+1
        em(ifo)=inu
        if (ifo<nsam) goto 102
end subroutine
!-------------------------------------
subroutine ObjFun(lam,mu2n,mu4n,Obf)
        implicit none
        integer, parameter :: MyR = selected_real_kind(15,100)
        real (kind=MyR), intent(in) :: lam(2),mu2n,mu4n
        real (kind=MyR), intent(out) :: Obf
        real (kind=MyR) :: One=1._MyR,Fiv=5._MyR,Thr=3._MyR
        real (kind=MyR) :: ws1,ws2,ws3,ws4,ws5,mu2a,mu4a,Fp(2)
        ws2=lam(2)+One
        ws1=Thr/lam(1)
        Call beta(ws1,ws2,ws3)
        ws1=One/lam(1)
        Call beta(ws1,ws2,ws4)
        mu2a=ws3/ws4
        ws1=Fiv/lam(1)
        Call beta(ws1,ws2,ws5)
        mu4a=ws5/ws4
        Fp(1)=mu2a-mu2n
        Fp(2)=mu4a-mu4n
        Obf=abs(Fp(1))+abs(Fp(2))
        return
end subroutine
!*********************************