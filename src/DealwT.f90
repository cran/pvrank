subroutine gPerm(n,neq,nep,kdq,kdp,Tbq,Tbp,Taq,Tap,iprp,iprq,Teq,Tep,Tfq,Tfp,Tgq,Tgp,cp,cq,Hilo,rc,comp,index)
  	implicit none
  	integer, parameter  :: MyR = selected_real_kind(15)
  	integer, parameter :: MyI = selected_int_kind (8)  	
  	integer, intent(in) :: n,index
	integer(kind=MyI), intent(in) :: kdp,kdq
	integer(kind=MyI), intent(in) :: Taq(n),Tbq(n,n),Tap(n),Tbp(n,n),Tgq(n)
	integer(kind=MyI), intent(in) :: Teq(n,n),Tep(n,n),Tfq(n),Tfp(n),Tgp(n)
  	integer(kind=MyI), intent(in) :: nep,neq,iprq(n),iprp(n)
  	logical, intent(in) :: comp
  	real(kind=MyR), intent(in out) :: cq(n),cp(n)
  	real(kind=MyR), intent(in)     :: Hilo(2)
  	real(kind=MyR), intent(out)    :: rc
  	real(kind=MyR)   :: w1,w2,ra,rad,thet,sig,Cg,Cs,Ck,repl
	integer(kind=MyI) :: m,i,j,kj,kn,jj,kk,jk
	integer(kind=MyI) :: T2(2,2),T3(6,3),T4(24,4),T5(120,5),T6(720,6),T7(5040,7),T8(40320,8),T9(362880,9)
	integer(kind=MyI) :: ivec(n)
	integer(kind=MyI) :: jsubq,jsubp,jp,jq
	logical qfor,qind	
    call Base(T2,T3,T4,T5,T6,T7,T8,T9)
    qfor = .false.; qind = .true.
    kn=mod(n,2);Cg=n**2-kn;Cs=6./(n**3-n);Ck=2./(n**2-n)
  	ra=0.5*(Hilo(1)+Hilo(2));rad=Hilo(1)-Hilo(2);w1=0.;w2=0.
  	repl=neq*nep
  	do jq = 1, neq
  		if (neq.gt.1) then
    		jsubq = jq
    		call simdo(qind,qfor,iprq,kdq,jsubq,ivec)
    		do i=1,kdq
    	    	kk=ivec(i);kj=Tfq(i);jk=Tgq(i)
    			select case (kj)
    			case (2)    				
    				do jj=1,2
    					m=T2(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (3)
    				do jj=1,3
    					m=T3(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (4)
    				do jj=1,4
    					m=T4(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (5)
    				do jj=1,5
    					m=T5(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (6)
    				do jj=1,6
    					m=T6(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (7)
    				do jj=1,7
    					m=T7(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (8)
    				do jj=1,8
    					m=T8(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (9)
    				do jj=1,9
    					m=T9(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			end select	
    		enddo
    	endif
    	do jp = 1,nep
    		if (nep.gt.1) then
    			jsubp = jp
    			call simdo(qind,qfor,iprp,kdp,jsubp,ivec)
    			do i=1,kdp
    				kk=ivec(i);kj=Tfp(i);jk=Tgp(i)
    				select case (kj)   				
    				case (2)
    					do jj=1,2
    						m=T2(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (3)
    					do jj=1,3
    						m=T3(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (4)
    					do jj=1,4
    						m=T4(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (5)
    					do jj=1,5
    						m=T5(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (6)
    					do jj=1,6
    						m=T6(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (7)
    					do jj=1,7
    						m=T7(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (8)
    					do jj=1,8
    						m=T8(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (9)
    					do jj=1,9
    						m=T9(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				end select	
    			enddo
    		endif
    		select case (index)   				
    			case (1)
    		    	rc=ra
    				if (comp) then
    					call Spearman(n,Cs,cq,cp,rc)
    				endif
    			case (2)
    		    	rc=ra
    				if(comp) then
    					call Kendall(n,Ck,cq,cp,rc)
    				endif
    			case (3)
    		    	rc=ra
    				if(comp) then 
    					call  Gini(n,Cg,cq,cp,rc)
    				endif
    			case (4)
    		    	rc=ra
    				if (comp) then
    					call r4(n,cq,cp,rc)
    				endif
    		end select
			if (comp) then 
				w1=w1+abs(rc-Hilo(2))
    			w2=w2+abs(rc-Hilo(1))
    		endif
    	enddo
  	enddo
  	if(comp) then
		rc=Hilo(2)/(w1+0.000001)+Hilo(1)/(w2+0.000001)
		rc=rc/(1./(w1+0.000001)+1./(w2+0.000001))
  	else
  		rc=ra
  	endif 
  	return	
end subroutine gPerm
!-------------------------------------------
subroutine wGini(n,q,p,repl,Hilo,index,rc)
  	implicit none
  	integer, parameter :: MyR = selected_real_kind(15)
  	integer, parameter :: MyI = selected_int_kind (8)  
  	integer, intent(in) :: n,index
  	real(kind=MyR), intent(in) :: repl
	real(kind=MyR), intent(in) :: q(n),p(n),Hilo(2)
  	real(kind=MyR), intent(out) :: rc
  	real(kind=MyR) :: cq(n),cp(n)
  	real(kind=MyR) :: ra,rad,w1,w2,thet,chk
	integer(kind=MyI) :: Taq(n),Tbq(n,n),Tcq(n),Tfq(n),Teq(n,n),Tgq(n)
	integer(kind=MyI) :: Tap(n),Tbp(n,n),Tcp(n),Tfp(n),Tep(n,n),Tgp(n)
	integer(kind=MyI) :: T2(2,2),T3(6,3),T4(24,4),T5(120,5),T6(720,6),T7(5040,7),T8(40320,8),T9(362880,9)
	integer(kind=MyI) :: used(n),Tdq(n),Tdp(n)
	integer(kind=MyI) :: kq,jq,jp,neq,nep,nepq,nq
	integer(kind=MyI) :: iprq(n),iprp(n)
  	integer(kind=MyI) :: i,j,kn,ntq,ntp,kdq,kdp
  	real(kind=MyR) 	 :: Cg,Cs,Ck,One,tm
	logical :: qfor,qind,comp
	nq=repl
  	ra=0.5*(Hilo(1)+Hilo(2));rad=Hilo(1)-Hilo(2)
  	comp=.true.;w1=0.;w2=0.
  	if(abs(rad).lt.0.01) then 
  		comp=.false.
  	endif
  	if(.not.comp) then
  		rc=ra
  		return
  	endif
 	call SerTies(n,q,ntq,kdq,Taq,Tbq,Teq,Tfq,Tgq,iprq,neq,cq)
    call SerTies(n,p,ntp,kdp,Tap,Tbp,Tep,Tfp,Tgp,iprp,nep,cp)
  	nepq=neq*nep
  	if(repl.ge.nepq) then
  		if(nepq.gt.0) then
  			call gPerm(n,neq,nep,kdq,kdp,Tbq,Tbp,Taq,Tap,iprp,iprq,Teq,Tep,Tfq,Tfp,Tgq,Tgp,cp,cq,Hilo,rc,comp,index)
  			return
  		else
  			return
  		endif
  	endif
	call Base(T2,T3,T4,T5,T6,T7,T8,T9)
    kn=mod(n,2);Cg=n**2-kn;Cs=6./(n**3-n);Ck=2./(n**2-n)
    do kq=1,nq
		call Oneper(n,kdq,kdp,nep,neq,T2,T3,T4,T5,T6,T7,T8,T9,Teq,Tep,Tbq,Tbp,Tgq,Tgp,Tfq,Tfp,iprq,iprp,cq,cp)
    	select case (index)   				
    		case (1)
    		    rc=ra
    			if(comp) call Spearman(n,Cs,cq,cp,rc)
    		case (2)
    		    rc=ra
    			if(comp) call Kendall(n,Ck,cq,cp,rc)
    		case (3)
    		    rc=ra
    			if(comp) call  Gini(n,Cg,cq,cp,rc)
    		case (4)
    		    rc=ra
    			if(comp) call r4(n,cq,cp,rc)
    	end select
		if (comp) then 
				w1=w1+abs(rc-Hilo(2))
    			w2=w2+abs(rc-Hilo(1))
    	endif
  	enddo
  	if(comp) then
			rc=Hilo(2)/(w1+0.000001)+Hilo(1)/(w2+0.000001)
			rc=rc/(1./(w1+0.000001)+1./(w2+0.000001))
  		else
  			rc=ra
  	endif 
  	return	
end subroutine wGini
!-----------------------------------
subroutine Gua(n,x,iy)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)  
  	integer, intent(in) :: n
	integer(kind=MyI), intent(in out) :: iy(n)
    real(kind=MyR), intent (in out)  :: x(n)
    integer :: i,j,m
    i=1
100	m=0
    do j=1,n
    	if(x(j).le.x(i)) m=m+1
    enddo
    iy(i)=m
    i=i+1
    if(i.le.n) goto 100		
    return
end subroutine Gua
!---------------------------------_________
subroutine Guar(n,arg1,arg2,ord)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
  	integer, intent(in) :: n
	integer(kind=MyI), intent(out)   :: ord(n)
    real(MyR), intent (in out)  :: arg1(n),arg2(n)
    real(MyR)    :: tmp
    real(MyR)    :: qs(n)
    integer(kind=MyI) :: i,j,k,itmp
    integer(kind=MyI) :: iar(n),iy(n),pos(n)
    integer :: m
    logical :: sorted
    do i=1,n
    	ord(i)=i
    enddo
    sorted=.false.;k=0
    do while(.not.sorted)
    	sorted=.true.;k=k+1
    	do j=1,n-k
    		if(arg1(j).gt.arg1(j+1)) then
    			tmp=arg1(j);arg1(j)=arg1(j+1);arg1(j+1)=tmp
    			tmp=arg2(j);arg2(j)=arg2(j+1);arg2(j+1)=tmp
    			itmp=ord(j);ord(j)=ord(j+1);ord(j+1)=itmp
    			sorted=.false.
    		endif
    	enddo
    enddo 
    i=1
    do while(i.lt.n)
    	m=1;iar(1)=i
    	do j=i+1,n
    		if(arg1(j).eq.arg1(i)) then
    			m=m+1;iar(m)=j
    		endif
    	enddo
		if(m.gt.1) then
			do k=1,m
				qs(k)=arg2(iar(k));pos(k)=ord(iar(k));iy(k)=k
			enddo
			call InSrtA(m,qs,iy)
			do k=1,m
				arg2(iar(k))=qs(k);ord(iar(k))=pos(iy(k))
			enddo
		endif
		i=i+m
	enddo	
	return
end subroutine Guar
subroutine GiUls(n,q,p,cqx,cpx,cqy,cpy)
	implicit none
  	integer, parameter :: MyR = selected_real_kind(15)
  	integer, parameter :: MyI = selected_int_kind(8)
  	integer, intent(in) :: n
	real(MyR), intent(in) :: q(n),p(n)
	real(MyR), intent(out):: cqx(n),cpx(n),cqy(n),cpy(n)
	integer(kind=MyI) :: ord(n),rank(n),iy(n)
	integer :: i,prime
	real(kind=MyR) :: arg1(n),arg2(n),xt(n),rky(n),rky1(n)
	real(kind=MyR) :: xrr(n),rx(n)
	do i=1,n
		arg1(i)=p(i);arg2(i)=q(i)
	enddo
	call Guar(n,arg1,arg2,ord)
	do i=1,n
		xt(i)=q(ord(i));rky(i)=i
		arg1(i)=xt(i);arg2(i)=rky(i)
	enddo
	call Guar(n,arg1,arg2,ord)
	do i=1,n
		rky1(i)=rky(ord(i))
	enddo
	do i=1,n
		arg1(i)=q(i);rank(i)=i
	enddo
	call Gua(n,arg1,rank)
	do i=1,n
		arg1(i)=rank(i)
	enddo
	call InSrtA(n,arg1,iy)
	prime=1
	call meanrank(n,arg1,rx,prime)
	do i=1,n
		xrr(i)=n+1 - arg1(rank(i));arg2(i)=p(i)
	enddo
	call Guar(n,arg2,xrr,ord)
	do i=1,n
		xt(i)=q(ord(i));rky(i)=i
		arg1(i)=xt(i);arg2(i)=n+1-i
	enddo
	call Guar(n,arg1,arg2,ord)
	do i=1,n
		cqx(i)=rky1(i);cpx(i)=i;cqy(i)=ord(i);cpy(i)=i
	enddo
	return
end subroutine GiUls
!----------------------------------------
subroutine GGH(n,q,p,Hilo,index)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
  	integer, intent(in) :: n,index
  	real(MyR), intent(out) :: Hilo(2)
  	real(MyR), intent(in) ::q(n),p(n)
  	real(MyR)     :: cqx(n),cpx(n),cqy(n),cpy(n)
  	real(MyR)     :: Cs,Ck,Cg
	integer :: i,kn
    call GiUls(n,q,p,cqx,cpx,cqy,cpy)
    select case (index)
    	case (1) 
    		Cs=6./(n**3-n)
			call Spearman(n,Cs,cqx,cpx,Hilo(2))
			call Spearman(n,Cs,cqy,cpy,Hilo(1))
		case (2) 
			Ck=2./(n**2-n)
			call Kendall(n,Ck,cqx,cpx,Hilo(2))
			call Kendall(n,Ck,cqy,cpy,Hilo(1))
		case (3)
		    kn=mod(n,2);Cg=n**2-kn
			call Gini(n,Cg,cqx,cpx,Hilo(2))
			call Gini(n,Cg,cqy,cpy,Hilo(1))
		case (4) 
			call r4(n,cqx,cpx,Hilo(2))
			call r4(n,cqy,cpy,Hilo(1))
	end select
	return
end subroutine GGH
!--------------------------------------
subroutine Midrank(n,q,p,prime,qx,px)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
	integer, intent(in) :: n,prime
	real (kind=MyR), intent(in)  :: q(n),p(n)
	real (kind=MyR), intent(out) :: qx(n),px(n)
	real (kind=MyR)	 :: cq(n),cp(n)
	integer(kind=MyI) :: iq(n),ip(n)
	integer :: i
	do i=1,n
		qx(i)=q(i);iq(i)=i
	enddo
	call InSrtA(n,qx,iq)
	do i=1,n
		px(i)=p(i);ip(i)=i
	enddo	
	call InSrtA(n,px,ip)
	call meanrank(n,qx,cq,prime)
	call meanrank(n,px,cp,prime)
	do i=1,n
		qx(iq(i))=cq(i);px(ip(i))=cp(i)
	enddo
	return
end subroutine
subroutine Nexper(n,A,mtc,nlast)
		implicit none
		integer, parameter :: MyI = selected_int_kind(8)
		integer, intent(in):: n
		logical, intent(in out) :: mtc
		integer(kind=MyI), intent(in out) :: A(n),nlast
		integer(kind=MyI) :: m,v,nf,j,j1,t,h,b,m1,h1,m2
      	save
      	if(n.le.0) then
        	 mtc=.false.
        	 return
      	endif
 10   if(n.eq.nlast) goto 20
 30   nlast=n
      	m=1;v=1;nf=1
     	 do j=1,n
         	nf=nf*j;A(j)=j
      	enddo
 40   mtc=m.ne.nf
      	return
 20   if(.not.mtc) goto 30
      	goto(70,80) v
 70   t=A(2)
      	A(2)=A(1);A(1)=t;v=2;m=m+1
      	goto 40
 80   h=3
      	m1=m/2
 90   b=mod(m1,h)
      	if(b.ne.0) goto 120
      	m1=m1/h;h=h+1
      	goto 90
 120  m1=n
      	h1=h-1
      	do 160 j=1,h1
        	m2=a(j)-a(h)
       		 if(m2.lt.0) m2=m2+n
       		 if(m2.ge.m1) goto 160
        	m1=m2
        	j1=j
 160    continue
        	t=A(h);A(h)=A(j1);A(j1)=t;v=1;m=m+1
      	return
end subroutine
Subroutine revers(itab, idim)
	implicit none
	integer, parameter :: MyI = selected_int_kind(8)
	integer(kind=MyI), intent(in) :: idim
	integer(kind=MyI), intent(in out):: itab(idim)
	integer(kind=MyI) :: i,k,ik,itemp,iter
	iter = idim / 2
	k = idim + 1
	do 10 i = 1, iter
		itemp = itab(i);ik = k - i
		itab(i) = itab(ik);itab(ik) = itemp
10 	continue
	return
end subroutine revers
!--------------------------------------------------------
subroutine simdo(qind,qfor,iprod,kdim,jsub,ivec)
	implicit none
	integer, parameter :: MyI = selected_int_kind(8)
	integer(kind=MyI), intent(in) :: kdim
	integer(kind=MyI), intent(in) :: iprod(kdim)
	integer(kind=MyI), intent(in out) ::jsub
	integer(kind=MyI), intent(out):: ivec(kdim)
	logical, intent(in) :: qind, qfor
	integer(kind=MyI) :: i,ij,ik,itempv,ifault	
	ifault = 0
	if (.not. qind) goto 12
	if (jsub .le. iprod(kdim)) goto 5
	ifault = 1
	return
5 	itempv = jsub - 1
	ij = kdim - 1
	if (ij.gt.0) then
		do i = 1, ij
			ik = kdim - i
			ivec(i)= itempv/iprod(ik)
			itempv= itempv-iprod(ik)*ivec(i)
			ivec(i)= ivec(i)+1
 		enddo
 	endif
	ivec(kdim) = itempv+1
	if (qfor) call revers(ivec, kdim)
	return
12 	if (ivec(1) .LE. iprod(1)) GOTO 14
	ifault = 2
	return
14 	do 15 i = 2, kdim
		if (ivec(i) .LE. iprod(i)/iprod(i-1)) GOTO 15
		ifault = 2
		return
15 	continue
	if (.not. qfor) call revers(ivec, kdim)
	jsub = ivec(1)
	do 20 i = 2, kdim
		jsub=jsub+(ivec(i)-1)*iprod(i-1)
20 	continue
	return
end subroutine
!-----------------------------------
subroutine Base(T2,T3,T4,T5,T6,T7,T8,T9)
  	implicit none
  	integer, parameter :: MyI = selected_int_kind(8)
  	integer(kind=MyI), intent(in out) :: T2(2,2),T3(6,3),T4(24,4),T5(120,5),T6(720,6)
  	integer(kind=MyI), intent(in out) :: T7(5040,7),T8(40320,8),T9(362880,9)
  	integer(kind=MyI) :: q5(5),q6(6),q7(7),q8(8),q9(9)
  	integer(kind=MyI) ::nlast,count1
  	integer :: i,n
  	logical :: mtc
  	data q5/1,2,3,4,5/;data q6/1,2,3,4,5,6/;data q7/1,2,3,4,5,6,7/
  	T2=reshape ((/1,2,2,1/),shape(T2))
  	T3=reshape ((/1,1,2,2,3,3,2,3,1,3,1,2,3,2,3,1,2,1/),shape(T3))
  	T4=reshape ((/1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,2,2,3,3,4,4,1,1,3,3,4,4,&
  	1,1,2,2,4,4,1,1,2,2,3,3,3,4,2,4,2,3,3,4,1,4,1,3,2,4,1,4,1,2,2,3,1,3,1,2,4,3,4,2,3,2,4,&
  	3,4,1,3,1,4,2,4,1,2,1,3,2,3,1,2,1/),shape(T4))
  	Mtc= .TRUE.;count1=0;nlast=0;n=5
	do while (Mtc)
    	call Nexper (n, q5, mtc, nlast)
		count1=count1+1
		do i=1,n
			T5(count1,i)=q5(i)
		enddo
	enddo
  	Mtc= .TRUE.;count1=0;nlast=0;n=6
	do while (Mtc)
    	call Nexper (n, q6, mtc, nlast)
		count1=count1+1
		do i=1,n
			T6(count1,i)=q6(i)
		enddo
	enddo
	Mtc= .TRUE.;count1=0;nlast=0;n=7
	do while (Mtc)
    	call Nexper (n, q7, mtc, nlast)
		count1=count1+1
		do i=1,n
			T7(count1,i)=q7(i)
		enddo
	enddo
	Mtc= .TRUE.;count1=0;nlast=0;n=8
	do while (Mtc)
    	call Nexper (n, q8, mtc, nlast)
		count1=count1+1
		do i=1,n
			T8(count1,i)=q8(i)
		enddo
	enddo
	Mtc= .TRUE.;count1=0;nlast=0;n=9
	do while (Mtc)
    	call Nexper (n, q9, mtc, nlast)
		count1=count1+1
		do i=1,n
			T9(count1,i)=q9(i)
		enddo
	enddo
	return
 end subroutine
subroutine Oneper(n,kdq,kdp,nep,neq,T2,T3,T4,T5,T6,T7,T8,T9,Teq,Tep,Tbq,Tbp,Tgq,Tgp,Tfq,Tfp,iprq,iprp,cq,cp)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
	integer, intent(in) :: n
	integer(kind=MyI), intent(in) :: kdp,kdq
	integer(kind=MyI), intent(in) :: T2(2,2),T3(6,3),T4(24,4),T5(120,5),T6(720,6)
	integer(kind=MyI), intent(in) :: T7(5040,7),T8(40320,8),T9(362880,9)
	integer(kind=MyI), intent(in) :: Tbq(n,n),Tbp(n,n),Tgq(n),Tgp(n)
	integer(kind=MyI), intent(in) :: Teq(n,n),Tep(n,n),Tfq(n),Tfp(n)
  	integer(kind=MyI), intent(in) :: nep,neq,iprq(n),iprp(n)
	real(kind=MyR), intent(in out) :: cq(n),cp(n)	
	integer(kind=MyI) :: ivecq(n),ivecp(n)
	real(kind=MyR) 	 :: unif
	integer(kind=MyI) :: m,i,j,kj,kn,jj,kk,jk,yn
	integer(kind=MyI) :: jsubq,jsubp,req,rep
	logical qfor,qind	
	qfor = .false.; qind = .true.	
101	CALL RANDOM_NUMBER(unif)
	jsubq=1+floor(unif*neq)
	if (jsubq.gt.neq) goto 101
102	CALL RANDOM_NUMBER(unif)
	jsubp=1+floor(unif*nep)
	if (jsubp.gt.nep) goto 102
	call simdo(qind,qfor,iprq,kdq,jsubq,ivecq)
	call simdo(qind,qfor,iprp,kdp,jsubp,ivecp)
    do i=1,kdq
    	kk=ivecq(i);kj=Tfq(i);jk=Tgq(i)
    	select case (kj)
    		case (2)    				
    			do jj=1,2
    				m=T2(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    			enddo
    		case (3)
    			do jj=1,3
    				m=T3(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    			enddo
    		case (4)
    			do jj=1,4
    				m=T4(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    			enddo
    		case (5)
    			do jj=1,5
    				m=T5(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    			enddo
    		case (6)
    			do jj=1,6
    				m=T6(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    			enddo
    		case (7)
    			do jj=1,7
    				m=T7(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    			enddo
    		case (8)
    			do jj=1,8
    				m=T8(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    			enddo
    		case (9)
    			do jj=1,9
    				m=T9(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    			enddo
    	end select	
    enddo
    do i=1,kdp
    	kk=ivecp(i);kj=Tfp(i);jk=Tgp(i)
    	select case (kj)   				
    		case (2)
    			do jj=1,2
    				m=T2(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    			enddo
    		case (3)
    			do jj=1,3
    				m=T3(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    			enddo
    		case (4)
    			do jj=1,4
    				m=T4(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    			enddo
    		case (5)
    			do jj=1,5
    				m=T5(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    			enddo
    		case (6)
    			do jj=1,6
    				m=T6(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    			enddo
    		case (7)
    			do jj=1,7
    				m=T7(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    			enddo
    		case (8)
    			do jj=1,8
    				m=T8(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    			enddo
    		case (9)
    			do jj=1,9
    				m=T9(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    			enddo
    	end select	
    enddo
    return
end subroutine Oneper
subroutine SamPerm(n,neq,nep,size,kdq,kdp,Tbq,Tbp,iprp,iprq,Teq,Tep,Tfq,Tfp,Tgq,Tgp,cp,cq,rc,index)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
  	integer, intent(in) :: n,index
  	real (kind=MyR), intent(in):: size
	integer(kind=MyI), intent(in) :: Tbq(n,n),Tbp(n,n),Tgq(n),iprq(n),iprp(n)
	integer(kind=MyI), intent(in) :: Teq(n,n),Tep(n,n),Tfq(n),Tfp(n),Tgp(n)
  	integer(kind=MyI), intent(in) :: kdp,kdq,nep,neq
  	real(kind=MyR), intent(in out) :: cq(n),cp(n)
  	real(kind=MyR), intent(out)    :: rc
  	real(kind=MyR) 				 :: Cg,mrc
	integer(kind=MyI) :: i,j,kn,jsubq,jsubp,jp,jq,kq,km
	integer(kind=MyI) :: T2(2,2),T3(6,3),T4(24,4),T5(120,5),T6(720,6),T7(5040,7),T8(40320,8),T9(362880,9)
	call Base(T2,T3,T4,T5,T6,T7,T8,T9)
    kn=mod(n,2);Cg=n**2-kn
	call seed_init1()
	mrc=0.;km=size
    do kq=1,km
		call Oneper(n,kdq,kdp,nep,neq,T2,T3,T4,T5,T6,T7,T8,T9,Teq,Tep,Tbq,Tbp,Tgq,Tgp,Tfq,Tfp,iprq,iprp,cq,cp)
    	if(index.eq.3) then 
    			call Gini(n,Cg,cq,cp,rc)
    	else
    			call r4(n,cq,cp,rc)
    	endif
    	mrc=mrc+rc
  	enddo
    rc=mrc/size
    return
end subroutine SamPerm
!--------------------------------
subroutine AllPerm(n,neq,nep,kdq,kdp,Tbq,Tbp,iprp,iprq,Teq,Tep,Tfq,Tfp,Tgq,Tgp,cp,cq,rc,index)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
  	integer, intent(in) :: n,index
	integer(kind=MyI), intent(in) :: kdp,kdq
	integer(kind=MyI), intent(in) :: Tbq(n,n),Tbp(n,n),Tgq(n)
	integer(kind=MyI), intent(in) :: Teq(n,n),Tep(n,n),Tfq(n),Tfp(n),Tgp(n)
  	integer(kind=MyI), intent(in) :: iprq(n),iprp(n)
  	integer(kind=MyI), intent(in) :: nep,neq
  	real(kind=MyR), intent(in out) :: cq(n),cp(n)
  	real(kind=MyR), intent(out)    :: rc
  	real(kind=MyR)   :: Cg,mrc
	integer(kind=MyI) :: m,i,j,kj,kn,jj,kk,jk
	integer(kind=MyI) :: T2(2,2),T3(6,3),T4(24,4),T5(120,5),T6(720,6),T7(5040,7),T8(40320,8),T9(362880,9)
	integer(kind=MyI) :: ivec(n)
	integer(kind=MyI) :: jsubq,jsubp,jp,jq,nepq
	logical qfor,qind	
    call Base(T2,T3,T4,T5,T6,T7,T8,T9)
    qfor = .false.; qind = .true.
    kn=mod(n,2);Cg=n**2-kn
    nepq=nep*neq
    mrc=0.
  	do jq = 1,neq
  		if (neq.gt.1) then
    		jsubq = jq
    		call simdo(qind,qfor,iprq,kdq,jsubq,ivec)
    		do i=1,kdq
    	    	kk=ivec(i);kj=Tfq(i);jk=Tgq(i)
    			select case (kj)
    			case (2)    				
    				do jj=1,2
    					m=T2(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (3)
    				do jj=1,3
    					m=T3(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (4)
    				do jj=1,4
    					m=T4(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (5)
    				do jj=1,5
    					m=T5(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (6)
    				do jj=1,6
    					m=T6(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (7)
    				do jj=1,7
    					m=T7(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (8)
    				do jj=1,8
    					m=T8(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			case (9)
    				do jj=1,9
    					m=T9(kk,jj);cq(Teq(jk,jj))=Tbq(jk,m)
    				enddo
    			end select	
    		enddo
    	endif
    	do jp = 1,nep
    		if(nep.gt.1) then
    			jsubp = jp
    			call simdo(qind,qfor,iprp,kdp,jsubp,ivec)
    			do i=1,kdp
    				kk=ivec(i);kj=Tfp(i);jk=Tgp(i)
    				select case (kj)   				
    				case (2)
    					do jj=1,2
    						m=T2(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (3)
    					do jj=1,3
    						m=T3(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (4)
    					do jj=1,4
    						m=T4(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (5)
    					do jj=1,5
    						m=T5(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (6)
    					do jj=1,6
    						m=T6(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (7)
    					do jj=1,7
    						m=T7(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (8)
    					do jj=1,8
    						m=T8(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				case (9)
    					do jj=1,9
    						m=T9(kk,jj);cp(Tep(jk,jj))=Tbp(jk,m)
    					enddo
    				end select	
    			enddo
    		endif
    		if(index.eq.3) then 
    			call Gini(n,Cg,cq,cp,rc)
    		else
    			call r4(n,cq,cp,rc)
    		endif
    		mrc=mrc+rc
    	enddo
  	enddo
	rc=mrc/nepq
	return
end subroutine AllPerm
subroutine seed_init1 
	   integer(kind=4) :: seed_size
       integer(kind=4),allocatable :: seed(:)
       call random_seed(size=seed_size) ! find out size of seed
       allocate(seed(seed_size))
       call random_seed(get=seed) ! get system generated seed
       seed=314196593 
       call random_seed(put=seed) ! set current seed
       call random_seed(get=seed) ! get current seed
       deallocate(seed)           ! safe
end subroutine seed_init1
!------------------------------------------
subroutine R4(n,x,y,rc)
		implicit none
		integer, parameter :: MyR = selected_real_kind(15)
		integer, parameter :: MyI = selected_int_kind(8)
		integer, intent(in) :: n
		real (kind=MyR), intent(in) :: x(n),y(n)
		real(kind=Myr), intent(out) :: rc
		integer(kind=MyI) :: q(n),qs(n),eta(n),etas(n)
		real(kind=Myr)   :: A(n,n),B1(n),B2(n),B3(n),B4(n)
		real(kind=Myr)   :: Mp,Ws1,Ws2,Ws,n1,n2,ix
		integer(kind=MyI) :: i,j,k,kn,k2x
		Ws=0;kn=mod(n,2);n2=n/2;k2x=floor(n2);n1=n+1
		do i=1,k2x
			ix=i;ws=ws+(n1-ix)/ix
		enddo
 		Mp=kn+2.*Ws; Mp=Mp**2-n**2
 		do i=1,n
 			B1(i)=x(i);B2(i)=y(i)
 		enddo
 		call InSrtB(n,B1,B2)
		do i=1,n
			eta(i)=B1(i);q(i)=B2(i);etas(i)=n1-i
			Ws1=i
			do j=1,n
				Ws2=j;A(i,j)=max(Ws1,Ws2)/min(Ws1,Ws2)
			enddo
		enddo
		do i=1,n
			qs(i)=n1-q(i)
		enddo
		do i=1,n
			B1(i)=A(eta(i),qs(i));B2(i)=A(etas(i),qs(i))
			B3(i)=A(eta(i),q(i));B4(i)=A(etas(i),q(i))
		enddo
		Ws1=0.;Ws2=0.
		do i=1,n
			do j=1,n
				Ws1=Ws1+B1(i)*B4(j);Ws2=Ws2+B2(i)*B3(j)
			enddo
		enddo
		rc=(Ws1-Ws2)/Mp
end subroutine R4	
!......................................................
subroutine Kendall(n,Ck,x,y,rc)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
  	integer, intent(in) :: n
  	real (kind=MyR), intent(in)  :: Ck
	real (kind=MyR), intent(in) :: x(n),y(n)
	real (kind=MyR), intent(out) :: rc
	real (kind=MyR) :: con,dis
	integer :: i,j
	con=0.;dis=0.
	do i=1,n-1
		do j=i+1,n
			if(x(j).gt.x(i) .and. y(j).gt.y(i)) then
				con=con+1;goto 100
			endif
			if(x(j).lt.x(i) .and. y(j).lt.y(i)) then
				con=con +1;goto 100
			endif
			if(x(j).gt.x(i) .and. y(j).lt.y(i)) then
				dis=dis+1;goto 100
			endif
			if(x(j).lt.x(i) .and. y(j).gt.y(i)) then
				dis=dis +1;goto 100
			endif
100	continue
		enddo
	enddo
	rc=(con-dis)*Ck
	return
end subroutine Kendall
subroutine Gini(n,Cg,x,y,rc)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
  	integer, intent(in) :: n
	real (kind=MyR), intent(in)  :: Cg
	real (kind=MyR), intent(in)  :: x(n),y(n)
	real (kind=MyR), intent(out) :: rc
	real (kind=MyR) :: G1,G2
	integer :: i,n1
	G1=0.;G2=0.;n1=n+1
	do i=1,n
		G1=G1+abs(n1-x(i)-y(i));G2=G2+abs(x(i)-y(i))
	enddo
	rc=(G1-G2)/Cg
	return
end subroutine Gini
!---------------------
subroutine Spearman(n,Cs,x,y,rc)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
  	integer, intent(in) :: n
	real (kind=MyR), intent(in)  :: Cs
	real (kind=MyR), intent(in)  :: x(n),y(n)
	real (kind=MyR), intent(out) :: rc
	real (kind=MyR) :: s
	integer :: i
	s=0.
	do i=1,n
		s=s+(x(i)-y(i))**2
	enddo
	rc=1.-Cs*s
	return
end subroutine Spearman
!------------------------------------------
subroutine Factorial(x,y)
    implicit none
    integer, parameter :: MyI = selected_int_kind(8)
    integer(kind=MyI), intent(in) :: x
    integer(kind=MyI), intent(out) :: y
    integer(kind=MyI) :: j
    y=1
	if (x.le.1) then
   			return
		else
			do j=1,x
   				y = y*j
   			enddo
	end if
end subroutine Factorial
subroutine InSrtA(n,x,iy)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
	integer, intent(in) :: n
	integer(kind=MyI), intent(in out) :: iy(n)
    real (kind=MyR), intent (in out)  :: x(n)
    real (kind=MyR)   :: temp
    integer :: i,j,k,itemp
    do i=1,n
    	iy(i)=i
    enddo
    do i = 2, n
    	if (x(i) .lt. x(i-1)) then
          do j = i-2, 1, -1
             if (x(i) .gt. x(j)) goto 71
          enddo
          j=0
71	temp=x(i);itemp=iy(i)
		do k=i,j+2,-1
			iy(k)=iy(k-1);x(k)=x(k-1)
		enddo
		x(j+1)=temp;iy(j+1)=itemp
        endif
    enddo   
    return
end subroutine InsrtA
subroutine InSrtB(n,x,y)
	implicit none
	integer, parameter  :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
	integer, intent(in) :: n
    real (kind=MyR), intent (in out) :: x(n),y(n)
    real (kind=MyR)   :: tempx,tempy
    integer :: i,j,k
    do i = 2, n
    	if (x(i) .lt. x(i-1)) then
          do j = i-2, 1, -1
             if (x(i) .gt. x(j)) goto 71
          enddo
        	j=0
71			tempx=x(i);tempy=y(i)
			do k=i,j+2,-1
				y(k)=y(k-1);x(k)=x(k-1)
			enddo
			x(j+1)=tempx;y(j+1)=tempy
        endif
    enddo 
    return
end subroutine InsrtB
!------------------------------
subroutine meanrank(n,x,rx,prime)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
	integer, intent(in) :: n, prime
	real (kind=MyR), intent(in)  :: x(n)
	real (kind=MyR), intent(out) :: rx(n)
	real (kind=MyR) :: ws,t0
	integer :: i,j,k,i1,t,m
	integer :: p(n)
	do i=1,n
		p(i)=i;rx(i)=0
	enddo
	i=1
	do while(i.lt.n)
		 if (x(i+1).ne.x(i)) then
		 	rx(i)=i;i=i+1
		 else 
		 	t=1;i1=i+1
		 	select case(prime)
		 		case(1)
		 			t0=i
		 		case(2)
		 			t0=i**2
		 	end select
		 	do j=i1,n
			   if (x(j).eq.x(i)) then
			   		select case(prime)
		 				case(1)
		 				t0=t0+p(j)
		 				case(2)
		 				t0=t0+p(j)**2
		 			end select
			        t=t+1
			    else
			    	exit
			    endif
			  enddo
			  select case(prime)
		 			case(1)
		 				ws=t0/t
		 			case(2)
		 				ws=(t0/t);ws=sqrt(ws)
		 			case(3)
		 				ws=i1-1 
		 			case(4)
		 				ws=i1+t-2
		 		end select
				m=i+t-1
				do k=i,m
					rx(k)=ws
				enddo
				i=m+1
			endif
	  enddo
	  if(rx(n).eq.0) rx(n)=n
	  return
end subroutine
!-----------------------------------
subroutine Myra(n,x,iy,AD)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
	integer, intent(in) :: n
	integer, intent(in out) :: iy(n)
    real(kind=MyR), intent (in out)  :: x(n)
    real(kind=MyR)   :: unif,div
    logical, intent(in) :: AD
    integer :: i,j,m,n1
    div=n*n;n1=n+1
    CALL RANDOM_NUMBER(unif)
    unif=unif/div
    do i=1,n		
		if(AD) then
			x(i)=x(i)+i*unif
		else
		 	x(i)=x(i)+(n1-i)*unif
		endif
	enddo
    i=1
100	m=0
    do j=1,n
    	if(x(j).le.x(i)) m=m+1
    enddo
    iy(i)=m
    i=i+1
    if(i.le.n) goto 100		
    return
end subroutine Myra
!=====================================================
subroutine SerTies(n,q,nt,n2,Ta,Tb,Te,Tf,Tg,ipr,ne,cq)
	implicit none
	integer, parameter :: MyR = selected_real_kind(15)
	integer, parameter :: MyI = selected_int_kind(8)
	integer, intent(in) :: n
	real (kind=MyR),  intent(in)   :: q(n)
	real (kind=MyR),  intent(out)  :: cq(n)
	integer(kind=MyI), intent(out)  :: ne
	integer(kind=MyI), intent(out)  :: nt,n2
	integer(kind=MyI), intent(out) :: Tb(n,n),Ta(n),Te(n,n),Tf(n),Tg(n)
	integer(kind=MyI), intent(out) :: ipr(n)
	integer(kind=MyI) :: ws,wt
	integer(kind=MyI) :: Td(n)
	integer :: i,j,m,mm,Sy(n),ord(n)
	real (kind=MyR) :: qs(n),qx(n)
	logical :: ksw,isw,jsw
	do i=1,n
		Sy(i)=0;qx(i)=q(i);ord(i)=i
	enddo
	ksw=.true.;call Myra(n,qx,ord,ksw)
	nt=1;qs(nt)=q(1);Ta(nt)=1;Sy(1)=1;Tb(nt,1)=ord(1)
	Te(nt,1)=1;mm=1
	ksw=.true.
	do while(ksw)
    	do j=1,n
    		if(Sy(j).ne.1) then
    			if(qs(nt).eq.q(j)) then
    				Ta(nt)=Ta(nt)+1;Sy(j)=1
    				Tb(nt,Ta(nt))=ord(j)
    				Te(nt,Ta(nt))=j
    			endif
    		endif
    	enddo
		m=0
		do i=1,n
			if(Sy(i).ne.1) then
				m=i
				exit
			endif
		enddo
		if(m.gt.0) then
			nt=nt+1
			qs(nt)=q(m);Ta(nt)=1
			isw=.false.
			do j=1,n
				if (Sy(j).eq.0) then
					Tb(nt,1)=ord(j);Te(nt,1)=j
					isw=.true.;Sy(j)=1
					exit
				endif
				if(isw) exit
			enddo
		  else
			ksw=.false.
		  endif
	enddo
	do i=1,n
    	cq(i)=0
    enddo
    do i=1,nt
    	if(Ta(i).eq.1) cq(Te(i,1))=Tb(i,1)
    enddo
    n2=0
	do i=1,nt
		if(Ta(i).gt.1) then
			ws=Ta(i);n2=n2+1
			call Factorial(ws,wt)
			Td(n2)=wt;Tf(n2)=Ta(i);Tg(n2)=i
		endif
	enddo
	if (n2.gt.0) then
		ipr(1) = Td(n2)
  		do i = 2, n2
    		ipr(i) = ipr(i-1) * Td(n2+1-i)
  		enddo
  		ne=ipr(n2)
  	else
  		n2=1;Td(1)=1;ne=1;ipr(1)=1
  	endif
  	return
end subroutine SerTies
!-------------------------------------
subroutine S_Wood(n,q,p,Taq,Tap,ntq,ntp,rc)
  	implicit none
  	integer, parameter :: MyR = selected_real_kind(15)
  	integer, parameter :: MyI = selected_int_kind(8)
	integer, intent(in) :: n
	real (kind=MyR),  intent(in)  :: q(n),p(n)
	integer(kind=MyI), intent(in) :: Taq(n),Tap(n)
	integer(kind=MyI), intent(in) :: ntp,ntq
  	real (kind=MyR),  intent(out) :: rc
  	real (kind=MyR)  :: cq(n),cp(n),qx(n),px(n)
  	real (kind=MyR)  :: sx,sz
  	integer(kind=MyI) :: iq(n),ip(n)  
  	integer :: i,prime	
  	do i=1,n
		qx(i)=q(i);iq(i)=i
	enddo
	call InSrtA(n,qx,iq)
	do i=1,n
		px(i)=p(i);ip(i)=i
	enddo	
	call InSrtA(n,px,ip)
	prime=1
	call meanrank(n,qx,cq,prime)
	call meanrank(n,px,cp,prime)
	do i=1,n
		qx(iq(i))=cq(i);px(ip(i))=cp(i)
	enddo
    sx=0
    do i=1,ntq
    	sx=sx+(Taq(i)**3-Taq(i))
    enddo
    do i=1,ntp
    	sx=sx+(Tap(i)**3-Tap(i))
    enddo
    sz=0
    do i=1,n
    	sz=sz+(qx(i)-px(i))**2
    enddo
    sx=sx/12.;rc=6.*(sz+sx)/(n**3-n);rc=1.-rc
    return
end subroutine S_Wood
!-------------------------------------------
subroutine Woodbury(n,q,p,size,index,rc,ifault)
  	implicit none
  	integer, parameter :: MyR = selected_real_kind(15)
  	integer, parameter :: MyI = selected_int_kind(8)
  	integer, intent(in) :: n,index
  	integer, intent(in out) :: ifault
  	real (kind=MyR), intent(in) :: q(n),p(n)
  	real (kind=MyR), intent(in):: size
  	real (kind=MyR), intent(out)    :: rc
  	real (kind=MyR)  :: cq(n),cp(n)
  	real (kind=MyR)  :: Ck
	integer(kind=MyI) :: Tbq(n,n),Tbp(n,n),iprp(n),iprq(n),Teq(n,n),Tep(n,n)
	integer(kind=MyI) :: Tfq(n),Tfp(n),Tgq(n),Tgp(n),Tap(n),Taq(n)
	integer(kind=MyI) :: neq,nep,n2q,n2p,ntq,ntp
	integer :: i,imaq,imap
	ifault=0
	if (index.eq.2) then
		Ck=2./n/(n-1)
		call Kendall(n,Ck,q,p,rc)
		return
	endif
	if (index.eq.1) then
		call SerTies(n,q,ntq,n2q,Taq,Tbq,Teq,Tfq,Tgq,iprq,neq,cq)
		call SerTies(n,p,ntp,n2p,Tap,Tbp,Tep,Tfp,Tgp,iprp,nep,cp)
		call S_Wood(n,q,p,Taq,Tap,ntq,ntp,rc)
		return
	endif	
	if ((index.eq.3) .or. (index.eq.4)) then
		call SerTies(n,q,ntq,n2q,Taq,Tbq,Teq,Tfq,Tgq,iprq,neq,cq)
		call SerTies(n,p,ntp,n2p,Tap,Tbp,Tep,Tfp,Tgp,iprp,nep,cp)
		imaq=0;imap=0
		if(n2q.gt.0) then
			do i=1,n2q
				if(imaq.lt.Taq(i)) imaq=Taq(i)
			enddo
		endif
		if(n2p.gt.0) then 
			do i=1,n2p
				if(imap.lt.Tap(i)) imap=Tap(i)
			enddo
		endif
		if (imaq.gt.9 .or. imap.gt.9) then
			ifault=1
			return
		endif
		if(neq*nep.gt.size) then
  			call SamPerm(n,neq,nep,size,n2q,n2p,Tbq,Tbp,iprp,iprq,Teq,Tep,Tfq,Tfp,Tgq,Tgp,cp,cq,rc,index)
  		else
  			call AllPerm(n,neq,nep,n2q,n2p,Tbq,Tbp,iprp,iprq,Teq,Tep,Tfq,Tfp,Tgq,Tgp,cp,cq,rc,index)
		endif
	endif
	return
end subroutine Woodbury
!===================================================================
subroutine DealwT(n,q,p,ities,index,size,repl,rc,Hilo,ifault)
  	implicit none
  	integer, parameter :: MyR = selected_real_kind(15)
  	integer, parameter :: MyI = selected_int_kind(8)
  	integer, intent(in) :: n,index,ities
  	real (kind=MyR), intent(in)	 :: q(n),p(n)
  	real (kind=MyR), intent(in)  :: size,repl
  	real (kind=MyR), intent(out) :: rc
  	real (kind=MyR), intent(out) :: Hilo(2)
  	integer, intent(out) 	 	 :: ifault
  	real (kind=MyR) :: Cs,Ck,Cg
  	real (kind=MyR) :: qx(n),px(n)
  	integer :: kn,prime
  	select case (ities)
    	case (1) 
  			call Woodbury(n,q,p,size,index,rc,ifault)
		case (2)
			call GGH(n,q,p,Hilo,index)
			rc=0.5*(Hilo(1)+Hilo(2))
		case (3)
			call GGH(n,q,p,Hilo,index)
			rc=0.5*(Hilo(1)+Hilo(2))
			call wGini(n,p,q,repl,Hilo,index,rc)
		case (4)
			prime=1
			call Midrank(n,q,p,prime,qx,px)
			select case (index)
				case(1)
					Cs=6./(n**3-n);call Spearman(n,Cs,qx,px,rc)
				case(2)
					Ck=2./n/(n-1);call Kendall(n,Ck,qx,px,rc)
				case(3)
  					kn=mod(n,2);Cg=(n**2-kn)/2.;call Gini(n,Cg,qx,px,rc)
  				case(4)
  					call R4(n,qx,px,rc)
  			end select
		case (5)
			prime=2
			call Midrank(n,q,p,prime,qx,px)
			select case (index)
				case(1)
					Cs=6./(n**3-n);call Spearman(n,Cs,qx,px,rc)
				case(2)
					Ck=2./n/(n-1);call Kendall(n,Ck,qx,px,rc)
				case(3)
  					kn=mod(n,2);Cg=(n**2-kn)/2.;call Gini(n,Cg,qx,px,rc)
  				case(4)
  					call R4(n,qx,px,rc)
  			end select
  		case (6)
  			select case (index)
				case(1)
					Cs=6./(n**3-n);call Spearman(n,Cs,q,p,rc)
				case(2)
					Ck=2./n/(n-1);call Kendall(n,Ck,q,p,rc)
				case(3)
  					kn=mod(n,2);Cg=(n**2-kn)/2.;call Gini(n,Cg,q,p,rc)
  				case(4)
  					call R4(n,q,p,rc)
  			end select
    end select	
    return
end subroutine DealwT