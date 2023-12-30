!
!=====================================================================*
! Shock tube problem
!=====================================================================*
    program shocktube
    implicit none
    integer::i,n
    integer,parameter::ndiv  =  1001
    real(8),parameter::rgasl =  288.696520833d0 !287.d0 ! Gas constant (Air)
!    real(8),parameter::rgash =  2077.d0 ! Gas constant (He)
    real(8),parameter::rgash =  288.696520833d0  !287.d0 ! Air, 2077.d0 ! Gas constant (He)
    real(8)::x0,xl,xr
    real(8),dimension(ndiv)::x,uvel,gamma,rgas,pres,dens,temp,sos,uvel_init,pres_init,dens_init,temp_init,sos_init
    real(8)::tlap,gammah,gammal,uvelh,uvell,presh,presl,densh,densl,temph,templ,sosh,sosl
    real(8)::qst,qstpre,falpha,fval,f,df
    real(8)::pres1,dens1,uvel1,pres2,dens2,uvel2,pres3,dens3,uvel3,pres4,dens4,uvel4,pres5,dens5,uvel5
    real(8)::sos4,vel2,velcd,mach
    character(len=25)::shocktubefile='profileshocktube.dat'
!
! Positions
    x0 =  0.d0
!    xr =  5.d-1
!    xl = -5.d-1
    xr =  250.d-3
    xl = -250.d-3
    do i=1,ndiv
      x(i) = xl + (xr-xl)*dble(i-1)/dble(ndiv-1)
    end do
!
! Parameters
    tlap  = 2.0d-4
!----------------------
!    gammah = 1.66d0
!    uvelh  = 0.d0
!    presh  = 55.7d6
!!    densh  = 6.26d0
!    temph  = 3687.d0
!----------------------
!    gammal = 1.40d0
!    uvell  = 0.d0
!    presl  = 1.d3
!!    densl  = 1.16d-2
!    templ  = 300.d0
!----------------------
!    uvelh = 0.d0
!    presh = 1.d0
!    densh = 1.d0
!    uvell = 0.d0
!    presl = 1.d-1
!    densl = 1.25d-1
!----------------------
    gammah = 1.40d0
    uvelh  = 0.d0
    presh  = 866.0895625d0
!    densh  = 1.d-2
    temph  = 300.d0
!----------------------
    gammal = 1.40d0
    uvell  = 0.d0
    presl  = 8.660895625d0
!    densl  = 1.d-4
    templ  = 300.d0
!----------------------
!
! Speed of sound and temperature
    densh = presh/(temph*rgash)
    densl = presl/(templ*rgasl)
    sosh  = dsqrt(gammah*presh/densh)
    sosl  = dsqrt(gammal*presl/densl)
!    temph = presh/(densh*rgas)
!    templ = presl/(densl*rgas)
!
! Initial conditions
    do i=1,ndiv
      if(x(i) <= 0.d0)then
        uvel(i)  = uvelh
        pres(i)  = presh
        dens(i)  = densh
        temp(i)  = temph
        sos(i)   = sosh
        gamma(i) = gammah
        rgas(i)  = rgash
      else
        uvel(i)  = uvell
        pres(i)  = presl
        dens(i)  = densl
        temp(i)  = templ
        sos(i)   = sosl
        gamma(i) = gammal
        rgas(i)  = rgasl
      end if
    end do
!
    do i=1,ndiv
      uvel_init(i)  = uvel(i)
      pres_init(i)  = pres(i)
      dens_init(i)  = dens(i)
      temp_init(i)  = temp(i)
      sos_init(i)   = sos(i)
    end do
!
! Solve the relation conditioning the shock tube problem
    falpha = 2.d0*gammah/(gammah-1.d0)
    qst    = 1.d0
    do n=1,100
      qstpre = qst
      fval = (gammal-1.d0)+(gammal+1.d0)*qst
      f  = dsqrt(2.d0/gammal)*(qst-1.d0)/dsqrt(fval)                            &
         - 2.d0/(gammah-1.d0)*sosh/sosl*(1.d0-(presl/presh*qst)**(1.d0/falpha))    &
         - (uvelh-uvell)/sosh
      df = dsqrt(2.d0/gammal)*(1.d0-0.5d0*(gammal+1.d0)*(qst-1.d0)/fval)/dsqrt(fval)    &
         - 2.d0/(gammah-1.d0)*sosh/sosl*(-1.d0/falpha*(presl/presh*qst)**(1.d0/falpha-1.d0)*presl/presh)
      qst = qst-f/df
      if((dabs(qst-qstpre)/qstpre) <= 1.d-10) exit
      if(n==100)then
        write(6,*)'There is no shock wave'
        stop
      end if
    end do
    write(6,*)qst
!
! shock tube distributions
! region 1 (shock wave - right-hand side wall)
    pres1 = presl
    dens1 = densl
    uvel1 = uvell
!
! region 2 (contact discontinuity - shock wave)
    pres2 = presl*qst
    dens2 = densl*((gammal-1.d0+(gammal+1.d0)*qst)/(gammal+1.d0+(gammal-1.d0)*qst))
    uvel2 = uvell+sosl*dsqrt(2.d0/gammal)*(qst-1.d0)/dsqrt(gammal-1.d0+(gammal+1.d0)*qst)
    vel2  = uvell+((qst-1.d0)*sosl**2.d0)/(gammal*(uvel2-uvell))
!
! region3 (x0 - contact discontinuity)
    pres3 = pres2
    dens3 = densh*(pres3/presh)**(1.d0/gammah)
    uvel3 = uvel2
    velcd = uvel2
!
! Distibution (Compress side)
    do i=1,ndiv
      if(x(i) >  x0+((gammah+1.d0)/2.d0*velcd-sosh-(gammah-1.d0)/2.d0*uvelh)*tlap .and.        &
         x(i) <= x0+velcd*tlap)then
        pres(i)  = pres3
        dens(i)  = dens3
        uvel(i)  = uvel3
        gamma(i) = gammah
        rgas(i)  = rgash
      elseif(x(i) >  x0+velcd*tlap .and.        &
             x(i) <= x0+vel2*tlap)then
        pres(i)  = pres2
        dens(i)  = dens2
        uvel(i)  = uvel2
        gamma(i) = gammal
        rgas(i)  = rgasl
      elseif(x(i) > x0+vel2*tlap)then
        pres(i)  = pres1
        dens(i)  = dens1
        uvel(i)  = uvel1
        gamma(i) = gammal        
        rgas(i)  = rgasl
      end if
      temp(i) = pres(i)/(dens(i)*rgas(i))
    end do
!
! (Expansion side)
    do i=1,ndiv
      if(x(i) >= x0+(uvelh-sosh)*tlap .and.        &
         x(i) <= x0+((gammah+1.d0)/2.d0*velcd-sosh-(gammah-1.d0)/2.d0*uvelh)*tlap)then
!
! region 4 (in expansion fan)
        uvel(i)  = 2.d0/(gammah+1.d0)*((x(i)-x0)/tlap+sosh+(gammah-1.d0)/2.d0*uvelh)
        sos4     = sosh-(gammah-1.d0)/2.d0*(uvel(i)-uvelh)
        pres(i)  = presh*(sos4/sosh)**falpha
        dens(i)  = densh*(pres(i)/presh)**(1.d0/gammah)
        gamma(i) = gammah        
        rgas(i)  = rgash
      elseif(x(i) < x0+(uvelh-sosh)*tlap)then
        pres(i)  = presh
        dens(i)  = densh
        uvel(i)  = uvelh        
        gamma(i) = gammah        
        rgas(i)  = rgash
      end if
      temp(i) = pres(i)/(dens(i)*rgas(i))
    end do
!
! Print out
     if(shocktubefile /= 'none')then    
      open(100,file=shocktubefile,status='unknown',form='formatted')
        write(100,"(a)")'variables = x,u,v,Ttr,p,rho,Mach,h'         
        write(100,*) ' zone T="Exact',tlap,'sec" , i= ',ndiv,' f=point '
        do i=1,ndiv
          sos(i) = dsqrt(gamma(i)*pres(i)/dens(i))
          mach = dabs(uvel(i))/sos(i)
          write(100,"(8(1pd18.10))")x(i),uvel(i),0.d0,temp(i),pres(i),dens(i),mach,0.d0
        end do
! Inital
        write(100,*) ' zone T="Exact_initial" , i= ',ndiv,' f=point '
        do i=1,ndiv
          mach = dabs(uvel_init(i))/sos_init(i)
          write(100,"(8(1pd18.10))")x(i),uvel_init(i),0.d0,temp_init(i),pres_init(i),dens_init(i),mach,0.d0
        end do
      close(100)
    end if
!
    stop
    end program
