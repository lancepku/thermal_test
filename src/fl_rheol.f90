
!  Rheology (Update stresses depending on rheology)
!  Calculate total finite strain and plastic strain  
    
subroutine fl_rheol
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dimension depl(4)
dimension s11p(4),s22p(4),s12p(4),s33p(4),s11v(4),s22v(4),s12v(4),s33v(4)
logical rh_sel
dimension sarc1(nz+1),sarc2(nz+1)!,dxinj(nz+1),rate_inject_d(nz+1)
!if( mod(nloop,10).eq.0 .OR. ireset.eq.1 ) then
!    rh_sel = .true.
!else
!    rh_sel = .false.
!endif
rh_sel = .true.

!XXX: irh==11, or irh>=11?
irh=irheol(mphase)
if(irh.eq.11) call init_visc

D_axial = 0.25*(cord(1,nx/2,2)+cord(2,nx/2,2)+cord(1,nx/2+1,2)+cord(2,nx/2+1,2))
!print *,'D',D_axial
!define the depth of AMC
Z1100 = 1
do j=1,nz
if (temp(j,nx/2).lt.1100.)then
Z1100=j+1
else if (temp(j,nx/2).ge.1100)then
Z1100 = Z1100
end if
end do
Zamc = abs(cord(int(Z1100),nx/2,2)+cord(int(Z1100)+1,nx/2+1,2)+cord(int(Z1100)+1,nx/2,2)+cord(int(Z1100),nx/2+1,2))+ D_axial
Z600 = 1
do j=1,nz
tempar = temp(j,nx/2)!+temp(j+1,nx/2)+temp(j,nx/2+1)+temp(j+1,nx/2+1))
!do j=1,nz
if (tempar.lt.600.)then
Z600=j+1
else if (tempar.ge.600)then
Z600 = Z600
end if
end do

HL=0.25*abs(cord(int(Z600),nx/2,2)+cord(int(Z600)+1,nx/2+1,2)+cord(int(Z600)+1,nx/2,2)+cord(int(Z600),nx/2+1,2))+D_axial
orgh1 = 3000 ! density between magma and rock
orgh2 = 20000 ! density between water and rock
pd = 1.e7 ! driving pressure
D = -D_axial ! make sure axial depth for valley is positive
!print *,D
HG1 = -3.*(orgh1*HL+orgh2*D)
HG2 = sqrt(9*((orgh1*HL+orgh2*D)**2)+24*orgh1*pd*HL)
!print *,HG2
HG= (HG1+HG2)/(2.*orgh1)
!if(iynts.eq.1) call init_temp
!print *,HL,HG
! Initial stress boundary condition
! Accretional Stresses
!write(*,*) rate_inject_0,rate_inject_1
!write(*,*) dA
rate_inject = 0.7*Vp_
!print *, rate_inject
!nelem_inject = nz-1
if (ny_inject.gt.0) then
         sarc1 = 0.
         sarc2 = 0. 
         if (ny_inject.eq.1) iinj = 1
         if (ny_inject.eq.2) iinj = nx/2 
!         write (*,*) nelem_inject
         nelem_inject = Z600+1
         !average dx for injection:
         dxinj = 0.
         do jinj = 1,nelem_inject
            iph=iphase(jinj,iinj)
            dxinj=dxinj+cord(jinj,iinj+2,1)-cord(jinj,iinj-1,1)
         enddo
         dxinj = dxinj/nelem_inject 
         ! Constants Elastic:
         poiss = 0.5*rl(iph)/(rl(iph)+rm(iph))
         young = rm(iph)*2.*(1.+poiss)   
         ! Additional Stress:
        do j = 1,nz-1
!          dxinj = cord(j,iinj+1,1)-cord(j,iinj,1)
          yc = abs(0.25 * (cord(j,iinj,2)+cord(j+1,iinj,2)+cord(j,iinj+1,2)+cord(j+1,iinj+1,2)))
          A = yc - HL
          AB = max(A,0.0)
          AC = HG
          rate_inject_a = rate_inject*(1-AB/AC)
!          print *,rate_inject_a
          rate_inject_d = 0.0
         sarc1 = -young/(1.-poiss*poiss)*rate_inject_d/dxinj*dt
!         sarc1(j)= -young*(1.-poiss)/((1.+poiss)*(1.-2.*poiss))*rate_inject/dxinj*dt
         sarc2(j)= sarc1(j)*poiss/(1.-poiss)
end do         
!write(*,*) sarc1(5),sarc2(5)
endif

irh_mark = 0

! max. deviatoric strain and area change of current time step
curr_devmax = devmax
curr_dvmax = dvmax

!$OMP Parallel Private(i,j,k,iph,irh,bulkm,rmu,coh,phi,psi, &
!$OMP                  stherm,hardn,vis, &
!$OMP                  de11,de22,de12,de33,dv, &
!$OMP                  s11p,s22p,s12p,s33p, &
!$OMP                  s11v,s22v,s12v,s33v, &
!$OMP                  depl,ipls,diss, &
!$OMP                  sII_plas,sII_visc, &
!$OMP                  quad_area,s0a,s0b,s0) &
!$OMP firstprivate(irh_mark)
!$OMP do schedule(guided) reduction(max: curr_devmax, curr_dvmax)
do 3 i = 1,nx-1
    do 3 j = 1,nz-1
        diss = 0.
        ! iphase (j,i) is number of a phase NOT a rheology
        iph = iphase(j,i)
        irh = irheol(iph)

!        if(ny_inject.gt.0.and.j.le.nelem_inject) then
!        if(i.eq.iinj.or.i.eq.iinj-1) irh_mark = 1
!        if(i.eq.iinj) irh = 3 
!        endif

        ! Elastic modules & viscosity & plastic properties
        bulkm = rl(iph) + 2.*rm(iph)/3.
        rmu   = rm(iph)

        ! Thermal stresses (alfa_v = 3.e-5 1/K)
        stherm = 0.
        if (istress_therm.gt.0) stherm = -alfa(iph)*bulkm*(temp(j,i)-temp0(j,i))


        ! Preparation of plastic properties
        if (irh.eq.6 .or. irh.ge.11) call pre_plast(i,j,coh,phi,psi,hardn)
              
        ! Re-evaluate viscosity
        if (irh.eq.3 .or. irh.eq.12) then 
            if( mod(nloop,ifreq_visc).eq.0 .OR. ireset.eq.1 ) visn(j,i) = Eff_visc(j,i)
            if (i.eq.1.or.i.eq.2) visn(j,i) = v_min
        endif
        vis = visn(j,i)

        ! Cycle by triangles
        do k = 1,4

            ! Incremental strains
            de11 = strainr(1,k,j,i)*dt
            de22 = strainr(2,k,j,i)*dt
            de12 = strainr(3,k,j,i)*dt
            de33 = 0.
            dv = dvol(j,i,k)
            s11p(k) = stress0(j,i,1,k) + stherm 
            s22p(k) = stress0(j,i,2,k) + stherm 
            if(ny_inject.gt.0.and.j.le.nelem_inject) then
                !XXX: iinj is un-init'd if ny_inject is not 1 or 2.
                if(i.ge.iinj-1.and. i.le.iinj+1) then
                    s11p(k) = stress0(j,i,1,k) + stherm +sarc1(j)
                    s22p(k) = stress0(j,i,2,k) + stherm +sarc2(j)
                    !!            irh = 1
                endif
            endif
            s12p(k) = stress0(j,i,3,k) 
            s33p(k) = stress0(j,i,4,k) + stherm
            s11v(k) = s11p(k)
            s22v(k) = s22p(k)
            s12v(k) = s12p(k)
            s33v(k) = s33p(k)
!!            if(abs(sarc11).gt.0.) write(*,*) i,j,sarc11,sarc22
            if (irh.eq.1) then
                ! elastic
                call elastic(bulkm,rmu,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de12)
                irheol_fl(j,i) = 0  
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)

            elseif (irh.eq.3) then
                ! viscous
                call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),de11,de22,de33,de12,dv,&
                     ndim,dt,curr_devmax,curr_dvmax)
                irheol_fl(j,i) = -1  
                stress0(j,i,1,k) = s11v(k)
                stress0(j,i,2,k) = s22v(k)
                stress0(j,i,3,k) = s12v(k)
                stress0(j,i,4,k) = s33v(k)

            elseif (irh.eq.6) then
                ! plastic
                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
                     ten_off,ndim,irh_mark)
                irheol_fl(j,i) = 1
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)

            elseif (irh.ge.11) then 
                ! Mixed rheology (Maxwell or plastic)
                if( rh_sel ) then
                    call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                        s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
                        ten_off,ndim,irh_mark)
                    call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),&
                        de11,de22,de33,de12,dv,&
                        ndim,dt,curr_devmax,curr_dvmax)
                else ! use previously defined rheology
                    if( irheol_fl(j,i) .eq. 1 ) then
                        call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                            s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
                            ten_off,ndim,irh_mark)
                        stress0(j,i,1,k) = s11p(k)
                        stress0(j,i,2,k) = s22p(k)
                        stress0(j,i,3,k) = s12p(k)
                        stress0(j,i,4,k) = s33p(k)
                    else  ! irheol_fl(j,i) = -1
                        call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),&
                            de11,de22,de33,de12,dv,&
                            ndim,dt,curr_devmax,curr_dvmax)
                        stress0(j,i,1,k) = s11v(k)
                        stress0(j,i,2,k) = s22v(k)
                        stress0(j,i,3,k) = s12v(k)
                        stress0(j,i,4,k) = s33v(k)
                    endif
                endif
            endif
        enddo

        if( irh.ge.11 .AND. rh_sel ) then
            ! deside - elasto-plastic or viscous deformation
            sII_plas = (s11p(1)+s11p(2)+s11p(3)+s11p(4)-s22p(1)-s22p(2)-s22p(3)-s22p(4))**2 &
                     + 4*(s12p(1)+s12p(2)+s12p(3)+s12p(4))**2

            sII_visc = (s11v(1)+s11v(2)+s11v(3)+s11v(4)-s22v(1)-s22v(2)-s22v(3)-s22v(4))**2 &
                     + 4*(s12v(1)+s12v(2)+s12v(3)+s12v(4))**2

            if (sII_plas .lt. sII_visc) then
                do k = 1, 4
                    stress0(j,i,1,k) = s11p(k)
                    stress0(j,i,2,k) = s22p(k)
                    stress0(j,i,3,k) = s12p(k)
                    stress0(j,i,4,k) = s33p(k)
                end do
                irheol_fl (j,i) = 1
            else 
                do k = 1, 4
                    stress0(j,i,1,k) = s11v(k)
                    stress0(j,i,2,k) = s22v(k)
                    stress0(j,i,3,k) = s12v(k)
                    stress0(j,i,4,k) = s33v(k)
                end do
                irheol_fl (j,i) = -1
            endif
        endif


        ! Averaging of isotropic stresses for pair of elements
        if (mix_stress .eq. 1 ) then
        
            ! For A and B couple:
            ! area(n,it) is INVERSE of "real" DOUBLE area (=1./det)
            quad_area = 1./(area(j,i,1)+area(j,i,2))
            s0a=0.5*(stress0(j,i,1,1)+stress0(j,i,2,1))
            s0b=0.5*(stress0(j,i,1,2)+stress0(j,i,2,2))
            s0=(s0a*area(j,i,2)+s0b*area(j,i,1))*quad_area
            stress0(j,i,1,1) = stress0(j,i,1,1) - s0a + s0
            stress0(j,i,2,1) = stress0(j,i,2,1) - s0a + s0
            stress0(j,i,1,2) = stress0(j,i,1,2) - s0b + s0
            stress0(j,i,2,2) = stress0(j,i,2,2) - s0b + s0

            ! For C and D couple:
            quad_area = 1./(area(j,i,3)+area(j,i,4))
            s0a=0.5*(stress0(j,i,1,3)+stress0(j,i,2,3))
            s0b=0.5*(stress0(j,i,1,4)+stress0(j,i,2,4))
            s0=(s0a*area(j,i,4)+s0b*area(j,i,3))*quad_area
            stress0(j,i,1,3) = stress0(j,i,1,3) - s0a + s0
            stress0(j,i,2,3) = stress0(j,i,2,3) - s0a + s0
            stress0(j,i,1,4) = stress0(j,i,1,4) - s0b + s0
            stress0(j,i,2,4) = stress0(j,i,2,4) - s0b + s0
        endif

        if (irh.eq.6 .or. irh.ge.11) then
            !  ACCUMULATED PLASTIC STRAIN
            ! Average the strain for pair of the triangles
            ! Note that area (n,it) is inverse of double area !!!!!
            aps(j,i) = aps(j,i) &
                 + 0.5*( depl(1)*area(j,i,2)+depl(2)*area(j,i,1) ) / (area(j,i,1)+area(j,i,2)) &
                 + 0.5*( depl(3)*area(j,i,4)+depl(4)*area(j,i,3) ) / (area(j,i,3)+area(j,i,4))
            if( aps(j,i) .lt. 0. ) aps(j,i) = 0.

            !	write(*,*) depl(1),depl(2),depl(3),depl(4),area(j,i,1),area(j,i,2),area(j,i,3),area(j,i,4)

        ! LINEAR HEALING OF THE PLASTIC STRAIN
            if (tau_heal .ne. 0.) &
!                 aps (j,i) = aps (j,i)/(1.+dt/tau_heal)
            !if (ny_inject.gt.0.and.i.eq.iinj) aps (j,i) = 0.
        ! HEALING OF THE PLASTIC STRAIN- now tau_heal is reference strain rate

      !      if (aps (j,i).gt.1.0) then
!                aps (j,i) = aps (j,i)*(1.0-dt*tau_heal)
       !     else
                aps (j,i) = aps (j,i)-dt*tau_heal
!            endif
!
            if (aps (j,i).lt.0.) aps (j,i) = 0.
            if (ny_inject.gt.0.and.i.le.iinj .and.i.ge.iinj) aps (j,i) = 0.00000  ! limited strain in Dike,Liu

        end if

        ! TOTAL FINITE STRAIN
        strain(j,i,1) = strain(j,i,1) + 0.25*dt*(strainr(1,1,j,i)+strainr(1,2,j,i)+strainr(1,3,j,i)+strainr(1,4,j,i))
        strain(j,i,2) = strain(j,i,2) + 0.25*dt*(strainr(2,1,j,i)+strainr(2,2,j,i)+strainr(2,3,j,i)+strainr(2,4,j,i))
        strain(j,i,3) = strain(j,i,3) + 0.25*dt*(strainr(3,1,j,i)+strainr(3,2,j,i)+strainr(3,3,j,i)+strainr(3,4,j,i))
! set shear heating rate equal to the plastic dissipation rate from rh_plastic.f90
        sshrheat(j,i) = diss
3 continue
!$OMP end do
!$OMP end parallel

    do i = 1,nx-1
        do j = 1,nz-1
            disstot = disstot + sshrheat(j,i)*dt
        enddo
    enddo

devmax = max(devmax, curr_devmax)
dvmax = max(dvmax, curr_dvmax)

return
end
