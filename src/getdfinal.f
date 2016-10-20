c------------------------------------------------------------------------
c----------------------- Two-Scale Damage Model--------------------------
c Inputs: Total global strain tensor at current time-step and previous time-step,
c History variables: Plastic strain, effective plastic strain, damage,
c                    stress state at last time-step,stored energy, backstress
c                    equivalent stress at last time-step, increment in effective_plastic_strain. 
c Outputs: Updated history variables and Damage !
c Author: Sagar Bhamare (1st December 2011)      
c------------------------------------------------------------------------      
      subroutine getdfinal(strain_tgo,strain_tg,strain_pl,p,D,ws,sig_eff
     .           ,backsig,del_po,sig_d,cmat,dmat)
      
      implicit real*8 (a-h,o-z)
      real*8 nu
      dimension strain_tg(6),strain_tgo(6),dstrain_tg(6)
      dimension strain_l(6),strain_el(6),strain_pl(6)
      dimension dstrain_tl(6),dstrain_el(6),dstrain_pl(6)
      dimension test_strain_l(6),test_strain_el(6)
      dimension dmat(6,6),cmat(6,6)
      dimension test_sig_eff(6),test_dev_sig_eff(6),backsig(6)
      dimension testsig(6)
      dimension dsig_g(6),Qs(6),dev_Qs(6),sigs(6),dev_sigs(6) 
      dimension sig_eff(6),dev_sig_eff(6)
      dimension sig_eff_p(3),sig_eff_n(3)
      
ccccccccccc Material Parameters cccccccccccccccccccccccc      
      E=1.97d5       ! MPa, Mg/mm^3, mm, sec
      nu=0.3d0
      dc=0.3d0
      h=0.2d0
      cy=1740.d0
      bigs=0.5d0      !S in Lemaitre's equation
      s=0.5d0         !s in Lemaitre's eq uation
      sig_f=180.0d0
      sig_u=577.00d0
      eps_D=0.08d0
!      pi=4.d0*asin(1.d0/sqrt(2.d0))
cccccc Define constitutive matrix cccccccccccccccccccccc      
      G = E/(2.d0*(1.d0+nu))
      dK = E/(3.d0*(1.d0-2.d0*nu)) ! using dK to replace K to avoid bugs
!      lambda = dK-((2.d0/3.d0)*G)

!      dmat(1,1)=lambda+2.d0*G;
!      dmat(1,2)=lambda;
!      dmat(1,3)=lambda;
!      dmat(1,4)=0.0;
!      dmat(1,5)=0.0;
!      dmat(1,6)=0.0;

!      dmat(2,1) = lambda;
!      dmat(2,2) = lambda+2.d0*G;
!      dmat(2,3) = lambda;
!      dmat(2,4) = 0.0;
!      dmat(2,5) = 0.0;
!      dmat(2,6) = 0.0;

!      dmat(3,1) = lambda;
!      dmat(3,2) = lambda;
!      dmat(3,3) = lambda+2.d0*G;
!      dmat(3,4) = 0.0;
!      dmat(3,5) = 0.0;
!      dmat(3,6) = 0.0;

!      dmat(4,1) = 0.0;
!      dmat(4,2) = 0.0;
!      dmat(4,3) = 0.0;
!      dmat(4,4) = G;
!      dmat(4,5) = 0.0;
!      dmat(4,6) = 0.0;

!      dmat(5,1) = 0.0;
!      dmat(5,2) = 0.0;
!      dmat(5,3) = 0.0;
!      dmat(5,4) = 0.0;
!      dmat(5,5) = G;
!      dmat(5,6) = 0.0;

!      dmat(6,1) = 0.0;
!      dmat(6,2) = 0.0;
!      dmat(6,3) = 0.0;
!      dmat(6,4) = 0.0;
!      dmat(6,5) = 0.0;
!      dmat(6,6) = G;
cccccc Define compliance matrix cccccccccccccccccccccc       
!      cmat(1,1)=1/E;
!      cmat(1,2)=-nu/E;
!      cmat(1,3)=-nu/E;
!      cmat(1,4)=0.0;
!      cmat(1,5)=0.0;
!      cmat(1,6)=0.0;

!      cmat(2,1) = -nu/E
!      cmat(2,2) = 1/E
!      cmat(2,3) = -nu/E
!      cmat(2,4) = 0.0
!      cmat(2,5) = 0.0
!      cmat(2,6) = 0.0

!      cmat(3,1) = -nu/E
!      cmat(3,2) = -nu/E
!      cmat(3,3) = 1/E
!      cmat(3,4) = 0.0;
!      cmat(3,5) = 0.0;
!      cmat(3,6) = 0.0;

!      cmat(4,1) = 0.0;
!      cmat(4,2) = 0.0;
!      cmat(4,3) = 0.0;
!      cmat(4,4) = 1/G;
!      cmat(4,5) = 0.0;
!      cmat(4,6) = 0.0;

!      cmat(5,1) = 0.0;
!      cmat(5,2) = 0.0;
!      cmat(5,3) = 0.0;
!      cmat(5,4) = 0.0;
!      cmat(5,5) = 1/G;
!      cmat(5,6) = 0.0;

!      cmat(6,1) = 0.0;
!      cmat(6,2) = 0.0;
!      cmat(6,3) = 0.0;
!      cmat(6,4) = 0.0;
!      cmat(6,5) = 0.0;
!      cmat(6,6) = 1/G;
cccccccccc Localization law parameters cccccccccccccccc
      a=(1.d0+nu)/(3.d0*(1.d0-nu))
      b=2.d0*(4.d0-5.d0*nu)/(15.d0*(1.d0-nu))     
cccccccccc Initialization of all arrays ccccccccccccccc
      do i=1,6
         strain_l(i)=0.d0      ! total local strain
         strain_el(i)=0.d0     ! total elastic strain
      enddo
c      
      del_p=0.d0
      del_D=0.d0  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccc Start the Timestep Loop ccccccccccccccccccccccc
 
c--------------Assign the Loading (only for stand-Alone code)--------
c Strain history is fully-reversed sine function--------------------
c     
c         
c------------Elastic Predictor stage---------------------------------
      C1=1.d0/(1.d0-b*D)                     ! eq(40) in two_scale_model report
      C2=(a-b)*D/((1.d0-b*D)*(1.d0-a*D))                                        ! Eqn. 148
      C3=b*(1.d0-D)/(1.d0-b*D)

c     assume elastic increment to calcualte current total elastic strain at microscale
c     test_strain_l= total test strain, and test_strain_el= elastic test strain  
c      
      test_strain_l=0.d0      ! total test local strain
      test_strain_el=0.d0     ! test local elastic strain   
      do j=1,6
         test_strain_l(j)=C1*strain_tg(j)+C3*strain_pl(j) ! eq(40)
      enddo   
c      
      trace_strain_tg=0.d0  ! determine trace for 2nd term in eq(40)
      do j=1,3
         trace_strain_tg=trace_strain_tg+strain_tg(j)   
      enddo
      trace_strain_tg=(C2*trace_strain_tg)/3.d0  ! eq(40)
c
      do j=1,3
        test_strain_l(j)=test_strain_l(j)+trace_strain_tg
      enddo 
c
      do j=1,6
         test_strain_el(j)=test_strain_l(j)-strain_pl(j) ! eq(41)
      enddo                          
c     calcualte test stress based on test elastic strain 

      test_sig_eff=0.d0                                ! eq(42)
      do j=1,6
         do k=1,6
           test_sig_eff(j)=test_sig_eff(j)+dmat(j,k)*test_strain_el(k)
         enddo
      enddo
c
      trace_test_sig_eff=0.d0
      do j=1,3
         trace_test_sig_eff=trace_test_sig_eff+test_sig_eff(j)   
      enddo
      trace_test_sig_eff=trace_test_sig_eff/3.d0 
c
      test_dev_sig_eff=0.d0
      test_dev_sig_eff(1)=test_sig_eff(1)-trace_test_sig_eff
      test_dev_sig_eff(2)=test_sig_eff(2)-trace_test_sig_eff
      test_dev_sig_eff(3)=test_sig_eff(3)-trace_test_sig_eff
      test_dev_sig_eff(4)=test_sig_eff(4)
      test_dev_sig_eff(5)=test_sig_eff(5)
      test_dev_sig_eff(6)=test_sig_eff(6)
c     
      testsig=0.d0 
      do j=1,6
         testsig(j)=test_dev_sig_eff(j)-backsig(j)
      enddo
c
      testsig_eq=1.5*(testsig(1)**2+testsig(2)**2+testsig(3)**2+
     .    2*(testsig(4)**2+testsig(5)**2+testsig(6)**2))
      testsig_eq=sqrt(testsig_eq)           
c
      if (testsig_eq .le. sig_f) then 
         sig_eff=test_sig_eff
         backsig=backsig
         D=D
         strain_l=test_strain_l        ! total strain at the end of increment
         strain_el=test_strain_el      ! total elastic -------//-------------
         strain_pl=strain_pl
         strain_tgo=strain_tg         ! strain_tgo is needed in next steps
         sig_d=testsig_eq             ! sig_d used for wd calculation when material becomes plastic
      else
c------------------ Plastic Corrector Phase -----------------------
c Form the Qs Matrix first
      dstrain_tg=0.d0                 ! increment in total global strain
      dstrain_tg=strain_tg-strain_tgo ! delta_Epsilon
c      
      dsig_g=0.d0
      do j=1,6
         do k=1,6
            dsig_g(j)=dsig_g(j)+dmat(j,k)*dstrain_tg(k)
         enddo     
      enddo

      do j=1,6
         dsig_g(j)=dsig_g(j)/(1.d0-b*D)                                         ! Eqn. 163
      enddo

c
      C4=dK*(a-b)*D/((1.d0-b*D)*(1.d0-a*D))                                      ! Eqn. 163, note that K here is incorrect
      trace_dstrain_tg=0.d0
      do j=1,3
         trace_dstrain_tg=trace_dstrain_tg+dstrain_tg(j)                        ! Eqn. 163
      enddo  ! Divison by 3 not needed, check eq(55)

c
      do j=1,3
         dsig_g(j)=dsig_g(j)+(C4*trace_dstrain_tg)                              ! Eqn. 163
      enddo

      Qs=0.d0
      do j=1,6
         Qs(j)=backsig(j)-sig_eff(j)-dsig_g(j)                                  ! Eqn. 163
      enddo

c Form the eta' in the delta_p expression
      eta=(cy*(1.d0-D))+(3.d0*G*((1.d0-b)/(1.d0-b*D)))                          ! Eqn. 162
c
      trace_Qs=0.d0
      do j=1,3
         trace_Qs=trace_Qs+Qs(j) 
      enddo            
      trace_Qs=trace_Qs/3.d0
c     
      dev_Qs=0.d0
      dev_Qs(1)=Qs(1)-trace_Qs
      dev_Qs(2)=Qs(2)-trace_Qs
      dev_Qs(3)=Qs(3)-trace_Qs
      dev_Qs(4)=Qs(4)
      dev_Qs(5)=Qs(5)
      dev_Qs(6)=Qs(6)
c
      Qs_eq=1.5*(dev_Qs(1)**2+dev_Qs(2)**2+dev_Qs(3)**2+
     .    2*(dev_Qs(4)**2+dev_Qs(5)**2+dev_Qs(6)**2))
      Qs_eq=sqrt(Qs_eq)
c
      trace_sigs=-trace_Qs
c
      del_p=(Qs_eq-sig_f)/eta                                                   ! Eqn. 170
      p=p+del_p
c 
      dev_sigs=0.d0   
      C5=1.d0+(eta*del_p/sig_f)                                                 ! Eqn. 171
      do j=1,6
         dev_sigs(j)=-dev_Qs(j)/c5                                              ! Eqn. 171
      enddo
c   
      sigs=0.d0                                                                 ! Eqn. 172
      sigs(1)=dev_sigs(1)+trace_sigs
      sigs(2)=dev_sigs(2)+trace_sigs
      sigs(3)=dev_sigs(3)+trace_sigs
      sigs(4)=dev_sigs(4)
      sigs(5)=dev_sigs(5)
      sigs(6)=dev_sigs(6)

c------------- Update the Variables --------------------
      dstrain_pl=0.d0
      do j=1,6
         dstrain_pl(j)=(3.d0*dev_sigs(j)*del_p/(2.d0*sig_f))                    ! Eqn. 173
      enddo
c
      strain_pl=strain_pl+dstrain_pl
c
      do j=1,6
         backsig(j)=backsig(j)+(2.d0*cy*(1.d0-D)*dstrain_pl(j)/3.d0)            ! Eqn. 174
      enddo            
c
c---------------------------------------------------------------------
c     calculate the effective stress at n+1 stage
      sig_eff=sigs+backsig                                                      ! Eqn. 175

c     sig_eff_eq is used for stroed energy calculation
      trace_sig_eff=0.d0
      do j=1,3
         trace_sig_eff=trace_sig_eff+sig_eff(j) 
      enddo  
      trace_sig_eff=trace_sig_eff/3.d0
c     
      dev_sig_eff=0.d0
      dev_sig_eff(1)=sig_eff(1)-trace_sig_eff
      dev_sig_eff(2)=sig_eff(2)-trace_sig_eff
      dev_sig_eff(3)=sig_eff(3)-trace_sig_eff
      dev_sig_eff(4)=sig_eff(4)
      dev_sig_eff(5)=sig_eff(5)
      dev_sig_eff(6)=sig_eff(6)
c
      sig_eff_eq=1.5*(dev_sig_eff(1)**2+dev_sig_eff(2)**2+
     .           dev_sig_eff(3)**2+2*(dev_sig_eff(4)**2+
     .           dev_sig_eff(5)**2+dev_sig_eff(6)**2))      
      sig_eff_eq=sqrt(sig_eff_eq)
c     calculate the total strain and use for plotting
      strain_el=0.d0                 
      do j=1,6
         do k=1,6
            strain_el(j)=strain_el(j)+cmat(j,k)*sig_eff(k)                      ! Eqn. 176
         enddo
      enddo
c
      strain_l=strain_el+strain_pl         
c---------------- update the damage--------------------------------------
c----------------- Calculate stored energy density to predict damage initiation---
      wd=(sig_u-sig_f)*eps_D                                                    ! Eqn. 147
      sig_new=sig_eff_eq-sig_f
      sig_old=sig_d-sig_f

      if (sig_old .lt. 0) then
         sig_old=abs(sig_old)
      endif
      if (sig_new .lt. 0) then
         sig_new=abs(sig_new)
      endif

      ws=ws+(((sig_new*del_p)+(sig_old*del_po))/2.d0)                           
c aleternate formulation for ws
c      w_new=0.d0
c      w_new=abs(backsig(1)*dstrain_pl(1))+abs(backsig(2)*dstrain_pl(2))+
c     .  abs(backsig(3)*dstrain_pl(3))+2.d0*abs(backsig(4)*dstrain_pl(4))
c     .     +2.d0*abs(backsig(5)*dstrain_pl(5))+
c     .      2.d0*abs(backsig(6)*dstrain_pl(6))
c      ws=ws+((w_old+w_new)/2.d0)
c      w_old=w_new ! use this energy in the next step     
c      write(11,*) 'ws', ws, 'wd', wd
c      
      if (ws .ge. wd) then     
c ------- decompose effective stress tensor in postivie and negative part
      sig_eff_p=0.d0
      sig_eff_n=0.d0
      do j=1,3
         if (sig_eff(j) .lt. 0.d0)then
             sig_eff_n(j)=sig_eff(j)
         else
             sig_eff_p(j)=sig_eff(j)
         endif
      enddo      
c------------ evaluate Y and Damage Variable -------------------------
      Y_P=((1.d0+nu)/(2.d0*E))*((sig_eff_p(1)**2+sig_eff_p(2)**2+               ! Eqn. 178
     .    sig_eff_p(3)**2)+h*((1.d0-D)/(1.d0-h*D))**2*(sig_eff_n(1)**2+
     .    sig_eff_n(2)**2+sig_eff_n(3)**2))
c
      trace_sig_eff=0.d0
      do j=1,3
         trace_sig_eff=trace_sig_eff+sig_eff(j) 
      enddo            
      if (trace_sig_eff .ge. 0.d0) then                                         ! Eqn. 178
         Y_n=(nu/(2.d0*E))*(trace_sig_eff**2)
      else
         Y_n=(nu/(2.d0*E))*(h*((1.d0-D)/(1.d0-h*D))**2*trace_sig_eff**2)
      endif
c      
      Y=Y_p-Y_n                                                                 ! Eqn. 178
c     
      del_D=(Y/bigs)**s*del_p                                                   ! Eqn. 177
      D=del_D+D
      else
      D=0.d0
      endif

c-------------- Update parameters at the end of current timestep
      strain_tgo=strain_tg  ! needed in next timestep for calculating del_tg 
      sig_d=sig_eff_eq      ! used for wd calculation in next step
      del_po=del_p          ! used for ws calculation in next step
      endif
      
ccccccccccccccccccc Endo of the Timestep Loop ccccccccccccccccccccc      
      return
      end
