      subroutine stress_eval(hnx,hny,disp,nstep,nmax,nndem,gstress,
     . nd,ns,nnd,gstrain,nele,nde,w,dt,igst,elesta,ngpstatus,
     . npl,xcnd,t1,iflag,nip,nipc)   
     
c This subroutine evaluates the stress at gauss points

      implicit real*8(a-h,o-z)
      dimension disp(12*nnd),hnx(ns,nd),hny(ns,nd),nndem(ns,nd),c(3,3)
      dimension gstrain(4*(nip+1)*nele,3),gstress(4*(nip+1)*nele,3)
      dimension bmat(3,8)
      dimension dstrain(8),ngtemp(ns),gton(ns,ns)
      dimension store(4*(nip+1)*nele,4)
      dimension nde(nele,ns)
      dimension hnex(6,nip+1),dispinter(2*nnd,nip+1),tt(nip+1)
      dimension ngpstatus(4*nele)
      integer   dnodes, ndo
      dimension dnodes(4)
c
      dimension xcnd(nnd,2),xcoord(nnd,2)
      character*30 filename 
      integer elesta(nele) 
      
c      
c      write(*,*) nstep,nmax,nd,ns,nnd,nele,w,dt,igst,nodeo
c      open(14,file='stress.dat',status='unknown')
      nch=4
c Stress= material_matrix*B_matrix*Displacement_vector      
c Define C matrix 
      youngs=1.97d5
      poisson=0.3d0
c     plane strain
c      c3(3)=youngs/2.d0/(1+posson)
c      c1(1)=c3(3)*2.d0*(1-posson)/(1-2.d0*posson)
c      c2(2)=c1(1)
c      c1(2)=c3(3)*2.d0*posson/(1-2.d0*posson)
c      c2(1)=c1(2)
c      c1(3)=0.0d0
c      c2(3)=0.0d0
c      c3(1)=0.0d0
c      c3(2)=0.0d0
c     plane stress
       cc=youngs/(1.0d0-poisson)/(1.0d0+poisson)
       c(1,1)=cc
       c(1,2)=cc*poisson
       c(1,3)=0.0d0
       c(2,1)=c(1,2)
       c(2,2)=c(1,1)
       c(2,3)=0.0d0
       c(3,1)=0.0d0
       c(3,2)=0.0d0
       c(3,3)=cc*(1.0d0-poisson)/2.0d0      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Interpolate the displacement in time direction For all the Nodes
c for given 'nex' points. 
c 'dispinter' = Spatial Nodes*2 (x and y), nex columns 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
c       t1=(nstep-1)*dt
       t2=t1+0.5d0*dt
       t3=t1+dt
       hnex=0.d0
       dispinter=0.d0
       tt=0.d0
       dstrain=0.d0
c       write(28,*) t1,t2,t3
       if (nstep.eq.nmax) then
          ien=nip+1
       else
          ien=nip
       endif
c       
       do i=1,ien
          tnex=t1+(i-1)*(dt/nip)
          hnex(1,i)=2.d0*(t2-tnex)*(t3-tnex)/dt**2
          hnex(2,i)=-4.d0*(t1-tnex)*(t3-tnex)/dt**2
          hnex(3,i)=2.d0*(t1-tnex)*(t2-tnex)/dt**2
          hnex(4,i)=hnex(1,i)*(dsin(w*tnex)-dsin(w*t1))
          hnex(5,i)=hnex(2,i)*(dsin(w*tnex)-dsin(w*t2))
          hnex(6,i)=hnex(3,i)*(dsin(w*tnex)-dsin(w*t3))
       enddo
c       
        do i=1,ien
           tt(i)=t1+(i-1)*(dt/nip)
c           print *, i, tt(i)
           do k=1,2*nnd
              do j=1,6
                 dispinter(k,i)=dispinter(k,i)+hnex(j,i)*
     &                       disp(2*(j-1)*nnd+k)
              enddo
           enddo
        enddo

!        open(100,FILE='disp',position='append')
!        do i = 1, ien
!            write(100,*) tt(i), dispinter(11,i), dispinter(20+nnd,i)
!        enddo
!        close(100)
                
     
c        write(18,*) 'ngfail', (ngfail(i),i=1,4)
c--------------------------------------------------------
c Evaluate strains and stresses at all the gauss points
c--------------------------------------------------------      
      gstrain=0.d0  ! intialize at each time step else superposition
      gstress=0.d0 
      do i=1,ien   ! Loop over the slabs 
       do j=1,igst  ! Loop over the all GPs in one slab  
c     
            if (ngpstatus(j) .eq. 0) then
                 do l=1,3
                    gstrain((i-1)*igst+j,l)=0.d0
                    gstress((i-1)*igst+j,l)=0.d0
                 enddo
C                 write(*,*) 'working'            
            else
            bmat=0.d0
c extract derivatives to form B Matrix         
            do l=1,2*nch ! loop over the number of columns of bmat
               if (l.lt.nch+1) then
                   bmat(1,l)=hnx(l,j)
                   bmat(3,l)=hny(l,j)
               elseif (l.gt.nch) then
                   bmat(2,l)=hny((l-nch),j)
                   bmat(3,l)=hnx((l-nch),j)
               endif
            enddo  ! End of loop over the number of columns of bmat
c Extract the node numbers for the current element related to gauss Point            
            do k=1,nch
               dnodes(k)=nndem(k,j)
            enddo

c Extract the displacements in x and y direction corresponding to the nodes above            
c dstrain is nothing but displacements not strains          
           do k=1,nch
              dstrain(k)=dispinter(dnodes(k),i)
              dstrain(k+nch)=dispinter(nnd+dnodes(k),i)
           enddo
c Calculate strain at the gauss point
           do l=1,3        ! loop over strain componenets
              do m=1,2*nch ! loop over the nodal displacements
              gstrain((i-1)*igst+j,l)=gstrain((i-1)*igst+j,l)+bmat(l,m)*
     &                                dstrain(m)
              enddo
           enddo
c Calculate stresses at the gauss point  
            do l=1,3        ! loop over stress componenets
              do m=1,3 ! loop over strain componenets
              gstress((i-1)*igst+j,l)=gstress((i-1)*igst+j,l)+c(l,m)*
     &                                gstrain((i-1)*igst+j,m)
              enddo
           enddo
c          write(*,*) 'working 1' 
c          

            endif
         enddo  ! End of loop on GPs in one slab
      enddo ! End of loop on the number of slabs for currnet dt
      !write(*,*) 'time =', time2 - time1 ! 0.92s / 1.24s
c---------------------------------------------------------------
c  Extrapolate stresses to the nodes 
c---------------------------------------------------------------
         
      gton(1,1)=1.d0+(sqrt(3.d0)/2.d0)
      gton(1,2)=-(1.d0)/(2.d0)
      gton(1,3)=1.d0-(sqrt(3.d0)/2.d0)
      gton(1,4)=-(1.d0)/(2.d0)
      
      gton(2,1)=-(1.d0)/(2.d0)
      gton(2,2)=1.d0+(sqrt(3.d0)/2.d0)
      gton(2,3)=-(1.d0)/(2.d0)
      gton(2,4)=1.d0-(sqrt(3.d0)/2.d0)
      
      gton(3,1)=1.d0-(sqrt(3.d0)/2.d0)
      gton(3,2)=-(1.d0)/(2.d0)
      gton(3,3)=1.d0+(sqrt(3.d0)/2.d0)
      gton(3,4)=-(1.d0)/(2.d0)
      
      gton(4,1)=-(1.d0)/(2.d0)
      gton(4,2)=1.d0-(sqrt(3.d0)/2.d0)
      gton(4,3)=-(1.d0)/(2.d0)
      gton(4,4)=1.d0+(sqrt(3.d0)/2.d0)
      
      do k=1,ien ! loop over number of slabs
         do i=1,nele    ! loop over number of elements heightwise
c ngtemp are the gauss point numbers for the given element         
         ngtemp(1)=(k-1)*igst+(4*i-3)  
         ngtemp(2)=(k-1)*igst+(4*i-2)
         ngtemp(3)=(k-1)*igst+(4*i-1)
         ngtemp(4)=(k-1)*igst+(4*i)
         
         do l=1,nch  ! loop for the nodes (rows)
            do m=1,3  ! loop for the stress components at node 'l'
               do j=1,4
c               print *, l, m, j, (k-1)*nnd+nde(i,l)
c               print *, gton(l,j),gstress(ngtemp(j),m) 
                  store((k-1)*nnd+nde(i,l),m)=store((k-1)*nnd
     &            +nde(i,l),m)+gton(l,j)*gstress(ngtemp(j),m)
              enddo !j=4
            enddo ! m=3
            store((k-1)*nnd+nde(i,l),4)=store((k-1)*nnd+nde(i,l),4)+1
         enddo  ! l=nch   
         enddo ! end of loop over number of elements
      enddo ! end of loop over number of slabs          
      
c  average the stresses at each node depending on its weight stored
c  as the 4th column in the stress vector for each node

      do k=1,ien
         do i=1,nnd
           do m=1,3
               store((k-1)*nnd+i,m)=store((k-1)*nnd+i,m)/store((k-1)*
     &         nnd+i,4)
            enddo
         enddo
      enddo
c
c      do i=1,ien
c            write(51,370) tt(i),(store((i-1)*nnd+nodeo,j),j=1,3)
c      enddo
C      if (nstep.eq.1) then
C      do i=1,ien
C            write(16,370) tt(i),(gstress((i-1)*igst+2400,j),j=1,3),
C     .      (gstress((i-1)*igst+2399,j),j=1,3),
C     .      (gstress((i-1)*igst+2398,j),j=1,3),
C     .      (gstress((i-1)*igst+2397,j),j=1,3)
C      enddo
C      endif
 370  format(13(2x,E20.11))      

c--------------------------- Write the Tecplot Output -----------------
      if (nstep.eq.nmax) then   ! nwrite decides the number of fiels
          nwrite=npl+1
      else
          nwrite=npl
      endif
c
      nfile=(npl)*(nstep-1)+1;   ! index for writing files
      
c write 'nt-1' files for 'nstep <nmax' and 'nt' files for nstep=nmax

      if (iflag.eq.1) then
      do i=1,nwrite   ! Loop on the number of files 

         xcoord=0.d0; ! declare current xcoord=zero.

c  for x_current= xold+ disp((1:nnd)*slab_number)
c  for y_current= yold+disp(nt*nnd+(1:nnd)*slab_number)
c          if ((nstep.eq.nmax).and.(i.eq.nwrite)) then
c             nc=ien-1
c          else
c             nc=(i-1)*(nipc/npl)+1
c          endif
c  Modified version for writing the output files
          if ((nstep.eq.nmax).and.(i.eq.nwrite)) then
             nc=ien
          else
             nc=(i-1)*((nip+nipc)/npl) ! nipc here can lead to different phase
             if (i.eq.1) then
                 nc=1
             endif
          endif
          do j=1,nnd
             xcoord(j,1)=xcnd(j,1)!+1.d0*dispinter(j,nc) 
             xcoord(j,2)=xcnd(j,2)!+1.d0*dispinter(nnd+j,nc)
          enddo
         
C         write(*,*) 'nc', nc
         filename="def0000000.plt"
         write(filename(4:10),'(i7.7)') nfile
         open(20,FILE=filename,status='unknown')
         
         write(20,*)  ' VARIABLES = "X","Y","U","V","S11","S22","S12",
     &   "Mises"'
         write(20,'(1x,a,i6,a,i6,a)')
     &         'ZONE N =',nnd,', E =',nele-nef,
     &         ', F = FEPOINT, ET = QUADRILATERAL'

         do 11 k=1,nnd
            write(20,350) (xcoord(k,n),n=1,2),dispinter(k,nc),
     &      dispinter(nnd+k,nc),(store((nc-1)*nnd+k,n),n=1,3),
     &      sqrt(store((nc-1)*nnd+k,1)**2+store((nc-1)*nnd+k,2)**2
     &      -store((nc-1)*nnd+k,1)*store((nc-1)*nnd+k,2))
11       continue
         do 12 l=1,nele
            if (elesta(l).eq.0) then
C               write(*,*) 
            else   
                write(20,360) (nde(l,n),n=1,4)
            endif
 12      continue
         close(20)
         
         
 350     format(8(2x,E20.11)) 
 360     format(1x,4i10) 
         nfile=nfile+1
      enddo      ! End of loop for number of files
         nc = nip-nipc/4*3 + 1
c         print *,nc, tt(nc)
c        write output file for abaqus 
         OPEN(20,FILE='sent-304l-002.dat',position='append')
c        Time
c         DO nc = 1, ien
         WRITE(20, '(E20.11)') tt(nc)
c        Displacements
         DO k = 1, nnd
            WRITE(20,410) dispinter(k,nc), dispinter(nnd+k,nc)
         ENDDO
c        Strain components on integration points
         DO k = 1, igst
            WRITE(20,430) (gstrain((nc-1)*igst+k,n),n=1,3)
         ENDDO
c        Stress components on integration points
         DO k = 1, igst
            WRITE(20,430) (gstress((nc-1)*igst+k,n),n=1,3)
         ENDDO
c         ENDDO
         CLOSE(20)                  
410      format(2(E20.11,2x))
420      format(4i10,1x) 
430      format(3(E20.11,2x))
      endif  ! iflag condition
c      

      return 
      end 
