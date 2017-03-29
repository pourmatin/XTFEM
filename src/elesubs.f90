! Thu Mar 31 11:13:33 CDT 2016 
! Subroutines for elemental matrix operations
!


SUBROUTINE glbmtx(ne, nn, ndofn, nnse, nquads, ngpe, xconn, xcoord, &
    & kappa, cp, rho, wqs, ns, dns, ki, kj, kv, knz, mi, mj, mv, mnz)
! Assembly global matrices
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: ne, nn, ndofn, nnse, nquads, ngpe
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), DIMENSION(nn, ndofn), INTENT(IN) :: xcoord
REAL(KIND=REKIND), INTENT(IN) :: kappa, cp, rho
REAL(KIND=REKIND), DIMENSION(nquads), INTENT(IN) :: wqs
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, ndofn*nnse), INTENT(IN) :: ns
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, nnse), INTENT(IN) :: dns
INTEGER, INTENT(IN) :: knz
REAL(KIND=REKIND), DIMENSION(knz), INTENT(OUT) :: kv
INTEGER, DIMENSION(knz), INTENT(OUT) :: ki, kj
INTEGER, INTENT(IN) :: mnz
INTEGER, DIMENSION(mnz), INTENT(OUT) :: mi, mj
REAL(KIND=REKIND), DIMENSION(mnz), INTENT(OUT) :: mv
! ---Internal variables---
INTEGER :: i, j, ndofe
INTEGER, DIMENSION(nnse) :: nodes
INTEGER, DIMENSION(nnse*ndofn) :: dofs
REAL(KIND=REKIND), DIMENSION(nnse, ndofn) :: ncoorde
REAL(KIND=REKIND), DIMENSION(nnse*ndofn, nnse*ndofn) :: ke, me
ndofe = (nnse*ndofn)**2
!!$omp parallel cpm_threads(24)
!!$omp do private(nodes, ncoorde, dofs, ke, me)
DO i = 1, ne
    nodes = xconn(i,:)
    ncoorde = xcoord(nodes,:)
    DO j = 1, nnse
        dofs(j*3-2) = nodes(j)*3 - 2
        dofs(j*3-1) = nodes(j)*3 - 1
        dofs(j*3-0) = nodes(j)*3 - 0
    ENDDO
    CALL elemtx(ndofn, nnse, nquads, ngpe, ncoorde, rho, kappa, cp, wqs, ns, dns, ke, me)
    !CALL elenxyz(ndofn, nnse, nquads, ngpe, ncmp, ncoorde, sd, dns, dnxyz, &
    !    & Nx(:, (i-1)*ngpe+1:i*ngpe), Ny(:, (i-1)*ngpe+1:i*ngpe), &
    !    & Nz(:, (i-1)*ngpe+1:i*ngpe))
    !CALL elemmtx(ndofn, nnse, nquads, ngpe, ncoorde, rho, wqs, ns, dns, me)
    !CALL elekmtx(ndofn, nnse, nquads, ngpe, ncmp, ncoorde, sd, stiff, &
    !    & wqs, dns, dnxyz, ke)
    CALL assembly(ki((i-1)*ndofe+1:i*ndofe), kj((i-1)*ndofe+1:i*ndofe), & 
        & kv((i-1)*ndofe+1:i*ndofe), ke, dofs, nnse, ndofn)
    CALL assembly(mi((i-1)*ndofe+1:i*ndofe), mj((i-1)*ndofe+1:i*ndofe), &
        & mv((i-1)*ndofe+1:i*ndofe), me, dofs, nnse, ndofn)
ENDDO
!!$omp end do
!!$omp end parallel
RETURN
END SUBROUTINE glbmtx


SUBROUTINE elemtx(ndofn, nnse, nquads, ngpe, ncoorde, rho, kappa, cp, &
    & wqs, ns, dns, ke, me)
! Get spatial matrices
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: nnse, nquads, ngpe, ndofn
REAL(KIND=REKIND), INTENT(IN) :: rho, kappa, cp
REAL(KIND=REKIND), DIMENSION(nnse, 3), INTENT(IN) :: ncoorde
REAL(KIND=REKIND), DIMENSION(nquads), INTENT(IN) :: wqs
REAL(KIND=REKIND), DIMENSION(ngpe, nnse), INTENT(IN) :: ns
REAL(KIND=REKIND), DIMENSION(ngpe, nnse, 3), INTENT(IN) :: dns
REAL(KIND=REKIND), DIMENSION(nnse*ndofn, nnse*ndofn), INTENT(OUT) :: ke, me
! ---Internal variables---
INTEGER :: ii, jj, i, j, k, l, m, p
REAL(KIND=REKIND) :: jacdet
REAL(KIND=REKIND), DIMENSION(3) :: Niix, Njjx
REAL(KIND=REKIND), DIMENSION(ngpe, 3, 3) :: jac, gam
ke = 0._REKIND
me = 0._REKIND
jac = 0._REKIND
l = 0

DO i = 1, nquads
    DO j = 1, nquads
        DO k = 1, nquads
            l = l + 1
            DO m = 1, 3
                DO p = 1, 3
                    DO ii = 1, nnse
                        jac(l, m, p) = jac(l, m, p) + dns(l, ii, p) * ncoorde(ii, m)
                    ENDDO
                ENDDO
            ENDDO
            jacdet = determinant(jac(l, :, :))
            gam(l, :, :) = inverse(jac(l, :, :), 3)
            DO ii = 1, nnse
                DO jj = 1, nnse
                    me( ii, jj ) = me( ii, jj ) + rho*cp*ns(l, ii)*ns(l, jj)*jacdet*wqs(i)*wqs(j)*wqs(k)
                    DO m = 1, 3
                       Niix = 0
                       Njjx = 0
                        DO p = 1, 3
                            Niix(m) = Niix(m) + gam(l, p, m) * dns(l, ii, p)
                            Njjx(m) = Njjx(m) + gam(l, p, m) * dns(l, jj, p)
                        ENDDO
                        ke( ii, jj ) = ke( ii, jj ) + kappa*Niix(m)*Njjx(m)*jacdet*wqs(i)*wqs(j)*wqs(k)
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
RETURN
END SUBROUTINE elemtx


SUBROUTINE elemmtx(ndofn, nnse, nquads, ngpe, ncoorde, rho, wqs, ns, dns, me)
! Get spatial M matrix
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: ndofn, nnse, nquads, ngpe
REAL(KIND=REKIND), INTENT(IN) :: rho
REAL(KIND=REKIND), DIMENSION(nnse, ndofn), INTENT(IN) :: ncoorde
REAL(KIND=REKIND), DIMENSION(nquads), INTENT(IN) :: wqs
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, ndofn*nnse), INTENT(IN) :: ns
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, nnse), INTENT(IN) :: dns
REAL(KIND=REKIND), DIMENSION(nnse*ndofn, nnse*ndofn), INTENT(OUT) :: me
! ---Internal variables---
INTEGER :: i, j, k, l
REAL(KIND=REKIND) :: jacdet
REAL(KIND=REKIND), DIMENSION(nnse*ndofn, nnse*ndofn) :: met
REAL(KIND=REKIND), DIMENSION(ndofn, ndofn) :: jac
me = 0._REKIND
l = 0
DO i = 1, nquads
    DO j = 1, nquads
        DO k = 1, nquads
            l = l + 1
            jac = MATMUL(dns((l-1)*ndofn+1:l*ndofn, :), ncoorde)
            jacdet = determinant(jac)
            met = MATMUL(TRANSPOSE(ns((l-1)*ndofn+1:l*ndofn, :)), &
                &   (rho*ns((l-1)*ndofn+1:l*ndofn, :)))
            me = me + met*jacdet*wqs(i)*wqs(j)*wqs(k)
        ENDDO
    ENDDO
ENDDO
RETURN
END SUBROUTINE elemmtx


SUBROUTINE elekmtx(ndofn, nnse, nquads, ngpe, ncmp, ncoorde, sd, stiff, &
    & wqs, dns, dnxyz, ke)
! Get spatial K matrix
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: ndofn, nnse, nquads, ngpe, ncmp
REAL(KIND=REKIND), DIMENSION(nnse, ndofn), INTENT(IN) :: ncoorde
REAL(KIND=REKIND), DIMENSION(6,6), INTENT(IN) :: stiff
REAL(KIND=REKIND), DIMENSION(6,9), INTENT(IN) :: sd
REAL(KIND=REKIND), DIMENSION(nquads), INTENT(IN) :: wqs
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, nnse), INTENT(IN) :: dns
REAL(KIND=REKIND), DIMENSION(ndofn**2*ngpe, ndofn*nnse), INTENT(IN) :: dnxyz
REAL(KIND=REKIND), DIMENSION(nnse*ndofn, nnse*ndofn), INTENT(OUT) :: ke
! ---Internal variables---
INTEGER :: i, j, k, l
REAL(KIND=REKIND) :: jacdet
REAL(KIND=REKIND), DIMENSION(nnse*ndofn, nnse*ndofn) :: ket
REAL(KIND=REKIND), DIMENSION(ndofn, ndofn) :: jac, gam
REAL(KIND=REKIND), DIMENSION(ndofn**2, ndofn**2) :: gamexpand
REAL(KIND=REKIND), DIMENSION(ncmp, nnse*ndofn) :: bmat
ke = 0._REKIND
l = 0
DO i = 1, nquads
    DO j = 1, nquads
        DO k = 1, nquads
            l = l + 1
            jac = MATMUL(dns((l-1)*ndofn+1:l*ndofn, :), ncoorde)
            jacdet = determinant(jac)
            gam = inverse(jac, ndofn)
            gamexpand = 0._REKIND
            gamexpand(1:3,1:3) = gam
            gamexpand(4:6,4:6) = gam
            gamexpand(7:9,7:9) = gam
            bmat = MATMUL(MATMUL(sd, gamexpand), &
                &   dnxyz((l-1)*ndofn**2+1:l*ndofn**2, :))
            ket = MATMUL(TRANSPOSE(bmat),MATMUL(stiff,bmat))
            ke = ke + ket*jacdet*wqs(i)*wqs(j)*wqs(k)
        ENDDO
    ENDDO
ENDDO
RETURN
END SUBROUTINE elekmtx


SUBROUTINE elenxyz(ndofn, nnse, nquads, ngpe, ncmp, ncoorde, sd, dns, dnxyz, &
    & Nx, Ny, Nz)
! Get spatial Nx, Ny, Nz
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: ndofn, nnse, nquads, ngpe, ncmp
REAL(KIND=REKIND), DIMENSION(nnse, ndofn), INTENT(IN) :: ncoorde
REAL(KIND=REKIND), DIMENSION(6,9), INTENT(IN) :: sd
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, nnse), INTENT(IN) :: dns
REAL(KIND=REKIND), DIMENSION(ndofn**2*ngpe, ndofn*nnse), INTENT(IN) :: dnxyz
REAL(KIND=REKIND), DIMENSION(nnse, ngpe), INTENT(OUT) :: Nx, Ny, Nz
! ---Internal variables---
INTEGER :: i, j, k, l, m
REAL(KIND=REKIND) :: jacdet
REAL(KIND=REKIND), DIMENSION(ndofn, ndofn) :: jac, gam
REAL(KIND=REKIND), DIMENSION(ndofn**2, ndofn**2) :: gamexpand
REAL(KIND=REKIND), DIMENSION(ncmp, nnse*ndofn) :: bmat
l = 0
DO i = 1, nquads
    DO j = 1, nquads
        DO k = 1, nquads
            l = l + 1
            jac = MATMUL(dns((l-1)*ndofn+1:l*ndofn, :), ncoorde)
            jacdet = determinant(jac)
            gam = inverse(jac, ndofn)
            gamexpand = 0._REKIND
            gamexpand(1:3,1:3) = gam
            gamexpand(4:6,4:6) = gam
            gamexpand(7:9,7:9) = gam
            bmat = MATMUL(MATMUL(sd, gamexpand), &
                &   dnxyz((l-1)*ndofn**2+1:l*ndofn**2, :))
            DO m = 1, nnse
                Nx(m, l) = bmat(1, (m-1)*ndofn + 1)
                Ny(m, l) = bmat(2, (m-1)*ndofn + 2)
                Nz(m, l) = bmat(3, (m-1)*ndofn + 3)
            ENDDO
        ENDDO
    ENDDO
ENDDO
RETURN
END SUBROUTINE elenxyz


SUBROUTINE setshapefunc(nnsd, nnse, nquads, ngpe, wqs, ns, dns)
! Set spatial shape functions
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: nnsd, nnse, nquads, ngpe
REAL(KIND=REKIND), DIMENSION(nquads), INTENT(OUT) :: wqs
REAL(KIND=REKIND), DIMENSION(ngpe, nnse), INTENT(OUT) :: ns
REAL(KIND=REKIND), DIMENSION(ngpe, nnse, 3), INTENT(OUT) :: dns
! ---Internal variables---
INTEGER :: i, j, k, l
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: xiqs
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: tempns
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: tempdns
ALLOCATE(xiqs(nquads))
ALLOCATE(tempdns(nnse, 3))
ALLOCATE(tempns(nnse))
CALL gaussquadrature(nquads, xiqs, wqs)
l = 0 ! Initialize counter
DO i = 1, nquads ! Loop over quadrature points
    DO j = 1, nquads ! Loop over quadrature points
        DO k = 1, nquads ! Loop over quadrature points
            l = l + 1 ! Current Gauss point index
            CALL shapedx(tempdns,xiqs(i),xiqs(j),xiqs(k),nnsd)
            dns(l, :, :) = tempdns
            CALL shapex(tempns,xiqs(i),xiqs(j),xiqs(k),nnsd)
            ns(l, :) = tempns
        ENDDO
    ENDDO
ENDDO
DEALLOCATE(xiqs, tempdns, tempns)
RETURN
END SUBROUTINE setshapefunc


SUBROUTINE setsdstiff(kappa, cp, sd, stiff)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
REAL(KIND=REKIND), INTENT(IN) :: kappa, cp
REAL(KIND=REKIND), DIMENSION(6,6), INTENT(OUT) :: stiff
REAL(KIND=REKIND), DIMENSION(6,9), INTENT(OUT) :: sd
stiff = kappa/((1._REKIND+cp)*(1._REKIND-2._REKIND*cp))* &
    & RESHAPE([1._REKIND-cp, cp, cp, 0._REKIND, 0._REKIND, 0._REKIND, &
        & cp, 1._REKIND-cp, cp, 0._REKIND, 0._REKIND, 0._REKIND, &
        & cp, cp, 1._REKIND-cp, 0._REKIND, 0._REKIND, 0._REKIND, &
        & 0._REKIND,0._REKIND,0._REKIND,0.5_REKIND-cp,0._REKIND,0._REKIND, &
        & 0._REKIND,0._REKIND,0._REKIND,0._REKIND,0.5_REKIND-cp,0._REKIND, &
        & 0._REKIND,0._REKIND,0._REKIND,0._REKIND,0._REKIND,0.5_REKIND-cp], &
        & [6, 6])
sd = RESHAPE([1._REKIND,0._REKIND,0._REKIND,0._REKIND,0._REKIND,0._REKIND, &
            & 0._REKIND,0._REKIND,0._REKIND,1._REKIND,0._REKIND,0._REKIND, &
            & 0._REKIND,0._REKIND,0._REKIND,0._REKIND,0._REKIND,1._REKIND, &
            & 0._REKIND,0._REKIND,0._REKIND,1._REKIND,0._REKIND,0._REKIND, &
            & 0._REKIND,1._REKIND,0._REKIND,0._REKIND,0._REKIND,0._REKIND, &
            & 0._REKIND,0._REKIND,0._REKIND,0._REKIND,1._REKIND,0._REKIND, &
            & 0._REKIND,0._REKIND,0._REKIND,0._REKIND,0._REKIND,1._REKIND, &
            & 0._REKIND,0._REKIND,0._REKIND,0._REKIND,1._REKIND,0._REKIND, &
            & 0._REKIND,0._REKIND,1._REKIND,0._REKIND,0._REKIND,0._REKIND], &
            & [6, 9])
RETURN
END SUBROUTINE setsdstiff
