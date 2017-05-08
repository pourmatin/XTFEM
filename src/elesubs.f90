
SUBROUTINE elemtx_m(ndofn, nnse, nquads, ngpe, ncmp, ncoorde, rho, sd, stiff, &
    & wqs, ns, dns, dnxyz, Nx, Ny, Nz, ke, me)
! Get spatial matrices and Nx, Ny, Nz
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: ndofn, nnse, nquads, ngpe, ncmp
REAL(KIND=REKIND), INTENT(IN) :: rho
REAL(KIND=REKIND), DIMENSION(nnse, ndofn), INTENT(IN) :: ncoorde
REAL(KIND=REKIND), DIMENSION(6,6,8), INTENT(IN) :: stiff
REAL(KIND=REKIND), DIMENSION(6,9), INTENT(IN) :: sd
REAL(KIND=REKIND), DIMENSION(nquads), INTENT(IN) :: wqs
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, ndofn*nnse), INTENT(IN) :: ns
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, nnse), INTENT(IN) :: dns
REAL(KIND=REKIND), DIMENSION(ndofn**2*ngpe, ndofn*nnse), INTENT(IN) :: dnxyz
REAL(KIND=REKIND), DIMENSION(nnse, ngpe), INTENT(OUT) :: Nx, Ny, Nz
REAL(KIND=REKIND), DIMENSION(nnse*ndofn, nnse*ndofn), INTENT(OUT) :: ke, me
! ---Internal variables---
INTEGER :: i, j, k, l, m
REAL(KIND=REKIND) :: jacdet
REAL(KIND=REKIND), DIMENSION(ndofn, ndofn) :: jac, gam
REAL(KIND=REKIND), DIMENSION(ndofn**2, ndofn**2) :: gamexpand
REAL(KIND=REKIND), DIMENSION(ncmp, nnse*ndofn) :: bmat
REAL(KIND=REKIND), DIMENSION(nnse*ndofn, nnse*ndofn) :: ket, met
ke = 0.
me = 0.
l = 0
DO i = 1, nquads
    DO j = 1, nquads
        DO k = 1, nquads
            l = l + 1
            jac = MATMUL(dns((l-1)*ndofn+1:l*ndofn, :), ncoorde)
            jacdet = determinant(jac)
            gam = inverse(jac, ndofn)
            gamexpand = 0.
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
            ket = MATMUL(TRANSPOSE(bmat),MATMUL(stiff(:,:,l),bmat))
            ke = ke + ket*jacdet*wqs(i)*wqs(j)*wqs(k)
            met = MATMUL(TRANSPOSE(ns((l-1)*ndofn+1:l*ndofn, :)), &
                &   (rho*ns((l-1)*ndofn+1:l*ndofn, :)))
            me = me + met*jacdet*wqs(i)*wqs(j)*wqs(k)
        ENDDO
    ENDDO
ENDDO
RETURN
END SUBROUTINE elemtx_m

SUBROUTINE elemtx_t(ndofn, nnse, nquads, ngpe, ncoorde, rho, kappa, cp, &
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
REAL(KIND=REKIND), DIMENSION(3, 3) :: jac, gam
ke = 0.
me = 0.
l = 0

DO i = 1, nquads
    DO j = 1, nquads
        DO k = 1, nquads
            l = l + 1
            jac = 0.
            DO m = 1, 3
                DO p = 1, 3
                    DO ii = 1, nnse
                        jac(m, p) = jac(m, p) + dns(l, ii, p) * ncoorde(ii, m)
                    ENDDO
                ENDDO
            ENDDO
            jacdet = determinant(jac)
            gam = inverse(jac, 3)
            DO ii = 1, nnse
                DO jj = 1, nnse
                    me( ii, jj ) = me( ii, jj ) + rho*cp*ns(l, ii)*ns(l, jj)*jacdet*wqs(i)*wqs(j)*wqs(k)
                    DO m = 1, 3
                       Niix = 0
                       Njjx = 0
                        DO p = 1, 3
                            Niix(m) = Niix(m) + gam(p, m) * dns(l, ii, p)
                            Njjx(m) = Njjx(m) + gam(p, m) * dns(l, jj, p)
                        ENDDO
                        ke( ii, jj ) = ke( ii, jj ) + kappa*Niix(m)*Njjx(m)*jacdet*wqs(i)*wqs(j)*wqs(k)
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
RETURN
END SUBROUTINE elemtx_t


SUBROUTINE elevec_m(ndofn, nnse, nquads, ngpe, ncoorde, wqs, stress, ns, dns, ve)
! Get spatial matrices
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: nnse, nquads, ngpe, ndofn
REAL(KIND=REKIND), DIMENSION(nnse, 3), INTENT(IN) :: ncoorde
REAL(KIND=REKIND), DIMENSION(nquads), INTENT(IN) :: wqs
REAL(KIND=REKIND), DIMENSION(ndofn,ngpe), INTENT(IN) :: stress
REAL(KIND=REKIND), DIMENSION(ngpe, nnse), INTENT(IN) :: ns
REAL(KIND=REKIND), DIMENSION(ngpe, nnse, 3), INTENT(IN) :: dns
REAL(KIND=REKIND), DIMENSION(nnse*ndofn), INTENT(OUT) :: ve
! ---Internal variables---
INTEGER :: ii, jj, i, j, k, l, m, p
REAL(KIND=REKIND) :: jacdet, Niix
REAL(KIND=REKIND), DIMENSION(3, 3) :: jac, gam
l = 0
ve = 0.
DO i = 1, nquads
    DO j = 1, nquads
        DO k = 1, nquads
            l = l + 1
            jac = 0.
            DO m = 1, 3
                DO p = 1, 3
                    DO ii = 1, nnse
                        jac(m, p) = jac(m, p) + dns(l, ii, p) * ncoorde(ii, m)
                    ENDDO
                ENDDO
            ENDDO
            jacdet = determinant(jac)
            gam = inverse(jac, 3)
            DO ii = 1, nnse
                DO m = 1, 3
                    Niix = 0
                    DO p = 1, 3
                        Niix = Niix + gam(p, m) * dns(l, ii, p)
                    ENDDO
                    ve((ii-1)*ndofn+m) = ve((ii-1)*ndofn+m) + &
                            & Niix*stress(m,l)*jacdet*wqs(i)*wqs(j)*wqs(k)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
RETURN
END SUBROUTINE elevec_m


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
me = 0.
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
ke = 0.
l = 0
DO i = 1, nquads
    DO j = 1, nquads
        DO k = 1, nquads
            l = l + 1
            jac = MATMUL(dns((l-1)*ndofn+1:l*ndofn, :), ncoorde)
            jacdet = determinant(jac)
            gam = inverse(jac, ndofn)
            gamexpand = 0.
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
            gamexpand = 0.
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

SUBROUTINE setshapefunc_m(ndofn, nnsd, nnse, nquads, ngpe, wqs, ns, dns, dnxyz)
! Set spatial shape functions
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: ndofn, nnsd, nnse, nquads, ngpe
REAL(KIND=REKIND), DIMENSION(nquads), INTENT(OUT) :: wqs
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, ndofn*nnse), INTENT(OUT) :: ns
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, nnse), INTENT(OUT) :: dns
REAL(KIND=REKIND), DIMENSION(ndofn**2*ngpe, ndofn*nnse), INTENT(OUT) :: dnxyz
! ---Internal variables---
INTEGER :: i, j, k, l
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: xiqs
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: tempdns, tempdnxyz, tempns3
ALLOCATE(xiqs(nquads))
ALLOCATE(tempdns(ndofn, nnse))
ALLOCATE(tempdnxyz(ndofn**2, ndofn*nnse))
ALLOCATE(tempns3(ndofn, ndofn*nnse))
CALL gaussquadrature(nquads, xiqs, wqs)
l = 0 ! Initialize counter
DO i = 1, nquads ! Loop over quadrature points
    DO j = 1, nquads ! Loop over quadrature points
        DO k = 1, nquads ! Loop over quadrature points
            l = l + 1 ! Current Gauss point index
            CALL shape3dx(tempdns,tempdnxyz,xiqs(i),xiqs(j),xiqs(k),nnsd)
            dns((l-1)*ndofn+1:l*ndofn, :) = tempdns
            dnxyz((l-1)*ndofn**2+1:l*ndofn**2, :) = tempdnxyz
            CALL shape3x(tempns3,xiqs(i),xiqs(j),xiqs(k),nnsd) 
            ns((l-1)*ndofn+1:l*ndofn, :) = tempns3
        ENDDO
    ENDDO
ENDDO
DEALLOCATE(xiqs, tempdns, tempdnxyz, tempns3)
RETURN
END SUBROUTINE setshapefunc_m

SUBROUTINE setshapefunc_t(nnsd, nnse, nquads, ngpe, wqs, ns, dns)
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
END SUBROUTINE setshapefunc_t


SUBROUTINE setsdstiff(ym, nu, sd, stiff)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
REAL(KIND=REKIND), INTENT(IN) :: ym, nu
REAL(KIND=REKIND), DIMENSION(6,6), INTENT(OUT) :: stiff
REAL(KIND=REKIND), DIMENSION(6,9), INTENT(OUT) :: sd
stiff = ym/((1.+nu)*(1.-2.*nu))* &
    & RESHAPE([1.-nu, nu, nu, 0._REKIND, 0._REKIND, 0._REKIND, &
        & nu, 1.-nu, nu, 0._REKIND, 0._REKIND, 0._REKIND, &
        & nu, nu, 1.-nu, 0._REKIND, 0._REKIND, 0._REKIND, &
        & 0._REKIND,0._REKIND,0._REKIND,0.5-nu,0._REKIND,0._REKIND, &
        & 0._REKIND,0._REKIND,0._REKIND,0._REKIND,0.5-nu,0._REKIND, &
        & 0._REKIND,0._REKIND,0._REKIND,0._REKIND,0._REKIND,0.5-nu], &
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
