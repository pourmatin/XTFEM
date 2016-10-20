! Thu Mar 31 11:13:33 CDT 2016 
! Subroutines for elemental matrix operations
!


SUBROUTINE glbmtx(ne, nn, ndofn, nnse, nquads, ngpe, ncmp, xconn, xcoord, &
    & ym, nu, rho, wqs, ns, dns, dnxyz, sd, stiff, Nx, Ny, Nz, ki, kj, kv, &
    & knz, mi, mj, mv, mnz)
! Assembly global matrices
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: ne, nn, ndofn, nnse, nquads, ngpe, ncmp
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), DIMENSION(nn, ndofn), INTENT(IN) :: xcoord
REAL(KIND=REKIND), INTENT(IN) :: ym, nu, rho
REAL(KIND=REKIND), DIMENSION(nquads), INTENT(IN) :: wqs
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, ndofn*nnse), INTENT(IN) :: ns
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, nnse), INTENT(IN) :: dns
REAL(KIND=REKIND), DIMENSION(ndofn**2*ngpe, ndofn*nnse), INTENT(IN) :: dnxyz
REAL(KIND=REKIND), DIMENSION(6,9), INTENT(OUT) :: sd
REAL(KIND=REKIND), DIMENSION(6,6), INTENT(OUT) :: stiff
REAL(KIND=REKIND), DIMENSION(nnse, ngpe*ne), INTENT(OUT) :: Nx, Ny, Nz
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
CALL setsdstiff(ym, nu, sd, stiff)
!!$omp parallel num_threads(24)
!!$omp do private(nodes, ncoorde, dofs, ke, me)
DO i = 1, ne
    nodes = xconn(i,:)
    ncoorde = xcoord(nodes,:)
    DO j = 1, nnse
        dofs(j*3-2) = nodes(j)*3 - 2
        dofs(j*3-1) = nodes(j)*3 - 1
        dofs(j*3-0) = nodes(j)*3 - 0
    ENDDO
    CALL elemtx(ndofn, nnse, nquads, ngpe, ncmp, ncoorde, rho, sd, stiff, &
        & wqs, ns, dns, dnxyz, Nx(:, (i-1)*ngpe+1:i*ngpe), &
        & Ny(:, (i-1)*ngpe+1:i*ngpe), Nz(:, (i-1)*ngpe+1:i*ngpe), ke, me)
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


SUBROUTINE elemtx(ndofn, nnse, nquads, ngpe, ncmp, ncoorde, rho, sd, stiff, &
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
REAL(KIND=REKIND), DIMENSION(6,6), INTENT(IN) :: stiff
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
ke = 0._REKIND
me = 0._REKIND
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
            ket = MATMUL(TRANSPOSE(bmat),MATMUL(stiff,bmat))
            ke = ke + ket*jacdet*wqs(i)*wqs(j)*wqs(k)
            met = MATMUL(TRANSPOSE(ns((l-1)*ndofn+1:l*ndofn, :)), &
                &   (rho*ns((l-1)*ndofn+1:l*ndofn, :)))
            me = me + met*jacdet*wqs(i)*wqs(j)*wqs(k)
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


SUBROUTINE setshapefunc(ndofn, nnsd, nnse, nquads, ngpe, wqs, ns, dns, dnxyz)
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
END SUBROUTINE setshapefunc


SUBROUTINE setsdstiff(ym, nu, sd, stiff)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
REAL(KIND=REKIND), INTENT(IN) :: ym, nu
REAL(KIND=REKIND), DIMENSION(6,6), INTENT(OUT) :: stiff
REAL(KIND=REKIND), DIMENSION(6,9), INTENT(OUT) :: sd
stiff = ym/((1._REKIND+nu)*(1._REKIND-2._REKIND*nu))* &
    & RESHAPE([1._REKIND-nu, nu, nu, 0._REKIND, 0._REKIND, 0._REKIND, &
        & nu, 1._REKIND-nu, nu, 0._REKIND, 0._REKIND, 0._REKIND, &
        & nu, nu, 1._REKIND-nu, 0._REKIND, 0._REKIND, 0._REKIND, &
        & 0._REKIND,0._REKIND,0._REKIND,0.5_REKIND-nu,0._REKIND,0._REKIND, &
        & 0._REKIND,0._REKIND,0._REKIND,0._REKIND,0.5_REKIND-nu,0._REKIND, &
        & 0._REKIND,0._REKIND,0._REKIND,0._REKIND,0._REKIND,0.5_REKIND-nu], &
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
