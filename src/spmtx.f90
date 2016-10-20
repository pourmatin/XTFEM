! Timestamp: Sun May 03 15:30:24 CDT 2015
!_______________________________________________________________________________
! 
! Subourines for sparse matrix operations in space-time finite element method
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! 
! ---Subroutines---
!    Name               Description
! ------------------  ------------------------------------------------------
!  assembly             assemble global matrix as raw triplet format         #
!  stassembly           assemble space-time matrix                           &
!  csrapplybc           apply BCs to CSR matrix
!  tricsr               triplet raw assembled matrix -> standard CSR format  *
!  rm0csr               remove 0s in CSR matrix
!  spgemv               perform one of the matrix-vector operations, like gemv
!  cs_compress          convert triplet format to CSR format
!  cs_cumsum            computes the cumulative sum
!  cs_transpose         transpose a CSR matrix
!  cs_dupl              sum up duplicate entries in CSR matrix
!  cs_add               perform addtion of two CSR matrix
!  cs_scatter           scatter a sparse vector into a dense one
!  cs_permute           permutes a CSR matrix
! 
! ---User guide--- 
!  1: How to assemble a spatial global matrix?
!     Step 1, CALL #, assembling the global matrix element-by-element
!     Step 2, CALL *, converting the global triplet matrix to CSR matrix
!  2: How to assemble a space-time matrix?
!     CALL stassembly
!  3: How to apply BCs to a CSR matrix?
!     CALL csrapplybc
!_______________________________________________________________________________

SUBROUTINE assembly0(gi, gj, gv, elemat, dof, nele, npe, dpn, nz)
! Assembles elemental matrix into global matrix in triplet (i,j,v) format
! 
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! 
!  ---Entries---
!    Name     Type      Size                  Description
!  -------- --------- --------------------- --------------------------------
!   gi       integer   nele*(npe*dpn)^2      rows of triplet global matrix
!   gj       integer   nele*(npe*dpn)^2      cols of triplet global matrix
!   gv       real*8    nele*(npe*dpn)^2      vals of triplet global matrix
!   elemat   real*8    [npe*dpn]*[npe*dpn]   elemental matrix
!   dof      integer   npe*dpn               elemental dofs index list
!   nele     integer   1                     total number of elements
!   npe      integer   1                     nodes per element
!   dpn      integer   1                     DOFs per node
!   nz       integer   1                     number of nonzero entries
! 
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: nele, npe, dpn
INTEGER, DIMENSION(nele*(npe*dpn)**2), INTENT(INOUT) :: gi, gj
REAL(kind = 8), DIMENSION(nele*(npe*dpn)**2), INTENT(INOUT) :: gv
REAL(kind = 8), DIMENSION(npe*dpn, npe*dpn), INTENT(IN) :: elemat
INTEGER, DIMENSION(npe*dpn), INTENT(IN) :: dof
INTEGER, INTENT(INOUT) :: nz
! ---Internal variables---
INTEGER :: i, j
! Assembling global matrix
DO i = 1, npe*dpn ! rows
    DO j = 1, npe*dpn ! columns
        nz = nz + 1
        gi(nz) = dof(i)
        gj(nz) = dof(j)
        gv(nz) = elemat(i,j)
    ENDDO
ENDDO
RETURN
END SUBROUTINE assembly0


SUBROUTINE assembly(gi, gj, gv, elemat, dof, npe, dpn)
! Assembles elemental matrix into global matrix in triplet (i,j,v) format
! 
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Mar., 2016
! 
!  ---Entries---
!    Name     Type      Size                  Description
!  -------- --------- --------------------- --------------------------------
!   gi       integer   nele*(npe*dpn)^2      rows of triplet global matrix
!   gj       integer   nele*(npe*dpn)^2      cols of triplet global matrix
!   gv       real*8    nele*(npe*dpn)^2      vals of triplet global matrix
!   elemat   real*8    [npe*dpn]*[npe*dpn]   elemental matrix
!   dof      integer   npe*dpn               elemental dofs index list
!   npe      integer   1                     nodes per element
!   dpn      integer   1                     DOFs per node
! 
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: npe, dpn
INTEGER, DIMENSION((npe*dpn)**2), INTENT(INOUT) :: gi, gj
REAL(kind = 8), DIMENSION((npe*dpn)**2), INTENT(INOUT) :: gv
REAL(kind = 8), DIMENSION(npe*dpn, npe*dpn), INTENT(IN) :: elemat
INTEGER, DIMENSION(npe*dpn), INTENT(IN) :: dof
! ---Internal variables---
INTEGER :: i, j, k
k = 1
! Assembling global matrix
DO i = 1, npe*dpn ! rows
    DO j = 1, npe*dpn ! columns
        gi(k) = dof(i)
        gj(k) = dof(j)
        gv(k) = elemat(i,j)
        k = k + 1
    ENDDO
ENDDO
RETURN
END SUBROUTINE assembly


SUBROUTINE stassembly(ist, jst, st, stn, stnz, ik, jk, k, knz, im, jm, m, &
                &      mnz, kmn, ck, cm, cn)
! Assembles space-time matrix in CSR format with openmp
! 
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), May, 2015
! 
!  ---Entries---
!    Name     Type      Size         Description
!  -------- --------- ------------ -----------------------------------------
!   ist      integer   stn+1        rows of space-time matrix
!   jst      integer   stnz         cols of space-time matrix
!   st       real*8    stnz         vals of space-time matrix
!   stn      integer   1            order of space-time matrix
!   stnz     integer   1            number of nonzero entries in s-t matrix
!   ik       integer   kmn+1        rows of spatial K matrix
!   jk       integer   knz          cols of spatial K matrix
!   k        real*8    knz          vals of spatial K matrix
!   knz      integer   1            number of nonzero entries in K matrix
!   im       integer   kmn+1        rows of spatial M matrix
!   jm       integer   knz          cols of spatial M matrix
!   m        real*8    knz          vals of spatial M matrix
!   mnz      integer   1            number of nonzero entries in M matrix
!   kmn      integer   1            order of spatial K or M matrix
!   ck       real*8    [cn]*[cn]    coeffecients for K in s-t matrix
!   cm       real*8    [cn]*[cn]    coeffecients for M in s-t matrix
!   cn       integer   1            order of ck or cm
! 
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: mnz
INTEGER, INTENT(IN) :: knz
INTEGER, INTENT(IN) :: kmn 
INTEGER, INTENT(IN) :: stn
INTEGER, INTENT(INOUT) :: stnz
INTEGER, INTENT(IN) :: cn
INTEGER, DIMENSION(stn+1), INTENT(INOUT) :: ist
INTEGER, DIMENSION(stnz), INTENT(INOUT) :: jst
REAL(kind = 8), DIMENSION(stnz), INTENT(INOUT) :: st
INTEGER, DIMENSION(kmn+1), INTENT(IN) :: ik
INTEGER, DIMENSION(knz), INTENT(IN) :: jk
REAL(kind = 8), DIMENSION(knz), INTENT(IN) :: k
INTEGER, DIMENSION(kmn+1), INTENT(IN) :: im
INTEGER, DIMENSION(mnz), INTENT(IN) :: jm
REAL(kind = 8), DIMENSION(mnz), INTENT(IN) :: m
REAL(kind = 8), DIMENSION(cn, cn), INTENT(IN) :: ck, cm
! ---Internal variables---
INTEGER, DIMENSION(kmn+1) :: ic
INTEGER, DIMENSION(knz) :: jc
INTEGER, DIMENSION(stnz) :: ti, tj
REAL(kind = 8), DIMENSION(stnz) :: tv
REAL(kind = 8), DIMENSION(knz) :: c
INTEGER :: i, j, ii, jj, cnz, nz
nz = 0
DO i = 1, cn
    DO j = 1, cn    
        ic = ik
        jc = jk
        c = k 
        cnz = knz
        CALL cs_add(ic,jc,c,cnz, im,jm,m,mnz, kmn, ck(i,j), cm(i,j))
        DO ii = 1, kmn
            DO jj = ic(ii), ic(ii+1)-1
                nz = nz + 1
                ti(nz) = ii + (i-1)*kmn
                tj(nz) = jc(jj) + (j-1)*kmn
                tv(nz) = c(jj)
            ENDDO
        ENDDO
    ENDDO
ENDDO
CALL tricsr(ti, tj, tv, stn, nz)
ist(1:stn+1) = ti(1:stn+1)
jst(1:nz) = tj(1:nz)
st(1:nz) = tv(1:nz)
stnz = nz
RETURN
END SUBROUTINE stassembly


!SUBROUTINE stassembly(ist, jst, st, stn, stnz, ik, jk, k, knz, im, jm, m, &
!                &      mnz, kmn, ck, cm, cn)
! Assembles space-time matrix in CSR format
! 
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! 
!  ---Entries---
!    Name     Type      Size         Description
!  -------- --------- ------------ -----------------------------------------
!   ist      integer   stn+1        rows of space-time matrix
!   jst      integer   stnz         cols of space-time matrix
!   st       real*8    stnz         vals of space-time matrix
!   stn      integer   1            order of space-time matrix
!   stnz     integer   1            number of nonzero entries in s-t matrix
!   ik       integer   kmn+1        rows of spatial K matrix
!   jk       integer   knz          cols of spatial K matrix
!   k        real*8    knz          vals of spatial K matrix
!   knz      integer   1            number of nonzero entries in K matrix
!   im       integer   kmn+1        rows of spatial M matrix
!   jm       integer   knz          cols of spatial M matrix
!   m        real*8    knz          vals of spatial M matrix
!   mnz      integer   1            number of nonzero entries in M matrix
!   kmn      integer   1            order of spatial K or M matrix
!   ck       real*8    [cn]*[cn]    coeffecients for K in s-t matrix
!   cm       real*8    [cn]*[cn]    coeffecients for M in s-t matrix
!   cn       integer   1            order of ck or cm
! 
!IMPLICIT NONE
! ---External variables--- 
!INTEGER, DIMENSION(stn+1), INTENT(INOUT) :: ist
!INTEGER, DIMENSION(stnz), INTENT(INOUT) :: jst
!REAL(kind = 8), DIMENSION(stnz), INTENT(INOUT) :: st
!INTEGER, INTENT(IN) :: stn
!INTEGER, INTENT(INOUT) :: stnz
!INTEGER, DIMENSION(kmn+1), INTENT(IN) :: ik
!INTEGER, DIMENSION(knz), INTENT(IN) :: jk
!REAL(kind = 8), DIMENSION(knz), INTENT(IN) :: k
!INTEGER, INTENT(IN) :: knz
!INTEGER, DIMENSION(kmn+1), INTENT(IN) :: im
!INTEGER, DIMENSION(mnz), INTENT(IN) :: jm
!REAL(kind = 8), DIMENSION(mnz), INTENT(IN) :: m
!INTEGER, INTENT(IN) :: mnz
!INTEGER, INTENT(IN) :: kmn
!REAL(kind = 8), DIMENSION(cn, cn), INTENT(IN) :: ck, cm
!INTEGER, INTENT(IN) :: cn
! ---Internal variables---
!INTEGER, DIMENSION(stn+1) :: ic
!INTEGER, DIMENSION(stnz) :: jc
!REAL(kind = 8), DIMENSION(stnz) :: c
!INTEGER :: ii, jj, nz, rows, cols
!REAL(kind = 8) :: cf
! Assemble the space-time matrix
! First epxands the spatial matrix, then performs addition
!nz = stnz
!DO ii = 1, cn
!    DO jj = 1, cn
!        rows = (ii-1)*kmn
!        cols = (jj-1)*kmn
!        cf = ck(ii,jj)
!        CALL csrexpand(ic,jc,c,stn,nz, ik,jk,k,kmn,knz, rows, cols)
!        CALL cs_add(ist,jst,st,stnz, ic,jc,c,knz, stn, 1.0D+00, cf)
!        stnz = nz
!        cf = cm(ii,jj)
!        CALL csrexpand(ic,jc,c,stn,nz, im,jm,m,kmn,mnz, rows, cols)
!        CALL cs_add(ist,jst,st,stnz, ic,jc,c,mnz, stn, 1.0D+00, cf)
!        stnz = nz
!    ENDDO
!ENDDO
!stnz = ist(stn+1)-1
!RETURN
!END SUBROUTINE stassembly


SUBROUTINE csrapplybc(ia, ja, a, n, nz, bc)
! Applies boundary condtion to system matrix by 1-0 method
! 
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! 
!  ---Entries---
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------------
!   ia       integer   n+1      rows of CSR matrix
!   ja       integer   nz       cols of CSR matrix
!   a        real*8    nz       vals of CSR matrix
!   n        integer   1        order of CSR matrix
!   nz       integer   1        number of nonzero entries in CSR matrix
!   bc       logical   n        boundary conditions
! 
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(INOUT) :: nz 
INTEGER, DIMENSION(n+1), INTENT(INOUT) :: ia
INTEGER, DIMENSION(nz), INTENT(INOUT) :: ja
REAL(kind = 8), DIMENSION(nz), INTENT(INOUT) :: a
INTEGER, DIMENSION(n), INTENT(IN) :: bc
! ---Internal variables---
INTEGER :: i, j, k
DO i = 1, n
    DO k = ia(i), ia(i+1)-1
        j = ja(k)
        IF ((bc(i) .NE. 0).OR.(bc(j) .NE. 0)) THEN
            a(k) = 0.0D+00 ! clear rows and cols
            IF (i.EQ.j) THEN
                a(k) = 1.0D+00 ! set diagonal to one
            ENDIF
        ENDIF
    ENDDO
ENDDO
! Removing zero entries
CALL rm0csr(ia, ja, a, n, nz)
RETURN
END SUBROUTINE csrapplybc


!SUBROUTINE csrexpand(ia,ja,a,an,anz, ib,jb,b,bn,bnz, rows, cols)
! Expands and shifts a CSR matrix to a higher order matrix
! 
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! 
!  ---Entries---
!    Name     Type      Size         Description
!  -------- --------- ------------ -----------------------------------------
!   ia       integer   an+1         rows
!   ja       integer   anz          cols
!   a        real*8    anz          vals
!   an       integer   1            order
!   anz      integer   1            number of nonzero entries
!   ib       integer   bn+1         rows
!   jb       integer   bnz          cols
!   b        real*8    bnz          vals
!   bn       integer   1            order
!   bnz      integer   1            number of nonzero entries 
!   rows     integer   1            shift rows
!   cols     integer   1            shift cols
! 
!IMPLICIT NONE
! ---External variables--- 
!INTEGER, DIMENSION(an+1), INTENT(INOUT) :: ia
!INTEGER, DIMENSION(anz), INTENT(INOUT) :: ja
!REAL(kind = 8), DIMENSION(anz), INTENT(INOUT) :: a
!INTEGER, INTENT(IN) :: an
!INTEGER, INTENT(INOUT) :: anz
!INTEGER, DIMENSION(bn+1), INTENT(IN) :: ib
!INTEGER, DIMENSION(bnz), INTENT(IN) :: jb
!REAL(kind = 8), DIMENSION(bnz), INTENT(IN) :: b
!INTEGER, INTENT(IN) :: bn
!INTEGER, INTENT(IN) :: bnz
!INTEGER, INTENT(IN) :: rows, cols
! ---Internal variables---
!INTEGER :: ii, jj, nz
! Initialization
!ia = 1
!ja = 0
!a = 0.0D+00
! Copy data
!ia(rows+1:rows+bn) = ib(1:bn)
!ia(rows+bn+1:an+1) = ib(bn+1)
!ja(1:bnz) = jb(1:bnz) + cols
!a(1:bnz) = b(1:bnz)
!RETURN
!END SUBROUTINE csrexpand


SUBROUTINE tricsr(ti, tj, tv, n, nz, nti)
! Processes the triplet global matrix to standard CSR format by
! commpressing, sorting, summing up duplicate entries, removing zero entries
! 
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! 
!  ---Entries---
! 
!  On INPUT
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------
!   ti       integer   nti      rows of triplet matrix
!   tj       integer   nz       cols of triplet matrix
!   tv       real*8    nz       vals of triplet matrix
!   n        integer   1        order of the full matrix
!   nz       integer   1        length of the triplet matrix
!   nti      integer   1        actual size of ti, = MAX(nz, n+1)
! 
!  On OUTPUT
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------
!   ti       integer   n+1      row pointers of CSR                 *
!   tj       integer   nz       cols of CSR
!   tv       real*8    nz       vals of CSR
!   nz       integer   nz       number of nonzero entries
! 
!  * The actural Size(ti) = nz remains the same, but entries will be 
!    set to zero when index > n+1.  
! 
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: n, nti
INTEGER, INTENT(INOUT) :: nz
INTEGER, DIMENSION(nti), INTENT(INOUT) :: ti
INTEGER, DIMENSION(nz), INTENT(INOUT) :: tj
REAL(kind = 8), DIMENSION(nz), INTENT(INOUT) :: tv
! ---Internal variables---
! N/A
! Compressing to CSR format
 CALL cs_compress(ti, tj, tv, n, nz, nti)
! Sorting rows and columns by transposing twice
 CALL cs_transpose(ti(1:n+1), tj, tv, n, nz) ! 1st
 CALL cs_transpose(ti(1:n+1), tj, tv, n, nz) ! 2nd
! Summing up duplicate entries
 CALL cs_dupl(ti(1:n+1), tj, tv, n, nz)
! Removing zero entries
 CALL rm0csr(ti(1:n+1), tj, tv, n, nz)
RETURN
END SUBROUTINE tricsr


SUBROUTINE rm0csr(ia, ja, a, n, nz)
! Removes zero entries in a CSR matrix.
! 
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! 
!  ---Entries---
! 
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------
!   ia       integer   n+1      row pointers of CSR
!   ja       integer   nz       cols of CSR
!   a        real*8    nz       vals of CSR
!   n        integer   1        order of the full matrix
!   nz       integer   1        number of nonzero entries
! 
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(INOUT) :: nz
INTEGER, DIMENSION(n+1), INTENT(INOUT) :: ia
INTEGER, DIMENSION(nz), INTENT(INOUT) :: ja
REAL(kind = 8), DIMENSION(nz), INTENT(INOUT) :: a
! ---Internal variables---
INTEGER :: i, j, p, q
INTEGER, DIMENSION(n+1) :: rp
INTEGER, DIMENSION(nz) :: cols
REAL(kind = 8), DIMENSION(nz) :: vals
p = 1
rp(1) = 1
DO i = 1, n
    DO j = ia(i), ia(i+1) - 1
        IF (a(j)/=0.0D+00) THEN ! nonzero entries
            vals(p) = a(j)
            cols(p) = ja(j)
            p = p + 1
        ENDIF
    ENDDO
    rp(i+1) = p
ENDDO
nz = p - 1
ia = rp
ja(1:nz) = cols
a(1:nz) = vals
RETURN
END SUBROUTINE rm0csr


SUBROUTINE spgemv(ia, ja, a, n, nz, x, y, alpha, beta)
! Performs on of the matrix-vector operations, similiar to LAPACK/DGEMV
!   y := alpha*A*x + beta*y
! 
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! 
!  ---Entries---
! 
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------    
!   ia       integer   n+1      row pointers of CSR
!   ja       integer   nz       cols of CSR
!   a        real*8    nz       vals of CSR
!   n        integer   1        order of the full matrix
!   nz       integer   1        number of nonzero entries
!   x        real*8    n        
!   y        real*8    n   
!   alpha    real*8    1        
!   beta     real*8    1   
! 
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: n, nz
INTEGER, DIMENSION(n+1), INTENT(IN) :: ia
INTEGER, DIMENSION(nz), INTENT(IN) :: ja
REAL(kind = 8), DIMENSION(nz), INTENT(IN) :: a
REAL(kind = 8), DIMENSION(n), INTENT(IN) :: x
REAL(kind = 8), DIMENSION(n), INTENT(INOUT) :: y
REAL(kind = 8) :: alpha, beta
! ---Internal variables---
INTEGER :: i, k
DO i = 1, n
    DO k = ia(i), ia(i+1)-1
        y(i) = alpha*a(k)*x(ja(k)) + beta*y(i)
    ENDDO
ENDDO
RETURN
END SUBROUTINE spgemv


SUBROUTINE cs_compress(ti, tj, tv, n, nz, nti)
! Converts the triplet matrix into a CSR format
! 
! Translated by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! from its C++ version in CSparse (Dr. Tim Davis, TAMU), and
! changed from CSC format to CSR format
! 
!  ---Entries---
! 
!  On INPUT
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------
!   ti       integer   nti      rows of triplet matrix
!   tj       integer   nz       cols of triplet matrix
!   tv       real*8    nz       vals of triplet matrix
!   n        integer   1        order of the full matrix
!   nz       integer   1        length of the triplet matrix
!   nti      integer   1        actual size of ti, = MAX(nz, n+1)
! 
!  On OUTPUT
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------
!   ti       integer   n+1      row pointers of CSR                 *
!   tj       integer   nz       cols of CSR
!   tv       real*8    nz       vals of CSR
! 
!  * The actural Size(ti) = nz remains the same, but entries will be 
!    set to zero when index > n+1.  
! 
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: n, nz, nti
INTEGER, DIMENSION(nti), INTENT(INOUT) :: ti
INTEGER, DIMENSION(nz), INTENT(INOUT) :: tj
REAL(kind = 8), DIMENSION(nz), INTENT(INOUT) :: tv
! ---Internal variables---
INTEGER :: k, p
INTEGER, DIMENSION(n+1) :: ia ! row pointer
INTEGER, DIMENSION(nz) :: ja ! cols
REAL(kind = 8), DIMENSION(nz) :: a ! values
INTEGER, DIMENSION(n) :: w ! work array
! Initializing work array
w = 0
DO k = 1, nz
! Counting numbers of entries in a row
    w(ti(k)) = w(ti(k)) + 1
ENDDO
! Calculating row pointer
 CALL cs_cumsum(ia, w, n)
! Calculating cols and vals
DO k = 1, nz
    p = w(ti(k))
    w(ti(k)) = w(ti(k)) + 1
    ja(p) = tj(k)
    a(p) = tv(k)
ENDDO
! Restoring to ti, tj, tv
tj = ja
tv = a
ti(1:n+1) = ia
RETURN
END SUBROUTINE cs_compress


SUBROUTINE cs_cumsum(p, c, n)
! Computes the cumulative sum
! Sets p(i) equal to the sum of c(1) through c(i), and p(n+1) equal to
! sum of c(1) through c(n) plus 1. On output, c(1..n) is overwritten
! with p(1...n).
! 
! Translated by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! from its C++ version in CSparse (Dr. Tim Davis, TAMU)
! 
!  ---Entries---
! 
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------
!   p        integer   n+1      N/A
!   c        integer   n        N/A
!   n        integer   1        N/A
! 
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: n 
INTEGER, DIMENSION(n+1), INTENT(OUT) :: p
INTEGER, DIMENSION(n), INTENT(INOUT) :: c
! ---Internal variables---
INTEGER :: i, nz
nz = 1
DO i = 1, n
    p(i) = nz
    nz = nz + c(i)
    c(i) = p(i)
ENDDO
p(n+1) = nz
RETURN
END SUBROUTINE cs_cumsum


SUBROUTINE cs_transpose(ia, ja, a, n, nz)
! Transposes a CSR sparse matrix
! This is also a linear-time bucker sort algorithm, CSR matrix can be
! sorted by transposing it twice. 
! 
! Translated by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! from its C++ version in CSparse (Dr. Tim Davis, TAMU)
! 
!  ---Entries---
! 
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------
!   ia       integer   n+1      row pointers of CSR
!   ja       integer   nz       cols of CSR
!   a        real*8    nz       vals of CSR
!   n        integer   1        order of the full matrix
!   nz       integer   1        number of nonzero entries 
! 
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: n, nz 
INTEGER, DIMENSION(n+1), INTENT(INOUT) :: ia
INTEGER, DIMENSION(nz), INTENT(INOUT) :: ja
REAL(kind = 8), DIMENSION(nz), INTENT(INOUT) :: a
! ---Internal variables---
INTEGER :: i, j, k, p, q
INTEGER, DIMENSION(n+1) :: rp
INTEGER, DIMENSION(nz) :: cols
REAL(kind = 8), DIMENSION(nz) :: vals
INTEGER, DIMENSION(n) :: w ! work array
! Initializing work array
w = 0
DO k = 1, ia(n+1)-1
! Counting numbers of entries in a row
    w(ja(k)) = w(ja(k)) + 1
ENDDO
! Calculating row pointer
 CALL cs_cumsum(rp, w, n)
! Calculating cols and vals
DO j = 1, n
    DO p = ia(j), ia(j+1)-1
        q = w(ja(p))
        w(ja(p)) = w(ja(p)) + 1
        cols(q) = j
        vals(q) = a(p)
    ENDDO
ENDDO
! Restoring to ia, ja, a
ia = rp
ja = cols
a = vals
RETURN
END SUBROUTINE cs_transpose


SUBROUTINE cs_dupl(ia, ja, a, n, nz)
! Sums up duplicate entries in global matrix (CSR format)
! 
! Translated by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! from its C++ version in CSparse (Dr. Tim Davis, TAMU)
! 
!  ---Entries---
! 
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------
!   ia       integer   n+1      row pointers of CSR
!   ja       integer   nz       cols of CSR
!   a        real*8    nz       vals of CSR
!   n        integer   1        order of the full matrix
!   nz       integer   1        number of nonzero entries
! 
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(INOUT) :: nz
INTEGER, DIMENSION(n+1), INTENT(INOUT) :: ia
INTEGER, DIMENSION(nz), INTENT(INOUT) :: ja
REAL(kind = 8), DIMENSION(nz), INTENT(INOUT) :: a
! ---Internal variables---
INTEGER :: i, j, p, q, nnz
INTEGER, DIMENSION(n) :: w ! work array
! Initializing work array
w = 0
nnz = 1
DO j = 1, n
    q = nnz
    DO p = ia(j), ia(j+1)-1
        i = ja(p)
        IF (w(i)>=q) THEN
            a(w(i)) = a(w(i)) + a(p)
        ELSE
            w(i) = nnz
            ja(nnz) = i
            a(nnz) = a(p)
            nnz = nnz + 1
        ENDIF
    ENDDO
    ia(j) = q
ENDDO
ia(n+1) = nnz
nz = nnz
RETURN
END SUBROUTINE cs_dupl


SUBROUTINE cs_add(ia, ja, a, anz, ib, jb, b, bnz, n, alpha, beta)
! Performs matrix addtion C = alpha*A + beta*B in CSR format
! 
! Note: C is returned by A, thus A must have enough memory to store C
! 
! Translated and modified by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! from its C++ version in CSparse (Dr. Tim Davis, TAMU)
! 
!  ---Entries---
! 
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------
!   ia       integer   n+1      row pointers of A
!   ja       integer   anz      cols of A
!   a        real*8    anz      vals of A
!   anz      integer   1        number of nonzero entries in A
!   ib       integer   n+1      row pointers of B
!   jb       integer   bnz      cols of B
!   b        real*8    bnz      vals of B
!   bnz      integer   1        number of nonzero entries in B
!   n        integer   1        order of the full matrix
!   alpha    real*8    1        parameter
!   beta     real*8    1        parameter
! 
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(INOUT) :: anz
INTEGER, INTENT(IN) :: n, bnz
INTEGER, DIMENSION(n+1), INTENT(INOUT) :: ia
INTEGER, DIMENSION(n+1), INTENT(IN) :: ib
INTEGER, DIMENSION(anz), INTENT(INOUT) :: ja
INTEGER, DIMENSION(bnz), INTENT(IN) :: jb
REAL(kind = 8), DIMENSION(anz), INTENT(INOUT) :: a
REAL(kind = 8), DIMENSION(bnz), INTENT(IN) :: b
REAL(kind = 8), INTENT(IN) :: alpha, beta
! ---Internal variables---
INTEGER :: j, p, nz
INTEGER, DIMENSION(n+1) :: ic
INTEGER, DIMENSION(anz + bnz) :: jc
REAL(kind = 8), DIMENSION(anz + bnz) :: c
INTEGER, DIMENSION(n) :: w ! work array
REAL(kind = 8), DIMENSION(n) :: x
! Initializing work array
w = 0
nz = 1
DO j = 1, n  
    ic(j) = nz
    CALL cs_scatter(ia,ja,a,anz, j, alpha, w, x, j+1, jc,anz+bnz, n, nz)
    CALL cs_scatter(ib,jb,b,bnz, j, beta, w, x, j+1, jc,anz+bnz, n, nz)
    DO p = ic(j), nz-1
        c(p) = x(jc(p))
    ENDDO
ENDDO
ic(n+1) = nz
nz = nz - 1
! Since original cs_add does not return with sorted columns, thus
! first sorting rows and columns by transposing twice
CALL cs_transpose(ic, jc, c, n, nz) ! 1st
CALL cs_transpose(ic, jc, c, n, nz) ! 2nd
! then, summing up duplicate entries
CALL cs_dupl(ic, jc, c, n, nz)
! and, removing zero entries
CALL rm0csr(ic, jc, c, n, nz)
! finally, output
ia = ic
ja(1:nz) = jc(1:nz)
a(1:nz) = c(1:nz)
anz = nz
RETURN
END SUBROUTINE cs_add


SUBROUTINE cs_scatter(ia,ja,a,anz, j, beta, w, x, mark, jc,cnz, n, nz)
! Copys a sparse vector into a dense one
! 
! Translated by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! from its C++ version in CSparse (Dr. Tim Davis, TAMU)
! 
!  ---Entries---
! 
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------
!   ia       integer   n+1      row pointers of A
!   ja       integer   anz      cols of A
!   a        real*8    anz      vals of A
!   anz      integer   1        number of nonzero entries in A
!   j        integer   1        current column
!   beta     real*8    1        parameters
!   w        integer   n        work array
!   x        real*8    n        work array
!   mark     integer   1        N/A
!   jc       integer   cnz      cols of C
!   cnz      integer   1        number of nonzero entries in C
!   n        integer   1        order of the full matrix
!   nz       integer   1        number of nonzero entries
! 
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: anz, cnz, j, mark, n
INTEGER, INTENT(INOUT) :: nz
INTEGER, DIMENSION(n+1), INTENT(IN) :: ia
INTEGER, DIMENSION(anz), INTENT(IN) :: ja
INTEGER, DIMENSION(cnz), INTENT(INOUT) :: jc
REAL(kind = 8), DIMENSION(anz), INTENT(IN) :: a
INTEGER, DIMENSION(n), INTENT(INOUT) :: w
REAL(kind = 8), DIMENSION(n), INTENT(INOUT) :: x
REAL(kind = 8), INTENT(IN) :: beta
! ---Internal variables---
INTEGER :: i, p
DO p = ia(j), ia(j+1)-1
    i = ja(p)    
    IF (w(i) < mark) THEN
        w(i) = mark
        jc(nz) = i
        nz = nz + 1
        x(i) = beta*a(p)
    ELSE
        x(i) = x(i) + beta*a(p)
    ENDIF
ENDDO
RETURN
END SUBROUTINE cs_scatter


SUBROUTINE cs_permute(ia, ja, a, n, nz, p)
! Permutes a CSR matrix A as C = A(P,P), where P is the permutation vector
! 
! Translated by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2015
! from its C++ version in CSparse (Dr. Tim Davis, TAMU)
! 
!  ---Entries---
! 
!    Name     Type      Size     Description
!  -------- --------- -------- ---------------------------------------
!   ia       integer   n+1      row pointers
!   ja       integer   nz       cols
!   a        real*8    nz       vals
!   nz       integer   1        number of nonzero entries
!   n        integer   1        order of the full matrix
!   p        integer   n        permutation vector
! 
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: n, nz
INTEGER, DIMENSION(n+1), INTENT(INOUT) :: ia
INTEGER, DIMENSION(nz), INTENT(INOUT) :: ja
REAL(kind = 8), DIMENSION(nz), INTENT(INOUT) :: a
INTEGER, DIMENSION(n), INTENT(IN) :: p
! ---Internal variables---
INTEGER :: i, t, j, k, nnz
INTEGER, DIMENSION(n+1) :: ic
INTEGER, DIMENSION(nz) :: jc
REAL(kind = 8), DIMENSION(nz) :: c
INTEGER, DIMENSION(n) :: pinv
! Initialize pinv, the inverse permutation vector
DO i = 1, n
    pinv(p(i)) = i
ENDDO
! Perform permutation
nnz = 1
DO k = 1, n
    ic(k) = nnz
    j = p(k)
    DO t = ia(j), ia(j+1)-1
        c(nnz) = a(t)
        jc(nnz) = pinv(ja(t))
        nnz = nnz + 1
    ENDDO
ENDDO
ic(n+1) = nnz
! Output
ia = ic
ja = jc
a = c
RETURN
END SUBROUTINE cs_permute
