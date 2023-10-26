
C**** Lapack routines for decomposition of overlap matrix

      SUBROUTINE CHOLSK(S,C,N)
      use numbers

      implicit none
      integer n
      complex(imag) s(n,n)
      complex(imag) c(n,n)
C     local
      integer i, j, info

      do i = 1, n
         do j = 1, i
            c(j,i) = conjg(s(i,j))
         enddo
      enddo
      call zpotrf('u', n, c, n, info)
      if (info /= 0) then
         call errnic(0,info,'CHOLSK', 'Basis set Linear dependence')
      endif
      do i = 1, n-1
         do j = i+1, n
            c(j,i) = cmplx(0.0_float,0.0_float)
         enddo
      enddo
      call ztrtri('u', 'n', n, c, n, info)
      if (info /= 0) then
         call errnic(0,info,'CHOLSK', 'Basis set Linear dependence')
      endif
      END

      SUBROUTINE CHOLSK_GPU(S,C,N)
      use numbers
      use paral1_module,only:iam,nproc,ngpunode,nsocket

      implicit none
      integer n
      complex(imag) s(n,n)
      complex(imag) c(n,n)
C     local
      integer i, j, info
      integer mydev,iamnew,iammapped

      do i = 1, n
         do j = 1, i
            c(j,i) = conjg(s(i,j))
         enddo
      enddo

      iamnew=(nproc-iam)-1
      iammapped=iamnew/nsocket
      mydev=mod(iammapped,ngpunode)

!$acc set device_num(mydev)
!$acc data copy(c)
!$acc host_data use_device(c)
      call zpotrf('u', n, c, n, info)
!$acc end host_data
      if (info /= 0) then
         call errnic(0,info,'CHOLSK', 'Basis set Linear dependence')
      endif
      do i = 1, n-1
         do j = i+1, n
            c(j,i) = cmplx(0.0_float,0.0_float)
         enddo
      enddo
      call ztrtri('u', 'n', n, c, n, info)
      if (info /= 0) then
         call errnic(0,info,'CHOLSK', 'Basis set Linear dependence')
      endif
!$acc end data
      END

      SUBROUTINE RHOLSK(S,C,N)
      use numbers

      implicit none
      integer n
      real(float) s(n,n)
      real(float) c(n,n)
C     local
      integer i, j, info

      do i = 1, n
         do j = 1, i
            c(j,i) = s(i,j)
         enddo
      enddo
      call dpotrf('u', n, c, n, info)
      if (info /= 0) then
         call errnic(0,info,'RHOLSK', 'Basis set Linear dependence')
      endif
      do i = 1, n-1
         do j = i+1, n
            c(j,i) = 0.0_float
         enddo
      enddo
      call dtrtri('u', 'n', n, c, n, info)
      if (info /= 0) then
         call errnic(0,info,'RHOLSK', 'Basis set Linear dependence')
      endif
      END

      SUBROUTINE RHOLSK_GPU(S,C,N)
      use numbers
      use paral1_module,only:iam,nproc,ngpunode,nsocket

      implicit none
      integer n
      real(float) s(n,n)
      real(float) c(n,n)
C     local
      integer i, j, info
      integer mydev,iamnew,iammapped

      do i = 1, n
         do j = 1, i
            c(j,i) = s(i,j)
         enddo
      enddo

      iamnew=(nproc-iam)-1
      iammapped=iamnew/nsocket
      mydev=mod(iammapped,ngpunode)

!$acc set device_num(mydev)
!$acc data copy(c)
!$acc host_data use_device(c)
      call dpotrf('u', n, c, n, info)
!$acc end host_data
      if (info /= 0) then
         call errnic(0,info,'RHOLSK', 'Basis set Linear dependence')
      endif
      do i = 1, n-1
         do j = i+1, n
            c(j,i) = 0.0_float
         enddo
      enddo
      call dtrtri('u', 'n', n, c, n, info)
      if (info /= 0) then
         call errnic(0,info,'RHOLSK', 'Basis set Linear dependence')
      endif
!$acc end data
      END

C*** Example Lapack Divide and Conquer based Eigenvalue/vector routines.

      SUBROUTINE REIGM(A,EIG,N)
      use numbers

      implicit none
      integer, intent(in) :: n
      real(float) a(n*n), eig(n)
C     local vars
      integer info, lwork, ilwork

      lwork = 3 * n
      ilwork = n
      call reigen_val_dc(a, eig, n, lwork, ilwork, info)
      if (info /= 0) then
         call errnic(0, info, 'RDIAG', 'Eigenvalue problem')
      endif
      END

      subroutine reigen_val_dc(a, eig, n, m, im, info)
      use numbers

      implicit none
      integer, intent(in) :: n, m, im
      real(float) a(n*n), eig(n)
      integer info
C     local vars
      real(float), dimension(:), allocatable :: work!(m)
      integer, dimension(:), allocatable :: lwork!(im)
      Allocate(work(m))
      Allocate(lwork(im))
      call dsyevd('N', 'L', n, a, n, eig, work, m, lwork, im, info)
      DeAllocate(lwork)
      DeAllocate(work)
      end

      SUBROUTINE REIGN(A,R,EIG,N,IORD,KHOLD)
      use numbers

      implicit none
      integer, intent(in) :: n, iord, khold
      real(float) a(n,n), r(n,n), eig(n)
C     local vars
      integer info, k, lwork, ilwork

      if (n < 1024) then
         k = 10
      else if (n < 2048) then
         k = 11
      else if (n < 4096) then
         k = 12
      else
         k = 15
      endif
      lwork = 6 * n + 2 * n * k + 3 * n * n
      ilwork = 6 * n + 2
      if (iord /= 0) then
         call reigen_mult_dc(a, r, eig, n, lwork, ilwork, info)
      else
         call reigen_vec_dc(a, r, eig, n, lwork, ilwork, info)
      endif
      if (info /= 0) then
         call errnic(0, info, 'RDIAG', 'Eigenvector problem')
      endif
      END

      SUBROUTINE REIGN_GPU(A,R,EIG,N,IORD,KHOLD)
      use numbers

      implicit none
      integer, intent(in) :: n, iord, khold
      real(float) a(n,n), r(n,n), eig(n)
C     local vars
      integer info, k, lwork, ilwork

      if (n < 1024) then
         k = 10
      else if (n < 2048) then
         k = 11
      else if (n < 4096) then
         k = 12
      else
         k = 15
      endif
      lwork = 6 * n + 2 * n * k + 3 * n * n
      ilwork = 6 * n + 2
      if (iord /= 0) then
         call reigen_mult_dc_gpu(a, r, eig, n, lwork, ilwork, info)
      else
         call reigen_vec_dc_gpu(a, r, eig, n, lwork, ilwork, info)
      endif
      if (info /= 0) then
         call errnic(0, info, 'RDIAG', 'Eigenvector problem')
      endif
      END

      subroutine reigen_vec_dc(a, r, eig, n, m, im, info)
      use numbers

      implicit none
      integer, intent(in) :: n, m, im
      real(float) a(n,n), r(n,n), eig(n)
      integer info
C     local vars
!      real(float) work(m)
!      integer lwork(im)
      real(float), dimension(:), allocatable :: work
      integer, dimension(:), allocatable :: lwork
      Allocate(work(m))
      Allocate(lwork(im))

      call dsyevd('V', 'L', n, a, n, eig, work, m, lwork, im, info)
      if (info == 0) then
         call dcopy(n * n, a, 1, r, 1)
      endif
      DeAllocate(lwork)
      DeAllocate(work)
      end

      subroutine reigen_mult_dc(a, r, eig, n, m, im, info)
      use numbers

      implicit none
      integer, intent(in) :: n, m, im
      real(float) a(n,n), r(n,n), eig(n)
      integer info
C     local vars
      !real(float) temp(n, n), work(m)
      !integer lwork(im)
      real(float), dimension(:), allocatable :: work
      real(float), dimension(:,:), allocatable :: temp 
      integer, dimension(:), allocatable :: lwork
      Allocate(work(m))
      Allocate(lwork(im))

      call dsyevd('V', 'L', n, a, n, eig, work, m, lwork, im, info)
      DeAllocate(lwork)
      DeAllocate(work)
      if (info == 0) then
      Allocate(temp(n,n))
         call dcopy(n * n, r, 1, temp, 1)
         call dgemm('n', 'n', n, n, n, 1.0_float, temp, n, a, n,
     *        0.0_float, r, n)
      DeAllocate(temp)
      endif
      end

      subroutine reigen_vec_dc_gpu(a, r, eig, n, m, im, info)
      use numbers
      use cublas
      use paral1_module,only:iam,nproc,ngpunode,nsocket
      use memory_use      

      implicit none
      integer, intent(in) :: n, m, im
      real(float) a(n,n), r(n,n), eig(n)
      integer info
C     local vars
!      real(float) work(m)
!      integer lwork(im)
      real(float), dimension(:), allocatable :: work
      integer, dimension(:), allocatable :: lwork
      integer mydev,iamnew,iammapped

!      Allocate(work(m))
!      Allocate(lwork(im))
      call cryalloc(work,m,'reigen_vec_dc_gpu','work')
      call cryalloc(lwork,im,'reigen_vec_dc_gpu','lwork')
      
      iamnew=(nproc-iam)-1
      iammapped=iamnew/nsocket
      mydev=mod(iammapped,ngpunode)

!$acc set device_num(mydev)
!$acc data copy(a) 
!$acc data copyout(eig) create(work,lwork)
!$acc host_data use_device(a,eig,work,lwork) 
      call dsyevd('V', 'L', n, a, n, eig, work, m, lwork, im, info)
!$acc end host_data
!$acc end data
      if (info == 0) then
!$acc data copy(r)
!$acc host_data use_device(a,r)                
         call dcopy(n * n, a, 1, r, 1)
!$acc end host_data
!$acc end data         
      endif

!$acc end data      
!      DeAllocate(lwork)
!      DeAllocate(work)
      call crydealloc(lwork,'reigen_vec_dc_gpu','lwork')
      call crydealloc(work,'reigen_vec_dc_gpu','work')

!      CALL TIMVRS('reigen_vec_dc_gpu    ')
  
      end

      subroutine reigen_mult_dc_gpu(a, r, eig, n, m, im, info)
      use numbers
      use cublas 
      use paral1_module,only:iam,nproc,ngpunode,nsocket
      use memory_use      

      implicit none
      integer, intent(in) :: n, m, im
      real(float) a(n,n), r(n,n), eig(n)
      integer info
C     local vars
      !real(float) temp(n, n), work(m)
      !integer lwork(im)
      real(float), dimension(:), allocatable :: work
      real(float), dimension(:,:), allocatable :: temp 
      integer, dimension(:), allocatable :: lwork
      integer mydev,iamnew,iammapped

!      Allocate(work(m))
!      Allocate(lwork(im))
      call cryalloc(work,m,'reigen_mult_dc_gpu','work')
      call cryalloc(lwork,im,'reigen_mult_dc_gpu','lwork')

      iamnew=(nproc-iam)-1
      iammapped=iamnew/nsocket
      mydev=mod(iammapped,ngpunode)

!$acc set device_num(mydev)
!$acc data copy(a)
!$acc data copyout(eig) create(work,lwork)
!$acc host_data use_device(a,eig,work,lwork) 
      call dsyevd('V', 'L', n, a, n, eig, work, m, lwork, im, info)
!$acc end host_data
!$acc end data

      if (info == 0) then
!      Allocate(temp(n,n))
      call cryalloc(temp,n,n,'reigen_mult_dc_gpu','temp')

!$acc data copy (r) create(temp)
!$acc host_data use_device(a,r,temp)      
         call dcopy(n * n, r, 1, temp, 1)
         call dgemm('n', 'n', n, n, n, 1.0_float, temp, n, a, n,
     *        0.0_float, r, n)
!$acc end host_data
!$acc end data
!      DeAllocate(temp)
      call crydealloc(temp,'reigen_mult_dc_gpu','temp')
      endif
!$acc end data

!      DeAllocate(lwork)
!      DeAllocate(work)
      call crydealloc(lwork,'reigen_mult_dc_gpu','lwork')
      call crydealloc(work,'reigen_mult_dc_gpu','work')

!      CALL TIMVRS('reigen_mult_dc_gpu    ')
     
      end

      SUBROUTINE CEIGM(A,EIG,N)
      use numbers

      implicit none
      integer, intent(in) :: n
      complex(imag) a(n,n)
      real(float) eig(n)
C     local vars
      integer lwork, ilwork, rwork, info

C     Slightly dangerous to treat a as complex
      lwork = 2 * n
      rwork = n
      ilwork = n
      call ceigen_val_dc(a, eig, n, lwork, rwork, ilwork, info)
      if (info /= 0) then
         call errnic(0, info, 'CDIAG', 'Eigenvalue problem')
      endif
      END

      subroutine ceigen_val_dc(a, eig, n, m, rm, im, info)
      use numbers

      implicit none
      integer, intent(in) :: n, m, rm, im
      complex(imag) a(n,n)
      real(float) eig(n)
      integer info
C     local vars
!      complex(imag) work(m)
!      real(float) rwork(rm)
!      integer lwork(im)
      complex(imag),dimension(:), allocatable :: work
      real(float),dimension(:), allocatable :: rwork
      integer,dimension(:), allocatable :: lwork
      Allocate(work(m))
      Allocate(rwork(rm))
      Allocate(lwork(im))
      call zheevd('N', 'L', n, a, n, eig, work, m, rwork, rm, lwork,
     *     im, info)
      DeAllocate(lwork)
      DeAllocate(rwork)
      DeAllocate(work)
      end

      SUBROUTINE CEIGN(A,R,EIG,N,IORD,KHOLD)
      use numbers

      implicit none
      integer, intent(in) :: n, iord, khold
      complex(imag) a(n,n), r(n,n)
      real(float) eig(n)
C     local vars
      integer k, lwork, rwork, ilwork, info

C     Slightly dangerous to treat a as complex
      if (n < 1024) then
         k = 10
      else if (n < 2048) then
         k = 11
      else if (n < 4096) then
         k = 12
      else
         k = 15
      endif
      lwork = 2 * n + n * n
      rwork = 5 * n + 2 * n * k + 3 * n * n
      ilwork = 6 * n + 2
      if (iord /= 0) then
         call ceigen_mult_dc(a, r, eig, n, lwork, rwork, ilwork, info)
      else
         call ceigen_vec_dc(a, r, eig, n, lwork, rwork, ilwork, info)
      endif
      if (info /= 0) then
         call errnic(0, info, 'CDIAG', 'Eigenvector problem')
      endif
      END

      SUBROUTINE CEIGN_GPU(A,R,EIG,N,IORD,KHOLD)
      use numbers

      implicit none
      integer, intent(in) :: n, iord, khold
      complex(imag) a(n,n), r(n,n)
      real(float) eig(n)
C     local vars
      integer k, lwork, rwork, ilwork, info

C     Slightly dangerous to treat a as complex
      if (n < 1024) then
         k = 10
      else if (n < 2048) then
         k = 11
      else if (n < 4096) then
         k = 12
      else
         k = 15
      endif
      lwork = 2 * n + n * n
      rwork = 5 * n + 2 * n * k + 3 * n * n
      ilwork = 6 * n + 2
      if (iord /= 0) then
         call ceigen_mult_dc_gpu(a, r, eig, n, lwork, rwork, ilwork,
     * info)
      else
         call ceigen_vec_dc_gpu(a, r, eig, n, lwork, rwork, ilwork,
     * info)
      endif
      if (info /= 0) then
         call errnic(0, info, 'CDIAG', 'Eigenvector problem')
      endif
      END

      subroutine ceigen_vec_dc(a, r, eig, n, m, rm, im, info)
      use numbers

      implicit none
      integer, intent(in) :: n, m, rm, im
      complex(imag) a(n,n), r(n,n)
      real(float) eig(n)
      integer info
C     local vars
!      complex(imag) work(m)
!      real(float) rwork(rm)
!      integer lwork(im)
      complex(imag),dimension(:), allocatable :: work
      real(float),dimension(:), allocatable :: rwork
      integer,dimension(:), allocatable :: lwork
      Allocate(work(m))
      Allocate(rwork(rm))
      Allocate(lwork(im))
      call zheevd('V', 'L', n, a, n, eig, work, m, rwork, rm, lwork,
     *     im, info)
      if (info == 0) then
         call zcopy(n * n, a, 1, r, 1)
      endif
      DeAllocate(lwork)
      DeAllocate(rwork)
      DeAllocate(work)
      end

      subroutine ceigen_mult_dc(a, r, eig, n, m, rm, im, info)
      use numbers

      implicit none
      integer, intent(in) :: n, m, rm, im
      complex(imag) a(n,n), r(n,n)
      real(float) eig(n)
      integer info
C     local vars
!      complex(imag) temp(n, n), work(m)
!      real(float) rwork(rm)
!      integer lwork(im)
      complex(imag),dimension(:), allocatable :: work
      complex(imag),dimension(:,:), allocatable :: temp
      real(float),dimension(:), allocatable :: rwork
      integer,dimension(:), allocatable :: lwork
      complex(imag) alpha, beta
      Allocate(work(m))
      Allocate(rwork(rm))
      Allocate(lwork(im))

      call zheevd('V', 'L', n, a, n, eig, work, m, rwork, rm, lwork,
     *     im, info)
      DeAllocate(lwork)
      DeAllocate(rwork)
      DeAllocate(work)
      if (info == 0) then
      Allocate(temp(n,n))
         call zcopy(n * n, r, 1, temp, 1)
         alpha = (1.0_float, 0.0_float)
         beta = (0.0_float, 0.0_float)
         call zgemm('n', 'n', n, n, n, alpha, temp, n, a, n, beta, r, n)
      DeAllocate(temp)
      endif
      end

      subroutine ceigen_vec_dc_gpu(a, r, eig, n, m, rm, im, info)
      use numbers
      use cublas
      use paral1_module,only:iam,nproc,ngpunode,nsocket
      use memory_use     

      implicit none
      integer, intent(in) :: n, m, rm, im
      complex(imag) a(n,n), r(n,n)
      real(float) eig(n)
      integer info
C     local vars
!      complex(imag) work(m)
!      real(float) rwork(rm)
!      integer lwork(im)
      complex(imag),dimension(:), allocatable :: work
      real(float),dimension(:), allocatable :: rwork
      integer,dimension(:), allocatable :: lwork
      integer mydev,iammapped,iamnew

!      Allocate(work(m))
!      Allocate(rwork(rm))
!      Allocate(lwork(im))

      call Cryalloc(work,m,'ceigen_vec_dc_gpu','work')
      call Cryalloc(rwork,rm,'ceigen_vec_dc_gpu','rwork')
      call Cryalloc(lwork,im,'ceigen_vec_dc_gpu','lwork')

      iamnew=(nproc-iam)-1
      iammapped=iamnew/nsocket
      mydev=mod(iammapped,ngpunode)

!$acc set device_num(mydev)
!$acc data copy(a) 
!$acc data copyout(eig) create(work,rwork,lwork)
!$acc host_data use_device(a,eig,work,rwork,lwork)
      call zheevd('V', 'L', n, a, n, eig, work, m, rwork, rm, lwork,
     *     im, info)
!$acc end host_data
!$acc end data      
      if (info == 0) then
!$acc data copy(r)
!$acc host_data use_device(a,r)              
         call zcopy(n * n, a, 1, r, 1)
!$acc end host_data
!$acc end data         
      endif
!$acc end data

!      DeAllocate(lwork)
!      DeAllocate(rwork)
!      DeAllocate(work)

      call Crydealloc(lwork,'ceigen_vec_dc_gpu','lwork')
      call Crydealloc(rwork,'ceigen_vec_dc_gpu','rwork')
      call Crydealloc(work,'ceigen_vec_dc_gpu','work')  

!      CALL TIMVRS('ceigen_vec_dc_gpu    ')

      end

      subroutine ceigen_mult_dc_gpu(a, r, eig, n, m, rm, im, info)
      use numbers
      use cublas
      use paral1_module,only:iam,nproc,ngpunode,nsocket
      use memory_use       

      implicit none
      integer, intent(in) :: n, m, rm, im
      complex(imag) a(n,n), r(n,n)
      real(float) eig(n)
      integer info
C     local vars
!      complex(imag) temp(n, n), work(m)
!      real(float) rwork(rm)
!      integer lwork(im)
      complex(imag),dimension(:), allocatable :: work
      complex(imag),dimension(:,:), allocatable :: temp
      real(float),dimension(:), allocatable :: rwork
      integer,dimension(:), allocatable :: lwork
      complex(imag) alpha, beta
      integer mydev,iammapped,iamnew

!      Allocate(work(m))
!      Allocate(rwork(rm))
!      Allocate(lwork(im))

      call Cryalloc(work,m,'ceigen_mult_dc_gpu','work')
      call Cryalloc(rwork,rm,'ceigen_mult_dc_gpu','rwork')
      call Cryalloc(lwork,im,'ceigen_mult_dc_gpu','lwork')

      iamnew=(nproc-iam)-1
      iammapped=iamnew/nsocket
      mydev=mod(iammapped,ngpunode)

!$acc set device_num(mydev)
!$acc data copy(a)
!$acc data copyout(eig) create(work,rwork,lwork)
!$acc host_data use_device(a,eig,work,rwork,lwork)
      call zheevd('V', 'L', n, a, n, eig, work, m, rwork, rm, lwork,
     *     im, info)
!$acc end host_data 
!$acc end data

      if (info == 0) then

!      Allocate(temp(n,n))

      call Cryalloc(temp,n,n,'ceigen_mult_dc_gpu','temp')

!$acc data copy(r) create(temp)
!$acc host_data use_device(r,temp,a)
         call zcopy(n * n, r, 1, temp, 1)
         alpha = (1.0_float, 0.0_float)
         beta = (0.0_float, 0.0_float)
         call zgemm('n', 'n', n, n, n, alpha, temp, n, a, n, beta, r, n)
!$acc end host_data
!$acc end data 

!      DeAllocate(temp)
      call Crydealloc(temp,'ceigen_mult_dc_gpu','temp')

      endif

!$acc end data

!      DeAllocate(lwork)
!      DeAllocate(rwork)
!      DeAllocate(work)

      call Crydealloc(lwork,'ceigen_mult_dc_gpu','lwork')
      call Crydealloc(rwork,'ceigen_mult_dc_gpu','rwork')
      call Crydealloc(work,'ceigen_mult_dc_gpu','work')

!      CALL TIMVRS('ceigen_mult_dc_gpu    ')
      end

C**** BLAS based similarity transformations 

      subroutine zsimilarity(h, nh, u, nu, temp, nt, n)
      use numbers

      implicit none
      integer, intent(in) :: nh, nu, nt, n
      complex(imag), intent(in) :: u(nu,n)
      complex(imag) h(nh,n), temp(nt,n)
C     local vars
      integer, parameter :: blocksize = 100
      integer nn, nextra, nblock, nmax
      complex(imag), parameter :: alpha = (1.0_float, 0.0_float)
      complex(imag), parameter :: beta = (0.0_float, 0.0_float)

      call zgemm('c', 'n', n, n, n, alpha, u, nu, h, nh, beta, temp, nt)
C     only do lower triangular part of second mult for diag
      nextra = mod(n, blocksize)
      nblock = n - nextra
C     make last bit blocksize to 2 * blocksize - 1, to avoid nextra = 1
      if (nblock >= blocksize) then
         nblock = nblock - blocksize
         nextra = nextra + blocksize
      endif
      nmax = n
      do nn = 1, nblock, blocksize
         call zgemm('n', 'n', nmax, blocksize, n, alpha, temp(nn, 1),
     *        nt, u(1, nn), nu, beta, h(nn, nn), nh)
         nmax = nmax - blocksize
      enddo
      call zgemm('n', 'n', nmax, nmax, n, alpha, temp(nblock + 1, 1),
     *     nt, u(1, nblock + 1), nu, beta,
     *     h(nblock + 1, nblock + 1), nh)
      end

      subroutine zsimilarity_gpu(h, nh, u, nu, temp, nt, n)
      use numbers
      use cublas
      use paral1_module,only:iam,nproc,ngpunode,nsocket

      implicit none
      integer, intent(in) :: nh, nu, nt, n
      complex(imag), intent(in) :: u(nu,n)
      complex(imag) h(nh,n), temp(nt,n)
C     local vars
      integer, parameter :: blocksize = 100
      integer nn, nextra, nblock, nmax
      complex(imag), parameter :: alpha = (1.0_float, 0.0_float)
      complex(imag), parameter :: beta = (0.0_float, 0.0_float)
      integer mydev,iammapped,iamnew

      iamnew=(nproc-iam)-1
      iammapped=iamnew/nsocket
      mydev=mod(iammapped,ngpunode)

!$acc set device_num(mydev)
!$acc data copyin(u) copy(h,temp)
!$acc host_data use_device(u,h,temp)
      call zgemm('c', 'n', n, n, n, alpha, u, nu, h, nh, beta, temp, nt)
C     only do lower triangular part of second mult for diag
      nextra = mod(n, blocksize)
      nblock = n - nextra
C     make last bit blocksize to 2 * blocksize - 1, to avoid nextra = 1
      if (nblock >= blocksize) then
         nblock = nblock - blocksize
         nextra = nextra + blocksize
      endif
      nmax = n
      do nn = 1, nblock, blocksize
         call zgemm('n', 'n', nmax, blocksize, n, alpha, temp(nn, 1),
     *        nt, u(1, nn), nu, beta, h(nn, nn), nh)
         nmax = nmax - blocksize
      enddo
      call zgemm('n', 'n', nmax, nmax, n, alpha, temp(nblock + 1, 1),
     *     nt, u(1, nblock + 1), nu, beta,
     *     h(nblock + 1, nblock + 1), nh)
!$acc end host_data
!$acc end data

      end

      subroutine dsimilarity(h, nh, u, nu, temp, nt, n)
      use numbers

      implicit none
      integer, intent(in) :: nh, nu, nt, n
      real(float), intent(in) :: u(nu,n)
      real(float) h(nh,n), temp(nt,n)
C     local vars
      integer, parameter :: blocksize = 100
      integer i, j, nn, nextra, nblock, nmax

      call dgemm('t', 'n', n, n, n, 1.0_float, u, nu, h, nh,
     *     0.0_float, temp, nt)
C     only do lower triangular part of second mult for diag
      nextra = mod(n, blocksize)
      nblock = n - nextra
C     make last bit blocksize to 2 * blocksize - 1, to avoid nextra = 1
      if (nblock >= blocksize) then
         nblock = nblock - blocksize
         nextra = nextra + blocksize
      endif
      nmax = n
      do nn = 1, nblock, blocksize
         call dgemm('n', 'n', nmax, blocksize, n, 1.0_float,
     *        temp(nn, 1), nt, u(1, nn), nu, 0.0_float, h(nn, nn), nh)
         nmax = nmax - blocksize
      enddo
      call dgemm('n', 'n', nmax, nmax, n, 1.0_float,
     *     temp(nblock + 1, 1), nt, u(1, nblock + 1), nu,
     *     0.0_float, h(nblock + 1, nblock + 1), nh)
      end

      subroutine dsimilarity_gpu(h, nh, u, nu, temp, nt, n)
      use numbers
      use cublas
      use paral1_module,only:iam,nproc,ngpunode,nsocket

      implicit none
      integer, intent(in) :: nh, nu, nt, n
      real(float), intent(in) :: u(nu,n)
      real(float) h(nh,n), temp(nt,n)
C     local vars
      integer, parameter :: blocksize = 100
      integer i, j, nn, nextra, nblock, nmax
      integer mydev,iamnew,iammapped

      iamnew=(nproc-iam)-1
      iammapped=iamnew/nsocket
      mydev=mod(iammapped,ngpunode)

!$acc set device_num(mydev)
!$acc data copyin(u) copy(h,temp) 
!$acc host_data use_device(u,h,temp)
      call dgemm('t', 'n', n, n, n, 1.0_float, u, nu, h, nh,
     *     0.0_float, temp, nt)
C     only do lower triangular part of second mult for diag
      nextra = mod(n, blocksize)
      nblock = n - nextra
C     make last bit blocksize to 2 * blocksize - 1, to avoid nextra = 1
      if (nblock >= blocksize) then
         nblock = nblock - blocksize
         nextra = nextra + blocksize
      endif
      nmax = n
      do nn = 1, nblock, blocksize
         call dgemm('n', 'n', nmax, blocksize, n, 1.0_float,
     *        temp(nn, 1), nt, u(1, nn), nu, 0.0_float, h(nn, nn), nh)
         nmax = nmax - blocksize
      enddo
      call dgemm('n', 'n', nmax, nmax, n, 1.0_float,
     *     temp(nblock + 1, 1), nt, u(1, nblock + 1), nu,
     *     0.0_float, h(nblock + 1, nblock + 1), nh)
!$acc end host_data
!$acc end data

      end

      subroutine zfsimilarity(h, nh, u, nu, temp, nt, n)
      use numbers

      implicit none
      integer, intent(in) :: nh, nu, nt, n
      complex(imag), intent(in) :: u(nu,n)
      complex(imag) h(nh,n), temp(nt,n)
C     local vars
      complex(imag), parameter :: alpha = (1.0_float, 0.0_float)
      complex(imag), parameter :: beta = (0.0_float, 0.0_float)

      call zgemm('c', 'n', n, n, n, alpha, u, nu, h, nh, beta, temp, nt)
      call zgemm('n', 'n', n, n, n, alpha, temp, nt, u, nu, beta, h, nh)
      end

      subroutine dfsimilarity(h, nh, u, nu, temp, nt, n)
      use numbers

      implicit none
      integer, intent(in) :: nh, nu, nt, n
      real(float), intent(in) :: u(nu,n)
      real(float) h(nh,n), temp(nt,n)

      call dgemm('t', 'n', n, n, n, 1.0_float, u, nu, h, nh,
     *     0.0_float, temp, nt)
      call dgemm('n', 'n', n, n, n, 1.0_float, temp, nt, u, nu,
     *     0.0_float, h, nh)
      end

      FUNCTION ddot_internal_local(n,dx,incx,dy,incy)
      USE numbers
      IMPLICIT NONE
      INTEGER :: incx , incy , n
      REAL(float) :: DDOT_internal_local, ddot
      REAL(float) , DIMENSION(n) :: dx , dy
      INTENT (IN) dx , dy , incx , incy , n

      DDOT_internal_local = ddot(n, dx, incx, dy, incy)
      end

      SUBROUTINE diago (A,W,N)
      use numbers
      use memory_use
      implicit none
      integer, intent(in) :: n
      real(float) a(n,n), w(n)
C     local vars
      integer info, k, m, im
      real(float), dimension(:), allocatable :: work
      integer, dimension(:), allocatable :: lwork

      if (n < 1024) then
         k = 10
      else if (n < 2048) then
         k = 11
      else if (n < 4096) then
         k = 12
      else
         k = 15
      endif
      m = 6 * n + 2 * n * k + 3 * n * n
      im = 6 * n + 2
      call cryalloc(work,m,'DIAGO','real workspace')
      call cryalloc(lwork,im,'DIAGO','int workspace')
      call dsyevd('V', 'L', n, a, n, w, work, m, lwork, im, info)
      if (info /= 0) then
         call errnic(0, info, 'DIAGO', 'Eigenvector problem')
      endif
      call crydealloc(lwork,'DIAGO','int workspace')
      call crydealloc(work,'DIAGO','real workspace')
      end

      Subroutine diis_solve(matrix,vector,isz)
      Use numbers
      Use memory_use
      Implicit none
      Real(float), dimension(*) :: matrix,vector
      Integer, dimension(:), allocatable :: ipiv
      Integer :: isz,info
      Call Cryalloc(ipiv,isz,'diis_solve','ipiv')
      Call dgesv(isz,1,matrix,isz,ipiv,vector,isz,info)
      Call Crydealloc(ipiv,'diis_solve','ipiv')
      End Subroutine diis_solve



