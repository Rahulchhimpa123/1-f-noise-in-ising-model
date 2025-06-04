!nl = linear size of the system
!m = no. of iteration or monte carlo steps
!j3 = equilibration steps or transient


implicit real (kind=8)(a-h,o-z)
integer,parameter::n=18,m=2**(n)
integer,parameter::nl=2**3,nl1=nl**2,j3=10**2

        call randomseed()
        call ps1(m,nl,nl1,j3)

stop
end

! Subroutines
subroutine ps1(m,nl,nl1,j3)
implicit real(kind=8)(a-h,o-z)
integer,dimension(m)::nx
integer,dimension(nl1)::nz
integer,dimension(nl1,4)::nn    !nearest neighbors
real (kind=8),dimension(5)::p1
integer, dimension (:), allocatable :: nx1
nx=0
nz=0
nn=0
p1=0.d0

        call initial_conf(nz,nl1,0.5d0)
        call neighbors(nl,nl1,nn)
        tem = 2.269185
        call cost_exp(p1,tem)
        call single_spin_flip_metropolis(nl,nl1,nz,m,tem,j3,nn,p1,nx)
        nm = maxval(nx)
        nm1 = minval(nx)
        allocate(nx1(nm))
        nx1=0
        do i=1,m
            nx1(nx(i)) = nx1(nx(i))+1
        enddo
        do i=nm1,nm
            nf = nx1(i)
            if(nf>0)write(13,*)i,nf
        enddo
        !goto 10
        s1 = dfloat(sum(nx))/dfloat(m)
        s2 = 0.d0
        do i=1,m
           s2 = s2 + (dfloat(nx(i))-s1)**2
        enddo
        s2  = (s2/dfloat(m))**0.5
        print*,s1,s2
10      continue
        
return
end subroutine ps1

