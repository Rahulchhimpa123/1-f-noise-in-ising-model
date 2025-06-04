!program correlation time series
! m = evolution time, m1 = lag time


implicit real (kind=8)(a-h,o-z)
integer,parameter::m=2**(20),m1=100001,m2=10,k=10!realizations

        call random_seed
        nl = 2**6
        j3 = 10**2
        nl1=nl**2
        call p_corr_1(m,m1,m2,nl,nl1,j3,k)
        
                
stop
end


subroutine p_corr_1(m,m1,m2,nl,nl1,j3,k)
implicit real(kind=8)(a-h,o-z)
real (kind=8),dimension(m1)::sx2,z2
z2=0.d0
sx2=0.d0
! Set the desired number of threads
!$call omp_set_num_threads(k)  ! Set to the desired number of threads
!$OMP PARALLEL DO PRIVATE(i)
    do i=1,k
     !print*,i
     call lcorr(m,m1,m2,nl,nl1,j3,sx2,i)
     z2=z2+sx2
    enddo
!$OMP END PARALLEL DO
    z2=z2/dfloat(k)
    do i=2,m1    ! power sepctrum: f and S(f)
       write(11,*)i-1,z2(i)
    enddo
    print*,z2(1)
     
end subroutine p_corr_1


! Subroutines
! Local correlation
subroutine lcorr(m,m1,m2,nl,nl1,j3,sx2,i)
implicit real(kind=8)(a-h,o-z)
real (kind=8),dimension(m)::x
real (kind=8),dimension(m1)::c1,c2
real (kind=8), intent(inout),dimension(m1)::sx2
integer,dimension(nl1)::nz
integer,dimension(nl1,4)::nn    !nearest neighbors
real (kind=8),dimension(5)::p1
c2=0.d0
nz=0
nn=0
p1=0.d0
        call initial_conf(nz,nl1,0.55d0)
        call neighbors(nl,nl1,nn)
        tem = 2.269185
        call cost_exp(p1,tem)
        

        do k1=1,m2
	   if(i==1) print*,k1
	   x=0.d0
	   c1=0.d0
           tem = 2.269185
           call single_spin_flip_metropolis(nl,nl1,nz,m,tem,j3,nn,p1,x)
           call correlation(m,x,m1,c1)
           c2 = c2 + c1
        enddo
        c2 = c2/dfloat(m2)
        sx2 = c2
return
end subroutine lcorr


!Correlation of one realization
subroutine correlation(n,x,m1,c1)
implicit real(kind=8)(a-h,o-z)
real (kind=8), intent(in),dimension(n)::x
real (kind=8), dimension(n)::c
real (kind=8), intent(inout),dimension(m1)::c1
 c = 0.d0
        
        xav = sum(x)/dfloat(n)
        do j=1,m1
           a = 0.d0
           do i=1,n-m1+1
              a = a + (x(i+j-1)-xav)*(x(i)-xav)
           enddo
           c(j) = a
           c1(j) = c(j)!/c(1)
           !write(13,*)(j-1),c(j)/c(1)
        enddo
        
return
end subroutine correlation

! noisy signal

!2d single_spin_flip
subroutine single_spin_flip_metropolis(nl,nl1,nz,m,tem,j3,nn,p1,x)
implicit real(kind=8)(a-h,o-z)
integer,intent(inout),dimension(nl1)::nz
real (kind=8),intent(inout),dimension(m)::x
integer,dimension(nl1)::nz2,nz1
integer,intent(in),dimension(nl1,4)::nn
real (kind=8),intent(in),dimension(5)::p1
x=0.d0
nz2 = 0
nz1=0
        !kc = 0
        do it = 1,m+j3
           call conf_update1(nl,nl1,nz,tem,nc,ns,c,nz2,nn,p1) ! GD
           if(it>j3)then
           x(it-j3) =dfloat(nc)
           endif
        enddo


        
return
end subroutine single_spin_flip_metropolis


! Glauber dynamics algorithms single spin flip
subroutine conf_update1(nl,nl1,nz,tem,nc,ns,c,nz2,nn,p1)
implicit real(kind=8)(a-h,o-z)
integer,intent(inout),dimension(nl1)::nz,nz2
integer,dimension(nl1)::nz1
integer,intent(in),dimension(nl1,4)::nn
real (kind=8),intent(in),dimension(5)::p1
nz1=0

        nc = 0
        c = 0.d0
        do i=1,nl1
              call random_site(nl1,j)
              nf = nz(j)
              !call nn_pbc(nl,nl1,j11,nz,ns)
              ns=nz(nn(j,1))+nz(nn(j,2))+nz(nn(j,3))+nz(nn(j,4))
              cost = 2*nf*ns
              k=int((cost+8)/4)+1
              p = p1(k)
              call random_number(r)
              if(r<p)then
                nf = nf*(-1)
                nc = nc + 1
                nz1(j) = 1
                c = c+r
              endif
              nz(j) = nf
        enddo
        ns = sum(nz1)
        nz2 = nz1

return
end subroutine conf_update1



subroutine cost_exp(p1,tem)
implicit real(kind=8)(a-h,o-z)
real (kind=8),intent(inout),dimension(5)::p1
    
    cost=-8.d0
    p1(1)=1.d0/(1.d0+exp(cost/tem))
    
    cost=-4.d0
    p1(2)=1.d0/(1.d0+exp(cost/tem))
    
    cost=0.d0
    p1(3)=1.d0/(1.d0+exp(cost/tem))
    
    cost=4.d0
    p1(4)=1.d0/(1.d0+exp(cost/tem))
    
    cost=8.d0
    p1(5)=1.d0/(1.d0+exp(cost/tem))

return
endsubroutine cost_exp



subroutine neighbors(nl,nl1,nn)
implicit real(kind=8)(a-h,o-z)
integer,intent(inout),dimension(nl1,4)::nn

    do i=1,nl1
        nn(i,1) = mod(i,nl1)+1 !right neighbor
        nn(i,2) = mod(i-2+nl1,nl1)+1   !left neighbor
        nn(i,3) = mod(i-1+nl,nl1)+1    !bottom neighbor
        nn(i,4) = mod(i-1-nl+nl1,nl1)+1 !top neighbor
        
        if(mod(i-1,nl)==0) nn(i,2)=i+nl-1 !left edge
        if(mod(i,nl)==0) nn(i,1)=i-nl+1 !right edge
    enddo

return
endsubroutine neighbors


! energy calculation of a configuration
subroutine energy(nl,nl1,nz, ey,nn)
implicit real(kind=8)(a-h,o-z)
integer,intent(in),dimension(nl1)::nz
integer,intent(in),dimension(nl1,4)::nn

        ne = 0
        do i=1,nl1
           ns=nz(nn(i,1))+nz(nn(i,2))+nz(nn(i,3))+nz(nn(i,4))
           ne = ne - nz(i)*ns
        enddo
        ey = dfloat(ne)/2.d0


return
end subroutine energy



! one random site
subroutine random_site(nl1,j11)
implicit real(kind=8)(a-h,o-z)

       call random_number(u1)
       j11 = int(nl1*u1)+1

return
end subroutine random_site


subroutine initial_conf(nz,nl1,p)
implicit real(kind=8)(a-h,o-z)
integer,intent(inout),dimension(nl1)::nz

        do i=1,nl1
           call random_number(v)
           if(v<p)then
             nz(i) = 1
           else
             nz(i) = -1
           endif
        enddo
        
return
end subroutine initial_conf




! Generate Random Seed
subroutine randomseed()
implicit real(kind=8)(a-h,o-z)
integer :: values(1:8), k
integer,dimension(:),allocatable :: seed

        call date_and_time(values=values)
        call random_seed(size=k)
        allocate(seed(1:k))
        seed(:) = values(8)
        call random_seed(put=seed)
       
return
end subroutine randomseed

