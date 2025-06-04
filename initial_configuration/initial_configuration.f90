!nl = linear size of the system
!m = no. of iteration or monte carlo steps
!j3 = equilibration steps or transient


implicit real (kind=8)(a-h,o-z)
integer,parameter::n=15,m=2**(n)
integer,parameter::nl=2**7,nl1=nl**2,j3=10**2

        call random_seed
        call ps1(m,nl,nl1,j3)
        
                
stop
end


subroutine ps1(m,nl,nl1,j3)
implicit real(kind=8)(a-h,o-z)
integer,dimension(nl1)::nz
integer,dimension(nl1,4)::nn    !nearest neighbors
real (kind=8),dimension(5)::p1
nz=0
nn=0
        call initial_conf(nz,nl1,0.5d0)
        call neighbors(nl,nl1,nn)
        tem = 2.269185
        call cost_exp(p1,tem)
        call single_spin_flip_metropolis(nl,nl1,nz,m,tem,j3,nn,p1)

                
return
end subroutine ps1

!2d single_spin_flip
subroutine single_spin_flip_metropolis(nl,nl1,nz,m,tem,j3,nn,p1)
implicit real(kind=8)(a-h,o-z)
integer,intent(inout),dimension(nl1)::nz
real (kind=8),dimension(m)::x
integer,dimension(nl1)::nz2,nz1
integer,intent(in),dimension(nl1,4)::nn
real (kind=8),intent(in),dimension(5)::p1
real (kind=8),dimension(nl)::xc,xc1
x=0.d0
nz2 = 0
nz1=0
        !kc = 0
        do it = 1,m+j3
           call conf_update1(nl,nl1,nz,tem,nc,ns,c,nz2,nn,p1) ! GD
           if(it>j3)then
           nz1 = nz1 + nz2
           i0 = count(nz1>0)
           if(i0==nl1)then
             it0 = it-j3
             goto 10
           endif
           x(it-j3) =dfloat(nc)
           endif
        enddo
10      continue
        print*,it0
        call final_conf1(nl,nl1,nz)
        call final_conf(nl,nl1,nz2)

        
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


subroutine final_conf1(nl,nl1,nz)
implicit real(kind=8)(a-h,o-z)
integer,intent(inout),dimension(nl1)::nz

        do i=1,nl1
           ix = mod(i,nl)
           if(ix==0)ix = nl
           iy = int(dfloat(i-ix)/dfloat(nl)) + 1
           if(nz(i)==1)write(10,*)ix,iy
        enddo
        
return
end subroutine final_conf1

subroutine final_conf(nl,nl1,nz)
implicit real(kind=8)(a-h,o-z)
integer,intent(inout),dimension(nl1)::nz

        do i=1,nl1
           ix = mod(i,nl)
           if(ix==0)ix = nl
           iy = int(dfloat(i-ix)/dfloat(nl)) + 1
           if(nz(i)==1)write(11,*)ix,iy
        enddo
        
return
end subroutine final_conf



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

