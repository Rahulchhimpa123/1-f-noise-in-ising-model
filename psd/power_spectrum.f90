



!-----------------------------------------------------------------------------

!This is the Danielson and Lanczos implementation of the fast Fourier transform as described in Numerical Recipes, Press et al in section 12.2.  It has been tested by comparing with the original COOLEY-TUKEY TRANSFORM, which is a fortran 4 implementation of the same code. Transform(k) = sum(data(j)*exp(isign* 2pi*sqrt(-1)*(j-1)*(k-1)/nn)). Summed over all j and k from 1 to nn. Data is in a one-dimensional complex array (i.e., the real and imaginary parts are adjacent in storage, such as fortran IV places them) whose length nn=2**k, k.GE.0 (If necessary append zeroes to the data). isign is +1 or -1. If a -1 transform is followed by a +1 (or a +1 by a -1) the original data reappears, multiplied by nn. Transform values are returned in array data, replacing the input.
subroutine fft(data,nn,isign)
implicit real(kind=8)(a-h,o-z)
real (kind=8),intent(inout),dimension(2*nn)::data

! bit reversal task
        n = 2*nn                                                           
        j = 1                                                               
        do i=1,n,2                                                       
           if(j>i)then                                                 
             tempr = data(j)                                                
             tempi = data(j+1)                                              
             data(j) = data(i)                                              
             data(j+1) = data(i+1)                                          
             data(i) = tempr
             data(i+1) = tempi                                               
           endif                                                          
           m = n/2                                                          
1          if((m>=2).and.(j>m))then                                  
             j = j-m                                                        
             m = m/2                                                        
             goto 1                                                       
           endif                                                           
           j = j + m                                                          
        enddo             
                                               
! Here begins the Danielson-Lanczos section (outer loop executed Log2 (nn) times
                                                  
        mmax = 2                                                           
2       if(n>mmax)then                                                
          istep = 2*mmax  
          pi = 4.d0*datan(1.d0)                                                  
          theta = (2.d0*pi)/dfloat(isign*mmax)                          
          wpr = -2*(dsin(0.5d0*theta))**2                                    
          wpi = dsin(theta)                                                
          wr = 1.d0                                                            
          wi = 0.d0                                                           
          do m=1,mmax,2                                                  
             do i=m,n,istep                                                
                j = i + mmax                                                   
                tempr = wr*data(j) - wi*data(j+1)                              
                tempi = wr*data(j+1) + wi*data(j)                            
                data(j) = data(i) - tempr                                      
                data(j+1) = data(i+1) - tempi                                  
                data(i) = data(i) + tempr         
                data(i+1) = data(i+1) + tempi  
             enddo                                                        
             wtemp =  wr 
             wr = wr*wpr - wi*wpi + wr                 
             wi = wi*wpr + wtemp*wpi + wi        
          enddo                                                          
          mmax = istep                                                     
          goto 2                                                          
        endif                                                            
       
return
end subroutine fft 

!arranging real and imaginary parts in a single coulum vector, x2(i), in an adjacent manner i.e. odd entries === real part and even entries == imaginary part 
subroutine arrange(m,x1,x2) 
implicit real(kind=8)(a-h,o-z)
real (kind=8),intent(inout),dimension(m,2)::x1
real (kind=8),intent(inout),dimension(2*m)::x2 
 
        k = 0
        do i=1,2*m-1,2
           k = k + 1	! counter 
           x2(i) = x1(k,1) 
        enddo
        k = 0
        do i=2,2*m,2
           k = k + 1 
           x2(i) = x1(k,2) 
        enddo
                   
return
end subroutine arrange


!rearranging real and imaginary part
subroutine rearrange(m,x2,x1) 
implicit real(kind=8)(a-h,o-z)
real (kind=8),intent(inout),dimension(m,2)::x1
real (kind=8),intent(inout),dimension(2*m)::x2 
 
        k = 0
        do i=1,2*m-1,2
           k = k + 1 
           x1(k,1) = x2(i)
        enddo
        k = 0
        do i=2,2*m,2
           k = k + 1 
           x1(k,2) = x2(i) 
        enddo
               
return
end subroutine rearrange


