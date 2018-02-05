        subroutine ludcmp(a,n,indx,d,code)
        ! By J-P Moreau, Paris
        parameter(nmax=100,tiny=1.5D-16)
        real*8  amax, dum, sum, a(n,n), vv(nmax)
        integer code, d, indx(n)
       
        d=1; code=0
       
        do i=1,n
          amax=0.d0
          do j=1,n
            if (dabs(a(i,j)).gt.amax) amax=dabs(a(i,j))
          end do ! j loop
          if(amax.lt.tiny) then
            code = 1
            return
          end if
          vv(i) = 1.d0/amax
        end do ! i loop
       
        do j=1,n
          do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
              sum=sum-a(i,k)*a(k,j) 
            end do ! k loop
            a(i,j)=sum
          end do ! i loop
          amax=0.d0
          do i=j,n
            sum=a(i,j)
            do k=1,j-1
              sum=sum-a(i,k)*a(k,j) 
            enddo ! k loop
            a(i,j)=sum
            dum=vv(i)*dabs(sum)
            if(dum.ge.amax) then
              imax=i
              amax=dum
            end if
          end do ! i loop  
         
          if(j.ne.imax) then
            do k=1,n
              dum=a(imax,k)
              a(imax,k)=a(j,k)
              a(j,k)=dum
            end do ! k loop
            d=-d
            vv(imax)=vv(j)
          end if
        
          indx(j)=imax
          if(dabs(a(j,j))<tiny) a(j,j)=tiny
        
          if(j.ne.n) then
            dum=1.d0/a(j,j)
            do i=j+1,n
              a(i,j) = a(i,j)*dum
            end do ! i loop
          end if
        end do ! j loop
       
        return
        end subroutine ludcmp
       
        subroutine lubksb(a,n,indx,b)
        ! By J-P Moreau, Paris
        real*8  sum, a(n,n), b(n)
        integer indx(n)
       
        ii=0
       
        do i=1,n
          ll=indx(i)
          sum=b(ll)
          b(ll)=b(i)
          if(ii.ne.0) then
            do j=ii,i-1
              sum=sum-a(i,j)*b(j)
            end do ! j loop
          else if(sum.ne.0.d0) then
            ii=i
          end if
          b(i)=sum
        end do ! i loop
       
        do i=n,1,-1
          sum=b(i)
          if(i<n) then
            do j=i+1,n
              sum=sum-a(i,j)*b(j)
            end do ! j loop
          end if
          b(i)=sum/a(i,i)
        end do ! i loop
       
        return
        end subroutine lubksb
