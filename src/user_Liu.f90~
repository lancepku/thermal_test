subroutine ReadHydro()
include 'precision.inc'
include 'params.inc'

call AdvanceToNextInputLine( 4 ) ! add to the last line of the input file and param.inc !Liu 7/2017                                          \
                                                                                                                                               
read(4,*)if_hydro, nusselt, xmaxdepth,xmaxt,Vp_                                                                                                  

return
end
!*************************
function HydroDiff( j, i )
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
!common /user_Liu/if_hydro,nusselt,xmaxdepth,xmaxt

    iph = iphase(j,i)
    diff = conduct(iph)/den(iph)/cp(iph)
!   if(i.gt.ixh1t+1.and.i.lt.nx-ixh1t-1) then                                                                                                
 do  i= 1,nx-1
       do j= 1,nz-1
    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
    yc = 0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
!    if(aps(j,k).gt.xmaxstr) write(*,*) tmpr,yc,aps(j,k)                                                                                      
      if( tmpr.le.xmaxt .and.yc.ge.xmaxdepth) then
!    xx0 = 0.25*(cord(j,i,1)+cord(j+1,i,1)+cord(j,i+1,1)+cord(j+1,i+1,1))                                                                     
!    xx = 0.25*(cord(j,k,1)+cord(j+1,k,1)+cord(j,k+1,1)+cord(j+1,k+1,1))                                                                      
        HydroDiff = diff * nusselt
        write(*,*) HydroDiff,nusselt
      else
        HydroDiff = diff
      endif
      end do
   end do
return
end
