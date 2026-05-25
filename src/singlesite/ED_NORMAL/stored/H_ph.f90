  !Phononic hamiltonian H_ph = w0 b^+ b + A(bdg + b)
  do iph=1,DimPh
     htmp = w0_ph*(iph-1)
     call sp_insert_element(spH0_ph,htmp,iph,iph)
  enddo
  if(A_ph/=0.d0)then
     do iph=1,DimPh
        if(iph < DimPh) then !bdg = sum_n |n+1> sqrt(n+1)<n|
           htmp = A_ph*sqrt(dble(iph))
           call sp_insert_element(spH0_ph,htmp,iph+1,iph)
        endif
        if(iph > 1) then !bdg = sum_n |n+1> sqrt(n+1)<n|
           htmp = A_ph*sqrt(dble(iph-1))
           call sp_insert_element(spH0_ph,htmp,iph-1,iph)
        endif
     enddo
  end if
