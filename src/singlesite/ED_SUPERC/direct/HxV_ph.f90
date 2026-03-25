!HAMILTONIAN OF PHONONS: W0*(number of phonons)
  htmp = w0_ph*(iph - 1)
  Hv(i-MpiIshift) = Hv(i-MpiIshift) + htmp*vin(i)


  if(A_ph/=0.d0)then
    if(iph < DimPh) then !bdg = sum_n |n+1> sqrt(n+1)<n|
        htmp = A_ph*sqrt(dble(iph))
        i = i_el + (iph)*DimEl
        hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)

    endif
    if(iph > 1) then !bdg = sum_n |n+1> sqrt(n+1)<n|
        htmp = A_ph*sqrt(dble(iph-1))
        i = i_el + (iph-2)*DimEl
        hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)
    endif
  end if