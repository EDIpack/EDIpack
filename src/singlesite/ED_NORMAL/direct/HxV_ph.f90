!HAMILTONIAN OF PHONONS: W0*(number of phonons)
  do i=1,Nloc
     iph = (i-1)/(DimUp*DimDw) + 1
     htmp = one*w0_ph*(iph - 1)
     Hv(i) = Hv(i) + htmp*vin(i)
  enddo
