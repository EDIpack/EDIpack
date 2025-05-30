  !diagonal, spin conserving:
  do iorb=1,Norb
     do kp=1,Nbath
        ms=getBathStride(iorb,kp)
        !
        ! IMP UP <--> BATH UP
        if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (ib(iorb)==1) .AND. (ib(ms)==0) )then
           call c(iorb,m,k1,sg1)
           call cdg(ms,k1,k2,sg2)
           j_el = binary_search(Hsector%H(1)%map,k2)
           j    = j_el + (iph-1)*DimEl
           htmp = diag_hybr(1,iorb,kp)*sg1*sg2
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
        if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (ib(iorb)==0) .AND. (ib(ms)==1) )then
           call c(ms,m,k1,sg1)
           call cdg(iorb,k1,k2,sg2)
           j_el = binary_search(Hsector%H(1)%map,k2)
           j    = j_el + (iph-1)*DimEl
           htmp = diag_hybr(1,iorb,kp)*sg1*sg2
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
        !
        !IMP DW <--> BATH DW
        if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ib(iorb+Ns)==1) .AND. (ib(ms+Ns)==0) )then
           call c(iorb+Ns,m,k1,sg1)
           call cdg(ms+Ns,k1,k2,sg2)
           j_el = binary_search(Hsector%H(1)%map,k2)
           j    = j_el + (iph-1)*DimEl
           htmp = diag_hybr(Nspin,iorb,kp)*sg1*sg2
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
        if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ib(iorb+Ns)==0) .AND. (ib(ms+Ns)==1) )then
           call c(ms+Ns,m,k1,sg1)
           call cdg(iorb+Ns,k1,k2,sg2)
           j_el = binary_search(Hsector%H(1)%map,k2)
           j    = j_el + (iph-1)*DimEl
           htmp = diag_hybr(Nspin,iorb,kp)*sg1*sg2
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
     enddo
  enddo


