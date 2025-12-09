  !Diagonal Elements, i.e. local part
  htmp = zero
  !
  do iorb=1,Norb
    htmp = htmp + (impHloc(1,1,iorb,iorb) + mfHloc(1,1,iorb,iorb))*nup(iorb)
    htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*ndw(iorb)
    htmp = htmp + mfHloc(2,2,iorb,iorb) * ndw(iorb)!Needed because these terms come from the anticommutator
    htmp = htmp - xmu*(nup(iorb)+ndw(iorb))
  enddo
  !
  hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
  !
  !Off-diagonal elements, i.e. non-local part
  !1. same spin:
  do iorb=1,Norb
     do jorb=1,Norb
        !UP
       Jcondition = &
            ((impHloc(1,1,iorb,jorb)/=zero) .OR. &
             (mfHloc(1,1,iorb,jorb)/=zero)).AND. &
            (ib(jorb)==1)                  .AND. &
            (ib(iorb)==0)
        if (Jcondition) then
           call c(jorb,m,k1,sg1)
           call cdg(iorb,k1,k2,sg2)
           j_el = binary_search(Hsector%H(1)%map,k2)
           j    = j_el + (iph-1)*DimEl
           htmp = (impHloc(1,1,iorb,jorb)+mfHloc(1,1,iorb,jorb))*sg1*sg2
           !
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
           !
        endif
        !DW
       Jcondition = &
            ((impHloc(Nspin,Nspin,iorb,jorb)/=zero) .OR. &
             (mfHloc(2,2,iorb,jorb)/=zero)).AND. &
            (ib(jorb+Ns)==1)                       .AND. &
            (ib(iorb+Ns)==0)
        if (Jcondition) then
           call c(jorb+Ns,m,k1,sg1)
           call cdg(iorb+Ns,k1,k2,sg2)
           j_el = binary_search(Hsector%H(1)%map,k2)
           j    = j_el + (iph-1)*DimEl
           htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
           htmp = htmp + mfHloc(2,2,iorb,jorb)*sg1*sg2
           !
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
           !
        endif
     enddo
  enddo
  !
  !Evaluate: Fd . D = Fd . (C^+_{a,up}C^+_{b,dw} + C_{b,dw}C_{a,up})
  if(any(pair_field/=0d0).or.any(abs(impHloc_anomalous)/=0d0))then
     do iorb=1,Norb
       do jorb=1,Norb
          !
          Jcondition = (ib(iorb)==1) .AND. (ib(jorb+Ns)==1)
          if(Jcondition)then
             call c(iorb,m,k1,sg1)
             call c(jorb+Ns,k1,k2,sg2)
             j_el = binary_search(Hsector%H(1)%map,k2)
             j    = j_el + (iph-1)*DimEl
             htmp=one*impHloc_anomalous(1,1,iorb,jorb)*sg1*sg2
             if(iorb == jorb) htmp= htmp + one*pair_field(iorb)*sg1*sg2
             !
             hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
             !
          endif
          !
          Jcondition = (ib(iorb)==0) .AND. (ib(jorb+Ns)==0)
          if(Jcondition)then
             call cdg(jorb+Ns,m,k1,sg1)
             call cdg(iorb,k1,k2,sg2)
             j_el = binary_search(Hsector%H(1)%map,k2)
             j    = j_el + (iph-1)*DimEl
             htmp=one*conjg(impHloc_anomalous(1,1,iorb,jorb))*sg1*sg2
             if(iorb == jorb) htmp= htmp + one*pair_field(iorb)*sg1*sg2
             !
             hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
             !
          endif
       enddo
     enddo
  endif
