  !Phononic hamiltonian H_ph = w0 b^+ b + A*sqrt(2)* (b^+ + b)
  do i=MpiIstart,MpiIend
     i_el = mod(i-1,DimEl) + 1
     iph = (i-1)/DimEl + 1
     htmp = w0_ph*(iph-1)
     ! Phonon kinetic energy
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0,htmp,i,i)
     end select
     ! Displacement operator
     if(A_ph/=0.d0)then
        if(iph < DimPh) then !bdg = sum_n |n+1> sqrt(n+1)<n|
           j = i_el + (iph)*DimEl
           htmp = A_ph*sqrt(dble(iph))
           select case(MpiStatus)
           case (.true.)
              call sp_insert_element(MpiComm,spH0,htmp,i,j)
           case (.false.)
              call sp_insert_element(spH0,htmp,i,j)
           end select
        endif
        if(iph > 1) then !bdg = sum_n |n+1> sqrt(n+1)<n|
           j = i_el + (iph-2)*DimEl
           htmp = A_ph*sqrt(dble(iph-1))
           select case(MpiStatus)
           case (.true.)
              call sp_insert_element(MpiComm,spH0,htmp,i,j)
           case (.false.)
              call sp_insert_element(spH0,htmp,i,j)
           end select
        endif
     end if
  end do
