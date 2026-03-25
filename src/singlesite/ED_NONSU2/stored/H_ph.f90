!Phononic hamiltonian H_ph = w0 bdag b + A(bdag + b)
do i=MpiIstart,MpiIend
    i_el = mod(i-1,DimEl) + 1
    iph = (i-1)/DimEl + 1
    !
    ! w0 bdag b
    htmp = w0_ph*(iph-1)
    select case(MpiStatus)
        case (.true.)
            call sp_insert_element(MpiComm,spH0,htmp,i,i)
        case (.false.)
            call sp_insert_element(spH0,htmp,i,i)
    end select
    !
    if(A_ph/=0.d0)then
        if(iph < DimPh) then !bdg = sum_n |n+1> sqrt(n+1)<n|
            htmp = A_ph*sqrt(dble(iph))
            j = i_el + (iph)*DimEl
            call sp_insert_element(spH0_ph,htmp,i,j)
        endif
        if(iph > 1) then !bdg = sum_n |n+1> sqrt(n+1)<n|
            htmp = A_ph*sqrt(dble(iph-1))
            j = i_el + (iph-2)*DimEl
            call sp_insert_element(spH0_ph,htmp,i,j)
        endif
    end if
    !
 enddo
