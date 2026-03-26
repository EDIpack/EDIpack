subroutine ed_get_phon_site(self)
      real(8), dimension(3)          :: Self
      self(1) = ed_dens_ph
      self(2) = ed_X_ph
      self(3) = ed_X2_ph
end subroutine ed_get_phon_site