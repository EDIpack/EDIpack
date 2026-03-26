subroutine ed_get_phon_site(self)
      real(8), dimension(3)          :: Self
      self(1) = dens_ph
      self(2) = X_ph
      self(3) = X2_ph
end subroutine ed_get_phon_site