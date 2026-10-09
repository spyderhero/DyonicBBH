    ! Conformal Christoffel symbols of h_ij = tilde gamma_ij.
    cf1 = 0.0d0
    do a = 1, 3
      do b = 1, 3
        do c = b, 3
          cf1(a,b,c) = 0.5d0 * (d1_hh(a,b,c) + d1_hh(a,c,b) - d1_hh(b,c,a))
        end do
      end do
    end do
    cf1(:,2,1) = cf1(:,1,2)
    cf1(:,3,1) = cf1(:,1,3)
    cf1(:,3,2) = cf1(:,2,3)

    cf2 = 0.0d0
    do a = 1, 3
      do b = 1, 3
        do c = b, 3
          do m = 1, 3
            cf2(a,b,c) = cf2(a,b,c) + hu(a,m) * cf1(m,b,c)
          end do
        end do
      end do
    end do
    cf2(:,2,1) = cf2(:,1,2)
    cf2(:,3,1) = cf2(:,1,3)
    cf2(:,3,2) = cf2(:,2,3)

    ! Physical divergence D_i V^i for gamma_ij = chi^{-p} h_ij.
    ! D_i V^i = partial_i V^i + tildeGamma^i_{m i} V^m - (3p/2) V^i partial_i chi / chi.
    divE = 0.0d0
    divB = 0.0d0
    do a = 1, 3
      divE = divE + d1_lE(a,a) - 1.5d0*conf_fac_exponent*lE(a)*d1_ch(a)/ch
      divB = divB + d1_lB(a,a) - 1.5d0*conf_fac_exponent*lB(a)*d1_ch(a)/ch
      do m = 1, 3
        divE = divE + cf2(a,m,a) * lE(m)
        divB = divB + cf2(a,m,a) * lB(m)
      end do
    end do
