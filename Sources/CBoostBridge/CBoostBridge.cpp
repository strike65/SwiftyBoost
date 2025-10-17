//
//  Created by VT on 11.10.25.
//  Copyright © 2025 Volker Thieme.
//  Licensed under the same terms as the rest of the project.
//
//  This file composes the C ABI bridge from smaller implementation units.

#include "CBoostBridge.h"

#include "internal/bs_internal.hpp"

// Order of includes loosely mirrors public headers for clarity
#include "impl/constants.inc"
#include "impl/complex.inc"
#include "impl/gamma_error.inc"
#include "impl/beta_digamma_owens.inc"
#include "impl/exp_trig_airy.inc"
#include "impl/bessel_legendre_elliptic.inc"
#include "impl/lambert_hypergeom.inc"
#include "impl/numbers_factorials_polys_splines.inc"

