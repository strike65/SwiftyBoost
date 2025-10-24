//
//  Created by Volker Thieme 2025.
//  Copyright © 2025 Volker Thieme.
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
//

// Bring public headers into this TU (aggregate)
#include "include/CBoostBridge.h"

#include "internal/bs_internal.hpp"

// Order of includes loosely mirrors public headers for clarity
#include "impl/special_functions/bs_a_hyper.hxx"
#include "impl/special_functions/bs_airy.hxx"
#include "impl/special_functions/bs_bessel.hxx"
#include "impl/special_functions/bs_beta.hxx"
#include "impl/special_functions/bs_cardinal_b_splines.hxx"
#include "impl/special_functions/bs_chebyshev.hxx"
#include "impl/special_functions/bs_digamma_polygamma_zeta.hxx"
#include "impl/special_functions/bs_elliptic_carlson_forms.hxx"
#include "impl/special_functions/bs_elliptic_legendre_forms.hxx"
#include "impl/special_functions/bs_error.hxx"
#include "impl/special_functions/bs_exponential_and_helpers.hxx"
#include "impl/special_functions/bs_factorials_combinatorics.hxx"
#include "impl/special_functions/bs_gamma.hxx"
#include "impl/special_functions/bs_gegenbauer.hxx"
#include "impl/special_functions/bs_hypergeometric.hxx"
#include "impl/special_functions/bs_laguerre.hxx"
#include "impl/special_functions/bs_lambert_w.hxx"
#include "impl/special_functions/bs_legendre.hxx"
#include "impl/special_functions/bs_numbers.hxx"
#include "impl/special_functions/bs_owens_t.hxx"
#include "impl/special_functions/bs_spherical_harmonics.hxx"
#include "impl/special_functions/bs_trig_helpers.hxx"
// Complex helpers and functions (needed by several special/complex routines)
#include "impl/special_functions/complex.hxx"
// Cardinal sine/cosine helpers (real/complex)
#include "impl/special_functions/bs_sin_cardinal.hxx"

// Unified generic distribution vtable (runtime factory)
#include "impl/distributions/bs_generic_distribution.hxx"
#include "impl/distributions/bs_distribution_helpers.hxx"
