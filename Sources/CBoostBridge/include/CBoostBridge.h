//
//  Created by Volker Thieme 2025.
//  Copyright © 2025 Volker Thieme.
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to
//  deal in the Software without restriction, including without limitation the
//  rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
//  sell copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
//  IN THE SOFTWARE.
//

#pragma once

// Umbrella header for CBoostBridge: composed of focused public headers.

#include "bs_constants.h"

#include "complex.h"
#include "special/bs_error.h"
#include "special/bs_gamma.h"
#include "special/bs_beta.h"
#include "special/bs_digamma_polygamma_zeta.h"
#include "special/bs_owens_t.h"

#include "special/bs_exponential_and_helpers.h"
#include "special/bs_trig_helpers.h"
#include "special/bs_airy.h"
#include "special/bs_bessel.h"
#include "special/bs_legendre.h"
#include "special/bs_elliptic_legendre_forms.h"
#include "special/bs_elliptic_carlson_forms.h"
#include "special/bs_lambert_w.h"

#include "special/bs_hypergeometric.h"
#include "special/bs_numbers.h"
#include "special/bs_factorials_combinatorics.h"

#include "special/bs_laguerre.h"
#include "special/bs_chebyshev.h"
#include "special/bs_spherical_harmonics.h"
#include "special/bs_cardinal_b_splines.h"
#include "special/bs_gegenbauer.h"
#include "special/bs_sin_cardinal.h"
#include "special/bs_a_hyper.h"

// Generic vtable-based distribution factory
#include "distributions/bs_generic_distribution.h"
#include "distributions/bs_student_t_distribution.h"
