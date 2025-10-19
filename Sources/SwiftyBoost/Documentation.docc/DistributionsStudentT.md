# ``SwiftyBoost/Distribution/StudentT``

Student’s t distribution with ν degrees of freedom (ν > 0).

## Definition

The probability density function is

f(x; ν) = Γ((ν+1)/2) / (√(νπ) Γ(ν/2)) · (1 + x²/ν)^{−(ν+1)/2}.

The distribution is symmetric about zero. The mode and median are 0. The mean is 0 for ν > 1 (undefined otherwise). The variance is ν/(ν−2) for ν > 2, infinite for 1 < ν ≤ 2, and undefined for ν ≤ 1.

## Type and Precision

``Distribution/StudentT`` is generic over `BinaryFloatingPoint` and available in these typical specializations:

- ``Distribution/StudentT``<`Double`>
- ``Distribution/StudentT``<`Float`>
- ``Distribution/StudentT``<`Float80`> (x86_64 only; falls back to `Double` elsewhere)

Each instance constructs a Boost.Math `students_t_distribution` once and keeps an internal opaque handle, reused for all evaluations.

## API Overview

- Initialization: ``Distribution/StudentT/init(degreesOfFreedom:)``
- Density: ``Distribution/StudentT/pdf(_:)``
- Distribution functions: ``Distribution/StudentT/cdf(_:)`` and ``Distribution/StudentT/sf(_:)``
- Inverse functions: ``Distribution/StudentT/quantile(_:)`` and ``Distribution/StudentT/quantileComplement(_:)``
- Moments: ``Distribution/StudentT/mean``, ``Distribution/StudentT/variance``, ``Distribution/StudentT/mode``
- Planning helper: ``Distribution/StudentT/findDegreesOfFreedom(differenceFromMean:alpha:beta:sd:hint:)``

## Usage

```swift
import SwiftyBoost

let t = try Distribution.StudentT<Double>(degreesOfFreedom: 5)

let pdf0 = t.pdf(0.0)        // peak at 0
let cdf2 = t.cdf(2.0)
let q975 = t.quantile(0.975) // two-sided 95% ≈ ±q975

// Power planning (Boost helper):
let nu = Distribution.StudentT<Double>.findDegreesOfFreedom(
  differenceFromMean: 1.0, alpha: 0.05, beta: 0.2, sd: 1.0, hint: 5.0)
```

## Mathematical Notes

- Symmetric: mode = median = 0
- Mean: 0 for ν > 1; undefined for ν ≤ 1
- Variance: ν/(ν−2) for ν > 2; ∞ for 1 < ν ≤ 2; undefined for ν ≤ 1
- Support: x ∈ (−∞, ∞)

## Implementation Details

All operations delegate to Boost via the C bridge. The Swift type holds an opaque pointer to a `boost::math::students_t_distribution` in the matching precision. Calls route through stable `bs_` functions (e.g., `bs_student_t_pdf_h`, `bs_student_t_cdf_h`, `bs_student_t_quantile_h`). The planning utility forwards to Boost’s static method `students_t_distribution<>::find_degrees_of_freedom`.
