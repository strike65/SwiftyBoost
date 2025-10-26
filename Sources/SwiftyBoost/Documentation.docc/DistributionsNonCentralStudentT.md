#
``SwiftyBoost/Distribution/NonCentralStudentT``

Non-central Student’s t distribution with degrees of freedom `ν` and non-centrality parameter `λ`.

## Definition

The non-central t distribution generalises the Student’s t law to the case where the underlying normal component has non-zero mean. It is supported on `x ∈ (−∞, +∞)` with density

```
f(x; ν, λ) = \frac{ν^{ν/2}}{\sqrt{π} Γ(ν/2)} e^{-(λ^2)/2}
  \int_{0}^{∞} \frac{1}{t^{1/2}} e^{-t} (t + λ^2/2)^{(ν−1)/2}
  \exp\!\left(-\frac{(x \sqrt{ν} - λ)^2}{2(t + λ^2/2)}\right) dt,
```

where `Γ` denotes the gamma function. Parameters satisfy `ν > 0` and `λ` may be any finite real number.

## Type and Precision

``Distribution/NonCentralStudentT`` is generic over `BinaryFloatingPoint`. Concrete instantiations exist for `Float`, `Double`, and (on supported architectures) `Float80`. All distribution queries delegate to Boost.Math via the dynamic factory.

## API Overview

- Initialisation: ``Distribution/NonCentralStudentT/init(degreesOfFreedom:nonCentrality:)``
- Density & log-density: ``Distribution/NonCentralStudentT/pdf(_:)``,
  ``Distribution/NonCentralStudentT/logPdf(_:)``
- Distribution functions: ``Distribution/NonCentralStudentT/cdf(_:)``,
  ``Distribution/NonCentralStudentT/sf(_:)``
- Inverses: ``Distribution/NonCentralStudentT/quantile(_:)``,
  ``Distribution/NonCentralStudentT/quantileComplement(_:)``
- Summary statistics: ``Distribution/NonCentralStudentT/mean``,
  ``Distribution/NonCentralStudentT/variance``,
  ``Distribution/NonCentralStudentT/skewness``,
  ``Distribution/NonCentralStudentT/kurtosis``,
  ``Distribution/NonCentralStudentT/kurtosisExcess``,
  ``Distribution/NonCentralStudentT/mode``,
  ``Distribution/NonCentralStudentT/median``
- Hazards: ``Distribution/NonCentralStudentT/hazard(_:)``,
  ``Distribution/NonCentralStudentT/chf(_:)``
- Divergence: ``Distribution/NonCentralStudentT/klDivergence(relativeTo:options:)``

## Usage

```swift
import SwiftyBoost

let t = try Distribution.NonCentralStudentT<Double>(
    degreesOfFreedom: 8.0,
    nonCentrality: -1.25
)

let x: Double = 1.0
let density = try t.pdf(x)
let cdfValue = try t.cdf(x)
let quantile90 = try t.quantile(0.90)
```

## Mathematical Notes

- Support: `(-∞, +∞)`
- Mean exists for `ν > 1` and equals
  `E[X] = λ \sqrt{ν/2} Γ((ν - 1)/2) / Γ(ν/2)` (mirrored by Boost when defined)
- Variance exists for `ν > 2`
- Reduces to the central Student’s t distribution when `λ = 0`
- No lattice structure: ``Distribution/NonCentralStudentT/latticeStep`` and
  ``Distribution/NonCentralStudentT/latticeOrigin`` return `nil`
- Higher moments are provided by Boost.Math when they exist

## Implementation Details

- Swift performs parameter validation (positive ν, finite λ) before constructing the dynamic delegate.
- PDF/CDF, quantiles, hazards, and moment calculations are provided by Boost.Math.
- Entropy is presently unavailable from the backend and surfaces as `nil`.
