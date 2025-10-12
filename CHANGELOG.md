# Changelog

All notable changes to this project are tracked here, following the principles of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and adhering to [Semantic Versioning](https://semver.org/).

## Unreleased
- Corrected README usage examples to match public API naming (`errorFunction`, `besselJ(v:x:)`).
- Aligned supported platforms in README with `Package.swift` (macOS, iOS).
- Updated SPM installation URL to `github.com/strike65/SwiftyBoost`.
- Clarified DocC generation and output paths; added a link to the changelog and versioning policy.

## [0.1.0] - 2025-10-12
- Expanded Boost-backed special function coverage, including additional elliptic, hypergeometric, and Lambert W variants.
- Added the initial Swift Testing suite covering factorial and combinatorics APIs.
- Removed the experimental `ellint3` exposure pending stabilization.

## [0.0.1] - 2025-10-11
- Established the `SwiftyBoost` Swift package with the `CBoostBridge` interop target.
- Vendored Boost.Math headers under `extern/boost`.
- Published core Gamma, Beta, error, Legendre, and Bessel function wrappers with Swift-friendly APIs and error handling.

[0.1.0]: https://github.com/strike65/SwiftyBoost/releases/tag/0.1.0
[0.0.1]: https://github.com/strike65/SwiftyBoost/releases/tag/0.0.1
