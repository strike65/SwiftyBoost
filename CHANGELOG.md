# Changelog

All notable changes to this project are tracked here, following the principles of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and adhering to [Semantic Versioning](https://semver.org/).

## Unreleased
- Added a generic complex number type `Complex<T: BinaryFloatingPoint>`
  - Conforms to `Sendable`, `Equatable`, `Hashable`, `Codable`; includes `ComplexD`, `ComplexF`, and (x86_64) `ComplexX` typealiases.
  - Stable arithmetic: Smith’s algorithm for division, numerically robust magnitude, principal square root.
  - Mixed scalar operations, polar helpers (`fromPolar`, `phase`), and normalization.
  - Double/Float/(x86_64) Float80 specializations call Boost-backed bridge functions (`bs_cexp`, `bs_clog`, `bs_csin`, `bs_ccos`, `bs_ctan`, `bs_csinh`, `bs_ccosh`, `bs_ctanh`, `bs_catan`).
  - Added README usage section and examples.
- Corrected README usage examples to match public API naming (`errorFunction`, `besselJ(v:x:)`).
- Aligned supported platforms in README with `Package.swift` (macOS, iOS).
- Updated SPM installation URL to `github.com/strike65/SwiftyBoost`.
- Clarified DocC generation and output paths; added a link to the changelog and versioning policy.
- Added GitHub Actions workflow to build DocC and deploy `Docs/` to GitHub Pages on pushes to `main` (see `.github/workflows/docs.yml`).
- Added README badges: documentation link to `https://strike65.github.io/SwiftyBoost/` and Pages workflow status badge.
- Added direct “Browse online” link under Documentation section in README.
- Docs site root now redirects to the main DocC page for a cleaner entry URL (`https://strike65.github.io/SwiftyBoost` → `/documentation/swiftyboost/`); implemented by writing `Docs/index.html` during Pages deployment.

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
