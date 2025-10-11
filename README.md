# SwiftyBoost

A lightweight Swift wrapper around selected Boost.Math special functions via a small C bridge (CBoostBridge). It provides type-safe Swift APIs for common special functions (Gamma family, Beta family, error functions, digamma/polygamma, zeta, Airy, Bessel, Legendre, elliptic integrals in both Legendre and Carlson symmetric forms, Owen’s T, Lambert W, and various numerically stable helpers such as expm1/log1p/powm1).

- Swift-first APIs with generic overloads for BinaryFloatingPoint (Double-backed), plus specialized Float and Float80 (x86_64) overloads where available.
- Input validation with descriptive Swift errors.
- Backed by Boost.Math for numerical robustness.

## Requirements

- Swift 5.9+
- Apple platforms supported by SwiftPM
- For Float80 overloads: x86_64 architecture only

## Installation (Swift Package Manager)

Add the package to your `Package.swift`:

```swift
// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "YourApp",
    platforms: [
        .iOS(.v15), .macOS(.v12), .tvOS(.v15), .watchOS(.v8)
    ],
    dependencies: [
        .package(url: "https://github.com/your-org/SwiftyBoost.git", from: "0.1.0")
    ],
    targets: [
        .target(
            name: "YourApp",
            dependencies: [
                .product(name: "SwiftyBoost", package: "SwiftyBoost")
            ]
        )
    ]
)
