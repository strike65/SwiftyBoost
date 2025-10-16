// swift-tools-version: 6.1
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "SwiftyBoost",
    defaultLocalization: "en",
    platforms: [.macOS(.v13), .iOS(.v16)],
    products: [
        .library(
            name: "SwiftyBoost",
            targets: ["SwiftyBoost"]
        ),
    ],
    dependencies: [
        .package(url: "https://github.com/apple/swift-docc-plugin", from: "1.3.0")
    ],
    targets: [
        .target(
            name: "CBoostBridge",
            publicHeadersPath: "include",
            cSettings: [
                // Für C/ObjC-Header-Suche (falls nötig)
                .headerSearchPath("../../extern/boost")
            ],
            cxxSettings: [
                // Für C++-Header-Suche (Boost)
                .headerSearchPath("../../extern/boost")
            ],
            cxxLanguageStandard: .cxx17
        ),
        .target(
            name: "SwiftyBoost",
            dependencies: ["CBoostBridge"]
        ),
        .testTarget(
            name: "SwiftyBoostTests",
            dependencies: ["SwiftyBoost", "CBoostBridge"]
        ),
    ]
)
