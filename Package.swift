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
        .executable(
            name: "swiftyboost-demo",
            targets: ["SwiftyBoostDemo"]
        ),
    ],
    dependencies: [
        .package(url: "https://github.com/apple/swift-docc-plugin", from: "1.3.0"),
        .package(url: "https://github.com/apple/swift-numerics.git", .upToNextMajor(from: "1.1.1")),
    ],
    targets: [
        .target(
            name: "CBoostBridge",
            publicHeadersPath: "include",
            cSettings: [
                // Für C/ObjC-Header-Suche (falls nötig)
                .headerSearchPath("../../extern/boost/include")
            ],
            cxxSettings: [
                // Für C++-Header-Suche (Boost)
                .unsafeFlags(["-xc++", "-std=c++20"]),
                .headerSearchPath("../../extern/boost/include")
            ],
        ),
        .target(
            name: "SwiftyBoostPrelude",
            dependencies: [
                .product(name: "Numerics", package: "swift-numerics"),
                .product(name: "RealModule", package: "swift-numerics"),
                .product(name: "ComplexModule", package: "swift-numerics"),
                "CBoostBridge"
            ],
            path: "Sources/SwiftyBoostPrelude"
        ),
        .target(
            name: "SwiftyBoost",
            dependencies: [
                "SwiftyBoostPrelude",
                "CBoostBridge",
                .product(name: "ComplexModule", package: "swift-numerics")
            ],
            path: "Sources/SwiftyBoost"
        ),
        .executableTarget(
            name: "SwiftyBoostDemo",
            dependencies: [
                "SwiftyBoost",
                "SwiftyBoostPrelude",
                "CBoostBridge",
                .product(name: "ComplexModule", package: "swift-numerics")
            ],
            path: "Sources/SwiftyBoostDemo"
        ),
        .testTarget(
            name: "SwiftyBoostTests",
            dependencies: ["SwiftyBoost", "CBoostBridge"],
            exclude: [
                "TestHelpers/analyze_sample.m"
            ],
            resources: [
                .process("TestData")
            ]
        ),
    ],
    cxxLanguageStandard: .cxx20
)
