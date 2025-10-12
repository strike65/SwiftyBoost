// swift-tools-version: 6.2
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "SwiftyBoost",
    defaultLocalization: "en",
    platforms: [.macOS(.v13), .iOS(.v16)],
    products: [
        // Products define the executables and libraries a package produces, making them visible to other packages.
        .library(
            name: "SwiftyBoost",
            targets: ["SwiftyBoost"]
        ),
    ],
    dependencies: [
        // DocC plugin for command-line documentation generation
        .package(url: "https://github.com/apple/swift-docc-plugin", from: "1.3.0")
    ],
    targets: [
        // Targets are the basic building blocks of a package, defining a module or a test suite.
        // Targets can depend on other targets in this package and products from dependencies.
        .target(
                name: "CBoostBridge",
                publicHeadersPath: "include",
                cxxSettings: [
                    .headerSearchPath("../../extern/boost/libs/math/include"),
                    .unsafeFlags([
                        "-std=c++20"
                    ])
                ]
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
