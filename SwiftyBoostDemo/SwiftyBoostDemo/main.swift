//
//  Created by VT on 11.10.25.
//  Copyright © 2025 Volker Thieme. All rights reserved.
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

import Foundation
import SwiftyBoost
//import CBoostBridge

let g = try Distribution.Gamma<Double>(shape: 3.0, scale: 0.5)
let s = try Distribution.StudentT<Double>(degreesOfFreedom: 400.5)
let f = try Distribution.FisherF(degreesOfFreedom1: 100, degreesOfFreedom2: 200)
print("\(f.entropy)")
print("\(g.kurtosis)")
print("\(g.mean)")
print("\(g.median)")
print("\(g.mode)")
print("\(g.range)")
print("\(g.supportLowerBound)")
print("\(g.supportUpperBound)")
print("\(try s.cdf(1.0))")
print("\(s.kurtosis!)")
