# Third-Party Licenses

This document lists the licenses of all third-party code and libraries used in `pylimer-tools`.

The main `pylimer-tools` project is licensed under **GPL-3.0-or-later** (see [LICENSE](LICENSE)).

---

## Correlator Implementation

**Location**: `src/pylimer_tools_cpp/calc/Correlator.h` and `Correlator.cpp`

**Author**: Jorge Ramirez (2010)

**License**: MIT License

```
Copyright (c) 2010 Jorge Ramirez

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## Vendor Dependencies

### pybind11

**Location**: `vendor/pybind11/`

**Repository**: https://github.com/pybind/pybind11

**License**: BSD-3-Clause License

```
Copyright (c) 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>, All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

---

### Eigen

**Location**: Fetched via CMake FetchContent (see `vendor/eigen.CMakeLists.cmake`)

**Repository**: https://gitlab.com/libeigen/eigen

**Version**: 5.0.1

**License**: Mozilla Public License 2.0 (MPL-2.0)

This project builds Eigen with `EIGEN_MPL2_ONLY=ON`, ensuring only MPL2-licensed code is used.

Full license text: https://www.mozilla.org/en-US/MPL/2.0/

**Note**: Eigen also has components under other licenses (BSD, LGPL), but we use only the MPL2 subset.

---

### Spectra

**Location**: Fetched via CMake FetchContent (see `vendor/spectra.CMakeLists.cmake`)

**Repository**: https://github.com/yixuan/spectra

**Version**: v1.0.1

**License**: Mozilla Public License 2.0 (MPL-2.0)

Full license text: https://www.mozilla.org/en-US/MPL/2.0/

---

### cereal

**Location**: `vendor/cereal/`

**Repository**: https://github.com/USCiLab/cereal

**License**: BSD-3-Clause License

```
Copyright (c) 2013-2022, Randolph Voorhies, Shane Grant
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

---

### igraph

**Location**: `vendor/igraph/`

**Repository**: https://github.com/igraph/igraph

**License**: GNU General Public License v2.0 (GPL-2.0)

Full license text available at: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html

**Additional components in igraph**:
- GLPK: GPL-3.0
- AMD, COLAMD: BSD-3-Clause  
- MiniSat: MIT License
- PCG: Apache-2.0 or MIT
- Qhull: Custom permissive license
- Infomap: GPL-3.0

---

### NLopt

**Location**: `vendor/nlopt/`

**Repository**: https://github.com/stevengj/nlopt

**License**: Dual-licensed as LGPL-2.1-or-later OR MIT License

By default, NLopt includes code under LGPL-2.1. However, it can be built without the LGPL components (luksan directory), in which case it's MIT licensed.

**MIT License (for non-luksan portions)**:
```
Copyright (c) 2007-2024 Massachusetts Institute of Technology

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```

Full LGPL text: https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html

---

### CMake Modules

**Location**: `vendor/cmake-modules/`

**License**: Boost Software License 1.0

```
Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
```

---

## How to Comply with These Licenses

When distributing `pylimer-tools` or derivative works:

1. **Include this file** (`THIRD_PARTY_LICENSES.md`) with your distribution
2. **Include the main LICENSE file** (GPL-3.0-or-later)
3. **Provide source code** or offer to provide it (GPL requirement)
4. **Preserve copyright notices** in the source files
5. **Document any modifications** you make (GPL requirement)

For more information, see:
- GPL-3.0: https://www.gnu.org/licenses/gpl-3.0.html
- MPL-2.0: https://www.mozilla.org/en-US/MPL/2.0/FAQ/
- BSD/MIT: Generally only require copyright notice preservation
