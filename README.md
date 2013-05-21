pyTheta
=======

A non-linear system for computing a particular chemical reaction. It's powered by the [mpmath][1] 
library to allow arbitrary precision in the calculations.


### Build Instructions
 - For arbitrary precision arithmetic, mpmath depends on the [gmpy][2] library.
   Obtain it with `$ pip install gmpy`.
   (Note that you'll also need the [GMP][3] and [MPIR][4] libraries installed on your system.)
 - To install mpmath itself, run `$ pip install mpmath`.
 - To make sure mpmath is using gmpy under the hood, you can run this snippet on the Python interactive console:

```
>>> import mpmath.libmp
>>> mpmath.libmp.BACKEND
  'gmpy'
```

### Usage
On the Linux command line, run the generated program like so:

`$ python theta.py > results.txt`


### MIT License
Copyright (c) 2013 Philip Conrad.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

   [1]: http://code.google.com/p/mpmath/
   [2]: http://code.google.com/p/gmpy/
   [3]: http://gmplib.org/
   [4]: http://www.mpir.org/
