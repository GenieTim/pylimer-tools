# ComputedDoubleValues

### *class* pylimer_tools_cpp.ComputedDoubleValues(\*values)

Bases: [`IntEnum`](https://docs.python.org/3/library/enum.html#enum.IntEnum)

Floating point output quantities

### Attributes Summary

| [`GAMMA`](#pylimer_tools_cpp.ComputedDoubleValues.GAMMA)             |                                                      |
|----------------------------------------------------------------------|------------------------------------------------------|
| [`MAX_B`](#pylimer_tools_cpp.ComputedDoubleValues.MAX_B)             |                                                      |
| [`MEAN_B`](#pylimer_tools_cpp.ComputedDoubleValues.MEAN_B)           |                                                      |
| [`MSD`](#pylimer_tools_cpp.ComputedDoubleValues.MSD)                 |                                                      |
| [`PRESSURE`](#pylimer_tools_cpp.ComputedDoubleValues.PRESSURE)       |                                                      |
| [`RESIDUAL`](#pylimer_tools_cpp.ComputedDoubleValues.RESIDUAL)       |                                                      |
| [`STRESS_NXY`](#pylimer_tools_cpp.ComputedDoubleValues.STRESS_NXY)   |                                                      |
| [`STRESS_NXZ`](#pylimer_tools_cpp.ComputedDoubleValues.STRESS_NXZ)   |                                                      |
| [`STRESS_NYZ`](#pylimer_tools_cpp.ComputedDoubleValues.STRESS_NYZ)   |                                                      |
| [`STRESS_XX`](#pylimer_tools_cpp.ComputedDoubleValues.STRESS_XX)     |                                                      |
| [`STRESS_XY`](#pylimer_tools_cpp.ComputedDoubleValues.STRESS_XY)     |                                                      |
| [`STRESS_XZ`](#pylimer_tools_cpp.ComputedDoubleValues.STRESS_XZ)     |                                                      |
| [`STRESS_YY`](#pylimer_tools_cpp.ComputedDoubleValues.STRESS_YY)     |                                                      |
| [`STRESS_YZ`](#pylimer_tools_cpp.ComputedDoubleValues.STRESS_YZ)     |                                                      |
| [`STRESS_ZZ`](#pylimer_tools_cpp.ComputedDoubleValues.STRESS_ZZ)     |                                                      |
| [`TEMPERATURE`](#pylimer_tools_cpp.ComputedDoubleValues.TEMPERATURE) |                                                      |
| [`TIME`](#pylimer_tools_cpp.ComputedDoubleValues.TIME)               |                                                      |
| [`TIMESTEP`](#pylimer_tools_cpp.ComputedDoubleValues.TIMESTEP)       |                                                      |
| [`VOLUME`](#pylimer_tools_cpp.ComputedDoubleValues.VOLUME)           |                                                      |
| [`denominator`](#pylimer_tools_cpp.ComputedDoubleValues.denominator) | the denominator of a rational number in lowest terms |
| [`imag`](#pylimer_tools_cpp.ComputedDoubleValues.imag)               | the imaginary part of a complex number               |
| [`numerator`](#pylimer_tools_cpp.ComputedDoubleValues.numerator)     | the numerator of a rational number in lowest terms   |
| [`real`](#pylimer_tools_cpp.ComputedDoubleValues.real)               | the real part of a complex number                    |

### Methods Summary

| [`as_integer_ratio`](#pylimer_tools_cpp.ComputedDoubleValues.as_integer_ratio)()               | Return a pair of integers, whose ratio is equal to the original int.       |
|------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------|
| [`bit_count`](#pylimer_tools_cpp.ComputedDoubleValues.bit_count)()                             | Number of ones in the binary representation of the absolute value of self. |
| [`bit_length`](#pylimer_tools_cpp.ComputedDoubleValues.bit_length)()                           | Number of bits necessary to represent self in binary.                      |
| [`conjugate`](#pylimer_tools_cpp.ComputedDoubleValues.conjugate)                               | Returns self, the complex conjugate of any int.                            |
| [`from_bytes`](#pylimer_tools_cpp.ComputedDoubleValues.from_bytes)(bytes[, byteorder, signed]) | Return the integer represented by the given array of bytes.                |
| [`is_integer`](#pylimer_tools_cpp.ComputedDoubleValues.is_integer)()                           | Returns True.                                                              |
| [`to_bytes`](#pylimer_tools_cpp.ComputedDoubleValues.to_bytes)([length, byteorder, signed])    | Return an array of bytes representing an integer.                          |

### Attributes Documentation

#### GAMMA *= 14*

#### MAX_B *= 17*

#### MEAN_B *= 16*

#### MSD *= 18*

#### PRESSURE *= 3*

#### RESIDUAL *= 15*

#### STRESS_NXY *= 11*

#### STRESS_NXZ *= 13*

#### STRESS_NYZ *= 12*

#### STRESS_XX *= 5*

#### STRESS_XY *= 8*

#### STRESS_XZ *= 10*

#### STRESS_YY *= 6*

#### STRESS_YZ *= 9*

#### STRESS_ZZ *= 7*

#### TEMPERATURE *= 4*

#### TIME *= 1*

#### TIMESTEP *= 0*

#### VOLUME *= 2*

#### denominator

the denominator of a rational number in lowest terms

#### imag

the imaginary part of a complex number

#### numerator

the numerator of a rational number in lowest terms

#### real

the real part of a complex number

### Methods Documentation

#### as_integer_ratio()

Return a pair of integers, whose ratio is equal to the original int.

The ratio is in lowest terms and has a positive denominator.

```pycon
>>> (10).as_integer_ratio()
(10, 1)
>>> (-10).as_integer_ratio()
(-10, 1)
>>> (0).as_integer_ratio()
(0, 1)
```

#### bit_count()

Number of ones in the binary representation of the absolute value of self.

Also known as the population count.

```pycon
>>> bin(13)
'0b1101'
>>> (13).bit_count()
3
```

#### bit_length()

Number of bits necessary to represent self in binary.

```pycon
>>> bin(37)
'0b100101'
>>> (37).bit_length()
6
```

#### conjugate()

Returns self, the complex conjugate of any int.

#### *classmethod* from_bytes(bytes, byteorder='big', , signed=False)

Return the integer represented by the given array of bytes.

bytes
: Holds the array of bytes to convert.  The argument must either
  support the buffer protocol or be an iterable object producing bytes.
  Bytes and bytearray are examples of built-in objects that support the
  buffer protocol.

byteorder
: The byte order used to represent the integer.  If byteorder is ‘big’,
  the most significant byte is at the beginning of the byte array.  If
  byteorder is ‘little’, the most significant byte is at the end of the
  byte array.  To request the native byte order of the host system, use
  <br/>
  ```
  `
  ```
  <br/>
  sys.byteorder’ as the byte order value.  Default is to use ‘big’.

signed
: Indicates whether two’s complement is used to represent the integer.

#### is_integer()

Returns True. Exists for duck type compatibility with float.is_integer.

#### to_bytes(length=1, byteorder='big', , signed=False)

Return an array of bytes representing an integer.

length
: Length of bytes object to use.  An OverflowError is raised if the
  integer is not representable with the given number of bytes.  Default
  is length 1.

byteorder
: The byte order used to represent the integer.  If byteorder is ‘big’,
  the most significant byte is at the beginning of the byte array.  If
  byteorder is ‘little’, the most significant byte is at the end of the
  byte array.  To request the native byte order of the host system, use
  <br/>
  ```
  `
  ```
  <br/>
  sys.byteorder’ as the byte order value.  Default is to use ‘big’.

signed
: Determines whether two’s complement is used to represent the integer.
  If signed is False and a negative integer is given, an OverflowError
  is raised.
