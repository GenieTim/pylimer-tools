# ComputedIntValues

### *class* pylimer_tools_cpp.ComputedIntValues(\*values)

Bases: [`IntEnum`](https://docs.python.org/3/library/enum.html#enum.IntEnum)

Integer output quantities

### Attributes Summary

| [`NUM_ATOMS`](#pylimer_tools_cpp.ComputedIntValues.NUM_ATOMS)                 |                                                      |
|-------------------------------------------------------------------------------|------------------------------------------------------|
| [`NUM_BONDS`](#pylimer_tools_cpp.ComputedIntValues.NUM_BONDS)                 |                                                      |
| [`NUM_BONDS_TO_FORM`](#pylimer_tools_cpp.ComputedIntValues.NUM_BONDS_TO_FORM) |                                                      |
| [`NUM_EXTRA_ATOMS`](#pylimer_tools_cpp.ComputedIntValues.NUM_EXTRA_ATOMS)     |                                                      |
| [`NUM_EXTRA_BONDS`](#pylimer_tools_cpp.ComputedIntValues.NUM_EXTRA_BONDS)     |                                                      |
| [`NUM_RELOC`](#pylimer_tools_cpp.ComputedIntValues.NUM_RELOC)                 |                                                      |
| [`NUM_SHIFT`](#pylimer_tools_cpp.ComputedIntValues.NUM_SHIFT)                 |                                                      |
| [`STEP`](#pylimer_tools_cpp.ComputedIntValues.STEP)                           |                                                      |
| [`denominator`](#pylimer_tools_cpp.ComputedIntValues.denominator)             | the denominator of a rational number in lowest terms |
| [`imag`](#pylimer_tools_cpp.ComputedIntValues.imag)                           | the imaginary part of a complex number               |
| [`numerator`](#pylimer_tools_cpp.ComputedIntValues.numerator)                 | the numerator of a rational number in lowest terms   |
| [`real`](#pylimer_tools_cpp.ComputedIntValues.real)                           | the real part of a complex number                    |

### Methods Summary

| [`as_integer_ratio`](#pylimer_tools_cpp.ComputedIntValues.as_integer_ratio)()               | Return a pair of integers, whose ratio is equal to the original int.       |
|---------------------------------------------------------------------------------------------|----------------------------------------------------------------------------|
| [`bit_count`](#pylimer_tools_cpp.ComputedIntValues.bit_count)()                             | Number of ones in the binary representation of the absolute value of self. |
| [`bit_length`](#pylimer_tools_cpp.ComputedIntValues.bit_length)()                           | Number of bits necessary to represent self in binary.                      |
| [`conjugate`](#pylimer_tools_cpp.ComputedIntValues.conjugate)                               | Returns self, the complex conjugate of any int.                            |
| [`from_bytes`](#pylimer_tools_cpp.ComputedIntValues.from_bytes)(bytes[, byteorder, signed]) | Return the integer represented by the given array of bytes.                |
| [`is_integer`](#pylimer_tools_cpp.ComputedIntValues.is_integer)()                           | Returns True.                                                              |
| [`to_bytes`](#pylimer_tools_cpp.ComputedIntValues.to_bytes)([length, byteorder, signed])    | Return an array of bytes representing an integer.                          |

### Attributes Documentation

#### NUM_ATOMS *= 3*

#### NUM_BONDS *= 5*

#### NUM_BONDS_TO_FORM *= 7*

#### NUM_EXTRA_ATOMS *= 4*

#### NUM_EXTRA_BONDS *= 6*

#### NUM_RELOC *= 2*

#### NUM_SHIFT *= 1*

#### STEP *= 0*

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
