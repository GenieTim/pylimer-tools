# Changelog

## main

- Implement Python binding for the MCUniverseGenerator's `link_strands_callback`, i.e., allow the Python users to stear and stop the cross-linking procedure
- Switch from MEHPForceRelaxation to MEHPForceBalance2 internally of the MCUniverseGenerator to improve performance
- Allow MEHPForceRelaxation and MEHPForceBalance2 to compute soluble fractions even when the Universe is empty

## v0.3.1

This version mainly fixes some issues in the CI.

## v0.3.0

This has become the new baseline version.
There are only breaking changes, compared to previous versions, such that a CHANGELOG has become useless.
Sorry for that, please enjoy the many new features.
