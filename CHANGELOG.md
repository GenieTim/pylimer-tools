# Changelog

## main

- Enable switching between the unit types in the `generate_network` command
- (Possibly) fix compilation if Eigen3 is not installed

## v0.3.2

- Add methods to the `MCUniverseGenerator` to spawn star-like crosslinkers (`add_star_crosslinkers`) as well as link those to each other (`link_strands_to_strands_to_conversion`)
- Implement Python binding for the `MCUniverseGenerator`'s `link_strands_callback`, i.e., allow the Python users to stear and stop the cross-linking procedure
- Switch from `MEHPForceRelaxation` to `MEHPForceBalance2` internally of the `MCUniverseGenerator` to improve performance
- Allow `MEHPForceRelaxation` and `MEHPForceBalance2` to compute soluble fractions even when the Universe is empty

## v0.3.1

This version mainly fixes some issues in the CI.

## v0.3.0

This has become the new baseline version, since it's the first public version since the very first try to publish many years ago.
There are only breaking changes, compared to previous versions, such that a CHANGELOG has become useless.
Sorry for that, please enjoy the many new features by taking a look at the revised documentation.
