#!/usr/bin/env python
"""
Custom Crosslinking with Callback Functions
============================================

This example demonstrates how to use the :meth:`~pylimer_tools_cpp.MCUniverseGenerator.link_strands_callback`
method to create custom crosslinking procedures
as an alternative to the built-in linking methods,
:meth:`~pylimer_tools_cpp.MCUniverseGenerator.link_strands_to_conversion` and
:meth:`~pylimer_tools_cpp.MCUniverseGenerator.link_strands_to_soluble_fraction`.
The callback function allows users to control when the crosslinking process should stop based on
any criteria they can implement.

The callback function receives two parameters:
- The current MCUniverseGenerator instance
- The current step number

And should return one of the BackTrackStatus values:
- STOP: Stop the linking process
- TRACK_FORWARD: Continue linking forward
- TRACK_BACKWARD: Track backward in the linking process
"""

import matplotlib.pyplot as plt

from pylimer_tools_cpp import BackTrackStatus, MCUniverseGenerator, MEHPForceBalance2

# %%
# Example 1: Stop at Conversion OR Soluble Fraction
# -------------------------------------------------
#
# Using a callback allows for more complex logic, such as
# stopping the process based on multiple criteria as illustrated here.

generator = MCUniverseGenerator(25.0, 25.0, 25.0)
generator.set_seed(42)

# Add crosslinkers and strands
generator.add_crosslinkers(25, 4, 2)  # 25 crosslinkers, functionality 4
generator.add_strands(50, [10] * 50, 1)  # 50 strands, 10 beads each

target_conversion = 0.9  # Stop at 90% conversion
target_wsol = 0.2  # Stop if soluble fraction is below 20%

# record data throughout the process
conversions = []
wsol_values = []
steps = []


def conversion_callback(gen, step):
    """
    Callback function to monitor crosslinking progress.
    This function checks the current crosslinker conversion and soluble fraction,
    and decides whether to stop the linking process.

    .. todo:
      Implement the search in a binary search fashion to improve efficiency
    """
    assert isinstance(gen, MCUniverseGenerator), (
        "Callback must receive a MCUniverseGenerator instance"
    )
    fb2 = gen.get_force_balance2()
    assert isinstance(
        fb2, MEHPForceBalance2), "Force balance must be MEHPForceBalance2"
    fb2.run_force_relaxation()
    current_wsol = fb2.get_soluble_weight_fraction()
    current_conversion = gen.get_current_crosslinker_conversion()
    conversions.append(current_conversion)
    steps.append(step)
    wsol_values.append(current_wsol)

    if current_conversion >= target_conversion or current_wsol <= target_wsol:
        print(
            f"  Target conversion {
                target_conversion:.1%} or soluble fraction {
                target_wsol:.1%} reached at step {step}"
        )
        return BackTrackStatus.STOP
    else:
        return BackTrackStatus.TRACK_FORWARD


generator.link_strands_callback(conversion_callback, 1.0)

final_conversion = generator.get_current_crosslinker_conversion()
print(
    f"  Final conversion: {final_conversion:.3f} and w_sol: {wsol_values[-1]:.3f}")

# %%
# Visualizing the Results
# ----------------------
#
# Let's plot the evolution of crosslinker conversion
# and soluble fraction over the steps.

fig, ax = plt.subplots()

ax.plot(steps, [c * 100 for c in conversions], label="Conversion")
ax.plot(steps, [wsol * 100 for wsol in wsol_values], label="Soluble Fraction")
ax.set_xlabel("Step")
ax.set_ylabel("Percent [%]")
ax.set(xlim=(0, max(steps)), ylim=(0, 100))
ax.legend()

plt.show()
