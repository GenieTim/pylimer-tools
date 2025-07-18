#!/usr/bin/env python
"""
Monte Carlo Configuration
=========================

Fine-tune the Monte Carlo generation process:
"""

from pylimer_tools_cpp import MCUniverseGenerator

generator = MCUniverseGenerator(
    50.0, 50.0, 50.0
)  # Create a generator for a 50x50x50 simulation box

# Configure the main options of the generator
generator.set_seed(12345)  # Set random seed for reproducibility
generator.set_bead_distance(1.0)  # Set mean bead-to-bead distance

# Configure whether to generate loops or not
generator.config_primary_loop_probability(0.0)  # Primary loops prohibited
generator.config_secondary_loop_probability(0.0)  # Secondary loops prohibited

# Improve performance of sampling pairs of crosslinks for large systems
# Does not sample all pairs of crosslinks anymore,
# only 99.9994% of those that are selected, but is ca. 60 times faster
generator.use_zscore_max_distance(3.29)

# Configure MC parameters
# If you do not trust the Brownian Bridge process we use,
# you can add additional Monte Carlo steps
generator.config_nr_of_mc_steps(5000)  # More MC steps for better equilibration
