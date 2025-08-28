"""
Lightweight data loader for polymer properties.

This module provides functionality to load polymer properties data
without requiring pandas or openpyxl as runtime dependencies.
"""

import json
import os
from typing import Dict, List, Any


class PolymerData:
    """
    A lightweight container for polymer property data that mimics
    the basic functionality needed from pandas DataFrame/Series.
    """

    def __init__(self, data_dict: Dict[str, Any]):
        """Initialize with polymer property dictionary."""
        self._data = data_dict

    def __getattr__(self, name: str) -> Any:
        """Allow attribute-style access to properties."""
        if name in self._data:
            return self._data[name]
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )

    def __getitem__(self, key: str) -> Any:
        """Allow dictionary-style access to properties."""
        return self._data[key]

    def get(self, key: str, default: Any = None) -> Any:
        """Get property value with optional default."""
        return self._data.get(key, default)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return self._data.copy()


class PolymerDataFrame:
    """
    A lightweight container that mimics the basic DataFrame functionality
    needed for polymer properties data.
    """

    def __init__(self, polymers: List[Dict[str, Any]], columns: List[str]):
        """Initialize with list of polymer dictionaries and column names."""
        self._polymers = polymers
        self._columns = columns

    def itertuples(self, index: bool = True, name: str = "Pandas"):
        """
        Iterate over DataFrame rows as named tuples.

        This mimics pandas.DataFrame.itertuples() behavior.
        """
        from collections import namedtuple

        # Create namedtuple class with the polymer properties
        PolymerTuple = namedtuple(
            name, ["Index"] + self._columns if index else self._columns
        )

        for i, polymer in enumerate(self._polymers):
            values = [polymer.get(col) for col in self._columns]
            if index:
                yield PolymerTuple(i, *values)
            else:
                yield PolymerTuple(*values)

    def iterrows(self):
        """
        Iterate over DataFrame rows as (index, Series) pairs.

        This mimics pandas.DataFrame.iterrows() behavior.
        """
        for i, polymer in enumerate(self._polymers):
            yield i, PolymerData(polymer)

    def __getitem__(self, column: str) -> List[Any]:
        """Get a column as a list."""
        return [polymer.get(column) for polymer in self._polymers]

    @property
    def columns(self) -> List[str]:
        """Get column names."""
        return self._columns.copy()

    def unique(self, column: str) -> List[Any]:
        """Get unique values from a column."""
        values = self[column]
        return list(set(value for value in values if value is not None))


def load_everaers_et_al_data() -> PolymerDataFrame:
    """
    Load the Everaers et al. (2020) unit properties data from JSON.

    This replaces the pandas.read_excel() functionality with a lightweight
    JSON-based approach that doesn't require external dependencies.

    :return: PolymerDataFrame containing polymer properties
    :rtype: PolymerDataFrame
    """
    # Get the path to the JSON file
    data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
    json_path = os.path.join(data_dir, "everaers_et_al_unit_properties.json")

    if not os.path.exists(json_path):
        raise FileNotFoundError(
            f"Polymer properties data file not found at {json_path}. "
            f"Run 'python bin/convert-excel-to-json.py' to generate it from the Excel source."
        )

    try:
        with open(json_path, "r") as f:
            data = json.load(f)

        return PolymerDataFrame(data["polymers"], data["columns"])

    except (json.JSONDecodeError, KeyError) as e:
        raise RuntimeError(f"Failed to load polymer properties data: {e}")


def get_available_polymers() -> List[str]:
    """
    List all available polymers for which we have LJ unit conversions.

    :return: List of polymer names
    :rtype: List[str]
    """
    data = load_everaers_et_al_data()
    return data.unique("name")


def get_polymer_by_name(polymer_name: str) -> PolymerData:
    """
    Get polymer data by name with fuzzy matching.

    :param polymer_name: Name of the polymer to find
    :type polymer_name: str
    :return: Polymer data object
    :rtype: PolymerData
    :raises ValueError: If polymer not found
    """
    data = load_everaers_et_al_data()

    # Normalize the input name for comparison
    normalized_input = "".join(filter(str.isalnum, polymer_name)).lower()

    for polymer_dict in data._polymers:
        polymer_data_name = polymer_dict.get("name", "")
        normalized_name = "".join(
            filter(
                str.isalnum,
                str(polymer_data_name))).lower()

        if normalized_input == normalized_name:
            return PolymerData(polymer_dict)

    available = get_available_polymers()
    raise ValueError(
        f"Polymer '{polymer_name}' not found. Available polymers: {available}"
    )
