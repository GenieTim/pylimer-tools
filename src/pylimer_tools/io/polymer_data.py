"""
Lightweight data loader for polymer properties.

This module provides functionality to load polymer properties data
without requiring pandas or openpyxl as runtime dependencies.
"""

import json
import os
from typing import Dict, List, Any, Union


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

    def __getitem__(
        self, key: Union[int, str, List[str]]
    ) -> Union[List[PolymerData], "PolymerDataFrame"]:
        """
        Access columns by name or list of names.
        Access rows by integer index.
        """
        if isinstance(key, str):
            return [polymer.get(key) for polymer in self._polymers]
        elif isinstance(key, list):
            # Multiple column selection
            filtered_polymers = [
                {col: polymer.get(col) for col in key} for polymer in self._polymers
            ]
            return PolymerDataFrame(filtered_polymers, key)
        elif isinstance(key, int):
            # Row selection by index
            if 0 <= key < len(self._polymers):
                return PolymerData(self._polymers[key])
            else:
                raise IndexError("Row index out of range")
        else:
            raise TypeError(
                f"Key must be string, list of strings, or int, got {type(key)}"
            )

    @property
    def columns(self) -> List[str]:
        """Get column names."""
        return self._columns.copy()

    @property
    def index(self) -> List[int]:
        """Get index (row numbers)."""
        return list(range(len(self._polymers)))

    @property
    def shape(self) -> tuple:
        """Get shape as (rows, columns)."""
        return (len(self._polymers), len(self._columns))

    def head(self, n: int = 5) -> "PolymerDataFrame":
        """Get first n rows."""
        return PolymerDataFrame(self._polymers[:n], self._columns)

    def tail(self, n: int = 5) -> "PolymerDataFrame":
        """Get last n rows."""
        return PolymerDataFrame(self._polymers[-n:], self._columns)

    def query(self, expr: str) -> "PolymerDataFrame":
        """
        Simple query functionality. Only supports basic column comparisons.
        Example: df.query("name == 'PS'")
        """
        # Very basic implementation - in real use you might want a proper
        # parser
        filtered_polymers = []

        # Simple parsing for "column operator value" expressions
        expr = expr.strip()
        if "==" in expr:
            col, value = expr.split("==", 1)
            col = col.strip()
            value = value.strip().strip("'\"")

            for polymer in self._polymers:
                if str(polymer.get(col, "")) == value:
                    filtered_polymers.append(polymer)

        return PolymerDataFrame(filtered_polymers, self._columns)

    def unique(self, column: str) -> List[Any]:
        """Get unique values from a column."""
        values = self[column]
        return list(set(value for value in values if value is not None))

    def __len__(self) -> int:
        """Get number of polymers."""
        return len(self._polymers)

    def __repr__(self) -> str:
        """String representation similar to pandas DataFrame."""
        if len(self._polymers) == 0:
            return f"Empty PolymerDataFrame\nColumns: {self._columns}"

        # Show first few rows
        lines = []
        lines.append(
            "PolymerDataFrame ({} rows x {} columns)".format(
                self.shape[0], self.shape[1]
            )
        )
        lines.append("   " +
                     "  ".join(f"{col:>12}" for col in self._columns[:5]))

        for i, polymer in enumerate(self._polymers[:5]):
            values = []
            for col in self._columns[:5]:
                val = polymer.get(col, "")
                if isinstance(val, float):
                    values.append(f"{val:>12.3f}")
                else:
                    values.append(f"{str(val):>12}")
            lines.append(f"{i:2d} " + "  ".join(values))

        if len(self._polymers) > 5:
            lines.append("...")

        return "\n".join(lines)


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
