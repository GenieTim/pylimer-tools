import pandas as pd


def getTail(data, percentage=0.2, minN=25, maxPercentage=0.5):
    """
    Extract the last few entries of a list

    Arguments:
        - data (list|pd.DataFrame|pd.Series): the list to extract the last few entries from
        - percentage: the percentage of entries to extract
        - minN: the minimum number of entries to extract
        - maxPercentage: the maximum percentage of entries to extract

    Returns:
        A (list|pd.DataFrame|pd.Series) with at maximum maxPercentage, 
            at least minN entries (assuming the initial data is as large), 
            but ideally `percentage` many percentage of the last entries.
    """
    assert(percentage <= 1)
    assert(maxPercentage <= 1)
    tailN = int(min(max(min(minN, maxPercentage*len(data)),
                        percentage*len(data)), len(data)))
    if (isinstance(data, pd.DataFrame) or isinstance(data, pd.Series)):
        return data.tail(tailN)
    else:
        return data[-tailN:]
