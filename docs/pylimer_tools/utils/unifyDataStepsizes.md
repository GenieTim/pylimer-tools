Module pylimer_tools.utils.unifyDataStepsizes
=============================================

Functions
---------

    
`unifyDataStepsizes(data: pandas.core.frame.DataFrame, key: str, stepSize: int = None, maxExpectedStepSize: int = 100)`
:   Get a DataFrame where all data points have the same step between the values in column given by `key`
    NOTE: this function is rather unstable, as it has a few dirty assumptions, such as:
    - steps are modulo stepsize. Breaks e.g. with steps start with 1 and go up by stepSize.
    - ideal step-size is max step difference. Breaks e.g. if there is one big gap
    
    Args: 
        data: the DataFrame to unify the step-size for
        key: the column name indicating the column containing the step-nr.
        maxExpectedStepSize: use to get a warning if the computed step-size is larger