### defining a bunch of basic calcs ###

def fpd_to_ipy(fpd: float) -> float:
    """function to convert feet per day into inches
    per year. provide feet per day as float argument."""
    conversion = 12 * 365.25
    return fpd * conversion
