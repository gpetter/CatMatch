import catmatch
from astropy.table import Table
from astropy.coordinates import SkyCoord
import pandas as pd
import pytest

@pytest.mark.parametrize(
    ("join_type"),
    (
        [
            "1+2",
            "all_from_1",
            "1_not_2"
        ]
    )
)
@pytest.mark.parametrize(
    ("match_type"),
    (
        "best",
        "best_symmetric",
        "all"
    )
)
def test_match_catalogs(match_type, join_type):
    sep = 3.
    ages = Table.read("testdata/AGES.fits")
    cdwfs = Table.read("testdata/CDWFS.fits")
    topcat_result = Table.read(f"testdata/match_{int(sep)}arcsec_{match_type}_{join_type}.fits").to_pandas()

    actual_result = catmatch.match_catalogs(ages, cdwfs, sep, match_type=match_type, join_type=join_type).to_pandas()

    if (match_type == "best" or match_type == "all") and join_type != "1_not_2":
        topcat_result.drop("GroupID", axis=1, inplace=True)

    actual_result = actual_result.sort_values(["RA", "DEC"])
    topcat_result = topcat_result.sort_values(["RA", "DEC"])

    pd.testing.assert_frame_equal(
        actual_result.reset_index(drop=True),
        topcat_result.reset_index(drop=True),
        atol=0.05,
        check_dtype=False
    )