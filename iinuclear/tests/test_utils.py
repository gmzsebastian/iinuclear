import os
import pytest
import numpy as np
from iinuclear.utils import (get_tns_coords, get_tns_credentials, get_ztf_name,
                             get_ztf_coordinates, get_coordinates)


# This marker will skip tests on GitHub Actions.
skip_on_ci = pytest.mark.skipif(
    os.environ.get('CI_TESTING') == 'true' or
    os.environ.get('GITHUB_ACTIONS') == 'true',
    reason='Test skipped on CI environment'
)


def test_get_tns_coords():
    # Check that the function returns None for an object not found
    ra, dec, ztf_name = get_tns_coords('potato')
    assert ra is None
    assert dec is None
    assert ztf_name is None


def test_get_ztf_name():
    # Test with known coordinates of 2018hyz
    ra_true = 151.711964138
    dec_true = 1.69279894089
    ztf_name = get_ztf_name(ra_true, dec_true)
    assert ztf_name == "ZTF18acpdvos"

    # Test with invalid input types
    ztf_name = get_ztf_name("potato", 1.0)
    assert ztf_name is None


def test_get_ztf_coordinates():
    # Test with known ZTF object
    ras, decs = get_ztf_coordinates("ZTF18acpdvos")
    assert ras is not None
    assert decs is not None
    assert isinstance(ras, np.ndarray)
    assert isinstance(decs, np.ndarray)
    assert len(ras) > 0
    assert len(decs) > 0

    # Test some known coordinates (approximately)
    assert np.isclose(np.median(ras), 151.712, rtol=1e-3)
    assert np.isclose(np.median(decs), 1.693, rtol=1e-3)

    # Test with invalid ZTF name
    ras, decs = get_ztf_coordinates("potato")
    assert ras is None
    assert decs is None


def test_get_coordinates():
    # Test with IAU name
    ras, decs, ztf_name, iau_name = get_coordinates("2018hyz")
    assert ras is not None
    assert decs is not None
    assert isinstance(ras, np.ndarray)
    assert isinstance(decs, np.ndarray)
    assert ztf_name == "ZTF18acpdvos"
    assert iau_name == "2018hyz"
    assert np.isclose(np.median(ras), 151.712, rtol=1e-3)
    assert np.isclose(np.median(decs), 1.693, rtol=1e-3)

    # Test with ZTF name
    ras, decs, ztf_name, iau_name = get_coordinates("ZTF18acpdvos")
    assert ras is not None
    assert decs is not None
    assert isinstance(ras, np.ndarray)
    assert isinstance(decs, np.ndarray)
    assert ztf_name == "ZTF18acpdvos"
    assert iau_name is None
    assert np.isclose(np.median(ras), 151.712, rtol=1e-3)
    assert np.isclose(np.median(decs), 1.693, rtol=1e-3)

    # Test with coordinates
    ras, decs, ztf_name, iau_name = get_coordinates(151.711964138, 1.69279894089)
    assert ras is not None
    assert decs is not None
    assert isinstance(ras, np.ndarray)
    assert isinstance(decs, np.ndarray)
    assert ztf_name == "ZTF18acpdvos"
    assert iau_name is None
    assert np.isclose(np.median(ras), 151.712, rtol=1e-3)
    assert np.isclose(np.median(decs), 1.693, rtol=1e-3)

    # Test with invalid input
    with pytest.raises(ValueError):
        get_coordinates()  # No arguments
    with pytest.raises(ValueError):
        get_coordinates(1, 2, 3)  # Too many arguments

    # Test with invalid object name
    ras, decs, ztf_name, iau_name = get_coordinates("potato")
    assert ras is None
    assert decs is None
    assert ztf_name is None
    assert iau_name == "potato"


#########
# These functions do not work with CI/CD on GitHub Actions #####
#########
@skip_on_ci
def test_get_tns_credentials():
    # Call the function to retrieve credentials.
    api_key, tns_id, username = get_tns_credentials()

    # Check that username is correct
    assert username == "BTDG_Bot"


@skip_on_ci
def test_get_tns_coords_success():
    # Test a valid object that returns complete data.
    ra, dec, ztf_name = get_tns_coords('2018hyz')
    assert isinstance(ra, float)
    assert isinstance(dec, float)

    # Check that the ZTF name is correct
    expected_ztf = "ZTF18acpdvos"
    assert isinstance(ztf_name, str)
    assert ztf_name == expected_ztf, f"Expected ztf_name to be {expected_ztf}, but got {ztf_name}"

    ra_true = 151.711964138
    dec_true = 1.69279894089
    # Check that the coordinates are within 1 arcsecond of the true values
    assert abs(ra - ra_true) * 3600 <= 1
    assert abs(dec - dec_true) * 3600 <= 1


@skip_on_ci
def test_get_tns_coords_not_found():
    # Test an object that should not be found.
    ra, dec, ztf_name = get_tns_coords('potato')
    assert ra is None
    assert dec is None
    assert ztf_name is None


@skip_on_ci
def test_get_tns_coords_no_ztf():
    # Test an object that returns valid coordinates but no ZTF name.
    ra, dec, ztf_name = get_tns_coords('2016iet')
    assert ztf_name is None
    assert isinstance(ra, float)
    assert isinstance(dec, float)
