import os
import pytest
from iinuclear.utils import (get_tns_coords, get_tns_credentials)


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


def test_check_github_actions_env():
    """Test to check if GITHUB_ACTIONS is in environment variables"""
    is_github = 'GITHUB_ACTIONS' in os.environ
    github_value = os.environ.get('GITHUB_ACTIONS', 'not found')
    print(f"\nGITHUB_ACTIONS in env vars: {is_github}")
    print(f"GITHUB_ACTIONS value: {github_value}")
    
    # This assertion will pass regardless of environment
    # It's just to make the test valid
    assert True

def test_print_all_env_vars():
    """Test that prints all environment variables"""
    print("\nAll Environment Variables:")
    print("-" * 50)
    for key, value in sorted(os.environ.items()):
        print(f"{key}: {value}")
    print("-" * 50)
    
    # This assertion will pass regardless of environment
    # It's just to make the test valid
    assert True


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
