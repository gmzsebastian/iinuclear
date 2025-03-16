from iinuclear.utils import (get_tns_coords, get_tns_credentials)


def test_get_tns_coords():
    # Check that the function returns None for an object not found
    ra, dec, ztf_name = get_tns_coords('potato')
    assert ra is None
    assert dec is None
    assert ztf_name is None


# This function is deprecated because it does not work with CI/CD on GitHub Actions
if False:
    def test_get_tns_coords():
        # Call the function with a test object name
        ra, dec, ztf_name = get_tns_coords('2018hyz')

        # Check that the returned coordinates are correct
        assert isinstance(ra, float)
        assert isinstance(dec, float)

        # Check that the returned ztf_name is a string and matches the expected value.
        expected_ztf = "ZTF18acpdvos"
        assert isinstance(ztf_name, str)
        assert ztf_name == expected_ztf, f"Expected ztf_name to be {expected_ztf}, but got {ztf_name}"

        # True coordinates
        ra_true = 151.711964138
        dec_true = 1.69279894089

        # Check that the coordinates are within 1 arcsec
        assert abs(ra - ra_true) * 3600 <= 1
        assert abs(dec - dec_true) * 3600 <= 1

        # Check that the function returns None for an object not found
        ra, dec, ztf_name = get_tns_coords('potato')
        assert ra is None
        assert dec is None
        assert ztf_name is None

        # Check that querying 2016iet returns ztf_name None but real ra and dec
        ra, dec, ztf_name = get_tns_coords('2016iet')
        assert ztf_name is None
        assert isinstance(ra, float)
        assert isinstance(dec, float)
