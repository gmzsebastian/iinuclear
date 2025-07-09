from .utils import get_data, get_galaxy_center, check_nuclear, rice_separation, calc_separations
from .plots import plot_all
import matplotlib.pyplot as plt


def isit(*args, error=0.1, save_all=True, base_dir='.', add_sdss=True, add_ps1=True, plot=True, search_radius=120,
         coords_directory='coords'):
    """
    Main function to determine whether a transient is nuclear or not, the function takes either
    an IAU name, a ZTF name, or a set of coordinates. The function requires the transient to have
    a ZTF light curve with more than 1 detection, as well as an object within 1-arcmin in SDSS
    or PanSTARRS.

    Parameters
    -----------
    *args : str or float
        Either:
        - A single string containing an object name from IAU or ZTF (e.g., "2018hyz" or "ZTF18acpdvos")
        - Two floats representing RA and DEC in degrees
    save_all : bool, optional
        Save all data to disk (default: True)
    base_dir : str, optional
        Base directory for saving data (default: current directory)
    add_sdss : bool, optional
        Include SDSS data in the calculation (default: True)
    add_ps1 : bool, optional
        Include PS1 data in the calculation (default: True)
    plot : bool, optional
        Plot the output diagnosis plots?

    Returns
    -------
    sigma : float
        Significance of separation
    chi2_val : float
        Chi-square statistic
    p_val : float
        P-value for the hypothesis test (null hypothesis: ZTF positions are consistent with galaxy center)
    is_nuclear : bool
        True if the ZTF positions are consistent with the galaxy center, False otherwise
    mean_separation : float
        Mean separation calculated from Rice statistics
    upper_err : float
        Upper error on the mean separation
    lower_err : float
        Lower error on the mean separation
    upper_limit : float
        Upper limit on the separation at the specified confidence level
    """

    ras, decs, ztf_name, iau_name, catalog_result, used_catalog, image_data, image_header = get_data(*args, save_all=save_all, base_dir=base_dir,
                                                                                                     coords_directory=coords_directory, plot=plot)
    if len(catalog_result) > 0:
        # Format object name
        if iau_name is not None:
            object_name = iau_name
        elif ztf_name is not None:
            object_name = ztf_name

        # Get Galaxy Center
        ra_galaxy, dec_galaxy, mean_ra_offset, mean_dec_offset, error_arcsec = get_galaxy_center(catalog_result, used_catalog, object_name, error,
                                                                                                 add_sdss, add_ps1, search_radius)

        # Implement offsets
        ras_mod = ras - mean_ra_offset / 3600
        decs_mod = decs - mean_dec_offset / 3600

        if plot:
            (sigma, chi2_val, p_val, is_nuclear,
             mean_separation, upper_err,
             lower_err, upper_limit) = plot_all(image_data, image_header, ras_mod, decs_mod, ra_galaxy, dec_galaxy, error_arcsec,
                                                object_name=object_name, used_catalog=used_catalog, mean_ra_offset=mean_dec_offset,
                                                mean_dec_offset=mean_dec_offset)
            plt.savefig(f'{object_name}_iinuclear.pdf', bbox_inches='tight')
            plt.clf()
            plt.close('all')
        else:
            # Calculate metrics without plotting
            sigma, chi2_val, p_val, is_nuclear = check_nuclear(ras_mod, decs_mod, ra_galaxy, dec_galaxy, error_arcsec)
            separations = calc_separations(ras_mod, decs_mod, ra_galaxy, dec_galaxy)
            mean_separation, lower_err, upper_err, _, upper_limit, _ = rice_separation(separations, error_arcsec)
    else:
        print("No sources found in input catalog.")
        return None, None, None, None, None, None, None, None

    return sigma, chi2_val, p_val, is_nuclear, mean_separation, upper_err, lower_err, upper_limit
