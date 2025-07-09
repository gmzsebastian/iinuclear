from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, LinearStretch, LogStretch, SqrtStretch
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import warnings
from astropy.wcs import FITSFixedWarning
from .utils import rice_separation, check_nuclear, calc_separations
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'font.family': 'serif'})


def plot_image(image_data, image_header, ras, decs, radius_arcsec=2,
               scale='sqrt', figsize=(8, 8), ax=None, object_name=None,
               ra_galaxy=None, dec_galaxy=None, error_arcsec=None, used_catalog=None,
               mean_ra_offset=None, mean_dec_offset=None):
    """
    Plot PS1 image with ZTF positions overlaid in WCS coordinates.

    Parameters
    -----------
    image_data : numpy.ndarray
        Image data from the FITS file
    image_header : astropy.io.fits.Header
        FITS header with WCS information
    ras : numpy.ndarray
        Array of Right Ascension values in degrees
    decs : numpy.ndarray
        Array of Declination values in degrees
    radius_arcsec : float, optional
        Radius of the region to plot in arcseconds (default: 3)
    scale : str, optional
        Scaling for the image ('log', 'linear', 'sqrt', default: 'sqrt')
    figsize : tuple, optional
        Figure size in inches (default: (8,8))
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on (default: None)
    object_name : str, optional
        Name of the object to include in the title (default: None)
    ra_galaxy : float, optional
        Galaxy center RA in degrees
    dec_galaxy : float, optional
        Galaxy center Dec in degrees
    error_arcsec : float, optional
        Error circle radius in arcseconds
    used_catalog : str, optional
        Name of the catalog used for galaxy center determination (default: None)
    mean_ra_offset : float, optional
        Mean RA offset applied to ZTF positions in arcsec (default: None)
    mean_dec_offset : float, optional
        Mean Dec offset applied to ZTF positions in arcsec (default: None)

    Returns
    --------
    ax : matplotlib.axes.Axes
        The axes object containing the plot
    """

    # Create a WCS object from the header
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        wcs = WCS(image_header)

    # Create figure and axis with WCS projection if needed
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection=wcs)
        single = True
    else:
        # Use existing axis but ensure it has WCS projection
        if not hasattr(ax, 'wcs'):
            oldax = ax
            fig = ax.figure
            # Get the subplot position
            position = oldax.get_position()
            # Remove old axis
            oldax.remove()
            # Add new axis with WCS projection in same position
            ax = fig.add_axes(position, projection=wcs)
        single = False

    # Set up normalization based on scale parameter
    if scale == 'log':
        stretch = LogStretch()
    elif scale == 'sqrt':
        stretch = SqrtStretch()
    else:
        stretch = LinearStretch()

    # Convert the median RA/DEC to pixel coordinates using the WCS
    ra_center = np.nanmedian(ras)
    dec_center = np.nanmedian(decs)
    center_pix = wcs.wcs_world2pix(ra_center, dec_center, 0)
    center_x, center_y = center_pix

    # Calculate the radius in pixels
    plate_scale = 0.25
    pixel_radius = radius_arcsec / plate_scale

    # Limits
    min_x, max_x = center_x - pixel_radius, center_x + pixel_radius
    min_y, max_y = center_y - pixel_radius, center_y + pixel_radius

    # Calculate image statistics for scaling
    image_crop = image_data[
        max(int(np.floor(min_y)), 0): min(int(np.ceil(max_y)), image_data.shape[0]),
        max(int(np.floor(min_x)), 0): min(int(np.ceil(max_x)), image_data.shape[1])
    ]
    vmin = np.percentile(image_crop, 0.001)
    vmax = np.percentile(image_crop, 99.99)
    norm = ImageNormalize(image_crop, stretch=stretch, vmin=vmin, vmax=vmax)

    # Display the image data
    ax.imshow(image_data, origin='lower', cmap='gray', norm=norm)

    # Plot ZTF positions
    ax.scatter(ras, decs, transform=ax.get_transform('world'),
               s=50, edgecolor='b', facecolor='none', marker='o',
               label='ZTF detections', alpha=0.8)

    # Plot median position
    ax.scatter(ra_center, dec_center, transform=ax.get_transform('world'),
               s=100, color='cyan', marker='+', label='ZTF Center')

    # If galaxy center and error are provided, plot a circle around the galaxy center.
    if (ra_galaxy is not None) and (dec_galaxy is not None) and (error_arcsec is not None):
        # Convert galaxy center (world coordinates) to pixel coordinates.
        galaxy_pix = wcs.wcs_world2pix(ra_galaxy, dec_galaxy, 0)
        galaxy_x, galaxy_y = galaxy_pix
        # Convert error in arcsec to pixels
        error_pixel_radius = error_arcsec / plate_scale
        galaxy_circle = Circle((galaxy_x, galaxy_y), error_pixel_radius,
                               edgecolor='r', facecolor='none', lw=2,
                               transform=ax.get_transform('pixel'),
                               label='Galaxy Error', linestyle='--',
                               alpha=0.8)
        ax.add_patch(galaxy_circle)

        # Plot median position
        if mean_ra_offset is not None and mean_dec_offset is not None:
            if mean_dec_offset > 0 and mean_ra_offset > 0:
                offset_label = f'Galaxy Center\nOffset: RA {mean_ra_offset:.2f}", DEC {mean_dec_offset:.2f}"'
            else:
                offset_label = 'Galaxy Center'
        else:
            offset_label = 'Galaxy Center'
        ax.scatter(galaxy_x, galaxy_y, s=100, color='orange', marker='*', label=offset_label)

    # Label axes
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')

    # Invert RA axis
    ax.invert_xaxis()

    # Add legend
    ax.legend(loc='upper left')

    # Set the plot limits in pixel coordinates (x corresponds to RA axis, y to DEC)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)

    # Set the title equal to the object name if it exists
    if single:
        if object_name is not None:
            ax.set_title(object_name)
        else:
            ax.set_title(f'RA: {ra_center:.6f}, DEC: {dec_center:.6f}')
    else:
        ax.set_title(f'{used_catalog} Image')

    return ax


def plot_histogram(separations, error_arcsec, n_bins=15, ax=None,
                   confidence_level=0.95, separation_threshold=3.0,
                   object_name=None):
    """
    Plot histogram of separations between ZTF detections and galaxy center.
    Accounts for Rice distribution since separation is a positive quantity.

    Parameters
    -----------
    separations : numpy.ndarray
        Array of separations in arcseconds
    error_arcsec : float, optional
        Galaxy center uncertainty in arcseconds, will be plotted as vertical line
    n_bins : int, optional
        Number of bins for histogram (default: 15)
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on (default: None)
    confidence_level : float, optional
        Confidence level for upper limit (default: 0.95 for 95% confidence)
    separation_threshold : float, optional
        SNR threshold for considering a detection significant (default: 3.0)
    object_name : str, optional
        Name of the object to include in the title (default: None)

    Returns
    --------
    mean_separation : float
        Mean separation calculated from Rice statistics
    upper_err : float
        Upper error on the mean separation
    lower_err : float
        Lower error on the mean separation
    upper_limit : float
        Upper limit on the separation at the specified confidence level
    """
    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
        if object_name is not None:
            ax.set_title(object_name)
        else:
            ax.set_title('Separations')
    else:
        ax.set_title('Separations')

    # Calculate Rice statistics
    mean_separation, lower_err, upper_err, snr, upper_limit, sigma = rice_separation(separations, error_arcsec,
                                                                                     confidence_level, separation_threshold)

    # Create histogram
    ax.hist(separations, bins=n_bins, density=True, alpha=0.6, color='gray')

    # Add kernel density estimate
    kde = stats.gaussian_kde(separations)
    x_range = np.linspace(0, max(max(separations), 1.0), 100)
    kde_values = kde(x_range)

    ax.plot(x_range, kde_values, 'k-', lw=2, label=f'KDE: N = {len(separations)}')

    # Plot best fit Rice distribution
    rice_dist = stats.rice.pdf(x_range, mean_separation/sigma, scale=sigma)
    ax.plot(x_range, rice_dist, 'orange', lw=2, label='Rice Fit')

    # Plot Rayleigh statistics
    ax.axvline(mean_separation, color='red', linestyle='-', linewidth=1,
               label=(
                   rf'Mean = ${mean_separation:.2f}^{{+{upper_err:.2f}}}_'
                   rf'{{-{np.abs(lower_err):.2f}}}$' + r'$^{\prime\prime}$'
               ))
    ax.axvline(mean_separation - lower_err, color='red', linestyle='--', linewidth=1)
    ax.axvline(mean_separation + upper_err, color='red', linestyle='--', linewidth=1)

    # Plot galaxy center uncertainty if provided
    ax.axvline(error_arcsec, color='blue', linestyle='-', linewidth=1,
               label=rf'Galaxy $\sigma$ = {error_arcsec:.2f}' + r'$^{\prime\prime}$')

    # Plot upper limit
    ax.axvline(upper_limit, color='green', linestyle='-', linewidth=1,
               label=rf'{int(confidence_level*100)}% Limit = {upper_limit:.2f}'
               + r'$^{\prime\prime}$')

    # Add labels and legend
    ax.set_xlabel('Separation (arcsec)')
    ax.set_ylabel('Density')
    ax.set_xlim(0, 1.0)
    ax.legend(loc='upper right')

    return mean_separation, upper_err, lower_err, upper_limit


def plot_detections(ras, decs, ra_galaxy=None, dec_galaxy=None, error_arcsec=None,
                    radius_arcsec=None, figsize=(8, 8), ax=None, object_name=None):
    """
    Create a scatter plot of ZTF detections with density contours.

    Parameters
    -----------
    ras : numpy.ndarray
        Array of Right Ascension values in degrees
    decs : numpy.ndarray
        Array of Declination values in degrees
    ra_galaxy : float, optional
        Galaxy center RA in degrees
    dec_galaxy : float, optional
        Galaxy center Dec in degrees
    error_arcsec : float, optional
        Error circle radius in arcseconds
    radius_arcsec : float, optional
        Radius of the region to plot in arcseconds (default: None, will auto-scale)
    figsize : tuple, optional
        Figure size in inches (default: (8,8))
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on (default: None)
    object_name : str, optional
        Name of the object to include in the title (default: None)

    Returns
    --------
    ax : matplotlib.axes.Axes
        The axes object containing the plot
    """

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        if object_name is not None:
            ax.set_title(object_name)
        else:
            ax.set_title('Detections')
    else:
        ax.set_title('Detections')

    # Set up formatters to avoid scientific notation
    ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))

    # Center RA and DEC on plot
    delta_ra, delta_dec = calc_separations(ras, decs, ra_galaxy, dec_galaxy, separate=True)

    # Calculate the center and spread of detections
    deltara_center = np.nanmedian(delta_ra)
    deltadec_center = np.nanmedian(delta_dec)

    # Create a fine grid for density plotting
    ra_range = np.max(delta_ra) - np.min(delta_ra)
    dec_range = np.max(delta_dec) - np.min(delta_dec)
    margin = max(ra_range, dec_range) * 0.05  # 5% margin

    x_grid = np.linspace(np.min(delta_ra) - margin, np.max(delta_ra) + margin, 100)
    y_grid = np.linspace(np.min(delta_dec) - margin, np.max(delta_dec) + margin, 100)
    xx, yy = np.meshgrid(x_grid, y_grid)

    # Create density estimate
    positions = np.vstack([delta_ra, delta_dec])
    if len(delta_ra) < 3:
        print('Not enough detections to create a KDE')
    else:
        kernel = stats.gaussian_kde(positions)

        # Evaluate kernel on grid
        positions_grid = np.vstack([xx.ravel(), yy.ravel()])
        z = np.reshape(kernel(positions_grid).T, xx.shape)

        # Plot density contours
        ax.contourf(xx, yy, z, levels=20, cmap='Blues', alpha=0.4)

    # Plot ZTF positions
    ax.scatter(delta_ra, delta_dec, s=50, edgecolor='Blue', facecolor='none',
               marker='o', label='ZTF detections')

    # Plot median position of detections
    ax.scatter(deltara_center, deltadec_center, s=100, color='cyan',
               marker='+', label='ZTF center')

    # Calculate Covariance of ZTF detections
    Y = np.column_stack((delta_ra, delta_dec))
    star_mean = np.mean(Y, axis=0)
    cov_Y = np.cov(Y, rowvar=False, ddof=1)

    # Compute the ellipse parameters based on the covariance matrix
    vals, vecs = np.linalg.eigh(cov_Y)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    width, height = 2 * np.sqrt(vals)  # 2 standard deviations

    ellipse = Ellipse(xy=star_mean, width=width, height=height, angle=angle,
                      edgecolor='g', facecolor='none', linewidth=2, label='Covariance')
    ax.add_patch(ellipse)

    # If galaxy center and error are provided, plot them
    if (ra_galaxy is not None) and (dec_galaxy is not None):
        delta_galaxy_ra, delta_galaxy_dec = 0, 0

        if error_arcsec is not None:
            # Convert error from arcsec to degrees
            galaxy_circle = Circle((delta_galaxy_ra, delta_galaxy_dec), error_arcsec,
                                   edgecolor='r', facecolor='none',
                                   lw=2, linestyle='--',
                                   label='Galaxy Error')
            ax.add_patch(galaxy_circle)

        ax.scatter(delta_galaxy_ra, delta_galaxy_dec, s=150, color='orange',
                   marker='*', label='Galaxy center')

    # Set plot limits
    if radius_arcsec is not None:
        radius_deg = radius_arcsec / 3600
    else:
        radius_deg = max(ra_range, dec_range) / 2
    ax.set_xlim(radius_deg, -radius_deg)
    ax.set_ylim(-radius_deg, radius_deg)
    ax.set_aspect('equal', adjustable='box')

    # Format tick labels to show more decimal places
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))

    # Labels and title
    ax.set_xlabel(r'$\Delta$RA [arcsec]')
    ax.set_ylabel(r'$\Delta$Dec [arcsec]')

    # Add legend
    ax.legend(loc='best')

    return ax


def plot_pvalue_curve(ras, decs, ra_galaxy, dec_galaxy, error_arcsec,
                      figsize=(8, 6), ax=None, object_name=None):
    """
    Plot p-value as a function of galaxy position uncertainty.

    Parameters
    ----------
    ras : numpy.ndarray
        Array of Right Ascension values in degrees
    decs : numpy.ndarray
        Array of Declination values in degrees
    ra_galaxy : float
        Galaxy center RA in degrees
    dec_galaxy : float
        Galaxy center Dec in degrees
    error_arcsec : float
        Measured 1-sigma uncertainty in the galaxy center position in arcseconds
    figsize : tuple, optional
        Figure size in inches (default: (8,6))
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on (default: None)
    object_name : str, optional
        Name of the object to include in the title (default: None)

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes object containing the plot
    """

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        if object_name is not None:
            ax.set_title(object_name)
        else:
            ax.set_title('Statistical Significance')
    else:
        ax.set_title('Statistical Significance')

    # Create array of error values to test
    # Go from 0.1 to 5 times the measured error
    error_range = np.linspace(0.1 * error_arcsec, 1, 100)
    p_values = []

    # Calculate p-value for each error
    for err in error_range:
        sigma, chi2_val, p_value, is_nuclear = check_nuclear(ras, decs, ra_galaxy, dec_galaxy, err)
        p_values.append(p_value)

    # Calculate actual values for legend
    sigma, chi2_val, p_val, is_nuclear = check_nuclear(ras, decs, ra_galaxy, dec_galaxy, error_arcsec)

    # Plot p-value curve
    ax.plot(error_range, p_values, 'k-', lw=2,
            label=rf'$\sigma = {sigma:.1f}$, $p = {p_val:.3f}$')

    # Add horizontal line at p = 0.05
    ax.axhline(0.05, color='red', linestyle='--', alpha=0.5)

    # Add vertical line at measured error
    ax.axvline(error_arcsec, color='b', linestyle='--', alpha=0.5,
               label=rf'Galaxy $\sigma$ = {error_arcsec:.2f}"')
    ax.axhline(p_val, color='b', linestyle='--', alpha=0.5)

    ax.scatter(error_arcsec, p_val, color='b', marker='*', zorder=999,
               s=150)

    # Add shaded regions
    ax.fill_between(error_range, 0.05, 1, color='lightgreen', alpha=0.3,
                    label='Nuclear')
    ax.fill_between(error_range, 0, 0.05, color='salmon', alpha=0.3,
                    label='Not Nuclear')

    # Set axis labels and scales
    ax.set_xlabel('Galaxy Center Uncertainty (arcsec)')
    ax.set_ylabel('p-value')
    ax.set_yscale('log')

    # Set reasonable limits
    ax.set_xlim(0, 1.0)
    ax.set_ylim(1e-4, 1)

    # Add legend
    ax.legend(loc='lower right')

    return ax


def plot_all(image_data, image_header, ras, decs, ra_galaxy, dec_galaxy,
             error_arcsec, radius_arcsec=2, object_name=None, used_catalog=None,
             mean_ra_offset=None, mean_dec_offset=None,
             figsize=(13, 13)):
    """
    Create a 2x2 figure showing all analysis plots for a transient.

    Parameters
    ----------
    image_data : numpy.ndarray
        Image data from the FITS file
    image_header : astropy.io.fits.Header
        FITS header with WCS information
    ras : numpy.ndarray
        Array of Right Ascension values in degrees
    decs : numpy.ndarray
        Array of Declination values in degrees
    ra_galaxy : float
        Galaxy center RA in degrees
    dec_galaxy : float
        Galaxy center Dec in degrees
    error_arcsec : float
        Galaxy center uncertainty in arcseconds
    radius_arcsec : float, optional
        Radius of the region to plot in arcseconds (default: 2)
    object_name : str, optional
        Name of the object to include in the title
    used_catalog : str, optional
        Name of the catalog used for galaxy center determination (default: None)
    mean_ra_offset : float, optional
        Mean RA offset applied to ZTF positions (default: None)
    mean_dec_offset : float, optional
        Mean Dec offset applied to ZTF positions (default: None)
    figsize : tuple, optional
        Figure size in inches (default: (12,12))

    Returns
    -------
    sigma : float
        Significance of separation
    chi2_val : float
        Chi-square statistic
    p_value : float
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

    # Create figure and gridspec
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    # Plot 1: PS1 image with detections
    ax1 = fig.add_subplot(gs[0, 0])
    plot_image(image_data, image_header, ras, decs,
               radius_arcsec=radius_arcsec, ax=ax1,
               object_name=object_name,
               ra_galaxy=ra_galaxy, dec_galaxy=dec_galaxy,
               error_arcsec=error_arcsec, used_catalog=used_catalog,
               mean_ra_offset=mean_ra_offset, mean_dec_offset=mean_dec_offset)

    # Plot 2: ZTF detections with density
    ax2 = fig.add_subplot(gs[1, 0])
    plot_detections(ras, decs,
                    ra_galaxy=ra_galaxy, dec_galaxy=dec_galaxy,
                    error_arcsec=error_arcsec,
                    object_name=object_name, ax=ax2)

    # Plot 3: Separation histogram
    ax3 = fig.add_subplot(gs[0, 1])
    separations = calc_separations(ras, decs, ra_galaxy, dec_galaxy)
    mean_separation, upper_err, lower_err, upper_limit = plot_histogram(separations, error_arcsec, object_name=object_name,
                                                                        ax=ax3)

    # Plot 4: P-value curve
    ax4 = fig.add_subplot(gs[1, 1])
    plot_pvalue_curve(ras, decs, ra_galaxy, dec_galaxy,
                      error_arcsec, ax=ax4)

    # Calculate actual values for legend
    sigma, chi2_val, p_val, is_nuclear = check_nuclear(ras, decs, ra_galaxy, dec_galaxy, error_arcsec)

    # Add overall title if object name is provided
    if object_name:
        ra_center = np.nanmedian(ras)
        dec_center = np.nanmedian(decs)
        fig.suptitle(f'{object_name} (RA: {ra_center:.6f}, DEC: {dec_center:.6f}) = {is_nuclear}', fontsize=14, y=0.93)

    return sigma, chi2_val, p_val, is_nuclear, mean_separation, upper_err, lower_err, upper_limit
