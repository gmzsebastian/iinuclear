from .utils import get_data, get_galaxy_center
from .plots import plot_all
import matplotlib.pyplot as plt


def isit(*args, error=0.1):

    ras, decs, ztf_name, iau_name, catalog_result, image_data, image_header = get_data(*args)
    ra_galaxy, dec_galaxy, error_arcsec = get_galaxy_center(catalog_result, error)

    # Format object name
    if iau_name is not None:
        object_name = iau_name
    elif ztf_name is not None:
        object_name = ztf_name

    plot_all(image_data, image_header, ras, decs, ra_galaxy, dec_galaxy, error_arcsec,
             object_name=object_name)
    plt.savefig(f'{object_name}_iinuclear.pdf', bbox_inches='tight')
    plt.clf()
    plt.close('all')
