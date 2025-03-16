from astropy import table
import numpy as np
import os
import json
import pathlib
import warnings
from collections import OrderedDict
import requests
from alerce.core import Alerce

# Suppress insecure request warnings (for development purposes only)
warnings.filterwarnings("ignore", message="Unverified HTTPS request")


def get_tns_credentials():
    """
    Retrieve TNS credentials from environment variables if available,
    otherwise fall back to reading the tns_key.txt file.
    """
    api_key = os.environ.get("TNS_API_KEY")
    tns_id = os.environ.get("TNS_ID")
    username = os.environ.get("TNS_USERNAME")

    if api_key and tns_id and username:
        return api_key, tns_id, username

    # Fall back to the key file in the user's home directory.
    key_path = pathlib.Path.home() / 'tns_key.txt'
    try:
        with open(key_path, 'r') as key_file:
            lines = [line.strip() for line in key_file if line.strip()]
        if len(lines) < 3:
            raise ValueError("TNS key file is incomplete. It must contain API key, TNS ID, and username.")
        return lines[0], lines[1], lines[2]
    except Exception as e:
        raise Exception("Error retrieving TNS credentials: " + str(e))


def get_tns_coords(tns_name):
    """
    Retrieve the Right Ascension (RA) and Declination (DEC) in degrees for a transient
    from the Transient Name Server (TNS) based on its IAU name, along with the ZTF name if available.

    This function requires a TNS API key file located in the user's home directory named 'tns_key.txt'.
    The file should contain three lines:
      1. API key
      2. TNS ID
      3. Username

    Parameters
    -----------
    tns_name : str
        The name of the transient (e.g. "2018hyz"). If the name starts with "AT" or "AT_",
        those prefixes will be removed.

    Returns
    --------
    tuple
        A tuple (ra_deg, dec_deg, ztf_name) where ra_deg and dec_deg are floats representing the coordinates,
        and ztf_name is a string starting with "ZTF" if found in the internal_names field, or None if not.
        Returns (None, None, None) if the transient is not found or if an error occurs.
    """
    # Normalize tns_name: remove leading "AT_" or "AT" if present.
    if tns_name.startswith("AT_"):
        tns_name = tns_name[3:]
    elif tns_name.startswith("AT"):
        tns_name = tns_name[2:]

    # Locate and read the TNS key file
    try:
        api_key, tns_id, username = get_tns_credentials()
    except Exception as e:
        print("Error retrieving TNS credentials:", e)
        return None, None, None

    # Build the URL and the query payload
    base_url = "https://www.wis-tns.org/api/get"
    object_endpoint = f"{base_url}/object"
    query_data = OrderedDict([
        ("objname", tns_name),
        ("photometry", "0"),
        ("tns_id", tns_id),
        ("type", "user"),
        ("name", username)
    ])
    payload = {
        'api_key': (None, api_key),
        'data': (None, json.dumps(query_data))
    }
    headers = {
        'User-Agent': f'tns_marker{{"tns_id":{tns_id}, "type":"bot", "name":"{username}"}}'
    }

    try:
        print(f"Querying TNS for coordinates for object '{tns_name}'...")
        response = requests.post(object_endpoint, files=payload, headers=headers)
        response.raise_for_status()
        response_json = response.json()

        # Check if the response contains valid data
        if 'data' not in response_json:
            error_message = response_json.get('id_message', 'Unknown error from TNS')
            print("TNS error:", error_message)
            return None, None, None

        data = response_json['data']
        ra = data.get('radeg')
        dec = data.get('decdeg')
        if ra is None or dec is None:
            print(f"Coordinates not found in TNS response for object '{tns_name}'.")
            return None, None, None

        # Extract ztf_name from internal_names field if available.
        internal_names = data.get('internal_names', '')
        ztf_names = [name.strip() for name in internal_names.split(',') if name.strip().startswith("ZTF")]
        ztf_name = None
        if ztf_names:
            if len(ztf_names) > 1:
                print("Warning: Multiple ZTF names found. Using the first one:", ztf_names[0])
            ztf_name = ztf_names[0]

        print(f"Found coordinates for {tns_name} at RA={ra}, DEC={dec} with ZTF name={ztf_name}")
        return ra, dec, ztf_name

    except requests.RequestException as req_err:
        print("HTTP Request error while querying TNS:", req_err)
    except Exception as e:
        print("An unexpected error occurred while querying TNS:", e)

    return None, None, None


def get_ztf_name(ra, dec, acceptance_radius=3):
    """
    Query the Alerce database to find the ZTF name of an object at given coordinates.

    Parameters
    -----------
    ra : float
        Right Ascension in degrees
    dec : float
        Declination in degrees
    acceptance_radius : float, optional
        Search radius in arcseconds (default: 3)

    Returns
    --------
    ztf_name : str or None
        The ZTF object name if found, None if no object is found at the given coordinates.
    """
    try:
        # Initialize Alerce client
        client = Alerce()

        # Query for objects at the given coordinates
        objects = client.query_objects(ra=ra, dec=dec, radius=acceptance_radius)

        # Return the ZTF name if an object was found, None otherwise
        if len(objects) > 0:
            ztf_name = objects['oid'][0]
        else:
            ztf_name = None

        return ztf_name

    except Exception as e:
        print(f"Error querying Alerce: {str(e)}")
        return None


def get_ztf_coordinates(ztf_name):
    """
    Get the list of all RAs and DECs of all detections of a ZTF object
    given its name by querying the Alerce database.

    Parameters
    -----------
    ztf_name : str
        The ZTF object name to query

    Returns
    --------
    tuple or (None, None)
        A tuple (ras, decs) where ras and decs are numpy arrays containing the
        RAs and DECs in degrees for all detections of the object. Returns (None, None)
        if no detections are found or if an error occurs.
    """

    try:
        # Initialize Alerce client
        client = Alerce()

        # Query detections for the object
        detections = client.query_detections(ztf_name, format="pandas")

        # If we have detections, convert to astropy table and get coordinates
        if len(detections) > 0:
            det_table = table.Table.from_pandas(detections)['mjd', 'magpsf', 'sigmapsf', 'ra', 'dec', 'fid']

            # Get arrays of all RAs and DECs
            ras = np.array(det_table['ra'])
            decs = np.array(det_table['dec'])

            return ras, decs
        else:
            print(f"No detections found for {ztf_name}")
            return None, None

    except Exception as e:
        print(f"Error querying Alerce: {str(e)}")
        return None, None


def get_coordinates(*args):
    """
    Get the coordinates (RA, DEC) for an transients using various identification methods.
    Can accept either an object name (IAU or ZTF) or RA/DEC coordinates.

    Parameters
    -----------
    *args : str or float
        Either:
        - A single string containing an object name from IAU or ZTF (e.g., "2018hyz" or "ZTF18acpdvos")
        - Two floats representing RA and DEC in degrees

    Returns
    --------
    tuple
        A tuple containing:
        - ras: numpy array of Right Ascension values in degrees
        - decs: numpy array of Declination values in degrees
        - ztf_name: str, the ZTF identifier for the object
        - iau_name: str, the IAU name of the object
        Returns (None, None, None, None) if the object cannot be found or if an error occurs
    """
    ztf_name = None
    iau_name = None

    # Check input arguments
    if len(args) == 1:
        # Case 1: Single argument - object name
        object_name = args[0]
        if object_name.startswith("ZTF"):
            ztf_name = object_name
            ras, decs = get_ztf_coordinates(ztf_name)
            return ras, decs, ztf_name, iau_name
        else:
            iau_name = object_name
            ra, dec, ztf_name = get_tns_coords(iau_name)
            if ra is None:
                return None, None, ztf_name, iau_name

            # If TNS didn't return a ZTF name, try to find it
            if ztf_name is None:
                ztf_name = get_ztf_name(ra, dec)

            if ztf_name is not None:
                ras, decs = get_ztf_coordinates(ztf_name)
                return ras, decs, ztf_name, iau_name
            return None, None, ztf_name, iau_name

    # Case 2: Two arguments - RA, DEC coordinates
    elif len(args) == 2:
        ra_deg, dec_deg = args
        ztf_name = get_ztf_name(ra_deg, dec_deg)
        if ztf_name is not None:
            ras, decs = get_ztf_coordinates(ztf_name)
            return ras, decs, ztf_name, iau_name
        return None, None, ztf_name, iau_name

    else:
        raise ValueError("Must provide either object name or both RA and DEC in degrees.")
