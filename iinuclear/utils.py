import os
import json
import pathlib
import warnings
from collections import OrderedDict
import requests

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
