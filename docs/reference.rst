.. _usage:

Usage
=====

There are a number of different ways to use ``iinuclear``.

Getting Coordinates
-----------------

The main function ``get_coordinates`` allows you to retrieve coordinates for transient objects using different identifiers:

Using IAU/TNS Names
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from iinuclear.utils import get_coordinates

    # Query using IAU name
    ras, decs, ztf_name, iau_name = get_coordinates("2018hyz")

Using ZTF Names
~~~~~~~~~~~~~~

.. code-block:: python

    # Query using ZTF name
    ras, decs, ztf_name, iau_name = get_coordinates("ZTF18acpdvos")

Using Coordinates
~~~~~~~~~~~~~~~

.. code-block:: python

    # Query using RA, DEC in degrees
    ras, decs, ztf_name, iau_name = get_coordinates(151.711964138, 1.69279894089)

Return Values
~~~~~~~~~~~~

The function returns a tuple containing:

- ``ras``: numpy array of Right Ascension values in degrees for all detections
- ``decs``: numpy array of Declination values in degrees for all detections
- ``ztf_name``: string containing the ZTF identifier for the object
- ``iau_name``: string containing the IAU name of the object (if available)

If the object cannot be found or an error occurs, the function returns ``(None, None, None, None)``.

Individual Functions
-----------------

You can also use the individual component functions:

.. code-block:: python

    from iinuclear.utils import get_ztf_name, get_ztf_coordinates, get_tns_coords

    # Get ZTF name from coordinates
    ztf_name = get_ztf_name(ra=151.711964138, dec=1.69279894089)

    # Get all detection coordinates for a ZTF object
    ras, decs = get_ztf_coordinates("ZTF18acpdvos")

    # Get coordinates and ZTF name from TNS
    ra, dec, ztf_name = get_tns_coords("2018hyz")
