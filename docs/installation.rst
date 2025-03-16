.. _install:

Installation
============

The ``iinuclear`` package can be installed directly from the GitHub source.

From source
-----------

You can install the ``iinuclear`` package directly from the source on 
`GitHub <https://github.com/gmzsebastian/iinuclear>`_ by running:

.. code-block:: bash

    git clone https://github.com/gmzsebastian/iinuclear.git
    cd iinuclear
    python -m pip install -e .

Dependencies
------------

Some features require a TNS API key to obtain the coordinates of a transient given its IAU name.
There are two ways to provide these credentials:

1. **Using Environment Variables**

   Set the following three environment variables in your system:

   .. code-block:: bash

       export TNS_API_KEY=api_key
       export TNS_ID=tns_id
       export TNS_USERNAME=tns_username

2. **Using a Local Key File**

   Alternatively, create a file called ``tns_key.txt`` in your home directory with the following content:

   .. code-block:: bash

       api_key
       tns_id
       tns_username

To obtain an API key, register at the `TNS website <https://wis-tns.weizmann.ac.il>`_.
