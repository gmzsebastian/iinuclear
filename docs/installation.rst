.. _install:

Installation
============

The ``iinuclear`` package can be installed directly from the Github source.

From source
-----------

You can also install the ``iinuclear`` package directly from the source on `GitHub <https://github.com/gmzsebastian/iinuclear>`_.

.. code-block:: bash

    git clone https://github.com/gmzsebastian/iinuclear.git
    cd iinuclear
    python -m pip install -e .

Dependencies
------------

Some features require having a TNS API key, namely to obtain the coordinates of a transient given its IAU name. There are two ways to use this feature, one is
to include these three environment variables in your system:

.. code-block:: bash
    export TNS_API_KEY=api_key
    export TNS_ID=tns_id
    export TNS_USERNAME=tns_username

Alternatively, you can create a file called ``tns_key.txt`` and place it in your home directly, with the same three variables as above:

.. code-block:: bash
api_key
tns_id
tns_username

To obtain an API key, you can register at the `TNS website <https://wis-tns.weizmann.ac.il>`_.