.. _installation:

Installation
============

There are two ways to install ``iinuclear``: using pip or from source.

Using pip
---------

To install the latest release version::

    pip install iinuclear

From source
-----------

To install the latest development version directly from GitHub::

    git clone https://github.com/gmzsebastian/iinuclear.git
    cd iinuclear
    pip install -e .

TNS API Setup
-------------

Some features require access to the Transient Name Server (TNS). You'll need API credentials to query TNS data when using IAU names.

Getting TNS API Credentials
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Register for a TNS account at https://www.wis-tns.org/
2. Once logged in, go to "User Area" → "API tokens"
3. Create a new API key by registering a "BOT" account
4. Note down these three items:
   * Your API key
   * Your TNS ID (Bot ID)
   * Your TNS username (Bot name)

Setting Up Credentials
~~~~~~~~~~~~~~~~~~~~~~

There are two ways to configure your TNS credentials:

1. **Environment Variables (Recommended)**

   Add these lines to your ``~/.bashrc``, ``~/.zshrc``, or equivalent::

       export TNS_API_KEY="your_api_key_here"
       export TNS_ID="your_tns_id_here"
       export TNS_USERNAME="your_bot_name_here"

   Then source your shell configuration file::

       source ~/.bashrc  # or ~/.zshrc

2. **Configuration File**

   Create a file named ``tns_key.txt`` in your home directory::

       touch ~/tns_key.txt

   Add your credentials to this file in the following format::

       your_api_key_here
       your_tns_id_here
       your_bot_name_here

   Set appropriate file permissions::

       chmod 600 ~/tns_key.txt
