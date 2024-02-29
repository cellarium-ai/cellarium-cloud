Running Locally
===============
To run Cellarium Cloud services locally, you will need to have the following installed:

- Python 3.10
- `Local Postgres Database cluster <https://www.docker.com/blog/how-to-use-the-postgres-docker-official-image>`_
- Project Environment variables in a `src/casp/services/.env` file. :ref:`More info <Project Secrets>`.
- Project dependencies installed. ``pip install -r requirements.txt``
- `src` directory added to your ``PYTHONPATH`` environment variable. E.g. ``export PYTHONPATH=$PYTHONPATH:/path/to/cellarium-cloud/src``


Once you have the above installed, you can run the services locally using the following command:

.. code-block:: bash

    python src/casp/services/<service_name>/main.py

To check the API methods that exist and their documentation, you can visit `API docs page <http://localhost:8000/api/docs>`_

Most of the methods will reuqire you to be authenticated. To do this, you'd need to deploy Admin service:

.. code-block:: bash

    python src/casp/services/admin/server.py


Once the Admin service is running, you can visit `Admin Dashboard <http://localhost:5000>`_ to create a user and token.
