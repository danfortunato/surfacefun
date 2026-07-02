Installation
============

To install the Surfacefun MATLAB package, follow these steps.

Install with `mip <https://mip.sh>`_ (recommended)
--------------------------------------------------

The recommended way to install Surfacefun in MATLAB is via the `mip <https://mip.sh>`_
package manager:

.. code-block:: matlabsession

   >> mip install surfacefun

To use it in your current MATLAB session, run:

.. code-block:: matlabsession

   >> mip load surfacefun

This will automatically add Surfacefun to your MATLAB path.

Manual installation
-------------------

Alternatively, you can install Surfacefun manually by following these steps.

1. In a terminal, clone the source repository from GitHub:

   .. code-block:: console

      git clone --recurse-submodules https://github.com/danfortunato/surfacefun.git

2. In MATLAB, change to the directory where the repository was cloned, then
   change to the ``surfacefun`` directory:

   .. code-block:: matlabsession

      >> cd 'surfacefun'

3. Now add the source files to your MATLAB path by running:

   .. code-block:: matlabsession

      >> setup

4. (Optional) If you would like this path to persist between MATLAB sessions, you may save
   it by running:

   .. code-block:: matlabsession

      >> savepath
