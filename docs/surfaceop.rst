Partial differential equations
==============================

The Surfacefun package is capable of solving general variable-coefficient,
second-order, linear, elliptic partial differential equations on open and
closed surfaces. Consider such a PDE on a given surface :math:`\Gamma`,

.. math::

    \mathcal{L}_\Gamma u(\boldsymbol{x}) = f(\boldsymbol{x}), \qquad \boldsymbol{x} \in \Gamma,

where :math:`\mathcal{L}_\Gamma` is an elliptic partial differential operator of
the form

.. math::

    \mathcal{L}_\Gamma u(\boldsymbol{x}) = \sum_{i=1}^3 \sum_{j=1}^3
    a_{ij}(\boldsymbol{x}) \, \partial_i^\Gamma \partial_j^\Gamma u(\boldsymbol{x}) +
    \sum_{i=1}^3 b_i(\boldsymbol{x}) \,\partial_i^\Gamma u(\boldsymbol{x}) +
    c(\boldsymbol{x}) u(\boldsymbol{x}),

with smooth coefficients :math:`a_{ij}`, :math:`b_i`, and :math:`c`. Here, we
identify :math:`\partial_1^\Gamma`, :math:`\partial_2^\Gamma`,
:math:`\partial_3^\Gamma` with the Cartesian :math:`x`-, :math:`y`-, and
:math:`z`-components of the surface gradient, respectively.
If :math:`\Gamma` is not a closed surface, the PDE may also be
subject to boundary conditions, e.g., :math:`u(\boldsymbol{x}) = g(\boldsymbol{x})`
for :math:`\boldsymbol{x} \in \partial\Gamma` and some function :math:`g`.

Partial differential operators in Surfacefun are specified as MATLAB structs
with properties defining the coefficients on each term appearing in
:math:`\mathcal{L}_\Gamma`. The coefficients on the second-order terms are
specified by setting ``pdo.dxx``, ``pdo.dyy``, ``pdo.dzz``, ``pdo.dxy``,
``pdo.dyx``, ``pdo.dyz``, ``pdo.dzy``, ``pdo.dxz``, and ``pdo.dzx``.
Coefficients on the first-order terms can be set via ``pdo.dx``, ``pdo.dy``,
and ``pdo.dz``. The zeroth-order coefficient is specified via ``pdo.c``. Each
coefficient may be a constant, a function handle, or a ``surfacefun``. For
instance, the Laplace-Beltrami operator can be specified via

.. code-block:: matlab

    pdo = [];
    pdo.dxx = 1;
    pdo.dyy = 1;
    pdo.dzz = 1;

or more simply as

.. code-block:: matlab

    pdo = [];
    pdo.lap = 1;

This sets the coefficients on the terms :math:`\partial_x^\Gamma \partial_x^\Gamma`,
:math:`\partial_y^\Gamma \partial_y^\Gamma`, and
:math:`\partial_z^\Gamma \partial_z^\Gamma` to one and the rest to zero.

Let's solve a simple Laplace--Beltrami problem on the surface of the sphere.
Since the spherical harmonics are eigenfunctions of the Laplace--Beltrami
operator on the sphere, we can construct a test problem analytically:

.. code-block:: matlab

    % Make a sphere mesh
    p = 16;
    nref = 2;
    dom = surfacemesh.sphere(p+1, nref);

    % Construct the true solution and corresponding righthand side
    l = 3; m = 2;
    sol = spherefun.sphharm(l, m);
    sol = surfacefun(@(x,y,z) sol(x,y,z), dom);
    f = -l*(l+1)*sol;

    % Specify the partial differential operator to be Laplace-Beltrami
    pdo = [];
    pdo.lap = 1;

To solve the PDE, we create a ``surfaceop``. The ``surfaceop`` object
encapsulates the factorization and solution routines corresponding to the fast
direct solver in [1].

.. code-block:: matlab

    L = surfaceop(dom, pdo, f);

It turns out that the Laplace--Beltrami problem on a closed surface is rank
deficient (by one). However, it is uniquely solvable under the mean-zero
conditions

.. math::

    \int_\Gamma u = \int_\Gamma f = 0.

We can impose this condition by setting

.. code-block:: matlab

    L.rankdef = true;

Now we can solve the PDE:

.. code-block:: matlab

    u = L.solve();
    plot(u)

.. container:: output-image

    .. figure:: images/lb_sphere.png
        :width: 400px
        :align: center

Let's check the error:

.. code-block:: matlab

    norm(u-sol)

.. container:: output-text

    .. raw:: html

        <pre style="line-height: 1.4;">
        ans =

             6.209136661992651e-12
        </pre>

Surface PDEs on surfaces of arbitrary genus may be solved using ``surfaceop``.
For example, here is the solution to a variable-coefficient surface Helmholtz
equation on a genus-1 stellarator geometry:

.. code-block:: matlab

    % Construct the stellarator geometry
    p = 16; nu = 8; nv = 24;
    dom = surfacemesh.stellarator(p+1, nu, nv);

    % Variable-coefficient surface Helmholtz equation
    pdo = [];
    pdo.lap = 1;
    pdo.c = @(x,y,z) 300*(1-z);

    f = -1;
    L = surfaceop(dom, pdo, f);
    u = L.solve();
    plot(u), colorbar

.. container:: output-image

    .. figure:: images/helmholtz.png
        :width: 450px
        :align: center

Now let's solve a problem on an open surface. We'll create an open surface by
extracting a subset of the patches from a closed surface:

.. code-block:: matlab

    rng(0)
    p = 16;
    nref = 2;
    dom = surfacemesh.blob(p+1, nref);
    dom = surfacemesh(dom.x(1:16), dom.y(1:16), dom.z(1:16));
    plot(dom), view(-110, 30), camlight

.. container:: output-image

    .. figure:: images/open_surface.png
        :width: 250px
        :align: center

We construct a ``surfaceop`` on an open surface in the same way as on a closed
surface, except now the ``L.solve()`` method requires Dirichlet boundary data to
be passed as an argument:

.. code-block:: matlab

    % Specify the righthand side and Dirichlet boundary data
    f = -1;
    bc = 0;

    pdo = [];
    pdo.lap = 1;
    L = surfaceop(dom, pdo, f);
    u = L.solve(bc);

    plot(u), view(-110, 30), colorbar

.. container:: output-image

    .. figure:: images/open_sol.png
        :width: 375px
        :align: center

Modifying an existing ``surfaceop``
-----------------------------------

A ``surfaceop`` is a direct solver for the specified surface PDE. This means
that once a factorization of the problem is constructed, the factorization may
be reused to solve the same PDE with different righthand sides and boundary data
in a manner that is more efficient than creating a new ``surfaceop`` again and
again.

To this end, a ``surfaceop`` may be factorized before it is given any data. This
is performed implicitly when ``L.solve()`` is called, but may be performed
explicitly by calling ``L.build()``:

.. code-block:: matlab

    L = surfaceop(dom, pdo);
    L.build();

The ``surfaceop`` can now be given any righthand side and will efficiently
update its factorization accordingly:

.. code-block:: matlab

    L.rhs = @(x,y,z) sin(x.*y);
    u = L.solve(bc);
    plot(u), view(-110, 30), colorbar

.. container:: output-image

    .. figure:: images/open_sol2.png
        :width: 375px
        :align: center

The ``surfaceop`` object is agnostic to boundary data. This means that an
existing solver may be used with any Dirichlet boundary data by simply passing
it to ``L.solve()``:

.. code-block:: matlab

    bc = @(x,y,z) z;
    u = L.solve(bc);
    plot(u), view(-110, 30), colorbar

.. container:: output-image

    .. figure:: images/open_sol3.png
        :width: 375px
        :align: center