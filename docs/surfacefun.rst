Functions and vector fields
===========================

Once a surface mesh has been constructed, we may define scalar functions and
vector fields on the surface.

Scalar functions
----------------

The fundamental object which represents a scalar function on a surface is a
``surfacefun``. A ``surfacefun`` may be constructed from a function handle
representing a given function in Cartesian :math:`(x,y,z)` coordinates on a
given ``surfacemesh``:

.. code-block:: matlab

    f = surfacefun(@(x,y,z) cos(6*x).*y + exp(z), dom)

.. container:: output-text

    .. raw:: html

        <pre style="line-height: 1.4;">
        f = 

          surfacefun with properties:

            domain: [1×1 surfacemesh]
              vals: {96×1 cell}
        </pre>

Let's plot the function:

.. code-block:: matlab

    plot(f), hold on, plot(dom), colorbar

.. container:: output-image

    .. figure:: images/func.png
        :width: 350px
        :align: center

Many standard MATLAB arithmetic functions have been overloaded.

.. code-block:: matlab

    x = surfacefun(@(x,y,z) x, dom);
    g = abs(f + 10*x);
    plot(g)

.. container:: output-image

    .. figure:: images/func_arith.png
        :width: 350px
        :align: center

We can also visualize a ``surfacefun`` using a contour plot:

.. code-block:: matlab

    contour(f, linewidth=2)
    axis off

.. container:: output-image

    .. figure:: images/contour.png
        :width: 200px
        :align: center

We may numerically differentiate a function using the built-in ``diff`` or
``grad`` routines, which automatically take into account the on-surface metric.
For example:

.. code-block:: matlab

    [fx, fy, fz] = grad(f);
    subplot(131), plot(fx)
    subplot(132), plot(fy)
    subplot(133), plot(fz)

.. container:: output-image

    .. figure:: images/diff_func.png
        :width: 650px
        :align: center

Higher-order derivatives may be constructed by composing these operations. For
example, here is the surface Laplacian---or the Laplace--Beltrami
operator---applied to our function:

.. code-block:: matlab

    plot(lap(f))

.. container:: output-image

    .. figure:: images/func_lap.png
        :width: 350px
        :align: center

The definite integral of a function over the surface is given by:

.. code-block:: matlab

    integral(f)

.. container:: output-text

    .. raw:: html

        <pre style="line-height: 1.4;">
        ans =

          20.413449092485330
        </pre>

Similarly, the mean of the function is the integral of the function divided by
the surface area:

.. code-block:: matlab

    mean(f)

.. container:: output-text

    .. raw:: html

        <pre style="line-height: 1.4;">
        ans =

           1.111334042648337
        </pre>

Norms
~~~~~

The :math:`L^2` norm of a ``surfacefun`` may be computed via:

.. code-block:: matlab

    norm(f)

.. container:: output-text

    .. raw:: html

        <pre style="line-height: 1.4;">
        ans =

           5.947309239751656
        </pre>

Other norms are implemented as well. The :math:`L^\infty` norm is computed via:

.. code-block:: matlab

    norm(f, inf)

.. container:: output-text

    .. raw:: html

        <pre style="line-height: 1.4;">
        ans =

           3.229329881902320
        </pre>

Vector fields
-------------

The ``surfacefunv`` object represents a three-component vector field over a
``surfacemesh``. Each component is itself represented as a scalar
``surfacefun``.

Let's make quiver plot of the normal vectors over our surface. We'll plot 6
vectors per patch and scale their lengths by 0.2:

.. code-block:: matlab

    v = normal(dom);
    quiver(v, 0.2, 6)

.. container:: output-image

    .. figure:: images/vec_normals.png
        :width: 350px
        :align: center

The surface gradient of a ``surfacefun`` is a ``surfacefunv``:

.. code-block:: matlab

    grad(f)

.. container:: output-text

    .. raw:: html

        <pre style="line-height: 1.4;">
        ans = 

          surfacefunv with properties:

              components: {1×3 cell}
            isTransposed: 0
        </pre>

The gradient is tangent to the surface, as we can see from a quiver plot:

.. code-block:: matlab

    quiver(grad(f), 0.05, 6)

.. container:: output-image

    .. figure:: images/vec_grad.png
        :width: 350px
        :align: center

The surface divergence of the surface gradient is equal to the surface
Laplacian:

.. code-block:: matlab

    norm(div(grad(f)) - lap(f))

.. container:: output-text

    .. raw:: html

        <pre style="line-height: 1.4;">
        ans = 

              0
        </pre>

The mean curvature of a surface can be related its the normal vector field
via the surface divergence:

.. code-block:: matlab

    plot(div(v)/2)

.. container:: output-image

    .. figure:: images/vec_div.png
        :width: 400px
        :align: center

We can also take the surface curl of a ``surfacefunv``:

.. code-block:: matlab

    v = surfacefunv(@(x,y,z) cos(2*x), ...
                    @(x,y,z) sin(4*y,  ...
                    @(x,y,z) sin(3*z), dom);

.. container:: output-image

    .. figure:: images/vec_curl.png
        :width: 350px
        :align: center
