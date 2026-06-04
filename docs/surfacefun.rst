Functions and vector fields
===========================

Once a surface mesh has been constructed, we may define scalar functions and
vector fields on the surface.

.. note::

   The examples on this page run live in your browser through `numbl
   <https://numbl.org>`_. Click **▶ Edit & run** beneath any figure to open an
   editable copy of the complete, self-contained script, tweak it, and re-run
   it. Nothing is computed until you click; the first run downloads surfacefun
   (and Chebfun), and later runs reuse the cached copy.

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

.. raw:: html

    <numbl-embed lazy label="▶ Edit &amp; run this example" preparing-label="Installing surfacefun…">
    <iframe width="100%" height="560" frameborder="0"></iframe>
    <script type="text/plain" class="numbl-preamble">
    mip load --install flatironinstitute/flatironinstitute/surfacefun
    </script>
    <script type="text/plain" class="numbl-script">
    rng(0);
    p = 16;
    nref = 2;
    dom = surfacemesh.blob(p + 1, nref);

    f = surfacefun(@(x, y, z) cos(6*x) .* y + exp(z), dom);

    figure(1);
    plot(f), hold on, plot(dom), colorbar
    title('f = cos(6x) y + e^z');
    </script>
    </numbl-embed>

Many standard MATLAB arithmetic functions have been overloaded.

.. code-block:: matlab

    x = surfacefun(@(x,y,z) x, dom);
    g = abs(f + 2*x);
    plot(g), colorbar

.. container:: output-image

    .. figure:: images/func_arith.png
        :width: 350px
        :align: center

.. raw:: html

    <numbl-embed lazy label="▶ Edit &amp; run this example" preparing-label="Installing surfacefun…">
    <iframe width="100%" height="560" frameborder="0"></iframe>
    <script type="text/plain" class="numbl-preamble">
    mip load --install flatironinstitute/flatironinstitute/surfacefun
    </script>
    <script type="text/plain" class="numbl-script">
    rng(0);
    p = 16;
    nref = 2;
    dom = surfacemesh.blob(p + 1, nref);

    f = surfacefun(@(x, y, z) cos(6*x) .* y + exp(z), dom);
    x = surfacefun(@(x, y, z) x, dom);

    g = abs(f + 2*x);

    figure(1);
    plot(g), colorbar
    title('g = |f + 2x|');
    </script>
    </numbl-embed>

We can also visualize a ``surfacefun`` using a contour plot:

.. code-block:: matlab

    contour(f, linewidth=2)
    axis off

.. container:: output-image

    .. figure:: images/contour.png
        :width: 200px
        :align: center

.. raw:: html

    <numbl-embed lazy label="▶ Edit &amp; run this example" preparing-label="Installing surfacefun…">
    <iframe width="100%" height="560" frameborder="0"></iframe>
    <script type="text/plain" class="numbl-preamble">
    mip load --install flatironinstitute/flatironinstitute/surfacefun
    </script>
    <script type="text/plain" class="numbl-script">
    rng(0);
    p = 8;
    nref = 0;
    dom = surfacemesh.blob(p + 1, nref);

    f = surfacefun(@(x, y, z) cos(6*x) .* y + exp(z), dom);

    figure(1);
    contour(f, linewidth=2)
    axis off
    title('Contours of f');
    </script>
    </numbl-embed>

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

.. raw:: html

    <numbl-embed lazy label="▶ Edit &amp; run this example" preparing-label="Installing surfacefun…">
    <iframe width="100%" height="560" frameborder="0"></iframe>
    <script type="text/plain" class="numbl-preamble">
    mip load --install flatironinstitute/flatironinstitute/surfacefun
    </script>
    <script type="text/plain" class="numbl-script">
    rng(0);
    p = 16;
    nref = 2;
    dom = surfacemesh.blob(p + 1, nref);

    f = surfacefun(@(x, y, z) cos(6*x) .* y + exp(z), dom);

    [fx, fy, fz] = grad(f);

    figure(1);
    subplot(131), plot(fx), title('\partial_x f')
    subplot(132), plot(fy), title('\partial_y f')
    subplot(133), plot(fz), title('\partial_z f')
    </script>
    </numbl-embed>

Higher-order derivatives may be constructed by composing these operations. For
example, here is the surface Laplacian---or the Laplace--Beltrami
operator---applied to our function:

.. code-block:: matlab

    plot(lap(f)), colorbar

.. container:: output-image

    .. figure:: images/func_lap.png
        :width: 350px
        :align: center

.. raw:: html

    <numbl-embed lazy label="▶ Edit &amp; run this example" preparing-label="Installing surfacefun…">
    <iframe width="100%" height="560" frameborder="0"></iframe>
    <script type="text/plain" class="numbl-preamble">
    mip load --install flatironinstitute/flatironinstitute/surfacefun
    </script>
    <script type="text/plain" class="numbl-script">
    rng(0);
    p = 16;
    nref = 2;
    dom = surfacemesh.blob(p + 1, nref);

    f = surfacefun(@(x, y, z) cos(6*x) .* y + exp(z), dom);

    figure(1);
    plot(lap(f)), colorbar
    title('\Delta_\Gamma f');
    </script>
    </numbl-embed>

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

.. raw:: html

    <numbl-embed lazy label="▶ Edit &amp; run this example (integral, mean, norms)" preparing-label="Installing surfacefun…">
    <iframe width="100%" height="560" frameborder="0"></iframe>
    <script type="text/plain" class="numbl-preamble">
    mip load --install flatironinstitute/flatironinstitute/surfacefun
    </script>
    <script type="text/plain" class="numbl-script">
    rng(0);
    p = 16;
    nref = 2;
    dom = surfacemesh.blob(p + 1, nref);

    f = surfacefun(@(x, y, z) cos(6*x) .* y + exp(z), dom);

    fprintf('integral(f) = %.6f\n', integral(f));     % definite integral over the surface
    fprintf('mean(f)     = %.6f\n', mean(f));          % integral / surface area
    fprintf('L2 norm     = %.6f\n', norm(f));          % sqrt(integral of f^2)
    fprintf('Linf norm   = %.6f\n', norm(f, inf));     % max |f|
    </script>
    </numbl-embed>

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

.. raw:: html

    <numbl-embed lazy label="▶ Edit &amp; run this example" preparing-label="Installing surfacefun…">
    <iframe width="100%" height="560" frameborder="0"></iframe>
    <script type="text/plain" class="numbl-preamble">
    mip load --install flatironinstitute/flatironinstitute/surfacefun
    </script>
    <script type="text/plain" class="numbl-script">
    rng(0);
    p = 16;
    nref = 2;
    dom = surfacemesh.blob(p + 1, nref);

    v = normal(dom);

    figure(1);
    quiver(v, 0.2, 6)
    title('Surface normals');
    </script>
    </numbl-embed>

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

.. raw:: html

    <numbl-embed lazy label="▶ Edit &amp; run this example" preparing-label="Installing surfacefun…">
    <iframe width="100%" height="560" frameborder="0"></iframe>
    <script type="text/plain" class="numbl-preamble">
    mip load --install flatironinstitute/flatironinstitute/surfacefun
    </script>
    <script type="text/plain" class="numbl-script">
    rng(0);
    p = 16;
    nref = 2;
    dom = surfacemesh.blob(p + 1, nref);

    f = surfacefun(@(x, y, z) cos(6*x) .* y + exp(z), dom);

    figure(1);
    quiver(grad(f), 0.05, 6)
    title('Surface gradient of f (tangent to the surface)');

    % Identity: the divergence of the gradient is the Laplacian.
    fprintf('|| div(grad f) - lap(f) || = %.3e\n', norm(div(grad(f)) - lap(f)));
    </script>
    </numbl-embed>

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

The mean curvature of a surface can be related to its the normal vector field
via the surface divergence:

.. code-block:: matlab

    plot(div(v)/2)

.. container:: output-image

    .. figure:: images/vec_div.png
        :width: 400px
        :align: center

.. raw:: html

    <numbl-embed lazy label="▶ Edit &amp; run this example" preparing-label="Installing surfacefun…">
    <iframe width="100%" height="560" frameborder="0"></iframe>
    <script type="text/plain" class="numbl-preamble">
    mip load --install flatironinstitute/flatironinstitute/surfacefun
    </script>
    <script type="text/plain" class="numbl-script">
    rng(0);
    p = 16;
    nref = 2;
    dom = surfacemesh.blob(p + 1, nref);

    v = normal(dom);

    figure(1);
    plot(div(v)/2), colorbar
    title('Mean curvature  div(n)/2');
    </script>
    </numbl-embed>

We can also take the surface curl of a ``surfacefunv``:

.. code-block:: matlab

    v = surfacefunv(@(x,y,z) cos(2*x), ...
                    @(x,y,z) sin(4*y), ...
                    @(x,y,z) sin(3*z), dom);
    quiver(curl(v), 0.1, 6)

.. container:: output-image

    .. figure:: images/vec_curl.png
        :width: 350px
        :align: center

.. raw:: html

    <numbl-embed lazy label="▶ Edit &amp; run this example" preparing-label="Installing surfacefun…">
    <iframe width="100%" height="560" frameborder="0"></iframe>
    <script type="text/plain" class="numbl-preamble">
    mip load --install flatironinstitute/flatironinstitute/surfacefun
    </script>
    <script type="text/plain" class="numbl-script">
    rng(0);
    p = 16;
    nref = 2;
    dom = surfacemesh.blob(p + 1, nref);

    v = surfacefunv(@(x, y, z) cos(2*x), ...
                    @(x, y, z) sin(4*y), ...
                    @(x, y, z) sin(3*z), dom);

    figure(1);
    quiver(curl(v), 0.1, 6)
    title('Surface curl of a vector field');
    </script>
    </numbl-embed>
