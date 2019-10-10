# Problem Set 1

Due 1 Nov 2019, 5 p.m. Eastern.  Solutions can be in any language you like. I
will take a snapshot of the repository on GitHub at the deadline; I encourage
you to develop the code over many commits and to upload early and often, as if
you were working on a software project.  The entire assignment is worth 100
points.  There are "challenge" problems at the end of the assignment; completion
of these problems can boost your point total beyond 100 points, but the maximum
score I will record is 100.  Note that *all* problems require tests that
demonstrate the correctness of your code; these can be in the form of a program
that is executed and reports success/failures, in the from of a Jupyter notebook
with comments interspersed with executed tests / plots, or (preferably) both.

## Problem 1

Solving quadratic equations.  

- [ ] (5 points) Write a function that solves quadratic equations given their coefficients.
  If you implemented such a program as a `solve_quadratic` function in Python,
  the following should work:
  ```python
  a = 1.0
  b = 2.0
  c = 3.0
  x1, x2 = solve_quadratic(a, b, c)
  abs(a*x1*x1 + b*x1 + c), abs(a*x2*x2 + b*x2 + c) # Should both be close to zero
  ```

- [ ] (5 points) Write a test suite for your new function to ensure that it is correct.
  Your test suite should include some randomized tests.

- [ ] (5 points) The solution does not change if the coefficients are multiplied by a large
  number.  Ensure that your program gives the same answer as above for `a =
  1e20`, `b = 2e20`, and `c = 3e20`.  

- [ ] (5 points) How accurately does your function solve the quadratic with roots `x1 =
  100` and `x2 = 1/1000`?  Modify it to achieve more precision in cases like
  this ("cases like this" are cases where the product of the roots is much
  smaller than their sum).

## Problem 2

Numerical derivatives.  Recall that the simplest numerical derivative estimate
implements the limit definition faithfully (e.g., in Julia):
```julia
"""    numerical_D(f)

Returns a function that computes the forward-difference estimate of the
derivative of `f`.
"""
function numerical_D(f)
  function df(x)
    h = sqrt(eps(x)) # sqrt of the smallest number we can add to x and get something different
    (f(x+h) - f(x))/h
  end
  df
end
```

- [ ] (5 points) The above approximation is first-order accurate (truncation
  error scales as `h^2`).  Derive (or intuit?) and implement a second-order
  accurate approximation.  Write a test suite for this function that verifies
  that it correctly computes derivatives and that the error term is second
  order.

- [ ] (5 points) For the reasons discussed in class, "single step" derivative
  computations like this are dominated by roundoff error well before they
  achieve accuracy close to machine precision.  Write a function that
  accumulates several (it is up to you how many) such approximations to the
  derivative, and then extrapolates to `h = 0`.  This technique, mentioned
  several times in class, is called [Richardson
  extrapolation](https://en.wikipedia.org/wiki/Richardson_extrapolation), and it
  is *extraordinarily powerful*.  Don't forget the tests!

- [ ] (10 points) Generalize your code above to *choose* how many single-step
  approximations it will accumulate before terminating and returning the
  derivative.  You should request an accuracy requirement from the user
  (preferably, if your language allows it, with a sensible default), and try to
  meet it using successive "single step" computations followed by higher and
  higher order extrapolation.  Set a reasonable limit on the total number of
  attempts (i.e. the maximum order of the extrapolation) before bailing out with
  an error message of some kind.  

# Problem 3

Cosmology!  (Don't worry---it's just an excuse to do integrals and find roots.)  

- [ ] (5 points) In class we discussed the trapezoidal rule for integrals.
  Using the fact that the error term for this rule scales as the step size
  squared, show how you can combine two trapezoidal steps with size `h/2` and
  one trapezoidal step with size `h` to eliminate the second-order error.  You
  have just discovered Simpson's method.  Write a function that computes the
  definite integral of an arbitrary function to a user-specified accuracy.  Your
  function should compare Simpson's rule at stepsize `h` *and* the trapezoidal
  rule at the same stepsize (no need for new function evaluations!) to estimate
  the Simpson's error; if the error is too large, it should divide the interval
  into two and sum the integrals over the sub-divided interval.  Include a
  test-suite for your function that demonstrates its accuracy on several trial
  functions that can be integrated analytically.

- [ ] (7 points) Computing cosmological distances involves an integral over
  spacetime.  (If you're curious, see [Hogg
  (1999)](https://arxiv.org/abs/astro-ph/9905116) or consult a textbook such as
  [Peebles
  (1993)](https://press.princeton.edu/books/paperback/9780691019338/principles-of-physical-cosmology)
  by one of this years' Nobel winners.)  For example, the luminosity distance in
  a flat universe is given in terms of the Hubble constant, and matter density
  by `d_L = c/H0*(1+z)*d(z)` where `d(z)` is an integral of `1/E(Z)` from `Z=0`
  to `z`, with (Python syntax)

  ```python
  def E(z, Om):
    return np.sqrt((1-Om) + Om*(1+z)*(1+z)*(1+z))
  ```

  Write a function that computes `d_L(z)` to specified accuracy using your
  function from the previous part.  A test suite here is harder, but you can at
  least compare a plot of your function's output at sensible cosmological
  parameters to Hogg (1999).

- [ ] (8 points) The [astropy](https://www.astropy.org) library includes a helpful
  cosmological function called `z_at_value`.  It inverts the a cosmoligical
  relation like `d_L(z)` using a root finder.  For example, the following
  returns the redshift at which the cosmology from [Ade, et al.
  (2015)](http://dx.doi.org/10.1051/0004-6361/201525830) gives a luminosity
  distance of 4 Gpc.
  ```python
  cosmo.z_at_value(Planck15.luminosity_distance, 4*u.Gpc) # z = 0.65
  ```
  Using a bisection method, implement your own version of `z_at_value`.  Your
  test suite should include some randomized tests.  You can probably assume that
  the user is not interested in `z < 0` (truly an error) or `z > 100` (a
  specialized regime relevant only to a few), but you should verify this on
  input.

# Problem 4

Integrating ODEs.

- [ ] (15 points) Look up the coefficients for the "embedded" fourth and fifth order
  Runge-Kutta method of Fehlberg (e.g.
  [here](https://en.wikipedia.org/wiki/List_of_Runge–Kutta_methods#Embedded_methods)).
  Implement a stepper function for the method that returns an error estimate and
  a control loop much like what we did in class to evolve the solution to an ODE
  system over a finite time with bounded error.  Don't forget the test suite!

- [ ] (5 points) Now that you have a good, high-order integrator with good stepsize
  control, use your integrator to solve the
  [Lorentz](https://en.wikipedia.org/wiki/Lorenz_system) system of equations and
  show some pretty pictures of their chaotic orbits!

# Problem 5

Variational / symplectic integrators.

- [ ] (15 points) Implement a leapfrog integrator that is variational / symplectic, as
  described in class.  You can assume that you have a mechanical system
  described by the usual kinetic energy `0.5*m*v^2` and an arbitrary.  You can
  either formulate the problem as an approximation to the action integral with a
  linear-in-time trajectory and the trapezoid rule, or as a simple mapping
  integrator that splits the Hamiltonian into `H = T+V`.  Work in terms of
  velocities, not momenta, so that you do not need to specify the masses of the
  bodies; assume that the caller has supplied a function that consumes the
  positions of the bodies in the system and supplies the (vector) acceleration
  on each body.  A sample calling sequence (Julia syntax) would be
  ```julia
  x_new, v_new = leapfrog_step(h, x, v, a_func)
  ```
  A sensible test suite would be to apply this to the simple harmonic oscillator
  and verify that linear and angular momentum are conserved (to machine
  precision) and that the energy error is bounded for reasonable step sizes
  (much smaller than the fundamental oscillation period).

- [ ] (5 points) Apply your new integrator to explore the motion of stars in a galactic
  potential: solve the
  [Henon-Heiles](https://en.wikipedia.org/wiki/Hénon–Heiles_system) Hamiltonian
  with `lambda = 1` and various initial conditions.  Show some orbits.  Can you
  find some that are chaotic?

# Challenge Problem 1

(50 points) Implement either forward or backward autodiff in your favorite
language.  You should have a test suite that verifies your code works for basic
arithmetic, exponentials, and trig functions at least.

# Challenge Problem 2

(50 points) Apply a variational integrator to measuring the temperature of a
coupled spring system with weakly non-linear coupling, as described in [Lew, et
al. (2004)](http://doi.wiley.com/10.1002/nme.958) Figure 1.

# Challenge Problem 3

(50 points) Use the [JPL Ephemeris](https://ssd.jpl.nasa.gov/?ephemerides) to
set up the locations of the outer solar system giants (Jupiter, Saturn, Uranus,
and Neptune) on the date this assignment is due, and integrate them at least 1
Myr into the future.  You can use either the variational integrator from Problem
5 or your Runge-Kutta integrator from problem 4.
