#! /usr/bin/env python3
#
def beale_f ( x, n ):

#*****************************************************************************80
#
## BEALE_F evaluates the Beale function.
#
#  Discussion:
#
#    The function is the sum of the squares of three functions.
#
#    This function has a valley approaching the line X(2) = 1.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    E Beale,
#    On an Iterative Method for Finding a Local Minimum of a Function
#    of More than One Variable,
#    Technical Report 25, Statistical Techniques Research Group,
#    Princeton University, 1958.
#
#    Richard Brent,
#    Algorithms for Finding Zeros and Extrema of Functions Without
#    Calculating Derivatives,
#    Stanford University Technical Report STAN-CS-71-198.
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  c1 = 1.5
  c2 = 2.25
  c3 = 2.625

  fx1 = c1 - x[0] * ( 1.0 - x[1]      )
  fx2 = c2 - x[0] * ( 1.0 - x[1] ** 2 )
  fx3 = c3 - x[0] * ( 1.0 - x[1] ** 3 )

  value = fx1 ** 2 + fx2 ** 2 + fx3 ** 2

  return value

def beale_test ( ):

#*****************************************************************************80
#
## BEALE_TEST calls PRAXIS for the Beale function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform
  
  n = 2

  print ( '' )
  print ( 'BEALE_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Beale function.' )

  t0 = 0.00001
  h0 = 0.25
  prin = 0

  x = np.array ( [ 0.1, 0.1 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( beale_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, beale_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( beale_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'BEALE_TEST:' )
  print ( '  Normal end of execution.' )
  return
 
def box_f ( x, n ):

#*****************************************************************************80
#
## BOX_F evaluates the Box function.
#
#  Discussion:
#
#    The function is formed by the sum of squares of 10 separate terms.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  import numpy as np

  value = 0.0

  for i in range ( 1, 11 ):

    c = - i / 10.0

    fx = np.exp ( c * x[0] ) - np.exp ( c * x[1] ) \
      - x[2] * ( np.exp ( c ) - np.exp ( 10.0 * c ) )

    value = value + fx ** 2

  return value

def box_test ( ):

#*****************************************************************************80
#
## BOX_TEST calls PRAXIS for the Box function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 3

  print ( '' )
  print ( 'BOX_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Box function.' )

  t0 = 0.00001
  h0 = 20.0
  prin = 0

  x = np.array ( [ 0.0, 10.0, 20.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( box_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, box_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( box_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'BOX_TEST:' )
  print ( '  Normal end of execution.' )
  return

def chebyquad_f ( x, n ):

#*****************************************************************************80
#
## CHEBYQUAD_F evaluates the Chebyquad function.
#
#  Discussion:
#
#    The function is formed by the sum of squares of N separate terms.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  import numpy as np

  fvec = np.zeros ( n )

  for j in range ( 0, n ):

    t1 = 1.0;
    t2 = 2.0 * x[j] - 1.0
    t = 2.0 * t2

    for i in range ( 0, n ):
      fvec[i] = fvec[i] + t2
      th = t * t2 - t1
      t1 = t2
      t2 = th

  for i in range ( 0, n ):
    fvec[i] = fvec[i] / n
    if ( ( i % 2 ) == 1 ):
      fvec[i] = fvec[i] + 1.0 / ( i * ( i + 2 ) )
#
#  Compute F.
#
  value = 0.0
  for i in range ( 0, n ):
    value = value + fvec[i] ** 2

  return value

def chebyquad_test ( ):

#*****************************************************************************80
#
## CHEBYQUAD_TEST calls PRAXIS for the Chebyquad function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 8

  print ( '' )
  print ( 'CHEBYQUAD_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Chebyquad function.' )

  t0 = 0.00001
  h0 = 0.1
  prin = 0

  x = np.zeros ( n )

  for i in range ( 0, n ):
    x[i] = float ( i + 1 ) / float ( n + 1 )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( chebyquad_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, chebyquad_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( chebyquad_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'CHEBYQUAD_TEST:' )
  print ( '  Normal end of execution.' )
  return

def cube_f ( x, n ):

#*****************************************************************************80
#
## CUBE_F evaluates the Cube function.
#
#  Discussion:
#
#    The function is the sum of the squares of two functions.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  fx1 = 10.0 * ( x[1] - x[0] ** 3 )
  fx2 = 1.0 - x[0]

  value = fx1 ** 2 + fx2 ** 2

  return value

def cube_test ( ):

#*****************************************************************************80
#
## CUBE_TEST calls PRAXIS for the Cube function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 2

  print ( '' )
  print ( 'CUBE_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Cube function.' )

  t0 = 0.00001
  h0 = 1.0
  prin = 0

  x = np.array ( [ -1.2, -1.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( cube_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, cube_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( cube_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'CUBE_TEST:' )
  print ( '  Normal end of execution.' )
  return

def helix_f ( x, n ):

#*****************************************************************************80
#
## HELIX_F evaluates the Helix function.
#
#  Discussion:
#
#    The function is the sum of the squares of three functions.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  import numpy as np

  r = np.linalg.norm ( x )

  if ( 0.0 <= x[0] ):
    theta = 0.5 * np.arctan2 ( x[1], x[0] ) / np.pi
  else:
    theta = 0.5 * ( np.arctan2 ( x[1], x[0] ) + np.pi ) / np.pi

  fx1 = 10.0 * ( x[2] - 10.0 * theta )
  fx2 = 10.0 * ( r - 1.0 )
  fx3 = x[2]

  value = fx1 ** 2 + fx2 ** 2 + fx3 ** 2

  return value

def helix_test ( ):

#*****************************************************************************80
#
## HELIX_TEST calls PRAXIS for the Fletcher-Powell Helix function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 3

  print ( '' )
  print ( 'HELIX_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Fletcher-Powell Helix function.' )

  t0 = 0.00001
  h0 = 1.0
  prin = 0

  x = np.array ( [ -1.0, 0.0, 0.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( helix_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, helix_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( helix_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'HELIX_TEST:' )
  print ( '  Normal end of execution.' )
  return

def hilbert_f ( x, n ):

#*****************************************************************************80
#
## HILBERT_F evaluates the Hilbert function.
#
#  Discussion:
#
#    The function is a positive definite quadratic function of the form
#
#      f(x) = x' A x
#
#    where A is the Hilbert matrix, A(I,J) = 1/(I+J-1).
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  value = 0.0

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      value = value + x[i] * x[j] / ( i + j + 1 )

  return value

def hilbert_test ( ):

#*****************************************************************************80
#
## HILBERT_TEST calls PRAXIS for the Hilbert function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 10

  print ( '' )
  print ( 'HILBERT_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Hilbert function.' )

  t0 = 0.00001
  h0 = 10.0
  prin = 0

  x = np.ones ( n )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( hilbert_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, hilbert_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( hilbert_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'HILBERT_TEST:' )
  print ( '  Normal end of execution.' )
  return

def powell3d_f ( x, n ):

#*****************************************************************************80
#
## POWELL3D_F evaluates the Powell 3D function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    M J D Powell,
#    An Efficient Method for Finding the Minimum of a Function of
#    Several Variables Without Calculating Derivatives,
#    Computer Journal,
#    Volume 7, Number 2, pages 155-162, 1964.
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  import numpy as np

  value = 3.0 - 1.0 / ( 1.0 + ( x[0] - x[1] ) ** 2 ) \
    - np.sin ( 0.5 * np.pi * x[1] * x[2] ) \
    - np.exp ( - ( ( x[0] - 2.0 * x[1] + x[2] ) / x[1] ) ** 2 )

  return value

def powell3d_test ( ):

#*****************************************************************************80
#
## POWELL3D_TEST calls PRAXIS for the Powell 3D function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 3

  print ( '' )
  print ( 'POWELL3D_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Powell 3D function.' )

  t0 = 0.00001
  h0 = 1.0
  prin = 0

  x = np.array ( [ 0.0, 1.0, 2.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( powell3d_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, powell3d_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( powell3d_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'POWELL3D_TEST:' )
  print ( '  Normal end of execution.' )
  return

def flin ( n, jsearch, l, f, x, nf, v, q0, q1, qd0, qd1, qa, qb, qc ):

#*****************************************************************************80
#
## FLIN is the function of one variable to be minimized by MINNY.
#
#  Discussion:
#
#    F(X) is a scalar function of a vector argument X.
#
#    A minimizer of F(X) is sought along a line or parabola.
#
#    This function has been modified, by removing the occurrence of a
#    common block, so that it looks more like a "normal" function that does
#    not rely so strongly on peculiarities of FORTRAN.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#  Parameters:
#
#    Input, integer N, the number of variables.
#
#    Input, integer JSEARCH, indicates the kind of search.
#    If J is a legal column index, linear search in direction of V(*,JSEARCH).
#    Otherwise, then the search is parabolic, based on X, Q0 and Q1.
#
#    Input, real L, is the parameter determining the particular
#    point at which F is to be evaluated.
#    For a linear search, L is the step size.
#    For a quadratic search, L is a parameter which specifies
#    a point in the plane of X, Q0 and Q1.
#
#    Input, real F ( X, N ), the function to be minimized.
#
#    Input, real X(N), the base point of the search.
#
#    Input/output, integer NF, the function evaluation counter.
#
#    Input, real V(N,N), a matrix whose columns constitute
#    search directions.
#
#    Input, real Q0(N), Q1(N), two auxiliary points used to
#    determine the plane when a quadratic search is performed.
#
#    Input, real QD0, QD1, values needed to compute the
#    coefficients QA, QB, QC.
#
#    Input/output, real QA, QB, QC, coefficients used to combine
#    Q0, X, and A1 if a quadratic search is used.  (Yes, technically
#    these are input quantities as well.)
#
#    Output, real VALUE, the value of the function at the
#    minimizing point.
#
  import numpy as np

  t = np.zeros ( n )
#
#  The search is linear.
#
  if ( 0 <= jsearch ):

    t[0:n] = x[0:n] + l * v[0:n,jsearch]
#
#  The search is along a parabolic space curve.
#
  else:

    qa =                 l * ( l - qd1 ) /       ( qd0 + qd1 ) / qd0
    qb = - ( l + qd0 ) *     ( l - qd1 ) / qd1                 / qd0
    qc =   ( l + qd0 ) * l               / qd1 / ( qd0 + qd1 )

    t[0:n] = qa * q0[0:n] + qb * x[0:n] + qc * q1[0:n]
#
#  The function evaluation counter NF is incremented.
#
  nf = nf + 1
#
#  Evaluate the function.
#
  value = f ( t, n )

  return value, nf, qa, qb, qc

def minfit ( n, tol, a ):

#*****************************************************************************80
#
## MINFIT computes the singular value decomposition of an N by N array.
#
#  Discussion:
#
#    This is an improved version of the EISPACK routine MINFIT
#    restricted to the case M = N and P = 0.
#
#    The singular values of the array A are returned in Q.  A is
#    overwritten with the orthogonal matrix V such that U * diag(Q) = A * V,
#    where U is another orthogonal matrix.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#    James Wilkinson, Christian Reinsch,
#    Handbook for Automatic Computation,
#    Volume II, Linear Algebra, Part 2,
#    Springer Verlag, 1971.
#
#    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, Yasuhiko Ikebe,
#    Virginia Klema, Cleve Moler,
#    Matrix Eigensystem Routines, EISPACK Guide,
#    Lecture Notes in Computer Science, Volume 6,
#    Springer Verlag, 1976,
#    ISBN13: 978-3540075462,
#    LC: QA193.M37.
#
#  Parameters:
#
#    Input, integer N, the order of the matrix A.
#
#    Input, real TOL, a tolerance which determines when a vector
#    (a column or part of a column of the matrix) may be considered
#    "essentially" equal to zero.
#
#    Input/output, real A(N,N).  On input, an N by N array whose
#    singular value decomposition is desired.  On output, the
#    SVD orthogonal matrix factor V.
#
#    Output, real Q(N), the singular values.
#
  import numpy as np

  kt_max = 30

  e = np.zeros ( n )
  q = np.zeros ( n )
#
#  Householder's reduction to bidiagonal form.
#
  if ( n == 1 ):
    q[0] = a[0,0]
    a[0,0] = 1.0
    return a, q

  g = 0.0
  x = 0.0

  for i in range ( 0, n ):

    e[i] = g
    l = i + 1

    s = 0.0
    for i1 in range ( i, n ):
      s = s + a[i1,i] ** 2

    g = 0.0

    if ( tol <= s ):

      f = a[i,i]

      g = np.sqrt ( s )
      if ( 0.0 <= f ):
        g = - g

      h = f * g - s
      a[i,i] = f - g

      for j in range ( l, n ):

        f = 0.0
        for i1 in range ( i, n ):
          f = f + a[i1,i] * a[i1,j]
        f = f / h

        for i1 in range ( i, n ):
          a[i1,j] = a[i1,j] + f * a[i1,i]

    q[i] = g

    s = 0.0
    for j1 in range ( l, n ):
      s = s + a[i,j1] ** 2

    g = 0.0

    if ( tol <= s ):

      if ( i < n - 1 ):
        f = a[i,i+1]

      g = np.sqrt ( s )
      if ( 0.0 <= f ):
        g = - g

      h = f * g - s

      if ( i < n - 1 ):

        a[i,i+1] = f - g

        for j1 in range ( l, n ):
          e[j1] = a[i,j1] / h

        for j in range ( l, n ):

          s = 0.0
          for j1 in range ( l, n ):
            s = s + a[j,j1] * a[i,j1]

          for j1 in range ( l, n ):
            a[j,j1] = a[j,j1] + s * e[j1]

    y = abs ( q[i] ) + abs ( e[i] )

    x = max ( x, y )
#
#  Accumulation of right-hand transformations.
#
  a[n-1,n-1] = 1.0
  g = e[n-1]
  l = n - 1

  for i in range ( n - 2, -1, -1 ):

    if ( g != 0.0 ):

      h = a[i,i+1] * g

      for i1 in range ( l, n ):
        a[i1,i] = a[i,i1] / h

      for j in range ( l, n ):

        s = 0.0
        for j1 in range ( l, n ):
          s = s + a[i,j1] * a[j1,j]

        for i1 in range ( l, n ):
          a[i1,j] = a[i1,j] + s * a[i1,i]

    for j1 in range ( l, n ):
      a[i,j1] = 0.0

    for i1 in range ( l, n ):
      a[i1,i] = 0.0

    a[i,i] = 1.0

    g = e[i]

    l = i
#
#  Diagonalization of the bidiagonal form.
#
  epsx = r8_epsilon ( ) * x

  for k in range ( n - 1, -1, -1 ):

    kt = 0

    while ( True ):

      kt = kt + 1

      if ( kt_max < kt ):
        e[k] = 0.0
        print ( '' )
        print ( 'MINFIT - Fatal error!' )
        print ( '  The QR algorithm failed to converge.' )
        exit ( 'MINFIT - Fatal error!' )

      skip = False

      for l2 in range ( k, -1, -1 ):

        l = l2

        if ( abs ( e[l] ) <= epsx ):
          skip = True
          break

        if ( 0 < l ):
          if ( abs ( q[l-1] ) <= epsx ):
            break
#
#  Cancellation of E(L) if 1 < L.
#
      if ( not skip ):

        c = 0.0
        s = 1.0

        for i in range ( l, k + 1 ):

          f = s * e[i]
          e[i] = c * e[i]

          if ( abs ( f ) <= epsx ):
            break

          g = q[i]
#
#  q(i) = h = sqrt(g*g + f*f).
#
          h = r8_hypot ( f, g )

          q[i] = h

          if ( h == 0.0 ):
            g = 1.0
            h = 1.0

          c =   g / h
          s = - f / h
#
#  Test for convergence for this index K.
#
      z = q[k]

      if ( l == k ):
        if ( z < 0.0 ):
          q[k] = - z
          for i1 in range ( 0, n ):
            a[i1,k] = - a[i1,k]
        break
#
#  Shift from bottom 2*2 minor.
#
      x = q[l]
      y = q[k-1]
      g = e[k-1]
      h = e[k]
      f = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) / ( 2.0 * h * y )

      g = r8_hypot ( f, 1.0 )

      if ( f < 0.0 ):
        temp = f - g
      else:
        temp = f + g

      f = ( ( x - z ) * ( x + z ) + h * ( y / temp - h ) ) / x
#
#  Next QR transformation.
#
      c = 1.0
      s = 1.0

      for i in range ( l + 1, k + 1 ):

        g = e[i]
        y = q[i]
        h = s * g
        g = g * c

        z = r8_hypot ( f, h )

        e[i-1] = z

        if ( z == 0.0 ):
          f = 1.0
          z = 1.0

        c = f / z
        s = h / z
        f =   x * c + g * s
        g = - x * s + g * c
        h = y * s
        y = y * c

        for j in range ( 0, n ):
          x = a[j,i-1]
          z = a[j,i]
          a[j,i-1] = x * c + z * s
          a[j,i] = - x * s + z * c

        z = r8_hypot ( f, h )

        q[i-1] = z

        if ( z == 0.0 ):
          f = 1.0
          z = 1.0

        c = f / z
        s = h / z
        f =   c * g + s * y
        x = - s * g + c * y

      e[l] = 0.0
      e[k] = f
      q[k] = x

  return a, q

def minfit_test ( ):

#*****************************************************************************80
#
## MINFIT_TEST tests MINFIT, which is a sort of SVD computation.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 5

  print ( '' )
  print ( 'MINFIT_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  MINFIT computes part of the SVD of a matrix A.' )
  print ( '    SVD: A = U * D * V\'' )
  print ( '  MINFIT is given A, and returns the diagonal D' )
  print ( '  and the orthogonal matrix V.' )

  a = np.zeros ( [ n, n ] )

  for i in range ( 0, n ):
    a[i,i] = 2.0

  for i in range ( 0, n - 1 ):
    a[i,i+1] = -1.0

  for i in range ( 1, n ):
    a[i,i-1] = -1.0

  r8mat_print ( n, n, a, '  The matrix A:' )
#
#  Numpy's EPS function is not easy to find!
#
  eps = r8_epsilon ( )

  tol = np.sqrt ( eps )

  a, d = minfit ( n, tol, a )

  r8mat_print ( n, n, a, '  The vector V:' )

  r8vec_print ( n, d, '  The singular values D:' )
#
#  Because A is positive definite symmetric, the "missing" matrix V = U.
#
  print ( '' )
  print ( '  Because A is positive definite symmetric,' )
  print ( '  we can reconstruct it as A = V * D * V\'' )

  a2 = np.zeros ( [ n, n ] )

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      for k in range ( 0, n ):
        a2[i,j] = a2[i,j] + a[i,k] * d[k] * a[j,k]

  r8mat_print ( n, n, a2, '  The product A2 = V * D * V\'' )
#
#  Terminate.
#
  print ( '' )
  print ( 'MINFIT_TEST:' )
  print ( '  Normal end of execution.' )
  return

def minny ( n, jsearch, nits, d2, x1, f1, fk, f, x, t, h, v, q0, q1, \
  nl, nf, dmin, ldt, fx, qa, qb, qc, qd0, qd1 ):

#*****************************************************************************80
#
## MINNY minimizes a scalar function of N variables along a line.
#
#  Discussion:
#
#    MINNY minimizes F along the line from X in the direction V(*,J) unless
#    J is less than 1, when a quadratic search is made in the plane
#    defined by Q0, Q1 and X.
#
#    If FK = true, then F1 is FLIN(X1).  Otherwise X1 and F1 are ignored
#    on entry unless final FX is greater than F1.
#
#    This function was modified by removing the common blocks
#    and the use of labeled statements, 28 July 2016.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#  Parameters:
#
#    Input, integer N, the number of variables.
#
#    Input, integer JSEARCH, indicates the kind of search.
#    If JSEARCH is a legal column index, linear search in direction of V(*,J).
#    Otherwise, the search is parabolic, based on X, Q0 and Q1.
#
#    Input, integer NITS, the maximum number of times the interval
#    may be halved to retry the calculation.
#
#    Input/output, real D2, is either zero, or an approximation to
#    the value of (1/2) times the second derivative of F.
#
#    Input/output, real X1, on entry, an estimate of the
#    distance from X to the minimum along V(*,J), or, if J = 0, a curve.
#    On output, the distance between X and the minimizer that was found.
#
#    Input/output, real F1, ?
#
#    Input, logical FK if FK is TRUE, then on input F1 contains
#    the value FLIN(X1).
#
#    Input, real F ( X, N ), the function to be minimized.  
#
#    Input/output, real X(N), ?
#
#    Input, real T, ?
#
#    Input, real H, ?
#
#    Input, real V(N,N), a matrix whose columns are direction
#    vectors along which the function may be minimized.
#
#    Input, real Q0(N), an auxiliary point used to define
#    a curve through X.
#
#    Input, real Q1(N), an auxiliary point used to define
#    a curve through X.
#
#    Input/output, integer NL, the number of linear searches.
#
#    Input/output, integer NF, the number of function evaluations.
#
#    Input, real DMIN, an estimate for the smallest eigenvalue.
#
#    Input, real LDT, the length of the step.
#
#    Input/output, real FX, the value of F(X,N).
#
#    Input/output, real QA, QB, QC, ?
#
#    Input, real QD0, QD1, ?.
#
  import numpy as np

  machep = r8_epsilon ( )
  small = machep ** 2
  m2 = np.sqrt ( machep )
  m4 = np.sqrt ( m2 )
  sf1 = f1
  sx1 = x1
  k = 0
  xm = 0.0
  fm = fx
  f0 = fx
  dz = ( d2 < machep )
#
#  Find the step size.
#
  s = np.linalg.norm ( x )

  if ( dz ):
    temp = dmin
  else:
    temp = d2

  t2 = m4 * np.sqrt ( abs ( fx ) / temp + s * ldt ) + m2 * ldt
  s = m4 * s + t
  if ( dz and s < t2 ):
    t2 = s

  t2 = max ( t2, small )
  t2 = min ( t2, 0.01 * h )

  if ( fk and f1 <= fm ):
    xm = x1
    fm = f1

  if ( ( not fk ) or abs ( x1 ) < t2 ):

    if ( 0.0 <= x1 ):
      temp = 1.0
    else:
      temp = - 1.0

    x1 = temp * t2

    f1, nf, qa, qb, qc = flin ( n, jsearch, x1, f, x, nf, v, \
      q0, q1, qd0, qd1, qa, qb, qc )

  if ( f1 <= fm ):
    xm = x1
    fm = f1
#
#  Evaluate FLIN at another point and estimate the second derivative.
#
  while ( True ):

    if ( dz ):

      if ( f1 <= f0 ):
        x2 = 2.0 * x1
      else:
        x2 = - x1

      f2, nf, qa, qb, qc = flin ( n, jsearch, x2, f, x, nf, v, \
        q0, q1, qd0, qd1, qa, qb, qc )

      if ( f2 <= fm ):
        xm = x2
        fm = f2

      d2 = ( x2 * ( f1 - f0 ) - x1 * ( f2 - f0 ) ) \
        / ( ( x1 * x2 ) * ( x1 - x2 ) )
#
#  Estimate the first derivative at 0.
#
    d1 = ( f1 - f0 ) / x1 - x1 * d2
    dz = True
#
#  Predict the minimum.
#
    if ( d2 <= small ):

      if ( 0.0 <= d1 ):
        x2 = - h
      else:
        x2 = h

    else:

      x2 = ( - 0.5 * d1 ) / d2

    if ( h < abs ( x2 ) ):

      if ( x2 <= 0.0 ):
        x2 = - h
      else:
        x2 = h
#
#  Evaluate F at the predicted minimum.
#
    ok = True

    while ( True ):

      f2, nf, qa, qb, qc = flin ( n, jsearch, x2, f, x, nf, v, \
        q0, q1, qd0, qd1, qa, qb, qc )

      if ( nits <= k or f2 <= f0 ):
        break

      k = k + 1

      if ( f0 < f1 and 0.0 < x1 * x2 ):
        ok = False
        break

      x2 = 0.5 * x2

    if ( ok ):
      break
#
#  Increment the one-dimensional search counter.
#
  nl = nl + 1

  if ( fm < f2 ):
    x2 = xm
  else:
    fm = f2
#
#  Get a new estimate of the second derivative.
#
  if ( small < abs ( x2 * ( x2 - x1 ) ) ):
    d2 = ( x2 * ( f1 - f0 ) - x1 * ( fm - f0 ) ) / ( ( x1 * x2 ) * ( x1 - x2 ) )
  else:
    if ( 0 < k ):
      d2 = 0.0

  d2 = max ( d2, small )

  x1 = x2
  fx = fm

  if ( sf1 < fx ):
    fx = sf1
    x1 = sx1
#
#  Update X for linear but not parabolic search.
#
  if ( 0 <= jsearch ):
    x[0:n] = x[0:n] + x1 * v[0:n,jsearch]

  return d2, x1, f1, x, nl, nf, fx, qa, qb, qc

def praxis ( t0, h0, n, prin, x, f ):

#*****************************************************************************80
#
## PRAXIS seeks an N-dimensional minimizer X of a scalar function F(X).
#
#  Discussion:
#
#    PRAXIS returns the minimum of the function F(X,N) of N variables
#    using the principal axis method.  The gradient of the function is
#    not required.
#
#    The approximating quadratic form is
#
#      Q(x') = F(x,n) + (1/2) * (x'-x)' * A * (x'-x)
#
#    where X is the best estimate of the minimum and
#
#      A = inverse(V') * D * inverse(V)
#
#    V(*,*) is the matrix of search directions
#    D(*) is the array of second differences.
#
#    If F(X) has continuous second derivatives near X0, then A will tend
#    to the hessian of F at X0 as X approaches X0.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#  Parameters:
#
#    Input, real T0, is a tolerance.  PRAXIS attempts to return
#    praxis = f(x) such that if X0 is the true local minimum near X, then
#    norm ( x - x0 ) < T0 + sqrt ( EPSILON ( X ) ) * norm ( X ),
#    where EPSILON ( X ) is the machine precision for X.
#
#    Input, real H0, is the maximum step size.  H0 should be
#    set to about the maximum distance from the initial guess to the minimum.
#    If H0 is set too large or too small, the initial rate of
#    convergence may be slow.
#
#    Input, integer N, the number of variables.
#
#    Input, integer PRIN, controls printing intermediate results.
#    0, nothing is printed.
#    1, F is printed after every n+1 or n+2 linear minimizations.
#       final X is printed, but intermediate X is printed only
#       if N is at most 4.
#    2, the scale factors and the principal values of the approximating
#       quadratic form are also printed.
#    3, X is also printed after every few linear minimizations.
#    4, the principal vectors of the approximating quadratic form are
#       also printed.
#
#    Input/output, real X(N), is an array containing on entry a
#    guess of the point of minimum, on return the estimated point of minimum.
#
#    Input, real F ( X, N ), the function to be minimized.
#
#    Output, real PRAXIS, the function value at the minimizer.
#
#  Local parameters:
#
#    Local, real DMIN, an estimate for the smallest eigenvalue.
#
#    Local, real FX, the value of F(X,N).
#
#    Local, logical ILLC, is TRUE if the system is ill-conditioned.
#
#    Local, real LDT, the length of the step.
#
#    Local, integer NF, the number of function evaluations.
#
#    Local, integer NL, the number of linear searches.
#
  import numpy as np
#
#  Initialization.
#
  machep = r8_epsilon ( )
  small = machep * machep
  vsmall = small * small
  large = 1.0 / small
  vlarge = 1.0 / vsmall
  m2 = np.sqrt ( machep )
  m4 = np.sqrt ( m2 )
#
#  Heuristic numbers:
#
#  If the axes may be badly scaled (which is to be avoided if
#  possible), then set SCBD = 10.  Otherwise set SCBD = 1.
#
#  If the problem is known to be ill-conditioned, initialize ILLC = true.
#
#  KTM is the number of iterations without improvement before the
#  algorithm terminates.  KTM = 4 is very cautious usually KTM = 1
#  is satisfactory.
#
  scbd = 1.0
  illc = False
  ktm = 1

  if ( illc ):
    ldfac = 0.1
  else:
    ldfac = 0.01

  kt = 0
  nl = 0
  nf = 1
  fx = f ( x, n )
  qf1 = fx
  t = small + abs ( t0 )
  t2 = t
  dmin = small
  h = h0
  h = max ( h, 100.0 * t )
  ldt = h
#
#  The initial set of search directions V is the identity matrix.
#
  v = np.zeros ( [ n, n ] )
  for i in range ( 0, n ):
    v[i,i] = 1.0

  d = np.zeros ( n )
  y = np.zeros ( n )
  z = np.zeros ( n )
  qa = 0.0
  qb = 0.0
  qc = 0.0
  qd0 = 0.0
  qd1 = 0.0
  q0 = x.copy ( )
  q1 = x.copy ( )

  if ( 0 < prin ):
    print2 ( n, x, prin, fx, nf, nl )
#
#  The main loop starts here.
#
  while ( True ):

    sf = d[0]
    d[0] = 0.0
#
#  Minimize along the first direction V(*,1).
#
    jsearch = 0
    nits = 2
    d2 = d[0]
    s = 0.0
    value = fx
    fk = False

    d2, s, value, x, nl, nf, fx, qa, qb, qc = minny ( n, jsearch, nits, \
      d2, s, value, fk, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, \
      fx, qa, qb, qc, qd0, qd1 )

    d[0] = d2

    if ( s <= 0.0 ):
      for i1 in range ( 0, n ):
        v[i1,0] = - v[i1,0]

    if ( sf <= 0.9 * d[0] or d[0] <= 0.9 * sf ):
      d[1:n] = 0.0
#
#  The inner loop starts here.
#
    for k in range ( 2, n + 1 ):

      y = x.copy ( )

      sf = fx

      if ( 0 < kt ):
        illc = True

      while ( True ):

        kl = k
        df = 0.0
#
#  A random step follows, to avoid resolution valleys.
#
        if ( illc ):

          for j in range ( 0, n ):
            r = np.random.rand ( 1 )
            s = ( 0.1 * ldt + t2 * 10.0 ** kt ) * ( r - 0.5 )
            z[j] = s
            x[0:n] = x[0:n] + s * v[0:n,j]

          fx = f ( x, n )
          nf = nf + 1
#
#  Minimize along the "non-conjugate" directions V(*,K),...,V(*,N).
#
        for k2 in range ( k, n + 1 ):

          sl = fx

          jsearch = k2 - 1
          nits = 2
          d2 = d[k2-1]
          s = 0.0
          value = fx
          fk = False

          d2, s, value, x, nl, nf, fx, qa, qb, qc = minny ( n, jsearch, nits, \
            d2, s, value, fk, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, \
            fx, qa, qb, qc, qd0, qd1 )

          d[k2-1] = d2

          if ( illc ):
            s = d[k2-1] * ( ( s + z[k2-1] ) ** 2 )
          else:
            s = sl - fx

          if ( df <= s ):
            df = s
            kl = k2
#
#  If there was not much improvement on the first try, set
#  ILLC = true and start the inner loop again.
#
        if ( illc ):
          break

        if ( abs ( 100.0 * machep * fx ) <= df ):
          break

        illc = True

      if ( k == 2 and 1 < prin ):
        r8vec_print ( n, d, '  The second difference array:' )
#
#  Minimize along the "conjugate" directions V(*,1),...,V(*,K-1).
#
      for k2 in range ( 1, k ):

        jsearch = k2 - 1
        nits = 2
        d2 = d[k2-1]
        s = 0.0
        value = fx
        fk = False

        d2, s, value, x, nl, nf, fx, qa, qb, qc = minny ( n, jsearch, nits, \
          d2, s, value, fk, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, \
          fx, qa, qb, qc, qd0, qd1 )

        d[k2-1] = d2

      f1 = fx
      fx = sf

      for i in range ( 0, n ):
        temp = x[i]
        x[i] = y[i]
        y[i] = temp - y[i]

      lds = np.linalg.norm ( y )
#
#  Discard direction V(*,kl).
#
#  If no random step was taken, V(*,KL) is the "non-conjugate"
#  direction along which the greatest improvement was made.
#
      if ( small < lds ):

        for j in range ( kl - 1, k - 1, -1 ):
          v[0:n,j] = v[0:n,j-1]
          d[j] = d[j-1]

        d[k-1] = 0.0

        v[0:n,k-1] = y[0:n] / lds
#
#  Minimize along the new "conjugate" direction V(*,k), which is
#  the normalized vector:  (new x) - (old x).
#
        jsearch = k - 1
        nits = 4
        d2 = d[k-1]
        value = f1
        fk = True

        d2, lds, value, x, nl, nf, fx, qa, qb, qc = minny ( n, jsearch, nits, \
          d2, lds, value, fk, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, \
          fx, qa, qb, qc, qd0, qd1 )

        d[k-1] = d2

        if ( lds <= 0.0 ):
          lds = - lds
          v[0:n,k-1] = - v[0:n,k-1]

      ldt = ldfac * ldt
      ldt = max ( ldt, lds )

      if ( 0 < prin ):
        print2 ( n, x, prin, fx, nf, nl )

      t2 = m2 * np.linalg.norm ( x ) + t
#
#  See whether the length of the step taken since starting the
#  inner loop exceeds half the tolerance.
#
      if ( 0.5 * t2 < ldt ):
        kt = - 1

      kt = kt + 1

      if ( ktm < kt ):

        if ( 0 < prin ):
          r8vec_print ( n, x, '  X:' )

        value = fx

        return value, x
#
#  The inner loop ends here.
#
#  Try quadratic extrapolation in case we are in a curved valley.
#
    x, q0, q1, nl, nf, fx, qf1, qa, qb, qc, qd0, qd1 = quad ( \
      n, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, fx, qf1, qa, qb, qc, qd0, qd1 )

    for j in range ( 0, n ):
      d[j] = 1.0 / np.sqrt ( d[j] )

    dn = max ( d )

    if ( 3 < prin ):
      r8mat_print ( n, n, v, '  The new direction vectors:' )

    for j in range ( 0, n ):
      v[0:n,j] = ( d[j] / dn ) * v[0:n,j]
#
#  Scale the axes to try to reduce the condition number.
#
    if ( 1.0 < scbd ):

      for i in range ( 0, n ):
        s = 0.0
        for j in range ( 0, n ):
          s = s + v[i,j] ** 2
        s = np.sqrt ( s )
        z[i] = max ( m4, s )

      s = min ( z )

      for i in range ( 0, n ):

        sl = s / z[i]
        z[i] = 1.0 / sl

        if ( scbd < z[i] ):
          sl = 1.0 / scbd
          z[i] = scbd

        v[i,0:n] = sl * v[i,0:n]
#
#  Calculate a new set of orthogonal directions before repeating
#  the main loop.
#
#  Transpose V for MINFIT:
#
    v = np.transpose ( v )
#
#  Call MINFIT to find the singular value decomposition of V.
#
#  This gives the principal values and principal directions of the
#  approximating quadratic form without squaring the condition number.
#
    v, d = minfit ( n, vsmall, v )
#
#  Unscale the axes.
#
    if ( 1.0 < scbd ):

      for i in range ( 0, n ):
        v[i,0:n] = z[i] * v[i,0:n]

      for j in range ( 0, n ):

        s = 0.0
        for i1 in range ( 0, n ):
          s = x + v[i1,j] ** 2
        s = sqrt ( s )

        d[j] = s * d[j]
        v[0:n,j] = v[0:n,j] / s

    for i in range ( 0, n ):

      dni = dn * d[i]

      if ( large < dni ):
        d[i] = vsmall
      elif ( dni < small ):
        d[i] = vlarge
      else:
        d[i] = 1.0 / dni ** 2
#
#  Sort the singular values and singular vectors.
#
    d, v = svsort ( n, d, v )
#
#  Determine the smallest eigenvalue.
#
    dmin = max ( d[n-1], small )
#
#  The ratio of the smallest to largest eigenvalue determines whether
#  the system is ill conditioned.
#
    if ( dmin < m2 * d[0] ):
      illc = True
    else:
      illc = False

    if ( 1 < prin ):

      if ( 1.0 < scbd ):
        r8vec_print ( n, z, '  The scale factors:' )

      r8vec_print ( n, d, '  Principal values of the quadratic form:' )

    if ( 3 < prin ):
      r8mat_print ( n, n, v, '  The principal axes:' )
#
#  The main loop ends here.
#
  if ( 0 < prin ):
    r8vec_print ( n, x, '  X:' )

  value = fx

  return value, x

def print2 ( n, x, prin, fx, nf, nl ):

#*****************************************************************************80
#
## PRINT2 prints certain data about the progress of the iteration.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#  Parameters:
#
#    Input, integer N, the number of variables.
#
#    Input, real X(N), the current estimate of the minimizer.
#
#    Input, integer PRIN, the user-specifed print level.
#    0, nothing is printed.
#    1, F is printed after every n+1 or n+2 linear minimizations.
#       final X is printed, but intermediate X is printed only
#       if N is at most 4.
#    2, the scale factors and the principal values of the approximating
#       quadratic form are also printed.
#    3, X is also printed after every few linear minimizations.
#    4, the principal vectors of the approximating quadratic form are
#       also printed.
#
#    Input, real FX, the smallest value of F(X) found so far.
#
#    Input, integer NF, the number of function evaluations.
#
#    Input, integer NL, the number of linear searches.
#
  print ( '' )
  print ( '  Linear searches      %d' % ( nl ) )
  print ( '  Function evaluations %d' % ( nf ) )
  print ( '  The function value FX = %g' % ( fx ) )

  if ( n <= 4 or 2 < prin ):
    r8vec_print ( n, x, '  X:' )

  return

def quad ( n, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, fx, qf1, qa, qb, \
  qc, qd0, qd1 ):

#*****************************************************************************80
#
## QUAD seeks to minimize the scalar function F along a particular curve.
#
#  Discussion:
#
#    The minimizer to be sought is required to lie on a curve defined
#    by Q0, Q1 and X.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#  Parameters:
#
#    Input, integer N, the number of variables.
#
#    Input, real F ( x, n ), the function.
#
#    Input/output, real X(N), ?
#
#    Input, real T, ?
#
#    Input, rea H, ?
#
#    Input, real V(N,N), the matrix of search directions.
#
#    Input/output, real Q0(N), an auxiliary point used to define
#    a curve through X.
#
#    Input/output, real Q1(N), an auxiliary point used to define
#    a curve through X.
#
#    Input/output, integer NL, the number of linear searches.
#
#    Input/output, integer NF, the number of function evaluations.
#
#    Input, real DMIN, an estimate for the smallest eigenvalue.
#
#    Input, real LDT, the length of the step.
#
#    Input/output, real FX, the value of F(X,N).
#
#    Input/output, real QF1, QA, QB, QC, QD0, QD1, ?
#
  import numpy as np

  temp = fx
  fx   = qf1
  qf1  = temp

  temp = x
  x = q1
  q1 = temp

  qd1 = np.linalg.norm ( x - q1 )

  l = qd1
  s = 0.0

  if ( qd0 <= 0.0 or qd1 <= 0.0 or nl < 3 * n * n ):

    fx = qf1
    qa = 0.0
    qb = 0.0
    qc = 1.0

  else:

    jsearch = -1
    nits = 2
    value = qf1
    fk = True

    s, l, value, x, nl, nf, fx, qa, qb, qc = minny ( n, jsearch, nits, \
      s, l, value, fk, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, \
      fx, qa, qb, qc, qd0, qd1 )

    qa =                 l * ( l - qd1 )       / ( qd0 + qd1 ) / qd0
    qb = - ( l + qd0 )     * ( l - qd1 ) / qd1                 / qd0
    qc =   ( l + qd0 ) * l               / qd1 / ( qd0 + qd1 )

  qd0 = qd1

  xnew = np.zeros ( n )
  xnew[0:n] = qa * q0[0:n] + qb * x[0:n] + qc * q1[0:n]

  q0[0:n] = x[0:n]
  x[0:n] = xnew[0:n]

  return x, q0, q1, nl, nf, fx, qf1, qa, qb, qc, qd0, qd1

def svsort ( n, d, v ):

#*****************************************************************************80
#
## SVSORT descending sorts singular values D and adjusts V.
#
#  Discussion:
#
#    A simple bubble sort is used on D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#  Parameters:
#
#    Input, integer N, the length of D, and the order of V.
#
#    Input/output, real D(N), the vector to be sorted.
#    On output, the entries of D are in descending order.
#
#    Input/output, real V(N,N), an N by N array to be adjusted
#    as D is sorted.  In particular, if the value that was in D(I) on input is
#    moved to D(J) on output, then the input column V(*,I) is moved to
#    the output column V(*,J).
#
  for j in range ( 0, n - 1 ):

    j3 = j
    for j2 in range ( j + 1, n ):
      if ( d[j3] < d[j2] ):
        j3 = j2

    t     = d[j]
    d[j]  = d[j3]
    d[j3] = t

    for i in range ( 0, n ):
      t       = v[i,j]
      v[i,j]  = v[i,j3]
      v[i,j3] = t

  return d, v

def svsort_test ( ):

#*****************************************************************************80
#
## SVSORT_TEST tests SVSORT, which sorts singular value information.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 5

  print ( '' )
  print ( 'SVSORT_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  SVSORT sorts a vector D, and the corresponding columns' )
  print ( '  of a matrix V.' )

  d = np.random.rand ( n )

  v = np.zeros ( [ n, n ] )

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      v[i,j] = 10 * ( i + 1 ) + ( j + 1 )

  print ( '' )
  print ( '  First row = entries of D.' )
  print ( '  Corresponding columns of V below.' )
  print ( '' )
  for j in range ( 0, n ):
    print ( '%14.6g' % ( d[j] ) ),
  print ( '' )
  print ( '' )
  for i in range ( 0, n ):
    for j in range ( 0, n ):
      print ( '%14.6g' % ( v[i,j] ) ),
    print ( '' )

  d, v = svsort ( n, d, v )

  print ( '' )
  print ( '  After sorting D and rearranging V:' )
  print ( '' )
  for j in range ( 0, n ):
    print ( '%14.6g' % ( d[j] ) ),
  print ( '' )
  print ( '' )
  for i in range ( 0, n ):
    for j in range ( 0, n ):
      print ( '%14.6g' % ( v[i,j] ) ),
    print ( '' )
#
#  Terminate.
#
  print ( '' )
  print ( 'SVSORT_TEST:' )
  print ( '  Normal end of execution.' )
  return

def praxis_test ( ):

#*****************************************************************************80
#
## PRAXIS_TEST tests the PRAXIS library.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 August 2016
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'PRAXIS_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Test the PRAXIS library.' )
#
#  Minimization tests.
#
  beale_test ( )
  box_test ( )
  chebyquad_test ( )
  cube_test ( )
  helix_test ( )
  hilbert_test ( )
  powell3d_test ( )
  rosenbrock_test ( )
  singular_test ( )
  tridiagonal_test ( )
  watson_test ( )
  wood_test ( )
#
#  Utility tests.
#
  minfit_test ( )
  svsort_test ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'PRAXIS_TEST:' )
  print ( '  Normal end of execution.' )
  return

def r8_epsilon ( ):

#*****************************************************************************80
#
## R8_EPSILON returns the R8 roundoff unit.
#
#  Discussion:
#
#    The roundoff unit is a number R which is a power of 2 with the 
#    property that, to the precision of the computer's arithmetic,
#      1 < 1 + R
#    but 
#      1 = ( 1 + R / 2 )
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 June 2013
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Output, real VALUE, the roundoff unit.
#
  value = 2.220446049250313E-016

  return value

def r8_epsilon_test ( ):

#*****************************************************************************80
#
## R8_EPSILON_TEST tests R8_EPSILON.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 September 2012
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'R8_EPSILON_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  R8_EPSILON produces the R8 roundoff unit.' )
  print ( '' )

  r = r8_epsilon ( )
  print ( '  R = R8_EPSILON()         = %e' % ( r ) )

  s = ( 1.0 + r ) - 1.0
  print ( '  ( 1 + R ) - 1            = %e' % ( s ) )

  s = ( 1.0 + ( r / 2.0 ) ) - 1.0
  print ( '  ( 1 + (R/2) ) - 1        = %e' % ( s ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'R8_EPSILON_TEST' )
  print ( '  Normal end of execution.' )
  return

def r8_hypot ( x, y ):

#*****************************************************************************80
#
## R8_HYPOT returns the value of sqrt ( X^2 + Y^2 ).
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    13 January 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real X, Y, the arguments.
#
#    Output, real VALUE, the value of sqrt ( X^2 + Y^2 ).
#
  import numpy as np

  if ( abs ( x ) < abs ( y ) ):
    a = abs ( y )
    b = abs ( x )
  else:
    a = abs ( x )
    b = abs ( y )
#
#  A contains the larger value.
#
  if ( a == 0.0 ):
    value = 0.0
  else:
    value = a * np.sqrt ( 1.0 + ( b / a ) ** 2 )

  return value

def r8_hypot_test ( ):

#*****************************************************************************80
#
## R8_HYPOT_TEST tests R8_HYPOT.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    13 January 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'R8_HYPOT_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  R8_HYPOT returns an accurate value for sqrt(A^2+B^2).' )
  print ( '' )
  print ( '             A          B          R8_HYPOT      sqrt(A^2+B^2)' )
  print ( '' )

  b = 2.0

  for i in range ( 0, 20 ):
    a = 1.0
    b = b / 2.0
    c = r8_hypot ( a, b )
    d = np.sqrt ( a ** 2 + b ** 2 )
    
    print ( '  %12g  %12g  %24.16g  %24.16g' % ( a, b, c, d ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'R8_HYPOT_TEST' )
  print ( '  Normal end of execution.' )
  return

def r8mat_print ( m, n, a, title ):

#*****************************************************************************80
#
## R8MAT_PRINT prints an R8MAT.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    31 August 2014
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer M, the number of rows in A.
#
#    Input, integer N, the number of columns in A.
#
#    Input, real A(M,N), the matrix.
#
#    Input, string TITLE, a title.
#
  r8mat_print_some ( m, n, a, 0, 0, m - 1, n - 1, title )

  return

def r8mat_print_test ( ):

#*****************************************************************************80
#
## R8MAT_PRINT_TEST tests R8MAT_PRINT.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    10 February 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'R8MAT_PRINT_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  R8MAT_PRINT prints an R8MAT.' )

  m = 4
  n = 6
  v = np.array ( [ \
    [ 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 ], 
    [ 21.0, 22.0, 23.0, 24.0, 25.0, 26.0 ], 
    [ 31.0, 32.0, 33.0, 34.0, 35.0, 36.0 ], 
    [ 41.0, 42.0, 43.0, 44.0, 45.0, 46.0 ] ], dtype = np.float64 )
  r8mat_print ( m, n, v, '  Here is an R8MAT:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'R8MAT_PRINT_TEST:' )
  print ( '  Normal end of execution.' )
  return

def r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title ):

#*****************************************************************************80
#
## R8MAT_PRINT_SOME prints out a portion of an R8MAT.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    10 February 2015
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer M, N, the number of rows and columns of the matrix.
#
#    Input, real A(M,N), an M by N matrix to be printed.
#
#    Input, integer ILO, JLO, the first row and column to print.
#
#    Input, integer IHI, JHI, the last row and column to print.
#
#    Input, string TITLE, a title.
#
  incx = 5

  print ( '' )
  print ( title )

  if ( m <= 0 or n <= 0 ):
    print ( '' )
    print ( '  (None)' )
    return

  for j2lo in range ( max ( jlo, 0 ), min ( jhi + 1, n ), incx ):

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )
    
    print ( '' )
    print ( '  Col: ', end = '' )

    for j in range ( j2lo, j2hi + 1 ):
      print ( '%7d       ' % ( j ), end = '' )

    print ( '' )
    print ( '  Row' )

    i2lo = max ( ilo, 0 )
    i2hi = min ( ihi, m )

    for i in range ( i2lo, i2hi + 1 ):

      print ( '%7d :' % ( i ), end = '' )
      
      for j in range ( j2lo, j2hi + 1 ):
        print ( '%12g  ' % ( a[i,j] ), end = '' )

      print ( '' )

  return

def r8mat_print_some_test ( ):

#*****************************************************************************80
#
## R8MAT_PRINT_SOME_TEST tests R8MAT_PRINT_SOME.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    31 October 2014
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'R8MAT_PRINT_SOME_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  R8MAT_PRINT_SOME prints some of an R8MAT.' )

  m = 4
  n = 6
  v = np.array ( [ \
    [ 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 ], 
    [ 21.0, 22.0, 23.0, 24.0, 25.0, 26.0 ], 
    [ 31.0, 32.0, 33.0, 34.0, 35.0, 36.0 ], 
    [ 41.0, 42.0, 43.0, 44.0, 45.0, 46.0 ] ], dtype = np.float64 )
  r8mat_print_some ( m, n, v, 0, 3, 2, 5, '  Here is an R8MAT:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'R8MAT_PRINT_SOME_TEST:' )
  print ( '  Normal end of execution.' )
  return

def r8vec_print ( n, a, title ):

#*****************************************************************************80
#
## R8VEC_PRINT prints an R8VEC.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    31 August 2014
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer N, the dimension of the vector.
#
#    Input, real A(N), the vector to be printed.
#
#    Input, string TITLE, a title.
#
  print ( '' )
  print ( title )
  print ( '' )
  for i in range ( 0, n ):
    print ( '%6d:  %12g' % ( i, a[i] ) )

def r8vec_print_test ( ):

#*****************************************************************************80
#
## R8VEC_PRINT_TEST tests R8VEC_PRINT.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 October 2014
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'R8VEC_PRINT_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  R8VEC_PRINT prints an R8VEC.' )

  n = 4
  v = np.array ( [ 123.456, 0.000005, -1.0E+06, 3.14159265 ], dtype = np.float64 )
  r8vec_print ( n, v, '  Here is an R8VEC:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'R8VEC_PRINT_TEST:' )
  print ( '  Normal end of execution.' )
  return

def rosenbrock_f ( x, n ):

#*****************************************************************************80
#
## ROSENBROCK_F evaluates the Rosenbrock function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  value = 0.0

  for j in range ( 0, n ):
    if ( ( j % 2 ) == 0 ):
      value = value + ( 1.0 - x[j] ) ** 2
    else:
      value = value + 100.0 * ( x[j] - x[j-1] ** 2 ) ** 2

  return value

def rosenbrock_test ( ):

#*****************************************************************************80
#
## ROSENBROCK_TEST calls PRAXIS for the Rosenbrock function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 2

  print ( '' )
  print ( 'ROSENBROCK_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Rosenbrock function.' )

  t0 = 0.00001
  h0 = 1.0
  prin = 0

  x = np.array ( [ -1.2, 1.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( rosenbrock_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, rosenbrock_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( rosenbrock_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'ROSENBROCK_TEST:' )
  print ( '  Normal end of execution.' )
  return

def singular_f ( x, n ):

#*****************************************************************************80
#
## SINGULAR_F evaluates the Powell Singular function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  value = 0.0

  for j in range ( 0, n, 4 ):

    if ( j + 1 <= n - 1 ):
      xjp1 = x[j+1]
    else:
      xjp1 = 0.0

    if ( j + 2 <= n - 1 ):
      xjp2 = x[j+2]
    else:
      xjp2 = 0.0

    if ( j + 3 <= n - 1 ):
      xjp3 = x[j+3]
    else:
      xjp3 = 0.0

    f1 = x[j] + 10.0 * xjp1

    if ( j + 1 <= n - 1 ):
      f2 = xjp2 - xjp3
    else:
      f2 = 0.0

    if ( j + 2 <= n - 1 ):
      f3 = xjp1 - 2.0 * xjp2
    else:
      f3 = 0.0

    if ( j + 3 <= n - 1 ):
      f4 = x[j] - xjp3
    else:
      f4 = 0.0

    value = value \
      +        f1 ** 2 \
      +  5.0 * f2 ** 2 \
      +        f3 ** 4 \
      + 10.0 * f4 ** 4

  return value

def singular_test ( ):

#*****************************************************************************80
#
## SINGULAR_TEST calls PRAXIS for the Powell Singular function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 4

  print ( '' )
  print ( 'SINGULAR_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Powell Singular function.' )

  t0 = 0.00001
  h0 = 1.0
  prin = 0

  x = np.array ( [ 3.0, -1.0, 0.0, 1.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( singular_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, singular_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( singular_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'SINGULAR_TEST:' )
  print ( '  Normal end of execution.' )
  return

def timestamp ( ):

#*****************************************************************************80
#
## TIMESTAMP prints the date as a timestamp.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    06 April 2013
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    None
#
  import time

  t = time.time ( )
  print ( time.ctime ( t ) )

  return None

def timestamp_test ( ):

#*****************************************************************************80
#
## TIMESTAMP_TEST tests TIMESTAMP.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    03 December 2014
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    None
#
  import platform

  print ( '' )
  print ( 'TIMESTAMP_TEST:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  TIMESTAMP prints a timestamp of the current date and time.' )
  print ( '' )

  timestamp ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'TIMESTAMP_TEST:' )
  print ( '  Normal end of execution.' )
  return

def tridiagonal_f ( x, n ):

#*****************************************************************************80
#
## TRIDIAGONAL_F evaluates the tridiagonal function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  value = x[0] ** 2 + 2.0 * sum ( x[1:n] ** 2 )

  for i in range ( 0, n - 1 ):
    value = value - 2.0 * x[i] * x[i+1]

  value = value - 2.0 * x[0]

  return value

def tridiagonal_test ( ):

#*****************************************************************************80
#
## TRIDIAGONAL_TEST calls PRAXIS for the Tridiagonal function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 4

  print ( '' )
  print ( 'TRIDIAGONAL_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Tridiagonal function.' )

  t0 = 0.00001
  h0 = 8.0
  prin = 0

  x = np.zeros ( n )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( tridiagonal_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, tridiagonal_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( tridiagonal_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'TRIDIAGONAL_TEST:' )
  print ( '  Normal end of execution.' )
  return

def watson_f ( x, n ):

#*****************************************************************************80
#
## WATSON_F evaluates the Watson function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  value = 0.0

  for i in range ( 0, 29 ):

    s1 = 0.0
    d = 1.0
    for j in range ( 1, n ):
      s1 = s1 + j * d * x[j]
      d = d * ( i + 1 ) / 29.0

    s2 = 0.0
    d = 1.0
    for j in range ( 0, n ):
      s2 = s2 + d * x[j]
      d = d * ( i + 1 ) / 29.0

    value = value + ( s1 - s2 * s2 - 1.0 ) ** 2

  value = value + x[0] ** 2 + ( x[1] - x[0] ** 2 - 1.0 ) ** 2
 
  return value

def watson_test ( ):

#*****************************************************************************80
#
## WATSON_TEST calls PRAXIS for the Watson function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 6

  print ( '' )
  print ( 'WATSON_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Watson function.' )

  t0 = 0.00001
  h0 = 1.0
  prin = 0

  x = np.zeros ( n )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( watson_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, watson_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( watson_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'WATSON_TEST:' )
  print ( '  Normal end of execution.' )
  return

def wood_f ( x, n ):

#*****************************************************************************80
#
## WOOD_F evaluates the Wood function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  f1 = x[1] - x[0] ** 2
  f2 = 1.0 - x[0]
  f3 = x[3] - x[2] ** 2
  f4 = 1.0 - x[2]
  f5 = x[1] + x[3] - 2.0
  f6 = x[1] - x[3]

  value = \
      100.0 * f1 ** 2 \
    +         f2 ** 2 \
    +  90.0 * f3 ** 2 \
    +         f4 ** 2 \
    +  10.0 * f5 ** 2 \
    +   0.1 * f6 ** 2

  return value

def wood_test ( ):

#*****************************************************************************80
#
## WOOD_TEST calls PRAXIS for the Wood function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 4

  print ( '' )
  print ( 'WOOD_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Wood function.' )

  t0 = 0.00001
  h0 = 10.0
  prin = 0

  x = np.array ( [ -3.0, -1.0, -3.0, -1.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( wood_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, wood_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( wood_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'WOOD_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  timestamp ( )
  praxis_test ( )
  timestamp ( )
 
