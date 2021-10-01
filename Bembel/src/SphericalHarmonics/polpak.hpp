#ifndef POLPAK_H_
#define POLPAK_H_

# include <cmath>
# include <complex>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

//# include "polpak.hpp"

namespace POLPAK{


    double r8_factorial ( int n )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    R8_FACTORIAL returns the factorial function as an R8.
    //
    //  Discussion:
    //
    //    Factorial ( N ) = Product ( 1 <= I <= N ) I
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    10 May 2003
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    Input, int N, the argument of the factorial function.
    //    0 <= N.
    //
    //    Output, double R8_FACTORIAL, the factorial of N.
    //
    {
      double fact;
      int i;
    //
    //  Check.
    //
      if ( n < 0 )
      {
        cerr << "\n";
        cerr << "R8_FACTORIAL - Fatal error!\n";
        cerr << "  N < 0.\n";
        exit ( 1 );
      }

      fact = 1.0;

      for ( i = 2; i <= n; i++ )
      {
        fact = fact * ( double ) i;
      }

      return fact;
    }

  //****************************************************************************80

  void legendre_associated ( int n, int m, double x, double cx[] )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    LEGENDRE_ASSOCIATED evaluates the associated Legendre functions.
  //
  //  Differential equation:
  //
  //    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
  //
  //  First terms:
  //
  //    M = 0  ( = Legendre polynomials of first kind P(N,X) )
  //
  //    P00 =    1
  //    P10 =    1 X
  //    P20 = (  3 X^2 -   1)/2
  //    P30 = (  5 X^3 -   3 X)/2
  //    P40 = ( 35 X^4 -  30 X^2 +   3)/8
  //    P50 = ( 63 X^5 -  70 X^3 +  15 X)/8
  //    P60 = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
  //    P70 = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
  //
  //    M = 1
  //
  //    P01 =   0
  //    P11 =   1 * SQRT(1-X*X)
  //    P21 =   3 * SQRT(1-X*X) * X
  //    P31 = 1.5 * SQRT(1-X*X) * (5*X*X-1)
  //    P41 = 2.5 * SQRT(1-X*X) * (7*X*X*X-3*X)
  //
  //    M = 2
  //
  //    P02 =   0
  //    P12 =   0
  //    P22 =   3 * (1-X*X)
  //    P32 =  15 * (1-X*X) * X
  //    P42 = 7.5 * (1-X*X) * (7*X*X-1)
  //
  //    M = 3
  //
  //    P03 =   0
  //    P13 =   0
  //    P23 =   0
  //    P33 =  15 * (1-X*X)^1.5
  //    P43 = 105 * (1-X*X)^1.5 * X
  //
  //    M = 4
  //
  //    P04 =   0
  //    P14 =   0
  //    P24 =   0
  //    P34 =   0
  //    P44 = 105 * (1-X*X)^2
  //
  //  Recursion:
  //
  //    if N < M:
  //      P(N,M) = 0
  //    if N = M:
  //      P(N,M) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
  //      all the odd integers less than or equal to N.
  //    if N = M+1:
  //      P(N,M) = X*(2*M+1)*P(M,M)
  //    if M+1 < N:
  //      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
  //
  //  Special values:
  //
  //    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
  //    function of the first kind equals the Legendre polynomial of the
  //    first kind.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    07 March 2005
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Milton Abramowitz, Irene Stegun,
  //    Handbook of Mathematical Functions,
  //    National Bureau of Standards, 1964,
  //    ISBN: 0-486-61272-4,
  //    LC: QA47.A34.
  //
  //  Parameters:
  //
  //    Input, int N, the maximum first index of the Legendre
  //    function, which must be at least 0.
  //
  //    Input, int M, the second index of the Legendre function,
  //    which must be at least 0, and no greater than N.
  //
  //    Input, double X, the point at which the function is to be
  //    evaluated.  X must satisfy -1 <= X <= 1.
  //
  //    Output, double CX[N+1], the values of the first N+1 function.
  //
  {
    double factor;
    int i;
    double somx2;

    if ( m < 0 )
    {
      cerr << "\n";
      cerr << "LEGENDRE_ASSOCIATED - Fatal error!\n";
      cerr << "  Input value of M is " << m << "\n";
      cerr << "  but M must be nonnegative.\n";
      exit ( 1 );
    }

    if ( n < m )
    {
      cerr << "\n";
      cerr << "LEGENDRE_ASSOCIATED - Fatal error!\n";
      cerr << "  Input value of M = " << m << "\n";
      cerr << "  Input value of N = " << n << "\n";
      cerr << "  but M must be less than or equal to N.\n";
      exit ( 1 );
    }

    if ( x < -1.0 )
    {
      cerr << "\n";
      cerr << "LEGENDRE_ASSOCIATED - Fatal error!\n";
      cerr << "  Input value of X = " << x << "\n";
      cerr << "  but X must be no less than -1.\n";
      exit ( 1 );
    }

    if ( 1.0 < x )
    {
      cerr << "\n";
      cerr << "LEGENDRE_ASSOCIATED - Fatal error!\n";
      cerr << "  Input value of X = " << x << "\n";
      cerr << "  but X must be no more than 1.\n";
      exit ( 1 );
    }

    for ( i = 0; i <= m-1; i++ )
    {
      cx[i] = 0.0;
    }
    cx[m] = 1.0;

    somx2 = sqrt ( 1.0 - x * x );

    factor = 1.0;
    for ( i = 1; i <= m; i++ )
    {
      cx[m] = -cx[m] * factor * somx2;
      factor = factor + 2.0;
    }

    if ( m == n )
    {
      return;
    }

    cx[m+1] = x * ( double ) ( 2 * m + 1 ) * cx[m];

    for ( i = m+2; i <= n; i++ )
    {
      cx[i] = ( ( double ) ( 2 * i     - 1 ) * x * cx[i-1]
              + ( double ) (   - i - m + 1 )     * cx[i-2] )
              / ( double ) (     i - m     );
    }

    return;
  }

  void legendre_associated_normalized ( int n, int m, double x, double cx[] )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    LEGENDRE_ASSOCIATED_NORMALIZED: normalized associated Legendre functions.
  //
  //  Discussion:
  //
  //    The unnormalized associated Legendre functions P_N^M(X) have
  //    the property that
  //
  //      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX
  //      = 2 * ( N + M )! / ( ( 2 * N + 1 ) * ( N - M )! )
  //
  //    By dividing the function by the square root of this term,
  //    the normalized associated Legendre functions have norm 1.
  //
  //    However, we plan to use these functions to build spherical
  //    harmonics, so we use a slightly different normalization factor of
  //
  //      sqrt ( ( ( 2 * N + 1 ) * ( N - M )! ) / ( 4 * pi * ( N + M )! ) )
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    07 March 2005
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Milton Abramowitz, Irene Stegun,
  //    Handbook of Mathematical Functions,
  //    National Bureau of Standards, 1964,
  //    ISBN: 0-486-61272-4,
  //    LC: QA47.A34.
  //
  //  Parameters:
  //
  //    Input, int N, the maximum first index of the Legendre
  //    function, which must be at least 0.
  //
  //    Input, int M, the second index of the Legendre function,
  //    which must be at least 0, and no greater than N.
  //
  //    Input, double X, the point at which the function is to be
  //    evaluated.  X must satisfy -1 <= X <= 1.
  //
  //    Output, double CX[N+1], the values of the first N+1 function.
  //
  {
    double factor;
    int i;
    int mm;
    const double r8_pi = 3.141592653589793;
    double somx2;

    if ( m < 0 )
    {
      cerr << "\n";
      cerr << "LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!\n";
      cerr << "  Input value of M is " << m << "\n";
      cerr << "  but M must be nonnegative.\n";
      exit ( 1 );
    }

    if ( n < m )
    {
      cerr << "\n";
      cerr << "LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!\n";
      cerr << "  Input value of M = " << m << "\n";
      cerr << "  Input value of N = " << n << "\n";
      cerr << "  but M must be less than or equal to N.\n";
      exit ( 1 );
    }

    if ( x < -1.0 )
    {
      cerr << "\n";
      cerr << "LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!\n";
      cerr << "  Input value of X = " << x << "\n";
      cerr << "  but X must be no less than -1.\n";
      exit ( 1 );
    }

    if ( 1.0 < x )
    {
      cerr << "\n";
      cerr << "LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!\n";
      cerr << "  Input value of X = " << x << "\n";
      cerr << "  but X must be no more than 1.\n";
      exit ( 1 );
    }

    for ( i = 0; i <= m-1; i++ )
    {
      cx[i] = 0.0;
    }
    cx[m] = 1.0;

    somx2 = sqrt ( 1.0 - x * x );

    factor = 1.0;
    for ( i = 1; i <= m; i++ )
    {
      cx[m] = -cx[m] * factor * somx2;
      factor = factor + 2.0;
    }

    if ( m+1 <= n )
    {
      cx[m+1] = x * ( double ) ( 2 * m + 1 ) * cx[m];
    }

    for ( i = m+2; i <= n; i++ )
    {
      cx[i] = ( ( double ) ( 2 * i     - 1 ) * x * cx[i-1]
              + ( double ) (   - i - m + 1 )     * cx[i-2] )
              / ( double ) (     i - m     );
    }
  //
  //  Normalization.
  //
    for ( mm = m; mm <= n; mm++ )
    {
      factor = sqrt ( ( ( double ) ( 2 * mm + 1 ) * r8_factorial ( mm - m ) )
        / ( 4.0 * r8_pi * r8_factorial ( mm + m ) ) );
      cx[mm] = cx[mm] * factor;
    }

    return;
  }

  void spherical_harmonic ( int l, int m, double theta, double phi,
    double c[], double s[] )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    SPHERICAL_HARMONIC evaluates spherical harmonic functions.
  //
  //  Discussion:
  //
  //    The spherical harmonic function Y(L,M,THETA,PHI,X) is the
  //    angular part of the solution to Laplace's equation in spherical
  //    coordinates.
  //
  //    Y(L,M,THETA,PHI,X) is related to the associated Legendre
  //    function as follows:
  //
  //      Y(L,M,THETA,PHI,X) = FACTOR * P(L,M,cos(THETA)) * exp ( i * M * PHI )
  //
  //    Here, FACTOR is a normalization factor:
  //
  //      FACTOR = sqrt ( ( 2 * L + 1 ) * ( L - M )! / ( 4 * PI * ( L + M )! ) )
  //
  //    In Mathematica, a spherical harmonic function can be evaluated by
  //
  //      SphericalHarmonicY [ l, m, theta, phi ]
  //
  //    Note that notational tradition in physics requires that THETA
  //    and PHI represent the reverse of what they would normally mean
  //    in mathematical notation; that is, THETA goes up and down, and
  //    PHI goes around.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    07 March 2005
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Milton Abramowitz, Irene Stegun,
  //    Handbook of Mathematical Functions,
  //    National Bureau of Standards, 1964,
  //    ISBN: 0-486-61272-4,
  //    LC: QA47.A34.
  //
  //    Eric Weisstein,
  //    CRC Concise Encyclopedia of Mathematics,
  //    CRC Press, 2002,
  //    Second edition,
  //    ISBN: 1584883472,
  //    LC: QA5.W45.
  //
  //    Stephen Wolfram,
  //    The Mathematica Book,
  //    Fourth Edition,
  //    Cambridge University Press, 1999,
  //    ISBN: 0-521-64314-7,
  //    LC: QA76.95.W65.
  //
  //  Parameters:
  //
  //    Input, int L, the first index of the spherical harmonic function.
  //    Normally, 0 <= L.
  //
  //    Input, int M, the second index of the spherical harmonic function.
  //    Normally, -L <= M <= L.
  //
  //    Input, double THETA, the polar angle, for which
  //    0 <= THETA <= PI.
  //
  //    Input, double PHI, the longitudinal angle, for which
  //    0 <= PHI <= 2*PI.
  //
  //    Output, double C[L+1], S[L+1], the real and imaginary
  //    parts of the functions Y(L,0:L,THETA,PHI).
  //

  {
    double angle;
    int i;
    int m_abs;
    double *plm;

    m_abs = abs ( m );

    plm = new double[l+1];

    legendre_associated_normalized ( l, m_abs, cos ( theta ), plm );

    angle = ( double ) ( m ) * phi;

    if ( 0 <= m )
    {
      for ( i = 0; i <= l; i++ )
      {
        c[i] = plm[i] * cos ( angle );
        s[i] = plm[i] * sin ( angle );
      }
    }
    else
    {
      for ( i = 0; i <= l; i++ )
      {
        c[i] = -plm[i] * cos ( angle );
        s[i] = -plm[i] * sin ( angle );
      }
    }

    delete [] plm;

    return;
  }

}



#endif
