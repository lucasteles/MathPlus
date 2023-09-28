// by Freya Holmér (https://github.com/FreyaHolmer/MathPlus)

using System;
using System.Runtime.CompilerServices;
using System.Numerics;

namespace MathPlus;

public struct Polynomial4D : IPolynomialCubic<Polynomial4D, Vector4>, IParamCurve3Diff<Vector4> {

	const MethodImplOptions INLINE = MethodImplOptions.AggressiveInlining;

	/// <inheritdoc cref="Polynomial.NaN"/>
	public static readonly Polynomial4D NaN = new Polynomial4D { x = Polynomial.NaN, y = Polynomial.NaN, z = Polynomial.NaN, w = Polynomial.NaN };

	public Polynomial x, y, z, w;

	public Polynomial4D( Polynomial x, Polynomial y, Polynomial z, Polynomial w ) => ( this.X,this.Y,this.Z,this.w ) = ( x, y, z, w );

	/// <inheritdoc cref="Polynomial(float,float,float,float)"/>
	public Polynomial4D( Vector4 c0, Vector4 c1, Vector4 c2, Vector4 c3 ) {
		this.X = new Polynomial( c0.X,c1.X,c2.X,c3.X );
		this.Y = new Polynomial( c0.Y,c1.Y,c2.Y,c3.Y );
		this.Z = new Polynomial( c0.Z,c1.Z,c2.Z,c3.Z );
		this.w = new Polynomial( c0.w, c1.w, c2.w, c3.w );
	}

	/// <inheritdoc cref="Polynomial(float,float,float)"/>
	public Polynomial4D( Vector4 c0, Vector4 c1, Vector4 c2 ) {
		this.X = new Polynomial( c0.X,c1.X,c2.X,0 );
		this.Y = new Polynomial( c0.Y,c1.Y,c2.Y,0 );
		this.Z = new Polynomial( c0.Z,c1.Z,c2.Z,0 );
		this.w = new Polynomial( c0.w, c1.w, c2.w, 0 );
	}

	/// <inheritdoc cref="Polynomial(float,float)"/>
	public Polynomial4D( Vector4 c0, Vector4 c1 ) {
		this.X = new Polynomial( c0.X,c1.X,0, 0 );
		this.Y = new Polynomial( c0.Y,c1.Y,0, 0 );
		this.Z = new Polynomial( c0.Z,c1.Z,0, 0 );
		this.w = new Polynomial( c0.w, c1.w, 0, 0 );
	}

	/// <inheritdoc cref="Polynomial(Matrix4x1)"/>
	public Polynomial4D( Vector4Matrix4x1 coefficients ) => ( x, y, z, w ) = ( new Polynomial( coefficients.X ), new Polynomial( coefficients.Y ), new Polynomial( coefficients.Z ), new Polynomial( coefficients.W ) );

	/// <inheritdoc cref="Polynomial(Matrix4x1)"/>
	public Polynomial4D( Vector4Matrix3x1 coefficients ) => ( x, y, z, w ) = ( new Polynomial( coefficients.X ), new Polynomial( coefficients.Y ), new Polynomial( coefficients.Z ), new Polynomial( coefficients.W ) );

	#region IPolynomialCubic

	public Vector4 C0 {
		[MethodImpl( INLINE )] get => new(x.c0, y.c0, z.c0, w.c0);
		[MethodImpl( INLINE )] set => ( x.c0, y.c0, z.c0, w.c0 ) = ( value.X,value.Y,value.Z,value.w );
	}
	public Vector4 C1 {
		[MethodImpl( INLINE )] get => new(x.c1, y.c1, z.c1, w.c1);
		[MethodImpl( INLINE )] set => ( x.c1, y.c1, z.c1, w.c1 ) = ( value.X,value.Y,value.Z,value.w );
	}
	public Vector4 C2 {
		[MethodImpl( INLINE )] get => new(x.c2, y.c2, z.c2, w.c2);
		[MethodImpl( INLINE )] set => ( x.c2, y.c2, z.c2, w.c2 ) = ( value.X,value.Y,value.Z,value.w );
	}
	public Vector4 C3 {
		[MethodImpl( INLINE )] get => new(x.c3, y.c3, z.c3, w.c3);
		[MethodImpl( INLINE )] set => ( x.c3, y.c3, z.c3, w.c3 ) = ( value.X,value.Y,value.Z,value.w );
	}

	public Polynomial this[ int i ] {
		get => i switch { 0 => x, 1 => y, 2 => z, 4 => w, _ => throw new IndexOutOfRangeException( "Polynomial4D component index has to be either 0, 1, 2, or 3" ) };
		set => _ = i switch { 0 => x = value, 1 => y = value, 2 => z = value, 3 => w = value, _ => throw new IndexOutOfRangeException() };
	}

	[MethodImpl( INLINE )] public Vector4 GetCoefficient( int degree ) =>
		degree switch {
			0 => C0,
			1 => C1,
			2 => C2,
			3 => C3,
			_ => throw new IndexOutOfRangeException( "Polynomial coefficient degree/index has to be between 0 and 3" )
		};

	[MethodImpl( INLINE )] public void SetCoefficient( int degree, Vector4 value ) {
		_ = degree switch {
			0 => C0 = value,
			1 => C1 = value,
			2 => C2 = value,
			3 => C3 = value,
			_ => throw new IndexOutOfRangeException( "Polynomial coefficient degree/index has to be between 0 and 3" )
		};
	}

	public Vector4 Eval( float t ) {
		var t2 = t * t;
		var t3 = t2 * t;
		return new Vector4(
			x.c3 * t3 + x.c2 * t2 + x.c1 * t + x.c0,
			y.c3 * t3 + y.c2 * t2 + y.c1 * t + y.c0,
			z.c3 * t3 + z.c2 * t2 + z.c1 * t + z.c0,
			w.c3 * t3 + w.c2 * t2 + w.c1 * t + w.c0
		);
	}

	[MethodImpl( INLINE )] public Vector4 Eval( float t, int n ) => Differentiate( n ).Eval( t );

	[MethodImpl( INLINE )] public Polynomial4D Differentiate( int n = 1 ) => new(x.Differentiate( n ), y.Differentiate( n ), z.Differentiate( n ), w.Differentiate( n ));

	public Polynomial4D ScaleParameterSpace( float factor ) {
		// ReSharper disable once CompareOfFloatsByEqualityOperator
		if( factor == 1f )
			return this;
		var factor2 = factor * factor;
		var factor3 = factor2 * factor;
		return new Polynomial4D(
			new Polynomial( x.c0, x.c1 / factor, x.c2 / factor2, x.c3 / factor3 ),
			new Polynomial( y.c0, y.c1 / factor, y.c2 / factor2, y.c3 / factor3 ),
			new Polynomial( z.c0, z.c1 / factor, z.c2 / factor2, z.c3 / factor3 ),
			new Polynomial( w.c0, w.c1 / factor, w.c2 / factor2, w.c3 / factor3 )
		);
	}

	public Polynomial4D Compose( float g0, float g1 ) => new(x.Compose( g0, g1 ), y.Compose( g0, g1 ), z.Compose( g0, g1 ), w.Compose( g0, g1 ));

	#endregion

	public static Polynomial4D FitCubicFrom0( float x1, float x2, float x3, Vector4 y0, Vector4 y1, Vector4 y2, Vector4 y3 ) {
		// precalcs
		var i12 = x2 - x1;
		var i13 = x3 - x1;
		var i23 = x3 - x2;
		var x1x2 = x1 * x2;
		var x1x3 = x1 * x3;
		var x2x3 = x2 * x3;
		var x1x2x3 = x1 * x2x3;
		var x0plusx1plusx2 = x1 + x2;
		var x0plusx1plusx3 = x1 + x3;
		var x2plusx3 = x2 + x3;
		var x1plusx2plusx3 = x1 + x2plusx3;
		var x1x2plusx1x3plusx2x3 = ( x1x2 + x1x3 + x2x3 );

		// scale factors
		var scl0 = y0 / -( x1 * x2 * x3 );
		var scl1 = y1 / +( x1 * i12 * i13 );
		var scl2 = y2 / -( x2 * i12 * i23 );
		var scl3 = y3 / +( x3 * i13 * i23 );

		// polynomial form
		Vector4 c0 = new(
			-( scl0.X * x1x2x3 ),
			-( scl0.Y * x1x2x3 ),
			-( scl0.Z * x1x2x3 ),
			-( scl0.w * x1x2x3 )
		);
		Vector4 c1 = new(
			scl0.X * x1x2plusx1x3plusx2x3 + scl1.X * x2x3 + scl2.X * x1x3 + scl3.X * x1x2,
			scl0.Y * x1x2plusx1x3plusx2x3 + scl1.Y * x2x3 + scl2.Y * x1x3 + scl3.Y * x1x2,
			scl0.Z * x1x2plusx1x3plusx2x3 + scl1.Z * x2x3 + scl2.Z * x1x3 + scl3.Z * x1x2,
			scl0.w * x1x2plusx1x3plusx2x3 + scl1.w * x2x3 + scl2.w * x1x3 + scl3.w * x1x2
		);
		Vector4 c2 = new(
			-( scl0.X * x1plusx2plusx3 + scl1.X * x2plusx3 + scl2.X * x0plusx1plusx3 + scl3.X * x0plusx1plusx2 ),
			-( scl0.Y * x1plusx2plusx3 + scl1.Y * x2plusx3 + scl2.Y * x0plusx1plusx3 + scl3.Y * x0plusx1plusx2 ),
			-( scl0.Z * x1plusx2plusx3 + scl1.Z * x2plusx3 + scl2.Z * x0plusx1plusx3 + scl3.Z * x0plusx1plusx2 ),
			-( scl0.w * x1plusx2plusx3 + scl1.w * x2plusx3 + scl2.w * x0plusx1plusx3 + scl3.w * x0plusx1plusx2 )
		);
		Vector4 c3 = new(
			scl0.X + scl1.X + scl2.X + scl3.x,
			scl0.Y + scl1.Y + scl2.Y + scl3.y,
			scl0.Z + scl1.Z + scl2.Z + scl3.z,
			scl0.w + scl1.w + scl2.w + scl3.w
		);

		return new Polynomial4D( c0, c1, c2, c3 );
	}

	/// <inheritdoc cref="Polynomial2D.GetBounds01"/>
	public (FloatRange x, FloatRange y, FloatRange z, FloatRange w) GetBounds01() => ( x.OutputRange01, y.OutputRange01, z.OutputRange01, w.OutputRange01 );

	/// <inheritdoc cref="Polynomial.Split01"/>
	public (Polynomial4D pre, Polynomial4D post) Split01( float u ) {
		( var xPre, var xPost ) = x.Split01( u );
		( var yPre, var yPost ) = y.Split01( u );
		( var zPre, var zPost ) = z.Split01( u );
		( var wPre, var wPost ) = w.Split01( u );
		return ( new Polynomial4D( xPre, yPre, zPre, wPre ), new Polynomial4D( xPost, yPost, zPost, wPost ) );
	}

	#region IParamCurve3Diff interface implementations

	public int Degree => Mathf.Max( x.Degree, y.Degree, z.Degree, w.Degree );
	public Vector4 EvalDerivative( float t ) => Differentiate().Eval( t );
	public Vector4 EvalSecondDerivative( float t ) => Differentiate( 2 ).Eval( t );
	public Vector4 EvalThirdDerivative( float t = 0 ) => Differentiate( 3 ).Eval( 0 );

	#endregion

	#region Project Point

	/// <inheritdoc cref="Polynomial2D.ProjectPoint(Vector2,int,int)"/>
	public Vector4 ProjectPoint( Vector4 point, int initialSubdivisions = 16, int refinementIterations = 4 ) => ProjectPoint( point, out _, initialSubdivisions, refinementIterations );

	struct PointProjectSample {
		public float t;
		public float distDeltaSq;
		public Vector4 f;
		public Vector4 fp;
	}

	static PointProjectSample[] pointProjectGuesses = { default, default, default };

	/// <inheritdoc cref="Polynomial2D.ProjectPoint(Vector2,out float,int,int)"/>
	public Vector4 ProjectPoint( Vector4 point, out float t, int initialSubdivisions = 16, int refinementIterations = 4 ) {
		// define a bezier relative to the test point
		var curve = this;
		curve.x.c0 -= point.x; // constant coefficient defines the start position
		curve.y.c0 -= point.y;
		curve.z.c0 -= point.z;
		var curveStart = curve.Eval( 0 );
		var curveEnd = curve.Eval( 1 );

		PointProjectSample SampleDistSqDelta( float tSmp ) {
			var s = new PointProjectSample { t = tSmp };
			( s.f, s.fp ) = ( curve.Eval( tSmp ), curve.EvalDerivative( tSmp ) );
			s.distDeltaSq = Vector4.Dot( s.f, s.fp );
			return s;
		}

		// find initial candidates
		var candidatesFound = 0;
		var prevSmp = SampleDistSqDelta( 0 );

		for( var i = 1; i < initialSubdivisions; i++ ) {
			var ti = i / ( initialSubdivisions - 1f );
			var smp = SampleDistSqDelta( ti );
			if( Mathf.SignAsInt( smp.distDeltaSq ) != Mathf.SignAsInt( prevSmp.distDeltaSq ) ) {
				pointProjectGuesses[candidatesFound++] = SampleDistSqDelta( ( prevSmp.t + smp.t ) / 2 );
				if( candidatesFound == 3 ) break; // no more than three possible candidates because of the polynomial degree
			}

			prevSmp = smp;
		}

		// refine each guess w. Newton-Raphson iterations
		void Refine( ref PointProjectSample smp ) {
			var fpp = curve.EvalSecondDerivative( smp.t );
			var tNew = smp.t - Vector4.Dot( smp.f, smp.fp ) / ( Vector4.Dot( smp.f, fpp ) + Vector4.Dot( smp.fp, smp.fp ) );
			smp = SampleDistSqDelta( tNew );
		}

		for( var p = 0; p < candidatesFound; p++ )
		for( var i = 0; i < refinementIterations; i++ )
			Refine( ref pointProjectGuesses[p] );

		// Now find closest. First include the endpoints
		float sqDist0 = curveStart.sqrMagnitude; // include endpoints
		float sqDist1 = curveEnd.sqrMagnitude;
		var firstClosest = sqDist0 < sqDist1;
		float tClosest = firstClosest ? 0 : 1;
		var ptClosest = ( firstClosest ? curveStart : curveEnd ) + point;
		var distSqClosest = firstClosest ? sqDist0 : sqDist1;

		// then check internal roots
		for( var i = 0; i < candidatesFound; i++ ) {
			float pSqmag = pointProjectGuesses[i].f.sqrMagnitude;
			if( pSqmag < distSqClosest ) {
				distSqClosest = pSqmag;
				tClosest = pointProjectGuesses[i].t;
				ptClosest = pointProjectGuesses[i].f + point;
			}
		}

		t = tClosest;
		return ptClosest;
	}

	#endregion

	#region Typecasting & Operators

	public static Polynomial4D operator /( Polynomial4D p, float v ) => new(p.C0 / v, p.C1 / v, p.C2 / v, p.C3 / v);
	public static Polynomial4D operator *( Polynomial4D p, float v ) => new(p.C0 * v, p.C1 * v, p.C2 * v, p.C3 * v);
	public static Polynomial4D operator *( float v, Polynomial4D p ) => p * v;

	public static explicit operator Vector4Matrix3x1( Polynomial4D poly ) => new(poly.C0, poly.C1, poly.C2);
	public static explicit operator Vector4Matrix4x1( Polynomial4D poly ) => new(poly.C0, poly.C1, poly.C2, poly.C3);
	public static explicit operator BezierQuad4D( Polynomial4D poly ) => poly.Degree < 3 ? new BezierQuad4D( CharMatrix.QuadraticBezierInverse * (Vector4Matrix3x1)poly ) : throw new InvalidCastException( "Cannot cast a cubic polynomial to a quadratic curve" );
	public static explicit operator BezierCubic4D( Polynomial4D poly ) => new(CharMatrix.CubicBezierInverse * (Vector4Matrix4x1)poly);
	public static explicit operator CatRomCubic4D( Polynomial4D poly ) => new(CharMatrix.CubicCatmullRomInverse * (Vector4Matrix4x1)poly);
	public static explicit operator HermiteCubic4D( Polynomial4D poly ) => new(CharMatrix.CubicHermiteInverse * (Vector4Matrix4x1)poly);
	public static explicit operator UBSCubic4D( Polynomial4D poly ) => new(CharMatrix.CubicUniformBsplineInverse * (Vector4Matrix4x1)poly);

	#endregion

}