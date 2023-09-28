// by Freya Holmér (https://github.com/FreyaHolmer/MathPlus)
// Do not manually edit - this file is generated by MathPlusCodegen.cs

using System;
using System.Runtime.CompilerServices;
using System.Numerics;

namespace MathPlus;

/// <summary>An optimized uniform 2D Quadratic bézier segment, with 3 control points</summary>
[Serializable] public struct BezierQuad2D : IParamSplineSegment<Polynomial2D,Vector2Matrix3x1> {

	const MethodImplOptions INLINE = MethodImplOptions.AggressiveInlining;

	Vector2Matrix3x1 pointMatrix;
	[NonSerialized] Polynomial2D curve;
	[NonSerialized] bool validCoefficients;

	/// <summary>Creates a uniform 2D Quadratic bézier segment, from 3 control points</summary>
	/// <param name="p0">The starting point of the curve</param>
	/// <param name="p1">The middle control point of the curve, sometimes called a tangent point</param>
	/// <param name="p2">The end point of the curve</param>
	public BezierQuad2D( Vector2 p0, Vector2 p1, Vector2 p2 ) : this(new Vector2Matrix3x1(p0, p1, p2)){}
	/// <summary>Creates a uniform 2D Quadratic bézier segment, from 3 control points</summary>
	/// <param name="pointMatrix">The matrix containing the control points of this spline</param>
	public BezierQuad2D( Vector2Matrix3x1 pointMatrix ) => (this.pointMatrix,curve,validCoefficients) = (pointMatrix,default,false);

	public Polynomial2D Curve {
		get {
			if( validCoefficients )
				return curve; // no need to update
			validCoefficients = true;
			return curve = new Polynomial2D(
				P0,
				2*(-P0+P1),
				P0-2*P1+P2
			);
		}
	}
	public Vector2Matrix3x1 PointMatrix {[MethodImpl( INLINE )] get => pointMatrix; [MethodImpl( INLINE )] set => _ = ( pointMatrix = value, validCoefficients = false ); }
	/// <summary>The starting point of the curve</summary>
	public Vector2 P0{ [MethodImpl( INLINE )] get => pointMatrix.m0; [MethodImpl( INLINE )] set => _ = ( pointMatrix.m0 = value, validCoefficients = false ); }
	/// <summary>The middle control point of the curve, sometimes called a tangent point</summary>
	public Vector2 P1{ [MethodImpl( INLINE )] get => pointMatrix.m1; [MethodImpl( INLINE )] set => _ = ( pointMatrix.m1 = value, validCoefficients = false ); }
	/// <summary>The end point of the curve</summary>
	public Vector2 P2{ [MethodImpl( INLINE )] get => pointMatrix.m2; [MethodImpl( INLINE )] set => _ = ( pointMatrix.m2 = value, validCoefficients = false ); }
	/// <summary>Get or set a control point position by index. Valid indices from 0 to 2</summary>
	public Vector2 this[ int i ] {
		get => i switch { 0 => P0, 1 => P1, 2 => P2, _ => throw new ArgumentOutOfRangeException( nameof(i), $"Index has to be in the 0 to 2 range, and I think {i} is outside that range you know" ) };
		set { switch( i ){ case 0: P0 = value; break; case 1: P1 = value; break; case 2: P2 = value; break; default: throw new ArgumentOutOfRangeException( nameof(i), $"Index has to be in the 0 to 2 range, and I think {i} is outside that range you know" ); }}
	}
	public static bool operator ==( BezierQuad2D a, BezierQuad2D b ) => a.pointMatrix == b.pointMatrix;
	public static bool operator !=( BezierQuad2D a, BezierQuad2D b ) => !( a == b );
	public bool Equals( BezierQuad2D other ) => P0.Equals( other.P0 ) && P1.Equals( other.P1 ) && P2.Equals( other.P2 );
	public override bool Equals( object obj ) => obj is BezierQuad2D other && pointMatrix.Equals( other.pointMatrix );
	public override int GetHashCode() => pointMatrix.GetHashCode();
	public override string ToString() => $"({pointMatrix.m0}, {pointMatrix.m1}, {pointMatrix.m2})";

	/// <summary>Returns a linear blend between two bézier curves</summary>
	/// <param name="a">The first spline segment</param>
	/// <param name="b">The second spline segment</param>
	/// <param name="t">A value from 0 to 1 to blend between <c>a</c> and <c>b</c></param>
	public static BezierQuad2D Lerp( BezierQuad2D a, BezierQuad2D b, float t ) =>
		new(
			Vector2.LerpUnclamped( a.P0, b.P0, t ),
			Vector2.LerpUnclamped( a.P1, b.P1, t ),
			Vector2.LerpUnclamped( a.P2, b.P2, t )
		);
	/// <summary>Splits this curve at the given t-value, into two curves that together form the exact same shape</summary>
	/// <param name="t">The t-value to split at</param>
	public (BezierQuad2D pre, BezierQuad2D post) Split( float t ) {
		var a = new Vector2(
			P0.X + ( P1.X - P0.X ) * t,
			P0.Y + ( P1.Y - P0.Y ) * t );
		var b = new Vector2(
			P1.X + ( P2.X - P1.X ) * t,
			P1.Y + ( P2.Y - P1.Y ) * t );
		var p = new Vector2(
			a.X + ( b.X - a.X ) * t,
			a.Y + ( b.Y - a.Y ) * t );
		return ( new BezierQuad2D( P0, a, p ), new BezierQuad2D( p, b, P2 ) );
	}
}