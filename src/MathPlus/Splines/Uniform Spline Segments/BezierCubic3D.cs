// by Freya Holmér (https://github.com/FreyaHolmer/MathPlus)
// Do not manually edit - this file is generated by MathPlusCodegen.cs

using System;
using System.Runtime.CompilerServices;
using System.Numerics;

namespace MathPlus;

/// <summary>An optimized uniform 3D Cubic bézier segment, with 4 control points</summary>
[Serializable] public struct BezierCubic3D : IParamSplineSegment<Polynomial3D,Vector3Matrix4x1> {

	const MethodImplOptions INLINE = MethodImplOptions.AggressiveInlining;

	Vector3Matrix4x1 pointMatrix;
	[NonSerialized] Polynomial3D curve;
	[NonSerialized] bool validCoefficients;

	/// <summary>Creates a uniform 3D Cubic bézier segment, from 4 control points</summary>
	/// <param name="p0">The starting point of the curve</param>
	/// <param name="p1">The second control point of the curve, sometimes called the start tangent point</param>
	/// <param name="p2">The third control point of the curve, sometimes called the end tangent point</param>
	/// <param name="p3">The end point of the curve</param>
	public BezierCubic3D( Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3 ) : this(new Vector3Matrix4x1(p0, p1, p2, p3)){}
	/// <summary>Creates a uniform 3D Cubic bézier segment, from 4 control points</summary>
	/// <param name="pointMatrix">The matrix containing the control points of this spline</param>
	public BezierCubic3D( Vector3Matrix4x1 pointMatrix ) => (this.pointMatrix,curve,validCoefficients) = (pointMatrix,default,false);

	public Polynomial3D Curve {
		get {
			if( validCoefficients )
				return curve; // no need to update
			validCoefficients = true;
			return curve = new Polynomial3D(
				P0,
				3*(-P0+P1),
				3*P0-6*P1+3*P2,
				-P0+3*P1-3*P2+P3
			);
		}
	}
	public Vector3Matrix4x1 PointMatrix {[MethodImpl( INLINE )] get => pointMatrix; [MethodImpl( INLINE )] set => _ = ( pointMatrix = value, validCoefficients = false ); }
	/// <summary>The starting point of the curve</summary>
	public Vector3 P0{ [MethodImpl( INLINE )] get => pointMatrix.m0; [MethodImpl( INLINE )] set => _ = ( pointMatrix.m0 = value, validCoefficients = false ); }
	/// <summary>The second control point of the curve, sometimes called the start tangent point</summary>
	public Vector3 P1{ [MethodImpl( INLINE )] get => pointMatrix.m1; [MethodImpl( INLINE )] set => _ = ( pointMatrix.m1 = value, validCoefficients = false ); }
	/// <summary>The third control point of the curve, sometimes called the end tangent point</summary>
	public Vector3 P2{ [MethodImpl( INLINE )] get => pointMatrix.m2; [MethodImpl( INLINE )] set => _ = ( pointMatrix.m2 = value, validCoefficients = false ); }
	/// <summary>The end point of the curve</summary>
	public Vector3 P3{ [MethodImpl( INLINE )] get => pointMatrix.m3; [MethodImpl( INLINE )] set => _ = ( pointMatrix.m3 = value, validCoefficients = false ); }
	/// <summary>Get or set a control point position by index. Valid indices from 0 to 3</summary>
	public Vector3 this[ int i ] {
		get => i switch { 0 => P0, 1 => P1, 2 => P2, 3 => P3, _ => throw new ArgumentOutOfRangeException( nameof(i), $"Index has to be in the 0 to 3 range, and I think {i} is outside that range you know" ) };
		set { switch( i ){ case 0: P0 = value; break; case 1: P1 = value; break; case 2: P2 = value; break; case 3: P3 = value; break; default: throw new ArgumentOutOfRangeException( nameof(i), $"Index has to be in the 0 to 3 range, and I think {i} is outside that range you know" ); }}
	}
	public static bool operator ==( BezierCubic3D a, BezierCubic3D b ) => a.pointMatrix == b.pointMatrix;
	public static bool operator !=( BezierCubic3D a, BezierCubic3D b ) => !( a == b );
	public bool Equals( BezierCubic3D other ) => P0.Equals( other.P0 ) && P1.Equals( other.P1 ) && P2.Equals( other.P2 ) && P3.Equals( other.P3 );
	public override bool Equals( object obj ) => obj is BezierCubic3D other && pointMatrix.Equals( other.pointMatrix );
	public override int GetHashCode() => pointMatrix.GetHashCode();
	public override string ToString() => $"({pointMatrix.m0}, {pointMatrix.m1}, {pointMatrix.m2}, {pointMatrix.m3})";

	/// <summary>Returns this curve flattened to 2D. Effectively setting z = 0</summary>
	/// <param name="curve3D">The 3D curve to flatten to the Z plane</param>
	public static explicit operator BezierCubic2D( BezierCubic3D curve3D ) => new BezierCubic2D( curve3D.P0, curve3D.P1, curve3D.P2, curve3D.P3 );
	public static explicit operator HermiteCubic3D( BezierCubic3D s ) =>
		new HermiteCubic3D(
			s.P0,
			3*(-s.P0+s.P1),
			s.P3,
			3*(-s.P2+s.P3)
		);
	public static explicit operator CatRomCubic3D( BezierCubic3D s ) =>
		new CatRomCubic3D(
			6*s.P0-6*s.P1+s.P3,
			s.P0,
			s.P3,
			s.P0-6*s.P2+6*s.P3
		);
	public static explicit operator UBSCubic3D( BezierCubic3D s ) =>
		new UBSCubic3D(
			6*s.P0-7*s.P1+2*s.P2,
			2*s.P1-s.P2,
			-s.P1+2*s.P2,
			2*s.P1-7*s.P2+6*s.P3
		);
	/// <summary>Returns a linear blend between two bézier curves</summary>
	/// <param name="a">The first spline segment</param>
	/// <param name="b">The second spline segment</param>
	/// <param name="t">A value from 0 to 1 to blend between <c>a</c> and <c>b</c></param>
	public static BezierCubic3D Lerp( BezierCubic3D a, BezierCubic3D b, float t ) =>
		new(
			Vector3.LerpUnclamped( a.P0, b.P0, t ),
			Vector3.LerpUnclamped( a.P1, b.P1, t ),
			Vector3.LerpUnclamped( a.P2, b.P2, t ),
			Vector3.LerpUnclamped( a.P3, b.P3, t )
		);

	/// <summary>Returns a linear blend between two bézier curves, where the tangent directions are spherically interpolated</summary>
	/// <param name="a">The first spline segment</param>
	/// <param name="b">The second spline segment</param>
	/// <param name="t">A value from 0 to 1 to blend between <c>a</c> and <c>b</c></param>
	public static BezierCubic3D Slerp( BezierCubic3D a, BezierCubic3D b, float t ) {
		Vector3 P0 = Vector3.LerpUnclamped( a.P0, b.P0, t );
		Vector3 P3 = Vector3.LerpUnclamped( a.P3, b.P3, t );
		return new BezierCubic3D(
			P0,
			P0 + Vector3.SlerpUnclamped( a.P1 - a.P0, b.P1 - b.P0, t ),
			P3 + Vector3.SlerpUnclamped( a.P2 - a.P3, b.P2 - b.P3, t ),
			P3
		);
	}
	/// <summary>Splits this curve at the given t-value, into two curves that together form the exact same shape</summary>
	/// <param name="t">The t-value to split at</param>
	public (BezierCubic3D pre, BezierCubic3D post) Split( float t ) {
		var a = new Vector3(
			P0.X + ( P1.X - P0.X ) * t,
			P0.Y + ( P1.Y - P0.Y ) * t,
			P0.Z + ( P1.Z - P0.Z ) * t );
		var b = new Vector3(
			P1.X + ( P2.X - P1.X ) * t,
			P1.Y + ( P2.Y - P1.Y ) * t,
			P1.Z + ( P2.Z - P1.Z ) * t );
		var c = new Vector3(
			P2.X + ( P3.X - P2.X ) * t,
			P2.Y + ( P3.Y - P2.Y ) * t,
			P2.Z + ( P3.Z - P2.Z ) * t );
		var d = new Vector3(
			a.X + ( b.X - a.X ) * t,
			a.Y + ( b.Y - a.Y ) * t,
			a.Z + ( b.Z - a.Z ) * t );
		var e = new Vector3(
			b.X + ( c.X - b.X ) * t,
			b.Y + ( c.Y - b.Y ) * t,
			b.Z + ( c.Z - b.Z ) * t );
		var p = new Vector3(
			d.X + ( e.X - d.X ) * t,
			d.Y + ( e.Y - d.Y ) * t,
			d.Z + ( e.Z - d.Z ) * t );
		return ( new BezierCubic3D( P0, a, d, p ), new BezierCubic3D( p, e, c, P3 ) );
	}
}