// by Freya Holmér (https://github.com/FreyaHolmer/Mathfs)
// Do not manually edit - this file is generated by MathfsCodegen.cs

using System;
using System.Runtime.CompilerServices;
using UnityEngine;

namespace Freya {

	/// <summary>An optimized uniform 2D Cubic hermite segment, with 4 control points</summary>
	[Serializable] public struct HermiteCubic2D : IParamCubicSplineSegment2D {

		const MethodImplOptions INLINE = MethodImplOptions.AggressiveInlining;

		/// <summary>Creates a uniform 2D Cubic hermite segment, from 4 control points</summary>
		/// <param name="p0">The starting point of the curve</param>
		/// <param name="v0">The rate of change (velocity) at the start of the curve</param>
		/// <param name="p1">The end point of the curve</param>
		/// <param name="v1">The rate of change (velocity) at the end of the curve</param>
		public HermiteCubic2D( Vector2 p0, Vector2 v0, Vector2 p1, Vector2 v1 ) {
			pointMatrix = new Vector2Matrix4x1( p0, v0, p1, v1 );
			validCoefficients = false;
			curve = default;
		}

		Polynomial2D curve;
		public Polynomial2D Curve {
			get {
				ReadyCoefficients();
				return curve;
			}
		}
		#region Control Points

		[SerializeField] Vector2Matrix4x1 pointMatrix;
		public Vector2Matrix4x1 PointMatrix => pointMatrix;

		/// <summary>The starting point of the curve</summary>
		public Vector2 P0 {
			[MethodImpl( INLINE )] get => pointMatrix.m0;
			[MethodImpl( INLINE )] set => _ = ( pointMatrix.m0 = value, validCoefficients = false );
		}

		/// <summary>The rate of change (velocity) at the start of the curve</summary>
		public Vector2 V0 {
			[MethodImpl( INLINE )] get => pointMatrix.m1;
			[MethodImpl( INLINE )] set => _ = ( pointMatrix.m1 = value, validCoefficients = false );
		}

		/// <summary>The end point of the curve</summary>
		public Vector2 P1 {
			[MethodImpl( INLINE )] get => pointMatrix.m2;
			[MethodImpl( INLINE )] set => _ = ( pointMatrix.m2 = value, validCoefficients = false );
		}

		/// <summary>The rate of change (velocity) at the end of the curve</summary>
		public Vector2 V1 {
			[MethodImpl( INLINE )] get => pointMatrix.m3;
			[MethodImpl( INLINE )] set => _ = ( pointMatrix.m3 = value, validCoefficients = false );
		}

		/// <summary>Get or set a control point position by index. Valid indices from 0 to 3</summary>
		public Vector2 this[ int i ] {
			get =>
				i switch {
					0 => P0,
					1 => V0,
					2 => P1,
					3 => V1,
					_ => throw new ArgumentOutOfRangeException( nameof(i), $"Index has to be in the 0 to 3 range, and I think {i} is outside that range you know" )
				};
			set {
				switch( i ) {
					case 0:
						P0 = value;
						break;
					case 1:
						V0 = value;
						break;
					case 2:
						P1 = value;
						break;
					case 3:
						V1 = value;
						break;
					default: throw new ArgumentOutOfRangeException( nameof(i), $"Index has to be in the 0 to 3 range, and I think {i} is outside that range you know" );
				}
			}
		}

		#endregion
		[NonSerialized] bool validCoefficients;

		[MethodImpl( INLINE )] void ReadyCoefficients() {
			if( validCoefficients )
				return; // no need to update
			validCoefficients = true;
			curve = new Polynomial2D(
				P0,
				V0,
				-3*P0-2*V0+3*P1-V1,
				2*P0+V0-2*P1+V1
			);
		}
		public static bool operator ==( HermiteCubic2D a, HermiteCubic2D b ) => a.pointMatrix == b.pointMatrix;
		public static bool operator !=( HermiteCubic2D a, HermiteCubic2D b ) => !( a == b );
		public bool Equals( HermiteCubic2D other ) => P0.Equals( other.P0 ) && V0.Equals( other.V0 ) && P1.Equals( other.P1 ) && V1.Equals( other.V1 );
		public override bool Equals( object obj ) => obj is HermiteCubic2D other && pointMatrix.Equals( other.pointMatrix );
		public override int GetHashCode() => pointMatrix.GetHashCode();
		public override string ToString() => $"({pointMatrix.m0}, {pointMatrix.m1}, {pointMatrix.m2}, {pointMatrix.m3})";

		/// <summary>Returns this spline segment in 3D, where z = 0</summary>
		/// <param name="curve2D">The 2D curve to cast to 3D</param>
		public static explicit operator HermiteCubic3D( HermiteCubic2D curve2D ) => new HermiteCubic3D( curve2D.P0, curve2D.V0, curve2D.P1, curve2D.V1 );
		public static explicit operator BezierCubic2D( HermiteCubic2D s ) =>
			new BezierCubic2D(
				s.P0,
				s.P0+(1/3f)*s.V0,
				s.P1-(1/3f)*s.V1,
				s.P1
			);
		public static explicit operator CatRomCubic2D( HermiteCubic2D s ) =>
			new CatRomCubic2D(
				-2*s.V0+s.P1,
				s.P0,
				s.P1,
				s.P0+2*s.V1
			);
		public static explicit operator UBSCubic2D( HermiteCubic2D s ) =>
			new UBSCubic2D(
				-s.P0-(7/3f)*s.V0+2*s.P1-(2/3f)*s.V1,
				2*s.P0+(2/3f)*s.V0-s.P1+(1/3f)*s.V1,
				-s.P0-(1/3f)*s.V0+2*s.P1-(2/3f)*s.V1,
				2*s.P0+(2/3f)*s.V0-s.P1+(7/3f)*s.V1
			);
		/// <summary>Returns a linear blend between two hermite curves</summary>
		/// <param name="a">The first spline segment</param>
		/// <param name="b">The second spline segment</param>
		/// <param name="t">A value from 0 to 1 to blend between <c>a</c> and <c>b</c></param>
		public static HermiteCubic2D Lerp( HermiteCubic2D a, HermiteCubic2D b, float t ) =>
			new(
				Vector2.LerpUnclamped( a.P0, b.P0, t ),
				Vector2.LerpUnclamped( a.V0, b.V0, t ),
				Vector2.LerpUnclamped( a.P1, b.P1, t ),
				Vector2.LerpUnclamped( a.V1, b.V1, t )
			);
	}
}
