using System;
using System.Numerics;

namespace MathPlus;

/// <summary>a 2D arc with support for straight lines</summary>
[Serializable]
public struct Arc2D {

	/// <summary>The starting point of the arc</summary>
	public Transform2D placement;
	/// <summary>The signed curvature of the arc, equal to 1/radius (0 = straight line, 1 = turning left, -1 = turning right)</summary>
	public float curvature;
	/// <summary>The length of the arc</summary>
	public float length;

	/// <summary>The radius of the circle traced by the arc. Returns infinity if this segment is linear, ie: if curvature is 0</summary>
	public float Radius => 1f / MathF.Abs( curvature );
	/// <summary>The center of the circle traced by the arc. Returns infinity if this segment is linear, ie: if curvature is 0</summary>
	public Vector2 CircleCenter => StartNormal / curvature;
	/// <summary>The normal direction at the start of the arc</summary>
	public Vector2 StartNormal => placement.AxisY;
	/// <summary>The tangent direction at the start of the arc</summary>
	public Vector2 StartTangent => placement.AxisX;
	/// <summary>The normal direction at the end of the arc</summary>
	public Vector2 EndNormal => GetNormal( length );
	/// <summary>The end point of the arc</summary>
	public Vector2 EndPoint => GetPosition( length );
	/// <summary>The signed angular span covered across the arc. This returns 0 if this segment is linear, ie: if curvature is 0</summary>
	public float AngularSpan => length * curvature; // s = ra  ▶  s = a/k  ▶  sk = a
	/// <summary>Whether or not this is a straight line rather than an arc, ie: if curvature is 0</summary>
	public bool IsStraight => Mathf.Approximately( curvature, 0 );

	/// <summary>Evaluates the position of this arc at the given arc length <c>s</c></summary>
	public Vector2 GetPosition( float s ) => Eval( s, nThDerivative: 0 );

	/// <summary>Evaluates the tangent direction of this arc at the given arc length <c>s</c></summary>
	public Vector2 GetTangent( float s ) => s == 0 ? StartTangent : Eval( s, nThDerivative: 1 ); // no need to normalize, it's already arc-length parameterized

	/// <summary>Evaluates the normal direction of this arc at the given arc length <c>s</c></summary>
	public Vector2 GetNormal( float s ) => Eval( s, nThDerivative: 1 ).Rotate90CCW(); // no need to normalize, it's already arc-length parameterized

	/// <summary>Evaluates the given derivative of this arc, by arc length <c>s</c></summary>
	public Vector2 Eval( float s, int nThDerivative = 0 ) {
		var ang = s * curvature;
		float x, y;

		switch( nThDerivative ) {
			case 0:
				x = s * Mathf.Sinc( ang );
				y = s * Mathf.Cosinc( ang );
				return placement.TransformPoint( x, y );
			case 1:
				x = MathF.Cos( ang );
				y = MathF.Sin( ang );
				break;
			case 2:
				x = -curvature * MathF.Sin( ang );
				y = +curvature * MathF.Cos( ang );
				break;
			case 3:
				var k2 = curvature * curvature;
				x = -k2 * MathF.Cos( ang );
				y = -k2 * MathF.Sin( ang );
				break;
			case 4:
				var k3 = curvature * curvature * curvature;
				x = +k3 * MathF.Sin( ang );
				y = -k3 * MathF.Cos( ang );
				break;
			case 5:
				var _k2 = curvature * curvature;
				var k4 = _k2 * _k2;
				x = k4 * MathF.Cos( ang );
				y = k4 * MathF.Sin( ang );
				break;
			default:
				// general form for n > 0
				var scale = MathF.Pow( curvature, nThDerivative - 1 );
				var xSgn = nThDerivative / 2 % 2 == 0 ? 1 : -1;
				var ySgn = ( nThDerivative - 1 ) / 2 % 2 == 0 ? 1 : -1;
				var even = nThDerivative % 2 == 0;
				x = xSgn * scale * ( even ? MathF.Sin( ang ) : MathF.Cos( ang ) );
				y = ySgn * scale * ( even ? MathF.Cos( ang ) : MathF.Sin( ang ) );
				break;
		}

		// space transformation
		return placement.TransformVector( x, y );
	}

}