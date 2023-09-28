// by Freya Holmér (https://github.com/FreyaHolmer/MathPlus)

using System.Collections.Generic;
using System.Numerics;

namespace MathPlus;

public class LagrangePolynomial2D {

	public List<Vector2> points = new List<Vector2>();
	public List<float> knots = null;
	public bool Uniform => knots == null;
	public FloatRange InternalKnotRange => Uniform ? ( 0, points.Count - 1 ) : ( knots[0], knots[^1] );

	public Vector2 Eval( float u ) {
		float l( int j ) {
			float prod = 1;
			for( var i = 0; i < points.Count; i++ ) {
				if( i == j )
					continue;
				if( Uniform )
					prod *= ( u - i ) / ( j - i );
				else
					prod *= Mathf.InverseLerp( knots[i], knots[j], u );
			}

			return prod;
		}

		var sum = Vector2.Zero;
		for( var j = 0; j < points.Count; j++ )
			sum += points[j] * l( j );

		return sum;
	}

}