// by Freya Holmér (https://github.com/FreyaHolmer/MathPlus)

using System;
using System.Runtime.CompilerServices;
using System.Numerics;
using static Freya.MathPlus;

namespace MathPlus;

/// <summary>A 2D rectangular box</summary>
[Serializable] public partial struct Box2D {

	const MethodImplOptions INLINE = MethodImplOptions.AggressiveInlining;

	/// <summary>The number of vertices in a rectangle</summary>
	public const int VERTEX_COUNT = 4;

	/// <summary>The center of this rectangle</summary>
	public Vector2 center;

	/// <summary>The extents of this rectangle (distance from the center to the edge) per axis</summary>
	public Vector2 extents;

}

/// <summary>A 3D rectangular box, also known as a cuboid</summary>
[Serializable] public partial struct Box3D {

	const MethodImplOptions INLINE = MethodImplOptions.AggressiveInlining;

	/// <summary>The number of vertices in a cuboid</summary>
	public const int VERTEX_COUNT = 8;

	/// <summary>The center of this cuboid</summary>
	public Vector3 center;

	/// <summary>The extents of this cuboid (distance from the center to the edge) per axis</summary>
	public Vector3 extents;

}

#region Size

public partial struct Box2D {
	/// <summary>The size of this box</summary>
	public Vector2 Size {
		[MethodImpl( INLINE )] get => extents * 2;
	}

	/// <summary>The width of this rectangle</summary>
	public float Width {
		[MethodImpl( INLINE )] get => extents.X * 2;
	}

	/// <summary>The height of this rectangle</summary>
	public float Height {
		[MethodImpl( INLINE )] get => extents.Y * 2;
	}
}

public partial struct Box3D {
	/// <inheritdoc cref="Box2D.Size"/>
	public Vector3 Size {
		[MethodImpl( INLINE )] get => extents * 2;
	}

	/// <inheritdoc cref="Box2D.Width"/>
	public float Width {
		[MethodImpl( INLINE )] get => extents.X * 2;
	}

	/// <inheritdoc cref="Box2D.Height"/>
	public float Height {
		[MethodImpl( INLINE )] get => extents.Y * 2;
	}

	/// <summary>The depth of this box</summary>
	public float Depth {
		[MethodImpl( INLINE )] get => extents.Z * 2;
	}
}

#endregion

#region Min/Max

public partial struct Box2D {
	/// <summary>The minimum coordinates inside the box, per axis</summary>
	public Vector2 Min {
		[MethodImpl( INLINE )] get => new Vector2( center.X - extents.X,center.Y - extents.Y );
	}
	/// <summary>The maximum coordinates inside the box, per axis</summary>
	public Vector2 Max {
		[MethodImpl( INLINE )] get => new Vector2( center.X + extents.X,center.Y + extents.Y );
	}
}

public partial struct Box3D {
	/// <inheritdoc cref="Box2D.Min"/>
	public Vector3 Min {
		[MethodImpl( INLINE )] get => new Vector3( center.X - extents.X,center.Y - extents.Y,center.Z - extents.Z );
	}
	/// <inheritdoc cref="Box2D.Max"/>
	public Vector3 Max {
		[MethodImpl( INLINE )] get => new Vector3( center.X + extents.X,center.Y + extents.Y,center.Z + extents.Z );
	}
}

#endregion

#region Vertices

public partial struct Box2D {
	/// <summary>Returns a vertex of this box by index</summary>
	/// <param name="index">The index of the vertex to retrieve</param>
	public Vector2 GetVertex( int index ) {
		switch( index ) {
			case 0:  return new Vector2( center.X - extents.X,center.Y - extents.Y );
			case 1:  return new Vector2( center.X - extents.X,center.Y + extents.Y );
			case 2:  return new Vector2( center.X + extents.X,center.Y - extents.Y );
			case 3:  return new Vector2( center.X + extents.X,center.Y + extents.Y );
			default: throw new ArgumentOutOfRangeException( nameof(index), $"Invalid index: {index}. Valid vertex indices range from 0 to {VERTEX_COUNT - 1}" );
		}
	}
}

public partial struct Box3D {
	/// <inheritdoc cref="Box2D.GetVertex"/>
	public Vector2 GetVertex( int index ) {
		switch( index ) {
			case 0:  return new Vector3( center.X - extents.X,center.Y - extents.Y,center.Z - extents.Z );
			case 1:  return new Vector3( center.X - extents.X,center.Y - extents.Y,center.Z + extents.Z );
			case 2:  return new Vector3( center.X - extents.X,center.Y + extents.Y,center.Z - extents.Z );
			case 3:  return new Vector3( center.X - extents.X,center.Y + extents.Y,center.Z + extents.Z );
			case 4:  return new Vector3( center.X + extents.X,center.Y - extents.Y,center.Z - extents.Z );
			case 5:  return new Vector3( center.X + extents.X,center.Y - extents.Y,center.Z + extents.Z );
			case 6:  return new Vector3( center.X + extents.X,center.Y + extents.Y,center.Z - extents.Z );
			case 7:  return new Vector3( center.X + extents.X,center.Y + extents.Y,center.Z + extents.Z );
			default: throw new ArgumentOutOfRangeException( nameof(index), $"Invalid index: {index}. Valid vertex indices range from 0 to {VERTEX_COUNT - 1}" );
		}
	}
}

#endregion

#region Volume

public partial struct Box2D {
	/// <summary>The area of this rectangle</summary>
	public float Area {
		[MethodImpl( INLINE )] get => ( extents.X * extents.Y ) * 4;
	}
}

public partial struct Box3D {
	/// <summary>The volume of this cuboid</summary>
	public float Volume {
		[MethodImpl( INLINE )] get => ( extents.X * extents.Y * extents.Z ) * 8;
	}
}

#endregion

#region Surface

public partial struct Box2D {
	/// <summary>The total perimeter length of this rectangle</summary>
	public float Perimeter {
		[MethodImpl( INLINE )] get => 4 * ( extents.X + extents.Y );
	}
}

public partial struct Box3D {
	/// <summary>The total surface area of this cuboid</summary>
	public float SurfaceArea {
		[MethodImpl( INLINE )] get => 4 * ( extents.Y * ( extents.X + extents.Z ) + extents.Z * extents.X );
	}
}

#endregion

#region Contains

public partial struct Box2D {
	/// <summary>Returns whether or not a point is inside this box</summary>
	/// <param name="point">The point to test if it's inside</param>
	[MethodImpl( INLINE )] public bool Contains( Vector2 point ) => Abs( point.X - center.X ) - extents.X <= 0 && Abs( point.Y - center.Y ) - extents.Y <= 0;
}

public partial struct Box3D {
	/// <inheritdoc cref="Box2D.Contains"/>
	[MethodImpl( INLINE )] public bool Contains( Vector3 point ) => Abs( point.X - center.X ) - extents.X <= 0 && Abs( point.Y - center.Y ) - extents.Y <= 0 && Abs( point.Z - center.Z ) - extents.Z <= 0;
}

#endregion

#region Encapsulate

public partial struct Box2D {
	/// <summary>Extends the boundary of this box to encapsulate a point</summary>
	/// <param name="point">The point to encapsulate</param>
	public void Encapsulate( Vector2 point ) {
		float minX = Min( center.X - extents.X,point.X );
		float minY = Min( center.Y - extents.Y,point.Y );
		float maxX = Max( center.X + extents.X,point.X );
		float maxY = Max( center.Y + extents.Y,point.Y );
		center.X = ( maxX + minX ) / 2;
		center.Y = ( maxY + minY ) / 2;
		extents.X = ( maxX - minX ) / 2;
		extents.Y = ( maxY - minY ) / 2;
	}
}

public partial struct Box3D {
	/// <inheritdoc cref="Box2D.Encapsulate"/>
	public void Encapsulate( Vector3 point ) {
		float minX = Min( center.X - extents.X,point.X );
		float minY = Min( center.Y - extents.Y,point.Y );
		float minZ = Min( center.Z - extents.Z,point.Z );
		float maxX = Max( center.X + extents.X,point.X );
		float maxY = Max( center.Y + extents.Y,point.Y );
		float maxZ = Max( center.Z + extents.Z,point.Z );
		center.X = ( maxX + minX ) / 2;
		center.Y = ( maxY + minY ) / 2;
		center.Z = ( maxZ + minZ ) / 2;
		extents.X = ( maxX - minX ) / 2;
		extents.Y = ( maxY - minY ) / 2;
		extents.Z = ( maxZ - minZ ) / 2;
	}
}

#endregion

#region Closest Corner

public partial struct Box2D {
	/// <summary>Returns the corner of this box closest to the given point</summary>
	/// <param name="point">The point to get the closest corner to</param>
	public Vector2 ClosestCorner( Vector2 point ) =>
		new Vector2(
			center.X + Sign( point.X - center.X ) * extents.x,
			center.Y + Sign( point.Y - center.Y ) * extents.y
		);
}

public partial struct Box3D {
	/// <inheritdoc cref="Box2D.ClosestCorner"/>
	public Vector3 ClosestCorner( Vector3 point ) =>
		new Vector3(
			center.X + Sign( point.X - center.X ) * extents.x,
			center.Y + Sign( point.Y - center.Y ) * extents.y,
			center.Z + Sign( point.Z - center.Z ) * extents.z
		);
}

#endregion

#region Closest point inside

public partial struct Box2D {
	/// <summary>Returns the point inside the box, closest to another point.
	/// Points already inside will return the same location</summary>
	/// <param name="point">The point to get the closest point to</param>
	[MethodImpl( INLINE )] public Vector2 ClosestPointInside( Vector2 point ) =>
		new Vector2(
			point.x.Clamp( center.X - extents.X,center.X + extents.X ),
			point.y.Clamp( center.Y - extents.Y,center.Y + extents.Y )
		);
}

public partial struct Box3D {
	/// <inheritdoc cref="Box2D.ClosestPointInside"/>
	[MethodImpl( INLINE )] public Vector3 ClosestPointInside( Vector3 point ) =>
		new Vector3(
			point.x.Clamp( center.X - extents.X,center.X + extents.X ),
			point.y.Clamp( center.Y - extents.Y,center.Y + extents.Y ),
			point.z.Clamp( center.Z - extents.Z,center.Z + extents.Z )
		);
}

#endregion

#region Closest point on bounds

public partial struct Box2D {
	/// <summary>Projects a point onto the boundary of this box. Points inside will be pushed out to the boundary</summary>
	/// <param name="point">The point to project onto the box boundary</param>
	[MethodImpl( INLINE )] public Vector2 ClosestPointOnBoundary( Vector2 point ) {
		var px = point.X - center.x;
		var py = point.Y - center.y;
		float ax = Abs( px );
		float ay = Abs( py );
		var dx = ax - extents.x;
		var dy = ay - extents.y;
		var caseX = dy <= dx;
		var caseY = caseX == false;
		return new Vector2(
			center.X + Sign( px ) * ( caseX ? extents.X : ax.AtMost( extents.X ) ),
			center.Y + Sign( py ) * ( caseY ? extents.Y : ay.AtMost( extents.Y ) )
		);
	}
}

public partial struct Box3D {
	/// <inheritdoc cref="Box2D.ClosestPointOnBoundary"/>
	[MethodImpl( INLINE )] public Vector3 ClosestPointOnBoundary( Vector3 point ) {
		var px = point.X - center.x;
		var py = point.Y - center.y;
		var pz = point.Z - center.z;
		float ax = Abs( px );
		float ay = Abs( py );
		float az = Abs( pz );
		var dx = ax - extents.x;
		var dy = ay - extents.y;
		var dz = az - extents.z;
		var caseX = dz <= dx && dy <= dx;
		var caseY = caseX == false && dx <= dy && dz <= dy;
		var caseZ = caseX == false && caseY == false;
		return new Vector3(
			center.X + Sign( px ) * ( caseX ? extents.X : ax.AtMost( extents.X ) ),
			center.Y + Sign( py ) * ( caseY ? extents.Y : ay.AtMost( extents.Y ) ),
			center.Z + Sign( pz ) * ( caseZ ? extents.Z : az.AtMost( extents.Z ) )
		);
	}
}

#endregion

#region Intersection Tests (2D only, for now)

public partial struct Box2D {
	/// <summary>Returns whether or not an infinite line intersects this box</summary>
	/// <param name="line">The line to see if it intersects</param>
	public bool Intersects( Line2D line ) => IntersectionTest.LineRectOverlap( extents, line.origin - center, line.dir );

	/// <summary>Returns the intersection points of this rectangle and a line</summary>
	/// <param name="line">The line to get intersection points of</param>
	public ResultsMax2<Vector2> Intersect( Line2D line ) => IntersectionTest.LinearRectPoints( center, extents, line.origin, line.dir );
}

#endregion