using System.Numerics;

namespace MathPlus;

/// <summary>An oriented 2D plane embedded in 3D space</summary>
public struct Plane2DIn3D
{
    public static readonly Plane2DIn3D XY = new(default, Vector3.right, Vector3.up);
    public static readonly Plane2DIn3D YZ = new(default, Vector3.up, Vector3.forward);
    public static readonly Plane2DIn3D ZX = new(default, Vector3.forward, Vector3.right);

    public Vector3 origin, axisX, axisY;

    /// <summary>Creates an oriented 2D plane embedded in 3D space</summary>
    /// <param name="origin">The origin of the plane</param>
    /// <param name="axisX">The x axis direction of the plane</param>
    /// <param name="axisY">The y axis direction of the plane</param>
    public Plane2DIn3D(Vector3 origin, Vector3 axisX, Vector3 axisY) =>
        (this.origin, this.axisX, this.axisY) =
        (origin, axisX.Normalized(), axisY.Normalized());

    /// <summary>Rotates this plane around the Y axis, setting the X axis,
    /// so that the given point <c>p</c> is in the plane where x > 0</summary>
    /// <param name="p">The point to include in the plane</param>
    /// <param name="pLocal">The included point in the 2D local space</param>
    public void RotateAroundYToInclude(Vector3 p, out Vector2 pLocal)
    {
        var pRel = p - origin;
        var yProj = Vector3.Dot(axisY, pRel);
        axisX = (pRel - axisY * yProj).Normalized();
        var xProj = Vector3.Dot(axisX, pRel);
        pLocal = new Vector2(xProj, yProj);
    }

    /// <summary>Transforms a local 2D point to a 3D world space point</summary>
    /// <param name="pt">The local space point to transform</param>
    public Vector3 TransformPoint(Vector2 pt)
    {
        return new( // unrolled for performance
            origin.X + axisX.X * pt.X + axisY.X * pt.Y,
            origin.Y + axisX.Y * pt.X + axisY.Y * pt.Y,
            origin.Z + axisX.Z * pt.X + axisY.Z * pt.Y
        );
    }

    /// <summary>Transforms a local 2D vector to a 3D world space vector, not taking position into account</summary>
    /// <param name="vec">The local space vector to transform</param>
    public Vector3 TransformVector(Vector2 vec)
    {
        return new( // unrolled for performance
            axisX.X * vec.X + axisY.X * vec.Y,
            axisX.Y * vec.X + axisY.Y * vec.Y,
            axisX.Z * vec.X + axisY.Z * vec.Y
        );
    }

    /// <summary>Transform a 3D world space point to a local 2D point</summary>
    /// <param name="pt">World space point</param>
    public Vector2 InverseTransformPoint(Vector3 pt)
    {
        var rx = pt.X - origin.x;
        var ry = pt.Y - origin.y;
        var rz = pt.Z - origin.z;
        return new(
            axisX.X * rx + axisX.Y * ry + axisX.Z * rz,
            axisY.X * rx + axisY.Y * ry + axisY.Z * rz
        );
    }

    /// <summary>Transform a 3D world space vector to a local 2D vector</summary>
    /// <param name="vec">World space vector</param>
    public Vector2 InverseTransformVector(Vector3 vec)
    {
        return new(
            axisX.X * vec.X + axisX.Y * vec.Y + axisX.Z * vec.z,
            axisY.X * vec.X + axisY.Y * vec.Y + axisY.Z * vec.z
        );
    }
}