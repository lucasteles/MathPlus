using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Serialization;

namespace MathPlus;

using System;
using System.Globalization;
using System.Runtime.CompilerServices;

// Representation of rays.
[Serializable, DataContract]
[StructLayout(LayoutKind.Sequential)]
public struct Ray : IFormattable, IEquatable<Ray>
{
    Vector3 origin;
    Vector3 direction;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Ray(Vector3 origin, Vector3 direction)
    {
        this.origin = origin;
        this.direction = direction.Normalized();
    }

    public Vector3 Origin
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => origin;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => origin = value;
    }

    public Vector3 Direction
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => direction;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => direction = value.Normalized();
    }

    public Vector3 GetPoint(float distance) => origin + direction * distance;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public override string ToString() => ToString(null, null);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public string ToString(string format) => ToString(format, null);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public string ToString(string? format, IFormatProvider? formatProvider)
    {
        if (string.IsNullOrEmpty(format)) format = "F2";

        formatProvider ??= CultureInfo.InvariantCulture.NumberFormat;
        return
            $"Origin: {origin.ToString(format, formatProvider)}, Dir: {direction.ToString(format, formatProvider)}";
    }

    public bool Equals(Ray other) =>
        origin.Equals(other.origin) && direction.Equals(other.direction);

    public override bool Equals(object? obj) => obj is Ray other && Equals(other);
    public override int GetHashCode() => HashCode.Combine(origin, direction);
    public static bool operator ==(Ray left, Ray right) => left.Equals(right);
    public static bool operator !=(Ray left, Ray right) => !(left == right);

    // adapted from http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-7-intersecting-simple-shapes/ray-box-intersection/
    /// <summary>
    /// Check if this <see cref="Ray"/> intersects the specified <see cref="Bounds"/>.
    /// </summary>
    public bool Intersects(Bounds box, out float distance)
    {
        const float epsilon = 1e-6f;

        float tMin = 0, tMax = 0;
        distance = 0;

        if (Math.Abs(Direction.X) < epsilon)
        {
            if (Origin.X < box.Min.X || Origin.X > box.Max.X)
                return false;
        }
        else
        {
            tMin = (box.Min.X - Origin.X) / Direction.X;
            tMax = (box.Max.X - Origin.X) / Direction.X;
            if (tMin > tMax) (tMin, tMax) = (tMax, tMin);
        }

        if (Math.Abs(Direction.Y) < epsilon)
        {
            if (Origin.Y < box.Min.Y || Origin.Y > box.Max.Y)
                return false;
        }
        else
        {
            var tMinY = (box.Min.Y - Origin.Y) / Direction.Y;
            var tMaxY = (box.Max.Y - Origin.Y) / Direction.Y;

            if (tMinY > tMaxY) (tMinY, tMaxY) = (tMaxY, tMinY);

            if (tMin > tMaxY || tMinY > tMax)
                return false;

            if (tMinY > tMin) tMin = tMinY;
            if (tMaxY < tMax) tMax = tMaxY;
        }

        if (Math.Abs(Direction.Z) < epsilon)
        {
            if (Origin.Z < box.Min.Z || Origin.Z > box.Max.Z)
                return false;
        }
        else
        {
            var tMinZ = (box.Min.Z - Origin.Z) / Direction.Z;
            var tMaxZ = (box.Max.Z - Origin.Z) / Direction.Z;

            if (tMinZ > tMaxZ) (tMinZ, tMaxZ) = (tMaxZ, tMinZ);

            if (tMin > tMaxZ || tMinZ > tMax)
                return false;

            if (tMinZ > tMin) tMin = tMinZ;
            if (tMaxZ < tMax) tMax = tMaxZ;
        }

        switch (tMin)
        {
            // having a positive tMax and a negative tMin means the ray is inside the box
            // we expect the intesection distance to be 0 in that case
            case < 0 when tMax > 0:
                return true;
            // a negative tMin means that the intersection point is behind the ray's origin
            // we discard these as not hitting the AABB
            case < 0:
                return false;
            default:
                distance = tMin;
                return true;
        }
    }
}