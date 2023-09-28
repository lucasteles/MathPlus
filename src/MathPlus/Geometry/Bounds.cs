using System.Numerics;

namespace MathPlus;

using System;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using scm = System.ComponentModel;

[Serializable]
[StructLayout(LayoutKind.Sequential)]
public struct Bounds : IEquatable<Bounds>, IFormattable
{
    public Vector3 Center;
    public Vector3 Extents;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Bounds(Vector3 center, Vector3 size)
    {
        Center = center;
        Extents = size * 0.5F;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public override int GetHashCode() => HashCode.Combine(Center, Extents);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public override bool Equals(object? other) => other is Bounds bounds && Equals(bounds);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Equals(Bounds other) =>
        Center.Equals(other.Center) && Extents.Equals(other.Extents);

    public Vector3 Size
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Extents * 2.0F;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => Extents = value * 0.5F;
    }

    public Vector3 Min
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Center - Extents;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => SetMinMax(value, Max);
    }

    public Vector3 Max
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Center + Extents;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => SetMinMax(Min, value);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool operator ==(Bounds a, Bounds b) =>
        a.Center == b.Center && a.Extents == b.Extents;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool operator !=(Bounds a, Bounds b) => !(a == b);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void SetMinMax(Vector3 min, Vector3 max)
    {
        Extents = (max - min) * 0.5F;
        Center = min + Extents;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Encapsulate(Vector3 point) =>
        SetMinMax(Vector3.Min(Min, point), Vector3.Max(Max, point));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Encapsulate(Bounds bounds)
    {
        Encapsulate(bounds.Center - bounds.Extents);
        Encapsulate(bounds.Center + bounds.Extents);
    }

    public void Expand(float amount)
    {
        amount *= .5f;
        Extents += new Vector3(amount, amount, amount);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Expand(Vector3 amount) => Extents += amount * .5f;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Intersects(Bounds bounds) =>
        Min.X <= bounds.Max.X && Max.X >= bounds.Min.X &&
        Min.Y <= bounds.Max.Y && Max.Y >= bounds.Min.Y &&
        Min.Z <= bounds.Max.Z && Max.Z >= bounds.Min.Z;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool IntersectRay(Ray ray) => ray.Intersects(this, out _);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool IntersectRay(Ray ray, out float distance) => ray.Intersects(this, out distance);

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
            $"Center: {Center.ToString(format, formatProvider)}, Extents: {Extents.ToString(format, formatProvider)}";
    }
}