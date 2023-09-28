using System.Numerics;

namespace MathPlus;

using System;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;


[Serializable]
[StructLayout(LayoutKind.Sequential)]
public struct Vector2Int : IEquatable<Vector2Int>, IFormattable
{
    public int X;
    public int Y;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Vector2Int(int x, int y)
    {
        X = x;
        Y = y;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Set(int x, int y)
    {
        X = x;
        Y = y;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public float Magnitude() => Mathf.Sqrt(X * X + Y * Y);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int SqrMagnitude() => X * X + Y * Y;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static float Distance(Vector2Int a, Vector2Int b)
    {
        float diffX = a.X - b.X;
        float diffY = a.Y - b.Y;

        return (float) Math.Sqrt(diffX * diffX + diffY * diffY);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int Min(Vector2Int lhs, Vector2Int rhs) =>
        new(Mathf.Min(lhs.X, rhs.X), Mathf.Min(lhs.Y, rhs.Y));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int Max(Vector2Int lhs, Vector2Int rhs) =>
        new(Mathf.Max(lhs.X, rhs.X), Mathf.Max(lhs.Y, rhs.Y));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int Scale(Vector2Int a, Vector2Int b) => new(a.X * b.X, a.Y * b.Y);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Scale(Vector2Int scale)
    {
        X *= scale.X;
        Y *= scale.Y;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Clamp(Vector2Int min, Vector2Int max)
    {
        X = Math.Max(min.X, X);
        X = Math.Min(max.X, X);
        Y = Math.Max(min.Y, Y);
        Y = Math.Min(max.Y, Y);
    }

    // Converts a Vector2Int to a [[Vector2]].
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static implicit operator Vector2(Vector2Int v) => new(v.X, v.Y);

    // Converts a Vector2Int to a [[Vector3Int]].
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static explicit operator Vector3Int(Vector2Int v)
    {
        return new Vector3Int(v.X, v.Y, 0);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int FloorToInt(Vector2 v) =>
        new(
            Mathf.FloorToInt(v.X),
            Mathf.FloorToInt(v.Y)
        );

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int CeilToInt(Vector2 v) =>
        new(
            Mathf.CeilToInt(v.X),
            Mathf.CeilToInt(v.Y)
        );

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int RoundToInt(Vector2 v) =>
        new(
            Mathf.RoundToInt(v.X),
            Mathf.RoundToInt(v.Y)
        );

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int operator -(Vector2Int v) => new(-v.X, -v.Y);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int operator +(Vector2Int a, Vector2Int b) => new(a.X + b.X, a.Y + b.Y);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int operator -(Vector2Int a, Vector2Int b) => new(a.X - b.X, a.Y - b.Y);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int operator *(Vector2Int a, Vector2Int b) => new(a.X * b.X, a.Y * b.Y);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int operator *(int a, Vector2Int b) => new(a * b.X, a * b.Y);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int operator *(Vector2Int a, int b) => new(a.X * b, a.Y * b);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2Int operator /(Vector2Int a, int b) => new(a.X / b, a.Y / b);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool operator ==(Vector2Int lhs, Vector2Int rhs) =>
        lhs.X == rhs.X && lhs.Y == rhs.Y;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool operator !=(Vector2Int lhs, Vector2Int rhs) => !(lhs == rhs);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public override bool Equals(object? other) =>
        other is Vector2Int vector2Int && Equals(vector2Int);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Equals(Vector2Int other) => X == other.X && Y == other.Y;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public override int GetHashCode() => HashCode.Combine(X, Y);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public override string ToString() => ToString(null, null);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public string ToString(string format) => ToString(format, null);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public string ToString(string? format, IFormatProvider? formatProvider)
    {
        formatProvider ??= CultureInfo.InvariantCulture.NumberFormat;
        return $"({X.ToString(format, formatProvider)}, {Y.ToString(format, formatProvider)})";
    }

    public static readonly Vector2Int Zero = new(0, 0);
    public static readonly Vector2Int One = new(1, 1);
    public static readonly Vector2Int Up = new(0, 1);
    public static readonly Vector2Int Down = new(0, -1);
    public static readonly Vector2Int Left = new(-1, 0);
    public static readonly Vector2Int Right = new(1, 0);
}