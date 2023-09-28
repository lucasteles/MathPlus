using System.Globalization;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.CompilerServices;

namespace MathPlus;

[Serializable]
[StructLayout(LayoutKind.Sequential)]
public struct Vector3Int : IEquatable<Vector3Int>, IFormattable
{
    public int X;

    public int Y;

    public int Z;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Vector3Int(int x, int y)
    {
        X = x;
        Y = y;
        Z = 0;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Vector3Int(int x, int y, int z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    // Set x, y and z components of an existing Vector.
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Set(int x, int y, int z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public float Magnitude() => Mathf.Sqrt(X * X + Y * Y + Z * Z);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int SqrMagnitude() => X * X + Y * Y + Z * Z;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static float Distance(Vector3Int a, Vector3Int b) => (a - b).Magnitude();

    // Returns a vector that is made from the smallest components of two vectors.
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int Min(Vector3Int lhs, Vector3Int rhs) =>
        new(Mathf.Min(lhs.X, rhs.X), Mathf.Min(lhs.Y, rhs.Y),
            Mathf.Min(lhs.Z, rhs.Z));

    // Returns a vector that is made from the largest components of two vectors.
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int Max(Vector3Int lhs, Vector3Int rhs) =>
        new(Mathf.Max(lhs.X, rhs.X), Mathf.Max(lhs.Y, rhs.Y),
            Mathf.Max(lhs.Z, rhs.Z));

    // Multiplies two vectors component-wise.
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int Scale(Vector3Int a, Vector3Int b) =>
        new(a.X * b.X, a.Y * b.Y, a.Z * b.Z);

    // Multiplies every component of this vector by the same component of /scale/.
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Scale(Vector3Int scale)
    {
        X *= scale.X;
        Y *= scale.Y;
        Z *= scale.Z;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Clamp(Vector3Int min, Vector3Int max)
    {
        X = Math.Max(min.X, X);
        X = Math.Min(max.X, X);
        Y = Math.Max(min.Y, Y);
        Y = Math.Min(max.Y, Y);
        Z = Math.Max(min.Z, Z);
        Z = Math.Min(max.Z, Z);
    }

    // Converts a Vector3Int to a [[Vector3]].
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static implicit operator Vector3(Vector3Int v)
    {
        return new Vector3(v.X, v.Y, v.Z);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static explicit operator Vector2Int(Vector3Int v) => new(v.X, v.Y);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int FloorToInt(Vector3 v) =>
        new(
            Mathf.FloorToInt(v.X),
            Mathf.FloorToInt(v.Y),
            Mathf.FloorToInt(v.Z)
        );

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int CeilToInt(Vector3 v) =>
        new(
            Mathf.CeilToInt(v.X),
            Mathf.CeilToInt(v.Y),
            Mathf.CeilToInt(v.Z)
        );

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int RoundToInt(Vector3 v,
        MidpointRounding midpointRounding = MidpointRounding.ToEven) =>
        new(
            Mathf.RoundToInt(v.X, midpointRounding),
            Mathf.RoundToInt(v.Y, midpointRounding),
            Mathf.RoundToInt(v.Z, midpointRounding)
        );

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int operator +(Vector3Int a, Vector3Int b) =>
        new(a.X + b.X, a.Y + b.Y, a.Z + b.Z);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int operator -(Vector3Int a, Vector3Int b) =>
        new(a.X - b.X, a.Y - b.Y, a.Z - b.Z);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int operator *(Vector3Int a, Vector3Int b) =>
        new(a.X * b.X, a.Y * b.Y, a.Z * b.Z);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int operator -(Vector3Int a) => new(-a.X, -a.Y, -a.Z);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int operator *(Vector3Int a, int b) => new(a.X * b, a.Y * b, a.Z * b);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int operator *(int a, Vector3Int b) => new(a * b.X, a * b.Y, a * b.Z);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector3Int operator /(Vector3Int a, int b) => new(a.X / b, a.Y / b, a.Z / b);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool operator ==(Vector3Int lhs, Vector3Int rhs) =>
        lhs.X == rhs.X && lhs.Y == rhs.Y && lhs.Z == rhs.Z;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool operator !=(Vector3Int lhs, Vector3Int rhs) => !(lhs == rhs);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public override bool Equals(object? other) =>
        other is Vector3Int vector3Int && Equals(vector3Int);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Equals(Vector3Int other) => this == other;

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
        return
            $"({X.ToString(format, formatProvider)}, {Y.ToString(format, formatProvider)}, {Z.ToString(format, formatProvider)})";
    }

    public static readonly Vector3Int Zero = new(0, 0, 0);
    public static readonly Vector3Int One = new(1, 1, 1);
    public static readonly Vector3Int Up = new(0, 1, 0);
    public static readonly Vector3Int Down = new(0, -1, 0);
    public static readonly Vector3Int Left = new(-1, 0, 0);
    public static readonly Vector3Int Right = new(1, 0, 0);
    public static readonly Vector3Int Forward = new(0, 0, 1);
    public static readonly Vector3Int Back = new(0, 0, -1);
}