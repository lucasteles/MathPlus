using System.Globalization;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Serialization;

namespace MathPlus;

[Serializable, DataContract]
[StructLayout(LayoutKind.Sequential)]
public struct Rect : IEquatable<Rect>, IFormattable
{
    float xMin;
    float yMin;
    float width;
    float height;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Rect(float x, float y, float width, float height)
    {
        xMin = x;
        yMin = y;
        this.width = width;
        this.height = height;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Rect(Vector2 position, Vector2 size)
    {
        xMin = position.X;
        yMin = position.Y;
        width = size.X;
        height = size.Y;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Rect(Rect source)
    {
        xMin = source.xMin;
        yMin = source.yMin;
        width = source.width;
        height = source.height;
    }

    public static Rect Zero => new(0.0f, 0.0f, 0.0f, 0.0f);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Rect MinMaxRect(float xmin, float ymin, float xmax, float ymax) =>
        new(xmin, ymin, xmax - xmin, ymax - ymin);

    public float X
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => xMin;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => xMin = value;
    }

    public float Y
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => yMin;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => yMin = value;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Set(float x, float y, float w, float h)
    {
        xMin = x;
        yMin = y;
        width = w;
        height = h;
    }

    public Vector2 Position
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(xMin, yMin);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set
        {
            xMin = value.X;
            yMin = value.Y;
        }
    }

    public Vector2 Center
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(X + width / 2f, Y + height / 2f);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set
        {
            xMin = value.X - width / 2f;
            yMin = value.Y - height / 2f;
        }
    }

    public Vector2 Min
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(XMin, YMin);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set
        {
            XMin = value.X;
            YMin = value.Y;
        }
    }

    public Vector2 Max
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(XMax, YMax);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set
        {
            XMax = value.X;
            YMax = value.Y;
        }
    }

    public float Width
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => width;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => width = value;
    }

    public float Height
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => height;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => height = value;
    }

    public Vector2 Size
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(width, height);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set
        {
            width = value.X;
            height = value.Y;
        }
    }

    public float XMin
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => xMin;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set
        {
            var old = XMax;
            xMin = value;
            width = old - xMin;
        }
    }

    public float YMin
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => yMin;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set
        {
            var old = YMax;
            yMin = value;
            height = old - yMin;
        }
    }

    public float XMax
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => width + xMin;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => width = value - xMin;
    }

    public float YMax
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => height + yMin;
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => height = value - yMin;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Contains(Vector2 point) =>
        point.X >= XMin && point.X < XMax && point.Y >= YMin && point.Y < YMax;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Contains(Vector3 point) =>
        point.X >= XMin && point.X < XMax && point.Y >= YMin && point.Y < YMax;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Contains(Vector3 point, bool allowInverse)
    {
        if (!allowInverse) return Contains(point);

        var xAxis = Width < 0f && point.X <= XMin && point.X > XMax
                    || Width >= 0f && point.X >= XMin && point.X < XMax;

        var yAxis = Height < 0f && point.Y <= YMin && point.Y > YMax
                    || Height >= 0f && point.Y >= YMin && point.Y < YMax;

        return xAxis && yAxis;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    static Rect OrderMinMax(Rect rect)
    {
        if (rect.XMin > rect.XMax)
            (rect.XMin, rect.XMax) = (rect.XMax, rect.XMin);

        if (rect.YMin > rect.YMax)
            (rect.YMin, rect.YMax) = (rect.YMax, rect.YMin);

        return rect;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Overlaps(Rect other) =>
        other.XMax > XMin &&
        other.XMin < XMax &&
        other.YMax > YMin &&
        other.YMin < YMax;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Overlaps(Rect other, bool allowInverse)
    {
        if (!allowInverse) return this.Overlaps(other);

        var ordered = OrderMinMax(this);
        other = OrderMinMax(other);

        return ordered.Overlaps(other);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2 NormalizedToPoint(Rect rectangle, Vector2 normalizedRectCoordinates) =>
        new(
            Mathf.Lerp(rectangle.X, rectangle.XMax, normalizedRectCoordinates.X),
            Mathf.Lerp(rectangle.Y, rectangle.YMax, normalizedRectCoordinates.Y)
        );

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector2 PointToNormalized(Rect rectangle, Vector2 point) =>
        new(
            Mathf.InverseLerp(rectangle.X, rectangle.XMax, point.X),
            Mathf.InverseLerp(rectangle.Y, rectangle.YMax, point.Y)
        );

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool operator !=(Rect lhs, Rect rhs) => !(lhs == rhs);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool operator ==(Rect lhs, Rect rhs) =>
        lhs.X == rhs.X && lhs.Y == rhs.Y && lhs.Width == rhs.Width &&
        lhs.Height == rhs.Height;

    public override int GetHashCode() => HashCode.Combine(X, Width, Y, Height);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public override bool Equals(object? other) => other is Rect rect && Equals(rect);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Equals(Rect other) =>
        X.Equals(other.X) && Y.Equals(other.Y) && Width.Equals(other.Width) &&
        Height.Equals(other.Height);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public override string ToString() => ToString(null, null);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public string ToString(string format) => ToString(format, null);

    public string ToString(string? format, IFormatProvider? formatProvider)
    {
        if (string.IsNullOrEmpty(format)) format = "F2";

        formatProvider ??= CultureInfo.InvariantCulture.NumberFormat;

        return
            $"(x:{X.ToString(format, formatProvider)}, y:{Y.ToString(format, formatProvider)}, width:{Width.ToString(format, formatProvider)}, height:{Height.ToString(format, formatProvider)})";
    }

    /// <summary>Expands the rectangle to encapsulate the point <c>p</c></summary>
    /// <param name="r">The rectangle to expand</param>
    /// <param name="p">The point to encapsulate</param>
    public Rect Encapsulate(Vector2 p)
    {
        Rect r = new(this);
        r.XMax = MathF.Max(r.XMax, p.X);
        r.XMin = MathF.Min(r.XMin, p.X);
        r.XMax = MathF.Max(r.XMax, p.Y);
        r.YMin = MathF.Min(r.YMin, p.Y);
        return r;
    }

    /// <summary>Interpolates a position within this rectangle, given a normalized position</summary>
    /// <param name="tPos">The normalized position within this rectangle</param>
    public Vector2 Lerp(Vector2 tPos) =>
        new(
            Mathf.Lerp(XMin, XMax, tPos.X),
            Mathf.Lerp(YMin, YMax, tPos.Y)
        );

    /// <summary>The x axis range of this rectangle</summary>
    public FloatRange RangeX() => new(XMin, XMax);

    /// <summary>The y axis range of this rectangle</summary>
    public FloatRange RangeY() => new(YMin, YMax);

    /// <summary>Places the center of this rectangle at its position,
    /// useful together with the constructor to define it by center instead of by corner</summary>
    public Rect ByCenter() => new(this) {Center = Position};
}