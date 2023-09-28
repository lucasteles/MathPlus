using System.Runtime.InteropServices;
using System.Runtime.Serialization;

namespace MathPlus;

/// <summary>A value range between two values</summary>
[Serializable, DataContract]
[StructLayout(LayoutKind.Sequential)]
public struct FloatRange : IEquatable<FloatRange>, IFormattable
{
    /// <summary>The unit interval of 0 to 1</summary>
    public static readonly FloatRange Unit = new(0, 1);

    /// <summary>The start of this range</summary>
    [DataMember]
    public float Start;

    /// <summary>The end of this range</summary>
    [DataMember]
    public float End;

    /// <summary>Creates a new value range</summary>
    /// <param name="a">The start of the range</param>
    /// <param name="b">The end of the range</param>
    public FloatRange(float a, float b) => (Start, End) = (a, b);

    /// <summary>The value at the center of this value range</summary>
    public float Center => (Start + End) / 2;

    /// <summary>The length/span of this value range</summary>
    public float Length => MathF.Abs(End - Start);

    /// <summary>The minimum value of this range</summary>
    public float Min => MathF.Min(Start, End);

    /// <summary>The maximum value of this range</summary>
    public float Max => MathF.Max(Start, End);

    /// <summary>The direction of this value range. Returns -1 if <c>b</c> is greater than <c>a</c>, otherwise returns 1</summary>
    public int Direction => End > Start ? 1 : -1;

    public void Deconstruct(out float start, out float end) =>
        (start, end) = (Start, End);

    /// <summary>Interpolates a value from <c>a</c> to <c>b</c>, based on a parameter <c>t</c></summary>
    /// <param name="t">The normalized interpolant from <c>a</c> to <c>b</c>. A value of 0 returns <c>a</c>, a value of 1 returns <c>b</c></param>
    public float Lerp(float t) => Mathf.Lerp(Start, End, t);

    /// <summary>Returns the normalized position of the input value <c>v</c> within this range</summary>
    /// <param name="v">The value to get the normalized position of</param>
    public float InverseLerp(float v) => Mathf.InverseLerp(Start, End, v);

    /// <summary>Returns whether or not this range contains the value <c>v</c> (inclusive)</summary>
    /// <param name="v">The value to see if it's inside</param>
    public bool Contains(float v) => v >= MathF.Min(Start, End) && v <= MathF.Max(Start, End);

    /// <summary>Returns whether or not this range contains the range <c>r</c></summary>
    /// <param name="r">The range to see if it's inside</param>
    public bool Contains(FloatRange r) => r.Min >= Min && r.Max <= Max;

    /// <summary>Remaps the input value from the <c>input</c> range to the <c>output</c> range</summary>
    /// <param name="value">The value to remap</param>
    /// <param name="input">The input range</param>
    /// <param name="output">The output range</param>
    public static float Remap(float value, FloatRange input, FloatRange output) =>
        output.Lerp(input.InverseLerp(value));

    /// <summary>Remaps a range from the <c>input</c> range to the <c>output</c> range</summary>
    /// <param name="value">The range to remap</param>
    /// <param name="input">The input range</param>
    /// <param name="output">The output range</param>
    public static FloatRange Remap(FloatRange value, FloatRange input, FloatRange output) =>
        new(Remap(value.Start, input, output), Remap(value.End, input, output));

    /// <summary>Returns whether or not this range overlaps another range</summary>
    /// <param name="other">The other range to test overlap with</param>
    public bool Overlaps(FloatRange other)
    {
        var separation = MathF.Abs(other.Center - Center);
        var rTotal = (other.Length + Length) / 2;
        return separation < rTotal;
    }

    /// <summary>Wraps/repeats the input value to stay within this range</summary>
    /// <param name="value">The value to wrap/repeat in this interval</param>
    public float Wrap(float value)
    {
        if (value >= Start && value < End) return value;
        return Start + Mathf.Repeat(value - Start, End - Start);
    }

    /// <summary>Clamps the input value to this range</summary>
    /// <param name="value">The value to clamp to this interval</param>
    public float Clamp(float value) => Mathf.Clamp(value, Min, Max);

    /// <summary>Expands the minimum or maximum value to contain the given <c>value</c></summary>
    /// <param name="value">The value to include</param>
    public FloatRange Encapsulate(float value) =>
        Direction switch
        {
            1 => new(Mathf.Min(Start, value),
                Mathf.Max(End, value)), // forward - a is min, b is max
            _ => new(Mathf.Min(End, value),
                Mathf.Max(Start, value)) // reversed - b is min, a is max
        };

    /// <summary>Expands the minimum or maximum value to contain the given <c>range</c></summary>
    /// <param name="range">The value range to include</param>
    public FloatRange Encapsulate(FloatRange range) =>
        (Direction switch
        {
            1 => new(Mathf.Min(Start, range.Start),
                Mathf.Max(End, range.End)), // forward - a is min, b is max
            _ => new(Mathf.Min(End, range.End),
                Mathf.Max(Start, range.Start)), // reversed - b is min, a is max
        });

    /// <summary>Returns a version of this range, scaled around its start value</summary>
    /// <param name="scale">The value to scale the range by</param>
    public FloatRange ScaleFromStart(float scale) =>
        new FloatRange(Start, Start + scale * (End - Start));

    /// <summary>Returns this range mirrored around a given value</summary>
    /// <param name="pivot">The value to mirror around</param>
    public FloatRange MirrorAround(float pivot) =>
        new(2 * pivot - Start, 2 * pivot - End);

    /// <summary>Returns a reversed version of this range, where a and b is swapped</summary>
    public FloatRange Reverse() => new(End, Start);

    public static FloatRange operator -(FloatRange range, float v) =>
        new(range.Start - v, range.End - v);

    public static FloatRange operator +(FloatRange range, float v) =>
        new(range.Start + v, range.End + v);

    public static FloatRange operator /(FloatRange range, int v) =>
        new(range.Start / v, range.End / v);

    public static FloatRange operator /(FloatRange range, float v) =>
        new(range.Start / v, range.End / v);

    public static FloatRange operator *(FloatRange range, int v) =>
        new(range.Start * v, range.End * v);

    public static FloatRange operator *(FloatRange range, float v) =>
        new(range.Start * v, range.End * v);

    public static explicit operator FloatRange((float, float) tuple) =>
        new(tuple.Item1, tuple.Item2);

    public static explicit operator FloatRange((int, int) tuple) =>
        new(tuple.Item1, tuple.Item2);

    public static implicit operator FloatRange(Range range) =>
        (FloatRange) range.ToExclusivePair();

    public static bool operator ==(FloatRange a, FloatRange b) =>
        a.Start.Approximately(b.Start) && a.End.Approximately(b.End);

    public static bool operator !=(FloatRange a, FloatRange b) => !(a == b);

    public bool Equals(FloatRange other) => Start.Equals(other.Start) && End.Equals(other.End);
    public override bool Equals(object? obj) => obj is FloatRange other && Equals(other);
    public override int GetHashCode() => HashCode.Combine(Start, End);

    public override string ToString() => ToString(null, null);

    public string ToString(string? format, IFormatProvider? formatProvider) =>
        $"[{Start.ToString(format, formatProvider)},{End.ToString(format, formatProvider)}]";

    public FloatEnumerable Enumerate(float step = 1) => new(this, step);

    public readonly struct FloatEnumerable
    {
        readonly FloatRange range;
        readonly float step;

        public FloatEnumerable(in FloatRange range, float step)
        {
            this.range = range;
            this.step = step;
        }
    }

    public struct FloatRangeEnumerator
    {
        readonly float step;
        readonly FloatRange range;
        float currValue;

        public FloatRangeEnumerator(in FloatRange range, float step)
        {
            this.step = step;
            this.range = range;
            currValue = range.Start - step;
        }

        public bool MoveNext()
        {
            currValue += step;
            return currValue <= range.End;
        }

        public float Current => currValue;
    }
}