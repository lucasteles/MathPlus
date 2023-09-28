using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Serialization;
using System.Text;

namespace MathPlus;

/// <summary>An integer range</summary>
[Serializable, DataContract]
[StructLayout(LayoutKind.Sequential)]
public readonly struct IntRange : IEquatable<IntRange>, IFormattable
{
    public static readonly IntRange Empty = new(0, 0);
    public static readonly IntRange Unit = new(0, 1);

    [DataMember]
    public readonly int Start;

    [DataMember]
    public readonly int Count;

    public int this[int i] => Start + i;

    public Vector<int> this[Range range]
    {
        get
        {
            Span<int> values = stackalloc int[Count];
            for (var i = 0; i <= Count; i++)
                values[i] = i + Start;
            return new(values[range]);
        }
    }

    /// <summary>The last integer in the range</summary>
    public int Last => Start + Count - 1;

    /// <summary>The distance from first to last integer. Equivalent to <c>count-1</c></summary>
    public int Distance => Count - 1;

    /// <summary>Creates a new integer range, given a start integer and how many integers to include in total</summary>
    /// <param name="start">The first integer</param>
    /// <param name="count">How many integers to include in the full range</param>
    public IntRange(int start, int count)
    {
        Start = start;
        Count = count;
    }

    public static explicit operator IntRange((int, int) tuple) =>
        new(tuple.Item1, tuple.Item2 - tuple.Item1);

    public static implicit operator IntRange(Range range) =>
        (IntRange) range.ToExclusivePair();

    /// <summary>Whether or not this range contains a given value (inclusive)</summary>
    /// <param name="value">The value to check if it's inside, or equal to the start or end</param>
    public bool Contains(int value) => value >= Start && value <= Last;

    /// <summary>Create an integer range from start to end (inclusive)</summary>
    /// <param name="first">The first integer</param>
    /// <param name="last">The last integer</param>
    public static IntRange FirstToLast(int first, int last) =>
        new(first, last - first + 1);

    static readonly StringBuilder toStrBuilder = new();

    public override string ToString() => ToString(null, null);

    public string ToString(string? format, IFormatProvider? formatProvider)
    {
        toStrBuilder.Clear();
        toStrBuilder.Append("{ ");
        var last = Last;
        for (var i = Start; i <= last; i++)
        {
            if (format is null)
                toStrBuilder.Append(i);
            else
                toStrBuilder.Append(i.ToString(format, formatProvider));

            if (i != last)
                toStrBuilder.Append(", ");
        }

        toStrBuilder.Append(" }");
        return toStrBuilder.ToString();
    }

    public bool Equals(IntRange other) => Start.Equals(other.Start) && Count.Equals(other.Count);
    public override bool Equals(object? obj) => obj is IntRange other && Equals(other);
    public override int GetHashCode() => HashCode.Combine(Start, Count);
    public static bool operator ==(IntRange left, IntRange right) => left.Equals(right);
    public static bool operator !=(IntRange left, IntRange right) => !(left == right);

    public IntRangeEnumerator GetEnumerator() => new(this);

    public struct IntRangeEnumerator
    {
        readonly IntRange intRange;
        int currValue;

        public IntRangeEnumerator(IntRange range) =>
            (intRange, currValue) = (range, range.Start - 1);

        public bool MoveNext() => ++currValue <= intRange.Last;
        public int Current => currValue;
    }
}