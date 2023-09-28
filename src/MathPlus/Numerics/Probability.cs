using System.Runtime.InteropServices;
using System.Runtime.Serialization;

namespace MathPlus;

/// <summary>A struct representing a probability (as a rational number)</summary>
[Serializable, DataContract]
[StructLayout(LayoutKind.Sequential)]
public struct Probability : IComparable<Probability>, IEquatable<Probability>, IFormattable
{
    public static readonly Rational Zero = new(0, 1);
    public static readonly Rational One = new(1, 1);

    /// <summary>The value of this probability</summary>
    [DataMember]
    public Rational Value;

    /// <summary>Creates a representation of probability using a rational number</summary>
    /// /// <param name="value">The probability value</param>
    public Probability(Rational value) => Value = value;

    /// <summary>Creates a representation of probability using a rational number</summary>
    /// <param name="num">The numerator of this probability</param>
    /// <param name="den">The denominator of this probability</param>
    public Probability(int num, int den) : this(new Rational(num, den))
    {
    }

    /// <summary>Randomly samples this probability, returning either true or false</summary>
    public bool Sample => Random.Shared.NextSingle(0, Value.Denominator) < Value.Numerator;

    public override string ToString() => ToString(null, null);

    public string ToString(string? format, IFormatProvider? formatProvider) =>
        $"{Value.ToString(format, formatProvider)} ({(float) Value * 100:#.#####}%)";

    // statics
    public static Probability operator &(Probability a, Probability b) =>
        new(a.Value * b.Value);

    public static Probability operator |(Probability a, Probability b) => !(!a & !b);

    public static Probability operator +(Probability a, Probability b) =>
        new(a.Value + b.Value);

    public static Probability operator !(Probability p) => new(1 - p.Value);

    // comparison operators
    public static bool operator ==(Probability a, Probability b) => a.Value == b.Value;
    public static bool operator !=(Probability a, Probability b) => a.Value != b.Value;
    public static bool operator <(Probability a, Probability b) => a.Value < b.Value;
    public static bool operator >(Probability a, Probability b) => a.Value > b.Value;
    public static bool operator <=(Probability a, Probability b) => a.Value <= b.Value;
    public static bool operator >=(Probability a, Probability b) => a.Value >= b.Value;

    // comparison functions
    public int CompareTo(Probability other) => Value.CompareTo(other.Value);
    public bool Equals(Probability other) => Value.Equals(other.Value);
    public override bool Equals(object? obj) => obj is Probability other && Equals(other);
    public override int GetHashCode() => Value.GetHashCode();
}