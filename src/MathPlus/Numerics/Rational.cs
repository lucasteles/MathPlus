// by Freya Holmér (https://github.com/FreyaHolmer/MathPlus)

using System.Runtime.InteropServices;
using System.Runtime.Serialization;

namespace MathPlus;

/// <summary>A struct representing exact rational numbers</summary>
[Serializable, DataContract]
[StructLayout(LayoutKind.Sequential)]
public readonly struct Rational : IEquatable<Rational>, IComparable<Rational>, IFormattable
{
    public static readonly Rational Zero = new(0, 1);
    public static readonly Rational One = new(1, 1);
    public static readonly Rational MaxValue = new(int.MaxValue, 1);
    public static readonly Rational MinValue = new(int.MinValue, 1);

    /// <summary>The numerator of this number</summary>
    [DataMember]
    public readonly int Numerator;

    /// <summary>The denominator of this number</summary>
    [DataMember]
    public readonly int Denominator;

    /// <summary>Creates an exact representation of a rational number</summary>
    /// <param name="num">The numerator of this number</param>
    /// <param name="den">The denominator of this number</param>
    public Rational(int num, int den)
    {
        switch (den)
        {
            case -1:
                (Numerator, Denominator) = (-num, -den);
                break;
            case 0: throw new DivideByZeroException("The denominator can't be 0");
            case 1:
                (Numerator, Denominator) = (num, den);
                break;
            default:
                if (num == 0)
                {
                    (Numerator, Denominator) = (0, 1);
                    break;
                }

                // ensure only the numerator carries the sign
                var sign = Mathf.Sign(den);
                Numerator = sign * num;
                Denominator = sign * den;

                if (Numerator is -1 or 1)
                    break; // no reduction needed

                // in this case, we have to try simplifying the expression
                var gcd = Mathf.Gcd(num, den);
                Numerator /= gcd;
                Denominator /= gcd;
                break;
        }
    }

    public void Deconstruct(out int num, out int den) =>
        (num, den) = (Numerator, Denominator);

    /// <summary>Returns the reciprocal of this number</summary>
    public Rational Reciprocal => new(Denominator, Numerator);

    public bool IsInteger => Denominator == 1;

    /// <summary>Returns the absolute value of this number</summary>
    public Rational Abs() => new(Numerator.Abs(), Denominator);

    /// <summary>Returns this number to the power of another integer <c>pow</c></summary>
    /// <param name="pow">The power to raise this number by</param>
    public Rational Pow(int pow) =>
        pow switch
        {
            <= -2 => Reciprocal.Pow(-pow),
            -1 => Reciprocal,
            0 => 1,
            1 => this,
            >= 2 => new Rational(Numerator.Pow(pow), Denominator.Pow(pow))
        };

    public string ToString(string? format = null, IFormatProvider? formatProvider = null) =>
        Denominator == 1
            ? Numerator.ToString(format, formatProvider)
            : $"{Numerator.ToString(format, formatProvider)}/{Denominator.ToString(format, formatProvider)}";

    // statics
    public static Rational Min(Rational a, Rational b) => a < b ? a : b;
    public static Rational Max(Rational a, Rational b) => a > b ? a : b;
    public static Rational Lerp(Rational a, Rational b, Rational t) => a + t * (b - a);
    public static Rational InverseLerp(Rational a, Rational b, Rational v) => (v - a) / (b - a);

    // type casting
    public static implicit operator Rational(int n) => new(n, 1);

    public static explicit operator int(Rational r) => r.IsInteger
        ? r.Numerator
        : throw new ArithmeticException($"Rational value {r} can't be cast to an integer");

    public static explicit operator float(Rational r) => (float) r.Numerator / r.Denominator;
    public static explicit operator double(Rational r) => (double) r.Numerator / r.Denominator;

    // unary operations
    public static Rational operator -(Rational r) => checked(new(-r.Numerator, r.Denominator));
    public static Rational operator +(Rational r) => r;

    // addition
    public static Rational operator +(Rational a, Rational b) =>
        checked(new(a.Numerator * b.Denominator + a.Denominator * b.Numerator,
            a.Denominator * b.Denominator));

    public static Rational operator +(Rational a, int b) =>
        checked(new(a.Numerator + a.Denominator * b, a.Denominator));

    public static Rational operator +(int a, Rational b) =>
        checked(new(a * b.Denominator + b.Numerator, b.Denominator));

    // subtraction
    public static Rational operator -(Rational a, Rational b) =>
        checked(new(a.Numerator * b.Denominator - a.Denominator * b.Numerator,
            a.Denominator * b.Denominator));

    public static Rational operator -(Rational a, int b) =>
        checked(new(a.Numerator - a.Denominator * b, a.Denominator));

    public static Rational operator -(int a, Rational b) =>
        checked(new(a * b.Denominator - b.Numerator, b.Denominator));

    // multiplication
    public static Rational operator *(Rational a, Rational b) =>
        checked(new(a.Numerator * b.Numerator, a.Denominator * b.Denominator));

    public static Rational operator *(Rational a, int b) =>
        checked(new(a.Numerator * b, a.Denominator));

    public static Rational operator *(int a, Rational b) =>
        checked(new(b.Numerator * a, b.Denominator));

    public static float operator *(Rational a, float b) => (a.Numerator * b) / a.Denominator;
    public static float operator *(float a, Rational b) => (b.Numerator * a) / b.Denominator;

    // division
    public static Rational operator /(Rational a, Rational b) =>
        checked(new(a.Numerator * b.Denominator, a.Denominator * b.Numerator));

    public static Rational operator /(Rational a, int b) =>
        checked(new(a.Numerator, a.Denominator * b));

    public static Rational operator /(int a, Rational b) =>
        checked(new(a * b.Denominator, b.Numerator));

    public static float operator /(Rational a, float b) => a.Numerator / (a.Denominator * b);
    public static float operator /(float a, Rational b) => (a * b.Denominator) / b.Numerator;

    // comparison operators
    public static bool operator ==(Rational a, Rational b) => a.CompareTo(b) == 0;
    public static bool operator !=(Rational a, Rational b) => a.CompareTo(b) != 0;
    public static bool operator <(Rational a, Rational b) => a.CompareTo(b) < 0;
    public static bool operator >(Rational a, Rational b) => a.CompareTo(b) > 0;
    public static bool operator <=(Rational a, Rational b) => a.CompareTo(b) <= 0;
    public static bool operator >=(Rational a, Rational b) => a.CompareTo(b) >= 0;

    // comparison functions
    public int CompareTo(Rational other) =>
        checked((Numerator * other.Denominator).CompareTo(Denominator * other.Numerator));

    public bool Equals(Rational other) =>
        Numerator == other.Numerator && Denominator == other.Denominator;

    public override bool Equals(object? obj) => obj is Rational other && Equals(other);
    public override int GetHashCode() => HashCode.Combine(Numerator, Denominator);
}