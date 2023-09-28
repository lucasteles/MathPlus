using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Serialization;

namespace MathPlus;

[Serializable, DataContract]
[StructLayout(LayoutKind.Sequential)]
public struct Pose : IEquatable<Pose>
{
    [DataMember]
    public Vector3 Position;

    [DataMember]
    public Quaternion Rotation;

    public Pose(Vector3 position, Quaternion rotation)
    {
        Position = position;
        Rotation = rotation;
    }

    public override string ToString() => $"({Position.ToString()}, {Rotation.ToString()})";

    public Pose GetTransformedBy(Pose lhs) =>
        new()
        {
            Position = lhs.Position + (lhs.Rotation.Mul(Position)),
            Rotation = lhs.Rotation * Rotation,
        };

    public Vector3 Forward => Rotation.Mul(Vector3V.Forward);

    public Vector3 Right => Rotation.Mul(Vector3V.Right);

    public Vector3 Up => Rotation.Mul(Vector3V.Up);

    public static Pose Identity { get; } = new(Vector3.Zero, Quaternion.Identity);

    public override bool Equals(object? obj) => obj is Pose pose && Equals(pose);

    public bool Equals(Pose other) =>
        Position == other.Position &&
        Rotation == other.Rotation;

    public override int GetHashCode() => HashCode.Combine(Position, Rotation);

    public static bool operator ==(Pose a, Pose b) => a.Equals(b);

    public static bool operator !=(Pose a, Pose b) => !(a == b);
}