// by Freya Holmér (https://github.com/FreyaHolmer/MathPlus)

using System.Runtime.CompilerServices;
using System.Numerics;

namespace MathPlus;

/// <summary>Various extensions for floats, vectors and colors</summary>
public static class MathPlusExtensions
{
    const MethodImplOptions Inline = MethodImplOptions.AggressiveInlining;

    #region Vector rotation and angles

    /// <summary>Returns the angle of this vector, in radians</summary>
    /// <param name="v">The vector to get the angle of. It does not have to be normalized</param>
    /// <seealso cref="Mathf.DirToAng"/>
    [MethodImpl(Inline)] public static float Angle(this Vector2 v) => MathF.Atan2(v.Y, v.X);

    /// <summary>Rotates the vector 90 degrees clockwise (negative Z axis rotation)</summary>
    [MethodImpl(Inline)] public static Vector2 Rotate90CW(this Vector2 v) =>
        new(v.Y, -v.X);

    /// <summary>Rotates the vector 90 degrees counter-clockwise (positive Z axis rotation)</summary>
    [MethodImpl(Inline)] public static Vector2 Rotate90CCW(this Vector2 v) =>
        new(-v.Y, v.X);

    /// <summary>Rotates the vector around <c>pivot</c> with the given angle (in radians)</summary>
    /// <param name="v">The vector to rotate</param>
    /// <param name="pivot">The point to rotate around</param>
    /// <param name="angRad">The angle to rotate by, in radians</param>
    [MethodImpl(Inline)]
    public static Vector2 RotateAround(this Vector2 v, Vector2 pivot, float angRad) =>
        Rotate(v - pivot, angRad) + pivot;

    /// <summary>Rotates the vector around <c>(0,0)</c> with the given angle (in radians)</summary>
    /// <param name="v">The vector to rotate</param>
    /// <param name="angRad">The angle to rotate by, in radians</param>
    public static Vector2 Rotate(this Vector2 v, float angRad)
    {
        var ca = MathF.Cos(angRad);
        var sa = MathF.Sin(angRad);
        return new Vector2(ca * v.X - sa * v.Y, sa * v.X + ca * v.Y);
    }

    public static Vector3 ToVector3(this in Vector2 v) => new(v.X, v.Y, 0);
    public static Vector3 ToVector3(this in Vector4 v) => new(v.X, v.Y, v.Z);
    public static Vector2 Normalized(this in Vector2 v) => Vector2.Normalize(v);
    public static Vector3 Normalized(this in Vector3 v) => Vector3.Normalize(v);
    public static Vector4 Normalized(this in Vector4 v) => Vector4.Normalize(v);


    /// <summary>Converts an angle in degrees to radians</summary>
    /// <param name="angDegrees">The angle, in degrees, to convert to radians</param>
    [MethodImpl(Inline)] public static float DegToRad(this float angDegrees) =>
        angDegrees * Mathf.Deg2Rad;

    /// <summary>Converts an angle in radians to degrees</summary>
    /// <param name="angRadians">The angle, in radians, to convert to degrees</param>
    [MethodImpl(Inline)] public static float RadToDeg(this float angRadians) =>
        angRadians * Mathf.Rad2Deg;

    /// <summary>Extracts the quaternion components into a Vector4</summary>
    /// <param name="q">The quaternion to get the components of</param>
    [MethodImpl(Inline)] public static Vector4 ToVector4(this Quaternion q) =>
        new(q.X, q.Y, q.Z, q.W);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Vector3 Axis, float Angle) ToAxisAngle(this Quaternion q)
    {
        // https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/
        // TODO: is this right?

        if (q.W > 1)
            q = Quaternion.Normalize(q);

        var angle = 2 * Mathf.Acos(q.W);

        var s = Mathf.Sqrt(1 - q.LengthSquared());
        Vector3 axis = s < 0.001
            ? new(q.X, q.Y, q.Z)
            : new(q.X / s, q.Y / s, q.Z / s);

        return (axis, angle);
    }

    /// <summary>Converts to a rotation vector (axis-angle where the angle is embedded in the magnitude, in radians)</summary>
    /// <param name="q">The quaternion to get the rotation vector of</param>
    [MethodImpl(Inline)] public static Vector3 ToRotationVector(this Quaternion q)
    {
        var (axis, angDeg) = q.ToAxisAngle();
        return axis * angDeg;
    }

    #endregion

    #region Swizzling

    /// <summary>Returns X and Y as a Vector2, equivalent to <c>new Vector2(v.X,v.y)</c></summary>
    [MethodImpl(Inline)] public static Vector2 ToXY(in this Vector2 v) => new(v.X, v.Y);

    /// <summary>Returns Y and X as a Vector2, equivalent to <c>new Vector2(v.y,v.X)</c></summary>
    [MethodImpl(Inline)] public static Vector2 ToYX(in this Vector2 v) => new(v.Y, v.X);


    /// <summary>Returns X and Z as a Vector2, equivalent to <c>new Vector2(v.X,v.z)</c></summary>
    [MethodImpl(Inline)] public static Vector2 ToXY(in this Vector3 v) => new(v.X, v.Y);

    [MethodImpl(Inline)] public static Vector2 ToXZ(in this Vector3 v) => new(v.X, v.Z);
    [MethodImpl(Inline)] public static Vector2 ToYX(in this Vector3 v) => new(v.Y, v.X);
    [MethodImpl(Inline)] public static Vector2 ToYZ(in this Vector3 v) => new(v.Y, v.Z);
    [MethodImpl(Inline)] public static Vector2 ToZX(in this Vector3 v) => new(v.Z, v.X);
    [MethodImpl(Inline)] public static Vector2 ToZY(in this Vector3 v) => new(v.Z, v.Y);


    /// <summary>Returns this vector as a Vector3, slotting X into X, and Y into Z, and the input value y into Y.
    /// Equivalent to <c>new Vector3(v.X,y,v.y)</c></summary>
    [MethodImpl(Inline)] public static Vector3 XZtoXYZ(this Vector2 v, float y = 0) =>
        new(v.X, y, v.Y);

    /// <summary>Returns this vector as a Vector3, slotting X into X, and Y into Y, and the input value z into Z.
    /// Equivalent to <c>new Vector3(v.X,v.y,z)</c></summary>
    [MethodImpl(Inline)] public static Vector3 XYtoXYZ(this Vector2 v, float z = 0) =>
        new(v.X, v.Y, z);

    /// <summary>Sets X to 0</summary>
    [MethodImpl(Inline)] public static Vector2 FlattenX(this Vector2 v) => v with
    {
        X = 0f,
    };

    /// <summary>Sets Y to 0</summary>
    [MethodImpl(Inline)] public static Vector2 FlattenY(this Vector2 v) => v with
    {
        Y = 0f,
    };

    /// <summary>Sets X to 0</summary>
    [MethodImpl(Inline)] public static Vector3 FlattenX(this Vector3 v) =>
        v with
        {
            X = 0f,
        };

    /// <summary>Sets Y to 0</summary>
    [MethodImpl(Inline)] public static Vector3 FlattenY(this Vector3 v) =>
        v with
        {
            Y = 0f
        };

    /// <summary>Sets Z to 0</summary>
    [MethodImpl(Inline)] public static Vector3 FlattenZ(this Vector3 v) =>
        v with
        {
            Z = 0f,
        };

    #endregion

    #region Vector directions & magnitudes

    /// <summary>Returns the chebyshev magnitude of this vector</summary>
    [MethodImpl(Inline)] public static float ChebyshevMagnitude(this Vector3 v) =>
        Mathf.Max(Abs(v.X), Abs(v.Y), Abs(v.Z));

    /// <summary>Returns the taxicab/rectilinear magnitude of this vector</summary>
    [MethodImpl(Inline)] public static float TaxicabMagnitude(this Vector3 v) =>
        Abs(v.X) + Abs(v.Y) + Abs(v.Z);

    /// <inheritdoc cref="ChebyshevMagnitude(Vector3)"/>
    [MethodImpl(Inline)] public static float ChebyshevMagnitude(this Vector2 v) =>
        Mathf.Max(Abs(v.X), Abs(v.Y));

    /// <inheritdoc cref="TaxicabMagnitude(Vector3)"/>
    [MethodImpl(Inline)] public static float TaxicabMagnitude(this Vector2 v) =>
        Abs(v.X) + Abs(v.Y);

    /// <summary>Returns a vector with the same direction, but with the given magnitude.
    /// Equivalent to <c>v.Normalized()*mag</c></summary>
    [MethodImpl(Inline)] public static Vector2 WithMagnitude(this Vector2 v, float mag) =>
        v.Normalized() * mag;

    /// <inheritdoc cref="WithMagnitude(Vector2,float)"/>
    [MethodImpl(Inline)] public static Vector3 WithMagnitude(this Vector3 v, float mag) =>
        v.Normalized() * mag;

    /// <summary>Returns a vector with the same direction, but extending the magnitude by the given amount</summary>
    [MethodImpl(Inline)]
    public static Vector2 AddMagnitude(this Vector2 v, float extraMagnitude) =>
        v * (1 + extraMagnitude / v.Length());

    /// <inheritdoc cref="AddMagnitude(Vector2,float)"/>
    [MethodImpl(Inline)]
    public static Vector3 AddMagnitude(this Vector3 v, float extraMagnitude) =>
        v * (1 + extraMagnitude / v.Length());

    /// <summary>Returns the vector going from one position to another, also known as the displacement.
    /// Equivalent to <c>target-v</c></summary>
    [MethodImpl(Inline)] public static Vector2 To(this Vector2 v, Vector2 target) => target - v;

    /// <summary>Returns the vector going from one position to another, also known as the displacement.
    /// Equivalent to <c>target-v</c></summary>
    [MethodImpl(Inline)] public static Vector3 To(this Vector3 v, Vector3 target) => target - v;

    /// <summary>Returns the normalized direction from this vector to the target.
    /// Equivalent to <c>(target-v).Normalized()</c> or <c>v.To(target).Normalized()</c></summary>
    [MethodImpl(Inline)] public static Vector2 DirTo(this Vector2 v, Vector2 target) =>
        (target - v).Normalized();

    /// <summary>Returns the normalized direction from this vector to the target.
    /// Equivalent to <c>(target-v).Normalized()</c> or <c>v.To(target).Normalized()</c></summary>
    [MethodImpl(Inline)] public static Vector3 DirTo(this Vector3 v, Vector3 target) =>
        (target - v).Normalized();

    /// <summary>Mirrors this vector around another point. Equivalent to rotating the vector 180° around the point</summary>
    /// <param name="p">The point to mirror</param>
    /// <param name="pivot">The point to mirror around</param>
    [MethodImpl(Inline)] public static Vector2 MirrorAround(this Vector2 p, Vector2 pivot) =>
        new(2 * pivot.X - p.X, 2 * pivot.Y - p.Y);

    /// <summary>Mirrors this vector around an x coordinate</summary>
    /// <param name="p">The point to mirror</param>
    /// <param name="xPivot">The x coordinate to mirror around</param>
    [MethodImpl(Inline)] public static Vector2 MirrorAroundX(this Vector2 p, float xPivot) =>
        new(2 * xPivot - p.X, p.Y);

    /// <summary>Mirrors this vector around a y coordinate</summary>
    /// <param name="p">The point to mirror</param>
    /// <param name="yPivot">The y coordinate to mirror around</param>
    [MethodImpl(Inline)] public static Vector2 MirrorAroundY(this Vector2 p, float yPivot) =>
        new(p.X, 2 * yPivot - p.Y);

    /// <inheritdoc cref="MirrorAroundX(Vector2,float)"/>
    [MethodImpl(Inline)] public static Vector3 MirrorAroundX(this Vector3 p, float xPivot) =>
        p with
        {
            X = 2 * xPivot - p.X,
        };

    /// <inheritdoc cref="MirrorAroundY(Vector2,float)"/>
    [MethodImpl(Inline)] public static Vector3 MirrorAroundY(this Vector3 p, float yPivot) =>
        p with
        {
            Y = 2 * yPivot - p.Y,
        };

    /// <summary>Mirrors this vector around a y coordinate</summary>
    /// <param name="p">The point to mirror</param>
    /// <param name="zPivot">The z coordinate to mirror around</param>
    [MethodImpl(Inline)] public static Vector3 MirrorAroundZ(this Vector3 p, float zPivot) =>
        p with
        {
            Z = 2 * zPivot - p.Z,
        };

    /// <inheritdoc cref="MirrorAround(Vector2,Vector2)"/>
    [MethodImpl(Inline)] public static Vector3 MirrorAround(this Vector3 p, Vector3 pivot) =>
        new(2 * pivot.X - p.X, 2 * pivot.Y - p.Y, 2 * pivot.Z - p.Z);

    /// <summary>Scale the point <c>p</c> around <c>pivot</c> by <c>scale</c></summary>
    /// <param name="p">The point to scale</param>
    /// <param name="pivot">The pivot to scale around</param>
    /// <param name="scale">The scale to scale by</param>
    [MethodImpl(Inline)]
    public static Vector2 ScaleAround(this Vector2 p, Vector2 pivot, Vector2 scale) =>
        new(pivot.X + (p.X - pivot.X) * scale.X, pivot.Y + (p.Y - pivot.Y) * scale.Y);

    /// <inheritdoc cref="ScaleAround(Vector2,Vector2,Vector2)"/>
    [MethodImpl(Inline)]
    public static Vector3 ScaleAround(this Vector3 p, Vector3 pivot, Vector3 scale) => new(
        pivot.X + (p.X - pivot.X) * scale.X, pivot.Y + (p.Y - pivot.Y) * scale.Y,
        pivot.Z + (p.Z - pivot.Z) * scale.Z);

    #endregion

    #region Quaternions

    /// <summary>Rotates 180° around the extrinsic pre-rotation X axis, sometimes this is interpreted as a world space rotation, as opposed to rotating around its own axes</summary>
    [MethodImpl(Inline)] public static Quaternion Rotate180AroundExtrX(this Quaternion q) =>
        new(q.W, -q.Z, q.Y, -q.X);

    /// <summary>Rotates 180° around the extrinsic pre-rotation Y axis, sometimes this is interpreted as a world space rotation, as opposed to rotating around its own axes</summary>
    [MethodImpl(Inline)] public static Quaternion Rotate180AroundExtrY(this Quaternion q) =>
        new(q.Z, q.W, -q.X, -q.Y);

    /// <summary>Rotates 180° around the extrinsic pre-rotation Z axis, sometimes this is interpreted as a world space rotation, as opposed to rotating around its own axes</summary>
    [MethodImpl(Inline)] public static Quaternion Rotate180AroundExtrZ(this Quaternion q) =>
        new(-q.Y, q.X, q.W, -q.Z);

    /// <summary>Rotates 180° around its local X axis</summary>
    [MethodImpl(Inline)] public static Quaternion Rotate180AroundSelfX(this Quaternion q) =>
        new(q.W, q.Z, -q.Y, -q.X);

    /// <summary>Rotates 180° around its local Y axis</summary>
    [MethodImpl(Inline)] public static Quaternion Rotate180AroundSelfY(this Quaternion q) =>
        new(-q.Z, q.W, q.X, -q.Y);

    /// <summary>Rotates 180° around its local Z axis</summary>
    [MethodImpl(Inline)] public static Quaternion Rotate180AroundSelfZ(this Quaternion q) =>
        new(q.Y, -q.X, q.W, -q.Z);

    /// <summary>Returns an 180° rotated version of this quaternion around the given axis</summary>
    /// <param name="q">The quaternion to rotate</param>
    /// <param name="axis">The axis to rotate around</param>
    /// <param name="space">The rotation space of the axis, if it should be intrinsic/self/local or extrinsic/"world"</param>
    public static Quaternion Rotate180Around(this Quaternion q, Axis axis,
        RotationSpace space = RotationSpace.Self) =>
        axis switch
        {
            Axis.X => space == RotationSpace.Self
                ? Rotate180AroundSelfX(q)
                : Rotate180AroundExtrX(q),
            Axis.Y => space == RotationSpace.Self
                ? Rotate180AroundSelfY(q)
                : Rotate180AroundExtrY(q),
            Axis.Z => space == RotationSpace.Self
                ? Rotate180AroundSelfZ(q)
                : Rotate180AroundExtrZ(q),
            _ => throw new ArgumentOutOfRangeException(nameof(axis),
                $"Invalid axis: {axis}. Expected 0, 1 or 2"),
        };

    /// <summary>Returns the quaternion rotated around the given axis by the given angle in radians</summary>
    /// <param name="q">The quaternion to rotate</param>
    /// <param name="axis">The axis to rotate around</param>
    /// <param name="angRad">The angle to rotate by (in radians)</param>
    /// <param name="space">The rotation space of the axis, if it should be intrinsic/self/local or extrinsic/"world"</param>
    public static Quaternion RotateAround(this Quaternion q, Axis axis, float angRad,
        RotationSpace space = RotationSpace.Self)
    {
        var aHalf = angRad / 2f;
        var c = MathF.Cos(aHalf);
        var s = MathF.Sin(aHalf);
        var xc = q.X * c;
        var yc = q.Y * c;
        var zc = q.Z * c;
        var wc = q.W * c;
        var xs = q.X * s;
        var ys = q.Y * s;
        var zs = q.Z * s;
        var ws = q.W * s;

        return space switch
        {
            RotationSpace.Self => axis switch
            {
                Axis.X => new(xc + ws, yc + zs, zc - ys, wc - xs),
                Axis.Y => new(xc - zs, yc + ws, zc + xs, wc - ys),
                Axis.Z => new(xc + ys, yc - xs, zc + ws, wc - zs),
                _ => throw new ArgumentOutOfRangeException(nameof(axis))
            },
            RotationSpace.Extrinsic => axis switch
            {
                Axis.X => new(xc + ws, yc - zs, zc + ys, wc - xs),
                Axis.Y => new(xc + zs, yc + ws, zc - xs, wc - ys),
                Axis.Z => new(xc - ys, yc + xs, zc + ws, wc - zs),
                _ => throw new ArgumentOutOfRangeException(nameof(axis))
            },
            _ => throw new ArgumentOutOfRangeException(nameof(space))
        };
    }

    /// <summary>Returns the quaternion rotated around the given axis by 90°</summary>
    /// <param name="q">The quaternion to rotate</param>
    /// <param name="axis">The axis to rotate around</param>
    /// <param name="space">The rotation space of the axis, if it should be intrinsic/self/local or extrinsic/"world"</param>
    public static Quaternion Rotate90Around(this Quaternion q,
        Axis axis, RotationSpace space = RotationSpace.Self)
    {
        const float v = Mathf.RSQRT2; // cos(90°/2) = sin(90°/2)
        var (x, y, z, w) = q;

        return space switch
        {
            RotationSpace.Self => axis switch
            {
                Axis.X => new Quaternion(v * (x + w), v * (y + z), v * (z - y), v * (w - x)),
                Axis.Y => new Quaternion(v * (x - z), v * (y + w), v * (z + x), v * (w - y)),
                Axis.Z => new Quaternion(v * (x + y), v * (y - x), v * (z + w), v * (w - z)),
                _ => throw new ArgumentOutOfRangeException(nameof(axis))
            },
            RotationSpace.Extrinsic => axis switch
            {
                Axis.X => new Quaternion(v * (x + w), v * (y - z), v * (z + y), v * (w - x)),
                Axis.Y => new Quaternion(v * (x + z), v * (y + w), v * (z - x), v * (w - y)),
                Axis.Z => new Quaternion(v * (x - y), v * (y + x), v * (z + w), v * (w - z)),
                _ => throw new ArgumentOutOfRangeException(nameof(axis))
            },
            _ => throw new ArgumentOutOfRangeException(nameof(space))
        };
    }

    /// <summary>Returns the quaternion rotated around the given axis by -90°</summary>
    /// <param name="q">The quaternion to rotate</param>
    /// <param name="axis">The axis to rotate around</param>
    /// <param name="space">The rotation space of the axis, if it should be intrinsic/self/local or extrinsic/"world"</param>
    public static Quaternion RotateNeg90Around(this Quaternion q, Axis axis,
        RotationSpace space = RotationSpace.Self)
    {
        const float v = Mathf.RSQRT2; // cos(90°/2) = sin(90°/2)
        var (x, y, z, w) = q;

        return space switch
        {
            RotationSpace.Self => axis switch
            {
                Axis.X => new Quaternion(v * (x - w), v * (y - z), v * (z + y), v * (w + x)),
                Axis.Y => new Quaternion(v * (x + z), v * (y - w), v * (z - x), v * (w + y)),
                Axis.Z => new Quaternion(v * (x - y), v * (y + x), v * (z - w), v * (w + z)),
                _ => throw new ArgumentOutOfRangeException(nameof(axis))
            },
            RotationSpace.Extrinsic => axis switch
            {
                Axis.X => new Quaternion(v * (x - w), v * (y + z), v * (z - y), v * (w + x)),
                Axis.Y => new Quaternion(v * (x - z), v * (y - w), v * (z + x), v * (w + y)),
                Axis.Z => new Quaternion(v * (x + y), v * (y - x), v * (z - w), v * (w + z)),
                _ => throw new ArgumentOutOfRangeException(nameof(axis))
            },
            _ => throw new ArgumentOutOfRangeException(nameof(space))
        };
    }

    /// <summary>Returns the given axis of this rotation (assumes this quaternion is normalized)</summary>
    public static Vector3 GetAxis(this Quaternion q, Axis axis) =>
        axis switch
        {
            Axis.X => q.Right(),
            Axis.Y => q.Up(),
            Axis.Z => q.Forward(),
            _ => throw new ArgumentOutOfRangeException(nameof(axis))
        };

    /// <summary>Returns the X axis of this rotation (assumes this quaternion is normalized)</summary>
    [MethodImpl(Inline)] public static Vector3 Right(this Quaternion q) =>
        new(
            q.X * q.X - q.Y * q.Y - q.Z * q.Z + q.W * q.W,
            2 * (q.X * q.Y + q.Z * q.W),
            2 * (q.X * q.Z - q.Y * q.W)
        );

    /// <summary>Returns the Y axis of this rotation (assumes this quaternion is normalized)</summary>
    [MethodImpl(Inline)] public static Vector3 Up(this Quaternion q) =>
        new(
            2 * (q.X * q.Y - q.Z * q.W),
            -q.X * q.X + q.Y * q.Y - q.Z * q.Z + q.W * q.W,
            2 * (q.X * q.W + q.Y * q.Z)
        );

    /// <summary>Returns the Z axis of this rotation (assumes this quaternion is normalized)</summary>
    [MethodImpl(Inline)] public static Vector3 Forward(this Quaternion q) =>
        new(
            2 * (q.X * q.Z + q.Y * q.W),
            2 * (q.Y * q.Z - q.X * q.W),
            -q.X * q.X - q.Y * q.Y + q.Z * q.Z + q.W * q.W
        );

    /// <summary>Converts this quaternion to a rotation matrix</summary>
    public static Matrix4x4 ToMatrix(this Quaternion q)
    {
        // you could just use Matrix4x4.Rotate( q ) but that's not as fun as doing this math myself
        var xx = q.X * q.X;
        var yy = q.Y * q.Y;
        var zz = q.Z * q.Z;
        var ww = q.W * q.W;
        var xy = q.X * q.Y;
        var yz = q.Y * q.Z;
        var zw = q.Z * q.W;
        var wx = q.W * q.X;
        var xz = q.X * q.Z;
        var yw = q.Y * q.W;

        return new()
        {
            M11 = xx - yy - zz + ww, // X
            M21 = 2 * (xy + zw),
            M31 = 2 * (xz - yw),
            M12 = 2 * (xy - zw), // Y
            M22 = -xx + yy - zz + ww,
            M32 = 2 * (wx + yz),
            M13 = 2 * (xz + yw), // Z
            M23 = 2 * (yz - wx),
            M33 = -xx - yy + zz + ww,
            M44 = 1,
        };
    }

    /// <summary>Returns the natural logarithm of a quaternion</summary>
    public static Quaternion Ln(this Quaternion q)
    {
        var vMagSq = (double) q.X * q.X + (double) q.Y * q.Y + (double) q.Z * q.Z;
        var vMag = Math.Sqrt(vMagSq);
        var qMag = Math.Sqrt(vMagSq + (double) q.W * q.W);
        var theta = Math.Atan2(vMag, q.W);
        var scV = vMag < 0.001 ? Mathf.SincRcp(theta) / qMag : theta / vMag;
        return new(
            (float) (scV * q.X),
            (float) (scV * q.Y),
            (float) (scV * q.Z),
            (float) Math.Log(qMag)
        );
    }

    /// <summary>Returns the natural logarithm of a unit quaternion</summary>
    public static Quaternion LnUnit(this Quaternion q)
    {
        var vMagSq = (double) q.X * q.X + (double) q.Y * q.Y + (double) q.Z * q.Z;
        var vMag = Math.Sqrt(vMagSq);
        var theta = Math.Atan2(vMag, q.W);
        var scV = vMag < 0.001 ? Mathf.SincRcp(theta) : theta / vMag;
        return new Quaternion(
            (float) (scV * q.X),
            (float) (scV * q.Y),
            (float) (scV * q.Z),
            0f
        );
    }

    /// <summary>Returns the natural exponent of a quaternion</summary>
    public static Quaternion Exp(this Quaternion q)
    {
        var vMag = Math.Sqrt((double) q.X * q.X + (double) q.Y * q.Y + (double) q.Z * q.Z);
        var sc = Math.Exp(q.W);
        var scV = sc * Mathf.Sinc(vMag);
        return new Quaternion((float) (scV * q.X), (float) (scV * q.Y), (float) (scV * q.Z),
            (float) (sc * Math.Cos(vMag)));
    }

    /// <summary>Returns the natural exponent of a pure imaginary quaternion</summary>
    public static Quaternion ExpPureIm(this Quaternion q)
    {
        var vMag = Math.Sqrt((double) q.X * q.X + (double) q.Y * q.Y + (double) q.Z * q.Z);
        var scV = Mathf.Sinc(vMag);
        return new Quaternion((float) (scV * q.X), (float) (scV * q.Y), (float) (scV * q.Z),
            (float) Math.Cos(vMag));
    }

    /// <summary>Returns the quaternion raised to a real power</summary>
    public static Quaternion Pow(this Quaternion q, float x)
    {
        return x switch
        {
            0f => new Quaternion(0, 0, 0, 1),
            1f => q,
            _ => q.Ln().Mul(x).Exp()
        };
    }

    /// <summary>Returns the unit quaternion raised to a real power</summary>
    public static Quaternion PowUnit(this Quaternion q, float x)
    {
        return x switch
        {
            0f => new Quaternion(0, 0, 0, 1),
            1f => q,
            _ => q.LnUnit().Mul(x).ExpPureIm()
        };
    }

    /// <summary>Returns the imaginary part of a quaternion as a vector</summary>
    public static Vector3 Imag(this Quaternion q) => new(q.X, q.Y, q.Z);

    /// <summary>Returns the squared magnitude of this quaternion</summary>
    public static float SqrMagnitude(this Quaternion q) => (float) ((double) q.W * q.W +
        (double) q.X * q.X + (double) q.Y * q.Y + (double) q.Z * q.Z);

    /// <summary>Returns the magnitude of this quaternion</summary>
    public static float Magnitude(this Quaternion q) => MathF.Sqrt(q.SqrMagnitude());

    /// <summary>Multiplies a quaternion by a scalar</summary>
    /// <param name="q">The quaternion to multiply</param>
    /// <param name="c">The scalar value to multiply with</param>
    public static Quaternion Mul(this Quaternion q, float c) =>
        Quaternion.Multiply(q, c);

    /// <summary>Multiplies a quaternion by a vector3</summary>
    /// <param name="q">The quaternion to multiply</param>
    /// <param name="v">The vector value to multiply with</param>
    public static Vector3 Mul(this Quaternion q, Vector3 v) =>
        Vector3.Transform(v, q);

    /// <summary>Adds a quaternion to an existing quaternion</summary>
    public static Quaternion Add(this Quaternion a, Quaternion b) =>
        new(a.X + b.X, a.Y + b.Y, a.Z + b.Z, a.W + b.W);

    /// <summary>Subtracts a quaternion from an existing quaternion</summary>
    public static Quaternion Sub(this Quaternion a, Quaternion b) =>
        new(a.X - b.X, a.Y - b.Y, a.Z - b.Z, a.W - b.W);

    /// <summary>The conjugate of a quaternion</summary>
    /// <param name="q">The quaternion to conjugate</param>
    public static Quaternion Conjugate(this Quaternion q) =>
        new(-q.X, -q.Y, -q.Z, q.W);

    /// <inheritdoc cref="Quaternion.Inverse(Quaternion)"/>
    public static Quaternion Inverse(this Quaternion q) => Quaternion.Inverse(q);

    /// <summary>The inverse of a unit quaternion, equivalent to the quaternion conjugate</summary>
    public static Quaternion InverseUnit(this Quaternion q)
    {
        return new Quaternion(-q.X, -q.Y, -q.Z, q.W);
    }

    /// <summary>The inverse of a pure imaginary unit quaternion, where w is assumed to be 0</summary>
    public static Quaternion InversePureIm(this Quaternion q)
    {
        var sqMag = q.X * q.X + q.Y * q.Y + q.Z * q.Z;
        return new Quaternion(-q.X / sqMag, -q.Y / sqMag, -q.Z / sqMag, 0);
    }

    /// <summary>Add to the magnitude of this quaternion</summary>
    public static Quaternion AddMagnitude(this Quaternion q, float amount) =>
        amount == 0f ? q : q.Mul(1 + amount / q.Magnitude());

    #endregion

    #region Simple float and int operations

    /// <summary>Returns whether or not two values are approximately equal.
    /// They are considered equal if they are within a <c>Epsilon*8</c> or <c>max(a,b)*0.000001f</c> range of each other</summary>
    [MethodImpl(Inline)] public static bool Approximately(this float a, float b) =>
        Mathf.Approximately(a, b);

    /// <summary>Returns true if v is between or equal to <c>min</c> &amp; <c>max</c></summary>
    /// <seealso cref="Between(float,float,float)"/>
    [MethodImpl(Inline)] public static bool Within(this float v, float min, float max) =>
        v >= min && v <= max;

    /// <summary>Returns true if v is between or equal to <c>min</c> &amp; <c>max</c></summary>
    /// <seealso cref="Between(int,int,int)"/>
    [MethodImpl(Inline)] public static bool Within(this int v, int min, int max) =>
        v >= min && v <= max;

    /// <summary>Returns true if v is between, but not equal to, <c>min</c> &amp; <c>max</c></summary>
    /// <seealso cref="Within(float,float,float)"/>
    [MethodImpl(Inline)] public static bool Between(this float v, float min, float max) =>
        v > min && v < max;

    /// <summary>Returns true if v is between, but not equal to, <c>min</c> &amp; <c>max</c></summary>
    /// <seealso cref="Within(int,int,int)"/>
    [MethodImpl(Inline)] public static bool Between(this int v, int min, int max) =>
        v > min && v < max;

    /// <summary>Clamps the value to be at least <c>min</c></summary>
    [MethodImpl(Inline)] public static float AtLeast(this float v, float min) =>
        v < min ? min : v;

    /// <summary>Clamps the value to be at least <c>min</c></summary>
    [MethodImpl(Inline)] public static int AtLeast(this int v, int min) => v < min ? min : v;

    /// <summary>Clamps the value to be at most <c>max</c></summary>
    [MethodImpl(Inline)] public static float AtMost(this float v, float max) =>
        v > max ? max : v;

    /// <summary>Clamps the value to be at most <c>max</c></summary>
    [MethodImpl(Inline)] public static int AtMost(this int v, int max) => v > max ? max : v;

    /// <summary>Squares the value. Equivalent to <c>v*v</c></summary>
    [MethodImpl(Inline)] public static float Square(this float v) => v * v;

    /// <summary>Cubes the value. Equivalent to <c>v*v*v</c></summary>
    [MethodImpl(Inline)] public static float Cube(this float v) => v * v * v;

    /// <summary>Squares the value. Equivalent to <c>v*v</c></summary>
    [MethodImpl(Inline)] public static int Square(this int v) => v * v;

    /// <summary>The next integer, modulo <c>length</c>. Behaves the way you want with negative values for stuff like array index access etc</summary>
    [MethodImpl(Inline)] public static int NextMod(this int value, int length) =>
        (value + 1).Mod(length);

    /// <summary>The previous integer, modulo <c>length</c>. Behaves the way you want with negative values for stuff like array index access etc</summary>
    [MethodImpl(Inline)] public static int PrevMod(this int value, int length) =>
        (value - 1).Mod(length);

    #endregion

    #region String extensions

    public static string ToValueTableString(this string[,] m)
    {
        var rowCount = m.GetLength(0);
        var colCount = m.GetLength(1);
        string[] r = new string[rowCount];
        for (var i = 0; i < rowCount; i++)
            r[i] = "";

        for (var c = 0; c < colCount; c++)
        {
            var endBit = c == colCount - 1 ? "" : ", ";

            var colWidth = 4; // min width
            string[] columnEntries = new string[rowCount];
            for (var row = 0; row < rowCount; row++)
            {
                var s = m[row, c].StartsWith('-') ? "" : " ";
                columnEntries[row] = $"{s}{m[row, c]}{endBit}";
                colWidth = Mathf.Max(colWidth, columnEntries[row].Length);
            }

            for (var row = 0; row < rowCount; row++)
            {
                r[row] += columnEntries[row].PadRight(colWidth, ' ');
            }
        }

        return string.Join('\n', r);
    }

    #endregion

    #region Matrix extensions

    public static Vector4 GetColumn(this Matrix4x4 v, int index) =>
        index switch
        {
            0 => new(v.M11, v.M21, v.M31, v.M41),
            1 => new(v.M12, v.M22, v.M32, v.M42),
            2 => new(v.M13, v.M23, v.M33, v.M43),
            3 => new(v.M14, v.M24, v.M34, v.M44),
            _ => throw new IndexOutOfRangeException("Invalid column index!"),
        };


    public static Vector4 GetRow(this Matrix4x4 v, int index) =>
        index switch
        {
            0 => new(v.M11, v.M12, v.M13, v.M14),
            1 => new(v.M21, v.M22, v.M23, v.M24),
            2 => new(v.M31, v.M32, v.M33, v.M34),
            3 => new(v.M41, v.M42, v.M43, v.M44),
            _ => throw new IndexOutOfRangeException("Invalid row index!"),
        };

    public static float AverageScale(this Matrix4x4 m) =>
    (
        m.GetColumn(0).ToVector3().Length() +
        m.GetColumn(1).ToVector3().Length() +
        m.GetColumn(2).ToVector3().Length()
    ) / 3;

    public static Matrix4x1 MultiplyColumnVector(this Matrix4x4 m, Matrix4x1 v) =>
        new(
            m.M11 * v.m0 + m.M12 * v.m1 + m.M13 * v.m2 + m.M14 * v.m3,
            m.M21 * v.m0 + m.M22 * v.m1 + m.M23 * v.m2 + m.M24 * v.m3,
            m.M31 * v.m0 + m.M32 * v.m1 + m.M33 * v.m2 + m.M34 * v.m3,
            m.M41 * v.m0 + m.M42 * v.m1 + m.M43 * v.m2 + m.M44 * v.m3
        );

    #endregion

    #region Extension method counterparts of the static Mathf functions - lots of boilerplate in here

    #region Math operations

    /// <inheritdoc cref="Mathf.Sqrt(float)"/>
    [MethodImpl(Inline)] public static float Sqrt(this float value) => MathF.Sqrt(value);

    /// <inheritdoc cref="Mathf.Sqrt(Vector2)"/>
    [MethodImpl(Inline)] public static Vector2 Sqrt(this Vector2 value) => Mathf.Sqrt(value);

    /// <inheritdoc cref="Mathf.Sqrt(Vector3)"/>
    [MethodImpl(Inline)] public static Vector3 Sqrt(this Vector3 value) => Mathf.Sqrt(value);

    /// <inheritdoc cref="Mathf.Sqrt(Vector4)"/>
    [MethodImpl(Inline)] public static Vector4 Sqrt(this Vector4 value) => Mathf.Sqrt(value);

    /// <inheritdoc cref="Mathf.Cbrt(float)"/>
    [MethodImpl(Inline)] public static float Cbrt(this float value) => MathF.Cbrt(value);

    /// <inheritdoc cref="Mathf.Pow(float, float)"/>
    [MethodImpl(Inline)] public static float Pow(this float value, float exponent) =>
        MathF.Pow(value, exponent);

    /// <summary>Calculates exact positive integer powers</summary>
    /// <param name="value"></param>
    /// <param name="pow">A positive integer power</param>
    [MethodImpl(Inline)] public static int Pow(this int value, int pow)
    {
        if (pow < 0)
            throw new ArithmeticException("int.Pow(int) doesn't support negative powers");
        checked
        {
            switch (pow)
            {
                case 0: return 1;
                case 1: return value;
                case 2: return value * value;
                case 3: return value * value * value;
                default:
                    if (value == 2)
                        return 1 << pow;
                    // from: https://stackoverflow.com/questions/383587/how-do-you-do-integer-exponentiation-in-c
                    var ret = 1;
                    while (pow != 0)
                    {
                        if ((pow & 1) == 1)
                            ret *= value;
                        value *= value;
                        pow >>= 1;
                    }

                    return ret;
            }
        }
    }

    #endregion

    #region Absolute Values

    /// <inheritdoc cref="Mathf.Abs(float)"/>
    [MethodImpl(Inline)] public static float Abs(this float value) => MathF.Abs(value);

    /// <inheritdoc cref="Mathf.Abs(int)"/>
    [MethodImpl(Inline)] public static int Abs(this int value) => Mathf.Abs(value);

    /// <inheritdoc cref="Mathf.Abs(Vector2)"/>
    [MethodImpl(Inline)] public static Vector2 Abs(this Vector2 v) => Mathf.Abs(v);

    /// <inheritdoc cref="Mathf.Abs(Vector3)"/>
    [MethodImpl(Inline)] public static Vector3 Abs(this Vector3 v) => Mathf.Abs(v);

    /// <inheritdoc cref="Mathf.Abs(Vector4)"/>
    [MethodImpl(Inline)] public static Vector4 Abs(this Vector4 v) => Mathf.Abs(v);

    #endregion

    #region Clamping

    /// <inheritdoc cref="Mathf.Clamp(float,float,float)"/>
    [MethodImpl(Inline)] public static float Clamp(this float value, float min, float max) =>
        Mathf.Clamp(value, min, max);

    /// <inheritdoc cref="Mathf.Clamp(Vector2,Vector2,Vector2)"/>
    [MethodImpl(Inline)]
    public static Vector2 Clamp(this Vector2 v, Vector2 min, Vector2 max) =>
        Mathf.Clamp(v, min, max);

    /// <inheritdoc cref="Mathf.Clamp(Vector3,Vector3,Vector3)"/>
    [MethodImpl(Inline)]
    public static Vector3 Clamp(this Vector3 v, Vector3 min, Vector3 max) =>
        Mathf.Clamp(v, min, max);

    /// <inheritdoc cref="Mathf.Clamp(Vector4,Vector4,Vector4)"/>
    [MethodImpl(Inline)]
    public static Vector4 Clamp(this Vector4 v, Vector4 min, Vector4 max) =>
        Mathf.Clamp(v, min, max);

    /// <inheritdoc cref="Mathf.Clamp(int,int,int)"/>
    [MethodImpl(Inline)] public static int Clamp(this int value, int min, int max) =>
        Mathf.Clamp(value, min, max);

    /// <inheritdoc cref="Mathf.Clamp01(float)"/>
    [MethodImpl(Inline)] public static float Clamp01(this float value) =>
        Mathf.Clamp01(value);

    /// <inheritdoc cref="Mathf.Clamp01(Vector2)"/>
    [MethodImpl(Inline)] public static Vector2 Clamp01(this Vector2 v) => Mathf.Clamp01(v);

    /// <inheritdoc cref="Mathf.Clamp01(Vector3)"/>
    [MethodImpl(Inline)] public static Vector3 Clamp01(this Vector3 v) => Mathf.Clamp01(v);

    /// <inheritdoc cref="Mathf.Clamp01(Vector4)"/>
    [MethodImpl(Inline)] public static Vector4 Clamp01(this Vector4 v) => Mathf.Clamp01(v);

    /// <inheritdoc cref="Mathf.ClampNeg1to1(float)"/>
    [MethodImpl(Inline)] public static float ClampNeg1To1(this float value) =>
        Mathf.ClampNeg1to1(value);

    /// <inheritdoc cref="Mathf.ClampNeg1to1(float)"/>
    [MethodImpl(Inline)] public static double ClampNeg1To1(this double value) =>
        Mathf.ClampNeg1to1(value);

    /// <inheritdoc cref="Mathf.ClampNeg1to1(Vector2)"/>
    [MethodImpl(Inline)] public static Vector2 ClampNeg1To1(this Vector2 v) =>
        Mathf.ClampNeg1to1(v);

    /// <inheritdoc cref="Mathf.ClampNeg1to1(Vector3)"/>
    [MethodImpl(Inline)] public static Vector3 ClampNeg1To1(this Vector3 v) =>
        Mathf.ClampNeg1to1(v);

    /// <inheritdoc cref="Mathf.ClampNeg1to1(Vector4)"/>
    [MethodImpl(Inline)] public static Vector4 ClampNeg1To1(this Vector4 v) =>
        Mathf.ClampNeg1to1(v);

    #endregion

    #region Min & Max

    /// <inheritdoc cref="Mathf.Min(Vector2)"/>
    [MethodImpl(Inline)] public static float Min(this Vector2 v) => Mathf.Min(v);

    /// <inheritdoc cref="Mathf.Min(Vector3)"/>
    [MethodImpl(Inline)] public static float Min(this Vector3 v) => Mathf.Min(v);

    /// <inheritdoc cref="Mathf.Min(Vector4)"/>
    [MethodImpl(Inline)] public static float Min(this Vector4 v) => Mathf.Min(v);

    /// <inheritdoc cref="Mathf.Max(Vector2)"/>
    [MethodImpl(Inline)] public static float Max(this Vector2 v) => Mathf.Max(v);

    /// <inheritdoc cref="Mathf.Max(Vector3)"/>
    [MethodImpl(Inline)] public static float Max(this Vector3 v) => Mathf.Max(v);

    /// <inheritdoc cref="Mathf.Max(Vector4)"/>
    [MethodImpl(Inline)] public static float Max(this Vector4 v) => Mathf.Max(v);

    #endregion

    #region Signs & Rounding

    /// <inheritdoc cref="Mathf.Sign(float)"/>
    [MethodImpl(Inline)] public static float Sign(this float value) => Mathf.Sign(value);

    /// <inheritdoc cref="Mathf.Sign(Vector2)"/>
    [MethodImpl(Inline)] public static Vector2 Sign(this Vector2 value) => Mathf.Sign(value);

    /// <inheritdoc cref="Mathf.Sign(Vector3)"/>
    [MethodImpl(Inline)] public static Vector3 Sign(this Vector3 value) => Mathf.Sign(value);

    /// <inheritdoc cref="Mathf.Sign(Vector4)"/>
    [MethodImpl(Inline)] public static Vector4 Sign(this Vector4 value) => Mathf.Sign(value);

    /// <inheritdoc cref="Mathf.Sign(int)"/>
    [MethodImpl(Inline)] public static int Sign(this int value) => Mathf.Sign(value);

    /// <inheritdoc cref="Mathf.SignAsInt(float)"/>
    [MethodImpl(Inline)] public static int SignAsInt(this float value) =>
        Mathf.SignAsInt(value);

    /// <inheritdoc cref="Mathf.SignWithZero(float,float)"/>
    [MethodImpl(Inline)]
    public static float SignWithZero(this float value, float zeroThreshold = 0.000001f) =>
        Mathf.SignWithZero(value, zeroThreshold);

    /// <inheritdoc cref="Mathf.SignWithZero(Vector2,float)"/>
    [MethodImpl(Inline)]
    public static Vector2 SignWithZero(this Vector2 value, float zeroThreshold = 0.000001f) =>
        Mathf.SignWithZero(value, zeroThreshold);

    /// <inheritdoc cref="Mathf.SignWithZero(Vector3,float)"/>
    [MethodImpl(Inline)]
    public static Vector3 SignWithZero(this Vector3 value, float zeroThreshold = 0.000001f) =>
        Mathf.SignWithZero(value, zeroThreshold);

    /// <inheritdoc cref="Mathf.SignWithZero(Vector4,float)"/>
    [MethodImpl(Inline)]
    public static Vector4 SignWithZero(this Vector4 value, float zeroThreshold = 0.000001f) =>
        Mathf.SignWithZero(value, zeroThreshold);

    /// <inheritdoc cref="Mathf.SignWithZero(int)"/>
    [MethodImpl(Inline)] public static int SignWithZero(this int value) =>
        Mathf.SignWithZero(value);

    /// <inheritdoc cref="Mathf.SignWithZeroAsInt(float,float)"/>
    [MethodImpl(Inline)]
    public static int SignWithZeroAsInt(this float value, float zeroThreshold = 0.000001f) =>
        Mathf.SignWithZeroAsInt(value, zeroThreshold);

    /// <inheritdoc cref="Mathf.Floor(float)"/>
    [MethodImpl(Inline)] public static float Floor(this float value) => Mathf.Floor(value);

    /// <inheritdoc cref="Mathf.Floor(Vector2)"/>
    [MethodImpl(Inline)] public static Vector2 Floor(this Vector2 value) =>
        Mathf.Floor(value);

    /// <inheritdoc cref="Mathf.Floor(Vector2)"/>
    [MethodImpl(Inline)] public static Vector3 Floor(this Vector3 value) =>
        Mathf.Floor(value);

    /// <inheritdoc cref="Mathf.Floor(Vector2)"/>
    [MethodImpl(Inline)] public static Vector4 Floor(this Vector4 value) =>
        Mathf.Floor(value);

    /// <inheritdoc cref="Mathf.FloorToInt(float)"/>
    [MethodImpl(Inline)] public static int FloorToInt(this float value) =>
        Mathf.FloorToInt(value);

    /// <inheritdoc cref="Mathf.FloorToInt(Vector2)"/>
    [MethodImpl(Inline)] public static Vector2Int FloorToInt(this Vector2 value) =>
        Mathf.FloorToInt(value);

    /// <inheritdoc cref="Mathf.FloorToInt(Vector2)"/>
    [MethodImpl(Inline)] public static Vector3Int FloorToInt(this Vector3 value) =>
        Mathf.FloorToInt(value);

    /// <inheritdoc cref="Mathf.Ceil(float)"/>
    [MethodImpl(Inline)] public static float Ceil(this float value) => Mathf.Ceil(value);

    /// <inheritdoc cref="Mathf.Ceil(Vector2)"/>
    [MethodImpl(Inline)] public static Vector2 Ceil(this Vector2 value) => Mathf.Ceil(value);

    /// <inheritdoc cref="Mathf.Ceil(Vector2)"/>
    [MethodImpl(Inline)] public static Vector3 Ceil(this Vector3 value) => Mathf.Ceil(value);

    /// <inheritdoc cref="Mathf.Ceil(Vector2)"/>
    [MethodImpl(Inline)] public static Vector4 Ceil(this Vector4 value) => Mathf.Ceil(value);

    /// <inheritdoc cref="Mathf.CeilToInt(float)"/>
    [MethodImpl(Inline)] public static int CeilToInt(this float value) =>
        Mathf.CeilToInt(value);

    /// <inheritdoc cref="Mathf.CeilToInt(Vector2)"/>
    [MethodImpl(Inline)] public static Vector2Int CeilToInt(this Vector2 value) =>
        Mathf.CeilToInt(value);

    /// <inheritdoc cref="Mathf.CeilToInt(Vector2)"/>
    [MethodImpl(Inline)] public static Vector3Int CeilToInt(this Vector3 value) =>
        Mathf.CeilToInt(value);

    /// <inheritdoc cref="Mathf.Round(float,MidpointRounding)"/>
    [MethodImpl(Inline)]
    public static float Round(this float value,
        MidpointRounding midpointRounding = MidpointRounding.ToEven) =>
        Mathf.Round(value, midpointRounding);

    /// <inheritdoc cref="Mathf.Round(Vector2,MidpointRounding)"/>
    [MethodImpl(Inline)]
    public static Vector2 Round(this Vector2 value,
        MidpointRounding midpointRounding = MidpointRounding.ToEven) =>
        Mathf.Round(value, midpointRounding);

    /// <inheritdoc cref="Mathf.Round(Vector2,MidpointRounding)"/>
    [MethodImpl(Inline)]
    public static Vector3 Round(this Vector3 value,
        MidpointRounding midpointRounding = MidpointRounding.ToEven) =>
        Mathf.Round(value, midpointRounding);

    /// <inheritdoc cref="Mathf.Round(Vector2,MidpointRounding)"/>
    [MethodImpl(Inline)]
    public static Vector4 Round(this Vector4 value,
        MidpointRounding midpointRounding = MidpointRounding.ToEven) =>
        Mathf.Round(value, midpointRounding);

    /// <inheritdoc cref="Mathf.Round(float,MidpointRounding)"/>
    [MethodImpl(Inline)]
    public static float Round(this float value, float snapInterval,
        MidpointRounding midpointRounding = MidpointRounding.ToEven) =>
        Mathf.Round(value, snapInterval, midpointRounding);

    /// <inheritdoc cref="Mathf.Round(Vector2,MidpointRounding)"/>
    [MethodImpl(Inline)]
    public static Vector2 Round(this Vector2 value, float snapInterval,
        MidpointRounding midpointRounding = MidpointRounding.ToEven) =>
        Mathf.Round(value, snapInterval, midpointRounding);

    /// <inheritdoc cref="Mathf.Round(Vector2,float,MidpointRounding)"/>
    [MethodImpl(Inline)]
    public static Vector3 Round(this Vector3 value, float snapInterval,
        MidpointRounding midpointRounding = MidpointRounding.ToEven) =>
        Mathf.Round(value, snapInterval, midpointRounding);

    /// <inheritdoc cref="Mathf.Round(Vector2,float,MidpointRounding)"/>
    [MethodImpl(Inline)]
    public static Vector4 Round(this Vector4 value, float snapInterval,
        MidpointRounding midpointRounding = MidpointRounding.ToEven) =>
        Mathf.Round(value, snapInterval, midpointRounding);

    /// <inheritdoc cref="Mathf.RoundToInt(float,MidpointRounding)"/>
    [MethodImpl(Inline)]
    public static int RoundToInt(this float value,
        MidpointRounding midpointRounding = MidpointRounding.ToEven) =>
        Mathf.RoundToInt(value, midpointRounding);

    /// <inheritdoc cref="Mathf.RoundToInt(Vector2,MidpointRounding)"/>
    [MethodImpl(Inline)]
    public static Vector2Int RoundToInt(this Vector2 value,
        MidpointRounding midpointRounding = MidpointRounding.ToEven) =>
        Mathf.RoundToInt(value, midpointRounding);

    /// <inheritdoc cref="Mathf.RoundToInt(Vector2,MidpointRounding)"/>
    [MethodImpl(Inline)]
    public static Vector3Int RoundToInt(this Vector3 value,
        MidpointRounding midpointRounding = MidpointRounding.ToEven) =>
        Mathf.RoundToInt(value, midpointRounding);

    #endregion

    #region Range Repeating

    /// <inheritdoc cref="Mathf.Frac(float)"/>
    [MethodImpl(Inline)] public static float Frac(this float x) => Mathf.Frac(x);

    /// <inheritdoc cref="Mathf.Frac(Vector2)"/>
    [MethodImpl(Inline)] public static Vector2 Frac(this Vector2 v) => Mathf.Frac(v);

    /// <inheritdoc cref="Mathf.Frac(Vector3)"/>
    [MethodImpl(Inline)] public static Vector3 Frac(this Vector3 v) => Mathf.Frac(v);

    /// <inheritdoc cref="Mathf.Frac(Vector4)"/>
    [MethodImpl(Inline)] public static Vector4 Frac(this Vector4 v) => Mathf.Frac(v);

    /// <inheritdoc cref="Mathf.Repeat(float,float)"/>
    [MethodImpl(Inline)] public static float Repeat(this float value, float length) =>
        Mathf.Repeat(value, length);

    /// <inheritdoc cref="Mathf.Mod(int,int)"/>
    [MethodImpl(Inline)] public static int Mod(this int value, int length) =>
        Mathf.Mod(value, length);

    #endregion

    #region Smoothing & Easing Curves

    /// <inheritdoc cref="Mathf.Smooth01(float)"/>
    [MethodImpl(Inline)] public static float Smooth01(this float x) => Mathf.Smooth01(x);

    /// <inheritdoc cref="Mathf.Smoother01(float)"/>
    [MethodImpl(Inline)] public static float Smoother01(this float x) => Mathf.Smoother01(x);

    /// <inheritdoc cref="Mathf.SmoothCos01(float)"/>
    [MethodImpl(Inline)] public static float SmoothCos01(this float x) =>
        Mathf.SmoothCos01(x);

    #endregion

    #region Value & Vector interpolation

    /// <inheritdoc cref="Mathf.Remap(float,float,float,float,float)"/>
    [MethodImpl(Inline)]
    public static float
        Remap(this float value, float iMin, float iMax, float oMin, float oMax) =>
        Mathf.Remap(iMin, iMax, oMin, oMax, value);

    /// <inheritdoc cref="Mathf.RemapClamped(float,float,float,float,float)"/>
    [MethodImpl(Inline)]
    public static float RemapClamped(this float value, float iMin, float iMax, float oMin,
        float oMax) => Mathf.RemapClamped(iMin, iMax, oMin, oMax, value);

    /// <inheritdoc cref="Mathf.Remap(FloatRange,FloatRange,float)"/>
    [MethodImpl(Inline)]
    public static float Remap(this float value, FloatRange inRange, FloatRange outRange) =>
        Mathf.Remap(inRange.Start, inRange.End, outRange.Start, outRange.End, value);

    /// <inheritdoc cref="Mathf.RemapClamped(FloatRange,FloatRange,float)"/>
    [MethodImpl(Inline)]
    public static float
        RemapClamped(this float value, FloatRange inRange, FloatRange outRange) =>
        Mathf.RemapClamped(inRange.Start, inRange.End, outRange.Start, outRange.End, value);

    /// <inheritdoc cref="Mathf.Remap(float,float,float,float,int)"/>
    [MethodImpl(Inline)]
    public static float Remap(this int value, float iMin, float iMax, float oMin, float oMax) =>
        Mathf.Remap(iMin, iMax, oMin, oMax, value);

    /// <inheritdoc cref="Mathf.RemapClamped(float,float,float,float,float)"/>
    [MethodImpl(Inline)]
    public static float RemapClamped(this int value, float iMin, float iMax, float oMin,
        float oMax) => Mathf.RemapClamped(iMin, iMax, oMin, oMax, value);

    /// <inheritdoc cref="Mathf.Lerp(float,float,float)"/>
    [MethodImpl(Inline)] public static float Lerp(this float t, float a, float b) =>
        Mathf.Lerp(a, b, t);

    /// <inheritdoc cref="Mathf.InverseLerp(float,float,float)"/>
    [MethodImpl(Inline)] public static float InverseLerp(this float value, float a, float b) =>
        Mathf.InverseLerp(a, b, value);

    /// <inheritdoc cref="Mathf.LerpClamped(float,float,float)"/>
    [MethodImpl(Inline)] public static float LerpClamped(this float t, float a, float b) =>
        Mathf.LerpClamped(a, b, t);

    /// <inheritdoc cref="Mathf.InverseLerpClamped(float,float,float)"/>
    [MethodImpl(Inline)]
    public static float InverseLerpClamped(this float value, float a, float b) =>
        Mathf.InverseLerpClamped(a, b, value);

    /// <inheritdoc cref="Mathf.Remap(Vector2,Vector2,Vector2,Vector2,Vector2)"/>
    [MethodImpl(Inline)]
    public static Vector2 Remap(this Vector2 v, Vector2 iMin, Vector2 iMax, Vector2 oMin,
        Vector2 oMax) => Mathf.Remap(iMin, iMax, oMin, oMax, v);

    /// <inheritdoc cref="Mathf.Remap(Vector3,Vector3,Vector3,Vector3,Vector3)"/>
    [MethodImpl(Inline)]
    public static Vector3 Remap(this Vector3 v, Vector3 iMin, Vector3 iMax, Vector3 oMin,
        Vector3 oMax) => Mathf.Remap(iMin, iMax, oMin, oMax, v);

    /// <inheritdoc cref="Mathf.Remap(Vector4,Vector4,Vector4,Vector4,Vector4)"/>
    [MethodImpl(Inline)]
    public static Vector4 Remap(this Vector4 v, Vector4 iMin, Vector4 iMax, Vector4 oMin,
        Vector4 oMax) => Mathf.Remap(iMin, iMax, oMin, oMax, v);

    /// <inheritdoc cref="Mathf.Eerp(float,float,float)"/>
    [MethodImpl(Inline)] public static float Eerp(this float t, float a, float b) =>
        Mathf.Eerp(a, b, t);

    /// <inheritdoc cref="Mathf.InverseEerp(float,float,float)"/>
    [MethodImpl(Inline)] public static float InverseEerp(this float v, float a, float b) =>
        Mathf.InverseEerp(a, b, v);

    #endregion

    #region Vector Math

    /// <inheritdoc cref="Mathf.GetDirAndMagnitude(Vector2)"/>
    [MethodImpl(Inline)]
    public static (Vector2 dir, float magnitude ) GetDirAndMagnitude(this Vector2 v) =>
        Mathf.GetDirAndMagnitude(v);

    /// <inheritdoc cref="Mathf.GetDirAndMagnitude(Vector3)"/>
    [MethodImpl(Inline)]
    public static (Vector3 dir, float magnitude ) GetDirAndMagnitude(this Vector3 v) =>
        Mathf.GetDirAndMagnitude(v);

    /// <inheritdoc cref="Mathf.ClampMagnitude(Vector2,float,float)"/>
    [MethodImpl(Inline)]
    public static Vector2 ClampMagnitude(this Vector2 v, float min, float max) =>
        Mathf.ClampMagnitude(v, min, max);

    /// <inheritdoc cref="Mathf.ClampMagnitude(Vector3,float,float)"/>
    [MethodImpl(Inline)]
    public static Vector3 ClampMagnitude(this Vector3 v, float min, float max) =>
        Mathf.ClampMagnitude(v, min, max);

    #endregion

    #endregion

    #region Random

    public static float Sign(this Random r) => r.NextSingle() > 0.5f ? 1f : -1f;

    public static float NextSingle(this Random r, float min, float max) =>
        r.NextSingle() * (max - min) + min;

    /// <summary>Returns a random float between -1 and 1</summary>
    public static float NextUnit(this Random r) =>
        r.NextSingle(-1f, 1f);

    public static double NextDouble(this Random r, double min, double max) =>
        r.NextDouble() * (max - min) + min;

    public static float Angle(this Random r) => r.NextSingle() * Mathf.TAU;

    /// <summary>Returns a random point on the unit circle</summary>
    public static Vector2 OnUnitCircle(this Random r) =>
        Mathf.AngToDir(r.Angle());

    //  TODO: Validate
    /// <summary>Returns a random point inside the unit circle</summary>
    public static Vector2 InUnitCircle(this Random r) =>
        r.OnUnitCircle().Normalized() * MathF.Max(r.NextSingle(), r.NextSingle());

    public static Vector2 InUnitSquare(this Random r) => new(r.NextSingle(), r.NextSingle());


    /// <summary>Returns a random point inside the unit cube. Values are between 0 to 1</summary>
    public static Vector3 InUnitCube(this Random r) =>
        new(r.NextSingle(), r.NextSingle(), r.NextSingle());

    //  TODO: how?
    // /// <summary>Returns a random point on the unit sphere</summary>
    // public static Vector3 OnUnitSphere => UnityRandom.onUnitSphere;


    //  TODO: how?
    // /// <summary>Returns a random point inside the unit sphere</summary>
    // public static Vector3 InUnitSphere => UnityRandom.insideUnitSphere;

    //  TODO: Validate
    /// <summary>Returns a random uniformly distributed rotation</summary>
    public static Quaternion Rotation(this Random r)
    {
        float norm, w, x, y, z;

        do
        {
            w = r.NextUnit();
            x = r.NextUnit();
            y = r.NextUnit();
            z = r.NextUnit();
            norm = w * w + x * x + y * y + z * z;
        } while (norm is > 1f or 0f);

        norm = MathF.Sqrt(norm);
        return new(
            w / norm,
            x / norm,
            y / norm,
            z / norm
        );
    }

    #endregion

    #region Deconstruction

    [MethodImpl(Inline)]
    public static void Deconstruct(this in Quaternion q,
        out float x, out float y, out float z, out float w) =>
        (x, y, z, w) = (q.X, q.Y, q.Z, q.W);


    public static void Deconstruct(this in Vector4 v,
        out float x, out float y, out float z, out float w) =>
        (x, y, z, w) = (v.X, v.Y, v.Z, v.W);

    [MethodImpl(Inline)]
    public static void Deconstruct(this in Vector3 v,
        out float x, out float y, out float z) =>
        (x, y, z) = (v.X, v.Y, v.Z);

    [MethodImpl(Inline)]
    public static void Deconstruct(this in Vector2 v,
        out float x, out float y) => (x, y) = (v.X, v.Y);

    #endregion

    /// <summary>
    /// Get a int pair from Range 
    /// ^ sets the value as exclusive
    /// </summary>
    public static (int, int) ToExclusivePair(this Range range)
    {
        var (start, end) = (range.Start.Value, range.End.Value);
        return (
            range.Start.IsFromEnd ? start + 1 : start,
            range.End.IsFromEnd ? end - 1 : end
        );
    }
}