using System.Runtime.CompilerServices;
using System.Numerics;

namespace MathPlus;

public static class VectorMathExt
{
    const MethodImplOptions Inline = MethodImplOptions.AggressiveInlining;

    [MethodImpl(Inline)]
    public static float SqMag<TV, TVm>(this TVm vm, TV v) where TVm : struct, IVectorMath<TV> =>
        vm.Dot(v, v);

    [MethodImpl(Inline)]
    public static float SqDist<TV, TVm>(this TVm vm, TV a, TV b)
        where TVm : struct, IVectorMath<TV> => vm.SqMag(vm.Sub(b, a));

    [MethodImpl(Inline)]
    public static TV VecProject<TV, TVm>(this TVm vm, TV p, TV to)
        where TVm : struct, IVectorMath<TV> => vm.Mul(to, vm.BasisProject(p, to));

    [MethodImpl(Inline)]
    public static TV VecReject<TV, TVm>(this TVm vm, TV p, TV to)
        where TVm : struct, IVectorMath<TV> => vm.Sub(p, vm.VecProject(p, to));

    [MethodImpl(Inline)]
    public static float BasisProject<TV, TVm>(this TVm vm, TV p, TV to)
        where TVm : struct, IVectorMath<TV> => vm.Dot(p, to) / vm.Dot(to, to);

    [MethodImpl(Inline)]
    public static TV VecFromLineToPoint<TV, TVm>(this TVm vm, TV o, TV n, TV p)
        where TVm : struct, IVectorMath<TV> => vm.VecReject(vm.Sub(p, o), n);

    [MethodImpl(Inline)]
    public static float SqDistFromPointToLine<TV, TVm>(this TVm vm, TV o, TV n, TV p)
        where TVm : struct, IVectorMath<TV> => vm.SqMag(vm.VecFromLineToPoint(o, n, p));

    [MethodImpl(Inline)]
    public static TV GetPointAlongLine<TV, TVm>(this TVm vm, TV o, TV n, float t)
        where TVm : struct, IVectorMath<TV> => vm.Add(o, vm.Mul(n, t));

    [MethodImpl(Inline)]
    public static float ProjPointToLineSegmentTValue<TV, TVm>(this TVm vm, TV a, TV b, TV p)
        where TVm : struct, IVectorMath<TV> =>
        Mathf.Clamp01(vm.ProjPointToLineTValue(a, vm.Sub(b, a), p));

    [MethodImpl(Inline)]
    public static float ProjPointToLineTValue<TV, TVm>(this TVm vm, TV o, TV n, TV p)
        where TVm : struct, IVectorMath<TV> => vm.BasisProject(vm.Sub(p, o), n);

    [MethodImpl(Inline)]
    public static TV ProjPointToLine<TV, TVm>(this TVm vm, TV o, TV n, TV p)
        where TVm : struct, IVectorMath<TV> => vm.Add(o, vm.VecProject(vm.Sub(p, o), n));

    public static (float tA, float tB) ClosestPointBetweenLinesTValues<TV, TVm>(this TVm vm,
        TV aOrigin, TV aDir, TV bOrigin, TV bDir) where TVm : struct, IVectorMath<TV>
    {
        // source: https://math.stackexchange.com/questions/2213165/find-shortest-distance-between-lines-in-3d
        var e = vm.Sub(aOrigin, bOrigin);
        var be = vm.Dot(aDir, e);
        var de = vm.Dot(bDir, e);
        var bd = vm.Dot(aDir, bDir);
        var b2 = vm.Dot(aDir, aDir);
        var d2 = vm.Dot(bDir, bDir);
        var a = -b2 * d2 + bd * bd;
        var s = (-b2 * de + be * bd) / a;
        var t = (d2 * be - de * bd) / a;
        return (t, s);
    }

    public static bool TryIntersectSphereAtOrigin<TV, TVm>(this TVm vm, TV o, TV n, float r,
        out (float tMin, float tMax) tValues) where TVm : struct, IVectorMath<TV>
    {
        var nn = vm.Dot(n, n);
        if (nn <= 0f)
        {
            // vector has zero length, there's no direction
            tValues = default;
            return false;
        }

        var oo = vm.Dot(o, o);
        var on = vm.Dot(o, n);

        // quadratic terms
        double a = nn;
        double b = 2 * on;
        double c = oo - r * r;

        // try root solving
        var discriminant = b * b - 4 * a * c;
        if (discriminant < 0)
        {
            // no root, line is outside the circle
            tValues = default;
            return false;
        }

        var sign = b < 0 ? -1 : 1;
        var u = -b - sign * Math.Sqrt(discriminant);

        var tA = (float) (u / (2 * a));
        var tB = (float) (2 * c / u);
        tValues = tA < tB ? (tA, tB) : (tB, tA);
        return true;
    }
}

public interface IVectorMath<TV>
{
    const MethodImplOptions Inline = MethodImplOptions.AggressiveInlining;
    [MethodImpl(Inline)] TV Add(TV a, TV b);
    [MethodImpl(Inline)] TV Sub(TV a, TV b);
    [MethodImpl(Inline)] TV Mul(TV v, float c);
    [MethodImpl(Inline)] TV Mul(float c, TV v);
    [MethodImpl(Inline)] TV Div(TV v, float c);
    [MethodImpl(Inline)] float Dot(TV a, TV b);
    [MethodImpl(Inline)] float Mag(TV v);
    [MethodImpl(Inline)] float Dist(TV a, TV b);
    [MethodImpl(Inline)] TV Normalize(TV v);
    [MethodImpl(Inline)] TV Lerp(TV a, TV b, float t);
}

public struct VectorMath1D : IVectorMath<float>
{
    const MethodImplOptions Inline = MethodImplOptions.AggressiveInlining;
    [MethodImpl(Inline)] public float Add(float a, float b) => a + b;
    [MethodImpl(Inline)] public float Sub(float a, float b) => a - b;
    [MethodImpl(Inline)] public float Mul(float v, float c) => v * c;
    [MethodImpl(Inline)] public float Div(float v, float c) => v / c;
    [MethodImpl(Inline)] public float Dot(float a, float b) => a * b;
    [MethodImpl(Inline)] public float Mag(float v) => MathF.Abs(v);
    [MethodImpl(Inline)] public float Dist(float a, float b) => MathF.Abs(b - a);
    [MethodImpl(Inline)] public float Normalize(float v) => v < 0 ? -1 : 1;

    // shared implementations
    [MethodImpl(Inline)] public float Lerp(float a, float b, float t) => (1f - t) * a + t * b;
}

public struct VectorMath2D : IVectorMath<Vector2>
{
    const MethodImplOptions Inline = MethodImplOptions.AggressiveInlining;
    [MethodImpl(Inline)] public Vector2 Add(Vector2 a, Vector2 b) => new(a.X + b.X, a.Y + b.Y);
    [MethodImpl(Inline)] public Vector2 Sub(Vector2 a, Vector2 b) => new(a.X - b.X, a.Y - b.Y);
    [MethodImpl(Inline)] public Vector2 Mul(Vector2 v, float c) => new(v.X * c, v.Y * c);
    [MethodImpl(Inline)] public Vector2 Mul(float c, Vector2 v) => new(v.X * c, v.Y * c);
    [MethodImpl(Inline)] public Vector2 Div(Vector2 v, float c) => new(v.X / c, v.Y / c);
    [MethodImpl(Inline)] public float Dot(Vector2 a, Vector2 b) => a.X * b.X + a.Y * b.Y;
    [MethodImpl(Inline)] public float Mag(Vector2 v) => MathF.Sqrt(Dot(v, v));
    [MethodImpl(Inline)] public float Dist(Vector2 a, Vector2 b) => Mag(Sub(b, a));
    [MethodImpl(Inline)] public Vector2 Normalize(Vector2 v) => Div(v, Mag(v));

    // shared implementations
    [MethodImpl(Inline)] public Vector2 Lerp(Vector2 a, Vector2 b, float t)
    {
        var omt = 1f - t;
        return new Vector2(omt * a.X + t * b.X, omt * a.Y + t * b.Y);
    }
}

public struct VectorMath3D : IVectorMath<Vector3>
{
    const MethodImplOptions Inline = MethodImplOptions.AggressiveInlining;

    [MethodImpl(Inline)] public Vector3 Add(Vector3 a, Vector3 b) =>
        new(a.X + b.X, a.Y + b.Y, a.Z + b.Z);

    [MethodImpl(Inline)] public Vector3 Sub(Vector3 a, Vector3 b) =>
        new(a.X - b.X, a.Y - b.Y, a.Z - b.Z);

    [MethodImpl(Inline)] public Vector3 Mul(Vector3 v, float c) =>
        new(v.X * c, v.Y * c, v.Z * c);

    [MethodImpl(Inline)] public Vector3 Mul(float c, Vector3 v) =>
        new(v.X * c, v.Y * c, v.Z * c);

    [MethodImpl(Inline)] public Vector3 Div(Vector3 v, float c) =>
        new(v.X / c, v.Y / c, v.Z / c);

    [MethodImpl(Inline)] public float Dot(Vector3 a, Vector3 b) =>
        a.X * b.X + a.Y * b.Y + a.Z * b.Z;

    [MethodImpl(Inline)] public float Mag(Vector3 v) => MathF.Sqrt(Dot(v, v));
    [MethodImpl(Inline)] public float Dist(Vector3 a, Vector3 b) => Mag(Sub(b, a));
    [MethodImpl(Inline)] public Vector3 Normalize(Vector3 v) => Div(v, Mag(v));

    // shared implementations
    [MethodImpl(Inline)] public Vector3 Lerp(Vector3 a, Vector3 b, float t)
    {
        var omt = 1f - t;
        return new Vector3(omt * a.X + t * b.X, omt * a.Y + t * b.Y, omt * a.Z + t * b.Z);
    }
}

public struct VectorMath4D : IVectorMath<Vector4>
{
    const MethodImplOptions Inline = MethodImplOptions.AggressiveInlining;

    [MethodImpl(Inline)] public Vector4 Add(Vector4 a, Vector4 b) =>
        new(a.X + b.X, a.Y + b.Y, a.Z + b.Z, a.W + b.W);

    [MethodImpl(Inline)] public Vector4 Sub(Vector4 a, Vector4 b) =>
        new(a.X - b.X, a.Y - b.Y, a.Z - b.Z, a.W - b.W);

    [MethodImpl(Inline)] public Vector4 Mul(Vector4 v, float c) =>
        new(v.X * c, v.Y * c, v.Z * c, v.W * c);

    [MethodImpl(Inline)] public Vector4 Mul(float c, Vector4 v) =>
        new(v.X * c, v.Y * c, v.Z * c, v.W * c);

    [MethodImpl(Inline)] public Vector4 Div(Vector4 v, float c) =>
        new(v.X / c, v.Y / c, v.Z / c, v.W / c);

    [MethodImpl(Inline)] public float Dot(Vector4 a, Vector4 b) =>
        a.X * b.X + a.Y * b.Y + a.Z * b.Z + a.W * b.W;

    [MethodImpl(Inline)] public float Mag(Vector4 v) => MathF.Sqrt(Dot(v, v));
    [MethodImpl(Inline)] public float Dist(Vector4 a, Vector4 b) => Mag(Sub(b, a));
    [MethodImpl(Inline)] public Vector4 Normalize(Vector4 v) => Div(v, Mag(v));

    // shared implementations
    [MethodImpl(Inline)] public Vector4 Lerp(Vector4 a, Vector4 b, float t)
    {
        var omt = 1f - t;
        return new Vector4(omt * a.X + t * b.X, omt * a.Y + t * b.Y, omt * a.Z + t * b.Z,
            omt * a.W + t * b.W);
    }
}

// todo: quaternions, as a treat. this is untested and unported basically lol
public struct VectorMathQuat : IVectorMath<Quaternion>
{
    const MethodImplOptions Inline = MethodImplOptions.AggressiveInlining;

    [MethodImpl(Inline)] public Quaternion Add(Quaternion a, Quaternion b) =>
        new(a.X + b.X, a.Y + b.Y, a.Z + b.Z, a.W + b.W);

    [MethodImpl(Inline)] public Quaternion Sub(Quaternion a, Quaternion b) =>
        new(a.X - b.X, a.Y - b.Y, a.Z - b.Z, a.W - b.W);

    [MethodImpl(Inline)] public Quaternion Mul(Quaternion v, float c) =>
        new(v.X * c, v.Y * c, v.Z * c, v.W * c);

    [MethodImpl(Inline)] public Quaternion Mul(float c, Quaternion v) =>
        new(v.X * c, v.Y * c, v.Z * c, v.W * c);

    [MethodImpl(Inline)] public Quaternion Div(Quaternion v, float c) =>
        new(v.X / c, v.Y / c, v.Z / c, v.W / c);

    [MethodImpl(Inline)] public float Dot(Quaternion a, Quaternion b) =>
        a.X * b.X + a.Y * b.Y + a.Z * b.Z + a.W * b.W;

    [MethodImpl(Inline)] public float Mag(Quaternion v) => MathF.Sqrt(Dot(v, v));
    [MethodImpl(Inline)] public float Dist(Quaternion a, Quaternion b) => Mathf.Angle(a, b);
    [MethodImpl(Inline)] public Quaternion Normalize(Quaternion v) => Div(v, Mag(v));

    // shared implementations
    [MethodImpl(Inline)] public Quaternion Lerp(Quaternion a, Quaternion b, float t)
    {
        var omt = 1f - t;
        return new Quaternion(omt * a.X + t * b.X, omt * a.Y + t * b.Y, omt * a.Z + t * b.Z,
            omt * a.W + t * b.W);
    }
}