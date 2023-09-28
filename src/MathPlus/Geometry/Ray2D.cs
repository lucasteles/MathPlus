using System.Runtime.CompilerServices;
using System.Numerics;

namespace MathPlus;

/// <summary>Similar to Unity's Ray2D, except this one allows you to not normalize the direction
/// which saves performance as well as allows you to work at different scales</summary>
[Serializable]
public struct Ray2D : ILinear2D
{
    const MethodImplOptions Inline = MethodImplOptions.AggressiveInlining;

    /// <summary>The origin of the ray</summary>
    public Vector2 Origin { [MethodImpl(Inline)] get; }

    /// <summary>The direction of the ray. Note: Ray2D allows non-normalized direction vectors</summary>
    public Vector2 Dir { [MethodImpl(Inline)] get; }

    /// <summary>Returns a normalized version of this ray. Normalized rays ensure t-values correspond to distance</summary>
    public Line2D Normalized => new Line2D(((ILinear2D) this).Origin, ((ILinear2D) this).Dir);

    /// <summary>Creates a 2D Ray. Note: direction does not have to be normalized, but if it is, the t-value will correspond to distance along the ray</summary>
    /// <param name="origin">The origin of the ray</param>
    /// <param name="dir">The direction of the ray. It does not have to be normalized, but if it is, the t-value when sampling will correspond to distance along the ray</param>
    public Ray2D(Vector2 origin, Vector2 dir) =>
        (Origin, Dir) = (origin, dir);

    /// <summary>Implicitly casts a Unity ray to a MathPlus ray</summary>
    /// <param name="ray">The ray to cast to a Unity ray</param>
    public static implicit operator Ray2D(Ray ray) =>
        new(ray.Origin.ToXY(), ray.Direction.ToXY());

    #region Internal interface stuff for generic line tests

    [MethodImpl(Inline)] bool ILinear2D.IsValidTValue(float t) => t >= 0;

    [MethodImpl(Inline)] float ILinear2D.ClampTValue(float t) => t < 0 ? 0 : t;

    #endregion
}