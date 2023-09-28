// by Freya Holmér (https://github.com/FreyaHolmer/MathPlus)

using System;
using System.Numerics;

namespace MathPlus;

public static class CharMatrix
{
    /// <summary>The characteristic matrix of a quadratic bézier curve</summary>
    public static readonly RationalMatrix3x3 QuadraticBezier = new(
        1, 0, 0,
        -2, 2, 0,
        1, -2, 1
    );

    /// <summary>The characteristic matrix of a cubic bézier curve</summary>
    public static readonly RationalMatrix4x4 CubicBezier = new(
        1, 0, 0, 0,
        -3, 3, 0, 0,
        3, -6, 3, 0,
        -1, 3, -3, 1
    );

    /// <summary>The characteristic matrix of a uniform cubic hermite curve</summary>
    public static readonly RationalMatrix4x4 CubicHermite = new(
        1, 0, 0, 0,
        0, 1, 0, 0,
        -3, -2, 3, -1,
        2, 1, -2, 1
    );

    public static readonly Polynomial[] CubicHermitePositionBasisFunctions =
    {
        GetBasisFunction(CubicHermite, 0),
        GetBasisFunction(CubicHermite, 2)
    };

    public static readonly Polynomial[] CubicHermiteVelocityBasisFunctions =
    {
        GetBasisFunction(CubicHermite, 1),
        GetBasisFunction(CubicHermite, 3)
    };

    /// <summary>The characteristic matrix of a uniform cubic catmull-rom curve</summary>
    public static readonly RationalMatrix4x4 CubicCatmullRom = new RationalMatrix4x4(
        0, 2, 0, 0,
        -1, 0, 1, 0,
        2, -5, 4, -1,
        -1, 3, -3, 1
    ) / 2;

    public static readonly Polynomial[] CubicCatmullRomBasisFunctions =
    {
        GetBasisFunction(CubicCatmullRom, 0),
        GetBasisFunction(CubicCatmullRom, 1),
        GetBasisFunction(CubicCatmullRom, 2),
        GetBasisFunction(CubicCatmullRom, 3)
    };

    /// <summary>The characteristic matrix of a uniform cubic B-spline curve</summary>
    public static readonly RationalMatrix4x4 CubicUniformBspline = new RationalMatrix4x4(
        1, 4, 1, 0,
        -3, 0, 3, 0,
        3, -6, 3, 0,
        -1, 3, -3, 1
    ) / 6;

    /// <summary>The inverse characteristic matrix of a quadratic bézier curve</summary>
    public static readonly RationalMatrix3x3 QuadraticBezierInverse = QuadraticBezier.Inverse;

    /// <summary>The inverse characteristic matrix of a cubic bézier curve</summary>
    public static readonly RationalMatrix4x4 CubicBezierInverse = CubicBezier.Inverse;

    /// <summary>The characteristic matrix of a uniform cubic hermite curve</summary>
    public static readonly RationalMatrix4x4 CubicHermiteInverse = CubicHermite.Inverse;

    /// <summary>The characteristic matrix of a uniform cubic catmull-rom curve</summary>
    public static readonly RationalMatrix4x4 CubicCatmullRomInverse = CubicCatmullRom.Inverse;

    /// <summary>The characteristic matrix of a uniform cubic B-spline curve</summary>
    public static readonly RationalMatrix4x4 CubicUniformBsplineInverse =
        CubicUniformBspline.Inverse;

    /// <summary>Returns the matrix to convert control points from one cubic spline to another, keeping the same curve intact</summary>
    /// <param name="from">The characteristic matrix of the spline to convert from</param>
    /// <param name="to">The characteristic matrix of the spline to convert from</param>
    public static RationalMatrix4x4 GetConversionMatrix(RationalMatrix4x4 from,
        RationalMatrix4x4 to) => to.Inverse * from;

    public static Matrix4x4 Create(float m00, float m01, float m02, float m03, float m10,
        float m11, float m12, float m13, float m20, float m21, float m22, float m23, float m30,
        float m31, float m32, float m33)
    {
        Matrix4x4 m;
        m.M11 = m00;
        m.M21 = m10;
        m.M31 = m20;
        m.M41 = m30;
        m.M12 = m01;
        m.M22 = m11;
        m.M32 = m21;
        m.M42 = m31;
        m.M13 = m02;
        m.M23 = m12;
        m.M33 = m22;
        m.M43 = m32;
        m.M14 = m03;
        m.M24 = m13;
        m.M34 = m23;
        m.M44 = m33;
        return m;
    }

    /// <summary>Returns the basis function (weight) for the given spline points by index <c>i</c>,
    /// equal to the t-matrix multiplied by the characteristic matrix</summary>
    /// <param name="c">The characteristic matrix to get the basis functions of</param>
    /// <param name="i">The point index to get the basis function of</param>
    public static Polynomial GetBasisFunction(RationalMatrix4x4 c, int i) =>
        i switch
        {
            0 => new Polynomial((float) c.m00, (float) c.m10, (float) c.m20, (float) c.m30),
            1 => new Polynomial((float) c.m01, (float) c.m11, (float) c.m21, (float) c.m31),
            2 => new Polynomial((float) c.m02, (float) c.m12, (float) c.m22, (float) c.m32),
            3 => new Polynomial((float) c.m03, (float) c.m13, (float) c.m23, (float) c.m33),
            _ => throw new IndexOutOfRangeException("Basis index needs to be between 0 and 3")
        };

    /// <inheritdoc cref="GetBasisFunction(RationalMatrix4x4,int)"/>
    public static Polynomial GetBasisFunction(Matrix4x4 c, int i) =>
        i switch
        {
            0 => new Polynomial(c.M11, c.M21, c.M31, c.M41),
            1 => new Polynomial(c.M12, c.M22, c.M32, c.M42),
            2 => new Polynomial(c.M13, c.M23, c.M33, c.M43),
            3 => new Polynomial(c.M14, c.M24, c.M34, c.M44),
            _ => throw new IndexOutOfRangeException("Basis index needs to be between 0 and 3")
        };
}