using System.Globalization;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Serialization;

namespace MathPlus;

/// <summary>A 3x3 matrix</summary>
/// 
[Serializable, DataContract]
[StructLayout(LayoutKind.Sequential)]
public readonly struct Matrix3X3
{
    public static readonly Matrix3X3 Identity = new(1, 0, 0, 0, 1, 0, 0, 0, 1);
    public static readonly Matrix3X3 Zero = new(0, 0, 0, 0, 0, 0, 0, 0, 0);

    [DataMember]
    public readonly float M11, M12, M13;

    [DataMember]
    public readonly float M21, M22, M23;

    [DataMember]
    public readonly float M31, M32, M33;

    public Matrix3X3(float m00, float m01, float m02, float m10, float m11, float m12,
        float m20, float m21, float m22)
    {
        (M11, M12, M13) = (m00, m01, m02);
        (M21, M22, M23) = (m10, m11, m12);
        (M31, M32, M33) = (m20, m21, m22);
    }

    public Matrix3X3(Vector3 col0, Vector3 col1, Vector3 col2)
    {
        M11 = col0.X;
        M21 = col0.Y;
        M31 = col0.Z;
        M12 = col1.X;
        M22 = col1.Y;
        M32 = col1.Z;
        M13 = col2.X;
        M23 = col2.Y;
        M33 = col2.Z;
    }

    public Matrix3X3(Matrix4x4 m)
    {
        (M11, M12, M13) = (m.M11, m.M12, m.M13);
        (M21, M22, M23) = (m.M21, m.M22, m.M23);
        (M31, M32, M33) = (m.M31, m.M32, m.M33);
    }

    public Matrix3X3(Quaternion q)
    {
        var m = Matrix4x4.CreateFromQuaternion(q);
        (M11, M12, M13) = (m.M11, m.M12, m.M13);
        (M21, M22, M23) = (m.M21, m.M22, m.M23);
        (M31, M32, M33) = (m.M31, m.M32, m.M33);
    }

    public static Matrix3X3 Scale(Vector3 s) =>
        new(
            s.X, 0, 0,
            0, s.Y, 0,
            0, 0, s.Z
        );

    public Matrix3X3 NormalizeColumns()
    {
        var lenX = Math.Sqrt(M11 * M11 + M21 * M21 + M31 * M31);
        var lenY = Math.Sqrt(M12 * M12 + M22 * M22 + M32 * M32);
        var lenZ = Math.Sqrt(M13 * M13 + M23 * M23 + M33 * M33);
        return new(
            (float) (M11 / lenX), (float) (M12 / lenY), (float) (M13 / lenZ),
            (float) (M21 / lenX), (float) (M22 / lenY), (float) (M23 / lenZ),
            (float) (M31 / lenX), (float) (M32 / lenY), (float) (M33 / lenZ));
    }

    public float this[int row, int column] =>
        (row, column) switch
        {
            (0, 0) => M11,
            (0, 1) => M12,
            (0, 2) => M13,
            (1, 0) => M21,
            (1, 1) => M22,
            (1, 2) => M23,
            (2, 0) => M31,
            (2, 1) => M32,
            (2, 2) => M33,
            _ => throw new IndexOutOfRangeException(
                $"Matrix row/column indices have to be from 0 to 2, got: ({row},{column})")
        };

    /// <summary>Returns the inverse of this matrix. Throws a division by zero exception if it's not invertible</summary>
    public Matrix3X3 Inverse
    {
        get
        {
            var a1212 = M22 * M33 - M23 * M32;
            var a0212 = M21 * M33 - M23 * M31;
            var a0112 = M21 * M32 - M22 * M31;
            var det = M11 * a1212 - M12 * a0212 + M13 * a0112;

            if (det == 0)
                throw new DivideByZeroException(
                    "The matrix is not invertible - its determinant is 0");

            return new Matrix3X3(
                a1212, M13 * M32 - M12 * M33, M12 * M23 - M13 * M22,
                -a0212, M11 * M33 - M13 * M31, M21 * M13 - M11 * M23,
                a0112, M31 * M12 - M11 * M32, M11 * M22 - M21 * M12
            ) / det;
        }
    }

    /// <summary>Returns the determinant of this matrix</summary>
    public float Determinant
    {
        get
        {
            var a1212 = M22 * M33 - M23 * M32;
            var a0212 = M21 * M33 - M23 * M31;
            var a0112 = M21 * M32 - M22 * M31;
            return M11 * a1212 - M12 * a0212 + M13 * a0112;
        }
    }

    public Matrix3X3 Transpose =>
        new(M11, M21, M31,
            M12, M22, M32,
            M13, M23, M33);

    public override string ToString() => ToStringMatrix().ToValueTableString();

    public string[,] ToStringMatrix() =>
        new[,]
        {
            {
                M11.ToString(CultureInfo.InvariantCulture),
                M12.ToString(CultureInfo.InvariantCulture),
                M13.ToString(CultureInfo.InvariantCulture),
            },
            {
                M21.ToString(CultureInfo.InvariantCulture),
                M22.ToString(CultureInfo.InvariantCulture),
                M23.ToString(CultureInfo.InvariantCulture),
            },
            {
                M31.ToString(CultureInfo.InvariantCulture),
                M32.ToString(CultureInfo.InvariantCulture),
                M33.ToString(CultureInfo.InvariantCulture),
            }
        };

    public static Matrix3X3 operator *(Matrix3X3 c, float v) =>
        new(c.M11 * v, c.M12 * v, c.M13 * v,
            c.M21 * v, c.M22 * v, c.M23 * v,
            c.M31 * v, c.M32 * v, c.M33 * v);

    public static Matrix3X3 operator /(Matrix3X3 c, float v) => c * (1f / v);

    public static Matrix3X3 operator *(Matrix3X3 a, Matrix3X3 b)
    {
        return new(
            GetEntry(0, 0), GetEntry(0, 1), GetEntry(0, 2),
            GetEntry(1, 0), GetEntry(1, 1), GetEntry(1, 2),
            GetEntry(2, 0), GetEntry(2, 1), GetEntry(2, 2)
        );

        float GetEntry(int r, int c) =>
            a[r, 0] * b[0, c] +
            a[r, 1] * b[1, c] +
            a[r, 2] * b[2, c];
    }

    public static Matrix3x1 operator *(Matrix3X3 c, Matrix3x1 m) =>
        new(m.m0 * c.M11 + m.m1 * c.M12 + m.m2 * c.M13,
            m.m0 * c.M21 + m.m1 * c.M22 + m.m2 * c.M23,
            m.m0 * c.M31 + m.m1 * c.M32 + m.m2 * c.M33);

    public static Vector3 operator *(Matrix3X3 c, Vector3 v) =>
        new(v.X * c.M11 + v.Y * c.M12 + v.Z * c.M13,
            v.X * c.M21 + v.Y * c.M22 + v.Z * c.M23,
            v.X * c.M31 + v.Y * c.M32 + v.Z * c.M33);

    public float AverageScale() =>
    (
        MathF.Sqrt(M11 * M11 + M21 * M21 + M31 * M31) +
        MathF.Sqrt(M12 * M12 + M22 * M22 + M32 * M32) +
        MathF.Sqrt(M13 * M13 + M23 * M23 + M33 * M33)
    ) / 3;
}