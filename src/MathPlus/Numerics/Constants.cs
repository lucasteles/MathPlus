using System.Numerics;
using System.Runtime.CompilerServices;

namespace MathPlus;

public static class Vector2V
{
    public static Vector2 Zero
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Vector2.Zero;
    }

    public static Vector2 One
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Vector2.One;
    }

    public static Vector2 Up
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Vector2.UnitY;
    }

    public static Vector2 Down
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(0.0f, -1.0f);
    }

    public static Vector2 Right
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Vector2.UnitX;
    }

    public static Vector2 Left
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(-1.0f, 0.0f);
    }

    public static Vector2 PositiveInfinity
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(float.PositiveInfinity, float.PositiveInfinity);
    }

    public static Vector2 NegativeInfinity
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(float.NegativeInfinity, float.NegativeInfinity);
    }
}

public static class Vector3V
{
    public static Vector3 Zero
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Vector3.Zero;
    }

    public static Vector3 One
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Vector3.One;
    }

    public static Vector3 Up
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Vector3.UnitY;
    }

    public static Vector3 Down
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(0f, -1f, 0f);
    }

    public static Vector3 Right
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Vector3.UnitX;
    }

    public static Vector3 Left
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(-1f, 0f, 0f);
    }

    public static Vector3 Forward
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Vector3.UnitZ;
    }

    public static Vector3 Back
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(0f, 0f, -1f);
    }

    public static Vector3 PositiveInfinity
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(float.PositiveInfinity, float.PositiveInfinity, float.PositiveInfinity);
    }

    public static Vector3 NegativeInfinity
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(float.NegativeInfinity, float.NegativeInfinity, float.NegativeInfinity);
    }
}