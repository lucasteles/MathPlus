﻿// by Freya Holmér (https://github.com/FreyaHolmer/MathPlus)

using System.Runtime.CompilerServices;
using System.Numerics;
using static MathPlus.Mathf;

namespace MathPlus;

/// <summary>Core intersection test functions.
/// Note: these are pretty esoteric, generally it's easier to use the instance methods in each shape,
/// such as <c>myLine.Intersect(otherThing)</c></summary>
public static partial class IntersectionTest
{
    // internal
    const float ParallelDeterminantThreshold = 0.00001f;
    const MethodImplOptions Inline = MethodImplOptions.AggressiveInlining;

    /// <summary>Returns whether or not these infinite lines intersect, and if they do, also returns the t-value for each infinite line</summary>
    /// <param name="aOrigin">First line origin</param>
    /// <param name="aDir">First line direction (does not have to be normalized)</param>
    /// <param name="bOrigin">Second line origin</param>
    /// <param name="bDir">Second line direction (does not have to be normalized)</param>
    /// <param name="tA">The t-value along the first line, where the intersection happened</param>
    /// <param name="tB">The t-value along the second line, where the intersection happened</param>
    public static bool LinearTValues(Vector2 aOrigin, Vector2 aDir, Vector2 bOrigin, Vector2 bDir,
        out float tA, out float tB)
    {
        var d = Determinant(aDir, bDir);
        if (Abs(d) < ParallelDeterminantThreshold)
        {
            tA = tB = default;
            return false;
        }

        var aToB = bOrigin - aOrigin;
        tA = Determinant(aToB, bDir) / d;
        tB = Determinant(aToB, aDir) / d;
        return true;
    }

    // based on https://stackoverflow.com/questions/1073336/circle-line-segment-collision-detection-algorithm
    /// <summary>Returns the intersections between an infinite line and a circle in the form of t-values along the line where the intersections lie. Or none, if there are none</summary>
    /// <param name="lineOrigin">Line origin</param>
    /// <param name="lineDir">Line direction (does not have to be normalized)</param>
    /// <param name="circleOrigin">Center or the circle</param>
    /// <param name="radius">Radius of the circle</param>
    public static ResultsMax2<float> LinearCircleTValues(Vector2 lineOrigin, Vector2 lineDir,
        Vector2 circleOrigin, float radius)
    {
        var circleToLineOrigin = lineOrigin - circleOrigin;
        var a = Vector2.Dot(lineDir, lineDir); // ray len sq
        var b = 2 * Vector2.Dot(circleToLineOrigin, lineDir);
        var c = Vector2.Dot(circleToLineOrigin, circleToLineOrigin) - radius.Square();
        var discriminant = b * b - 4 * a * c;
        if (discriminant > 0)
        {
            discriminant = Sqrt(discriminant);
            if (discriminant < 0.00001f) // line is tangent to the circle, one intersection
                return new(-b / (2 * a));
            var tA = (-b + discriminant) / (2 * a);
            var tB = (-b - discriminant) / (2 * a); // line has two intersections
            return new(tA, tB);
        }

        return default; // line doesn't hit it at all
    }

    /// <summary>Returns whether or not two circles intersect, and the two intersection points (if they exist)</summary>
    /// <param name="aPos">The position of the first circle</param>
    /// <param name="aRadius">The radius of the first circle</param>
    /// <param name="bPos">The position of the second circle</param>
    /// <param name="bRadius">The radius of the second circle</param>
    public static ResultsMax2<Vector2> CirclesIntersectionPoints(Vector2 aPos, float aRadius,
        Vector2 bPos, float bRadius)
    {
        var distSq = DistanceSquared(aPos, bPos);
        var dist = Sqrt(distSq);
        var differentPosition = dist > 0.00001f;
        var maxRad = Max(aRadius, bRadius);
        var minRad = Min(aRadius, bRadius);
        var ringsTouching = MathF.Abs(dist - maxRad) < minRad;

        if (ringsTouching && differentPosition)
        {
            var aRadSq = aRadius * aRadius;
            var bRadSq = bRadius * bRadius;
            var lateralOffset = (distSq - bRadSq + aRadSq) / (2 * dist);
            var normalOffset = (0.5f / dist) *
                               Sqrt(4 * distSq * aRadSq - (distSq - bRadSq + aRadSq).Square());
            var tangent = (bPos - aPos) / dist;
            var normal = tangent.Rotate90CCW();
            var chordCenter = aPos + tangent * lateralOffset;
            if (normalOffset < 0.00001f)
                return new(chordCenter); // double intersection at one point
            return new( // two intersections
                chordCenter + normal * normalOffset,
                chordCenter - normal * normalOffset
            );
        }

        return default; // no intersections
    }

    /// <summary>Returns whether or not two circles overlap</summary>
    /// <param name="aPos">The position of the first circle</param>
    /// <param name="aRadius">The radius of the first circle</param>
    /// <param name="bPos">The position of the second circle</param>
    /// <param name="bRadius">The radius of the second circle</param>
    public static bool CirclesOverlap(Vector2 aPos, float aRadius, Vector2 bPos, float bRadius)
    {
        var dist = Vector2.Distance(aPos, bPos);
        var maxRad = Max(aRadius, bRadius);
        var minRad = Min(aRadius, bRadius);
        return MathF.Abs(dist - maxRad) < minRad;
    }

    /// <summary>Returns whether or not a line passes through a box centered at (0,0)</summary>
    /// <param name="extents">Box extents/"radius" per axis</param>
    /// <param name="pt">A point in the line</param>
    /// <param name="dir">The direction of the line</param>
    public static bool LineRectOverlap(Vector2 extents, Vector2 pt, Vector2 dir)
    {
        var corner = new Vector2(extents.X, extents.Y * -Sign(dir.X * dir.Y));
        return SignAsInt(Determinant(dir, corner - pt)) !=
               SignAsInt(Determinant(dir, -corner - pt));
    }

    /// <summary>Returns the intersection points of a line passing through a box</summary>
    /// <param name="center">Center of the box</param>
    /// <param name="extents">Box extents/"radius" per axis</param>
    /// <param name="pt">A point in the line</param>
    /// <param name="dir">The direction of the line</param>
    public static ResultsMax2<Vector2> LinearRectPoints(Vector2 center, Vector2 extents, Vector2 pt,
        Vector2 dir)
    {
        const float flatThresh = 0.000001f;

        // place the line relative to the box
        pt.X -= center.X;
        pt.Y -= center.Y;

        // Vertical line
        if (dir.X.Abs() < flatThresh)
        {
            if (pt.X.Abs() <= extents.X) // inside - two intersections
                return new(
                    new(center.X + pt.X, center.Y - extents.Y),
                    new(center.X + pt.X, center.Y + extents.Y));

            return default; // outside the box
        }

        // Horizontal line
        if (dir.Y.Abs() < flatThresh)
        {
            if (pt.Y.Abs() <= extents.Y) // inside - two intersections
                return new(
                    new(center.X - extents.X, center.Y + pt.Y),
                    new(center.X + extents.X, center.Y + pt.Y));

            return default; // outside the box
        }


        // slope intercept form y = ax+b
        var a = dir.Y / dir.X;
        var b = pt.Y - pt.X * a;

        // y coords on vertical lines
        var xpy = a * extents.X + b; // x = extents.x
        var xny = -a * extents.X + b; // x = -extents.x
        // x coords on horizontal lines
        var ypx = (extents.Y - b) / a; // y = extents.y
        var ynx = (-extents.Y - b) / a; // y = -extents.y

        // validity checks
        var xp = Abs(xpy) <= extents.Y;
        var xn = Abs(xny) <= extents.Y;
        var yp = Abs(ypx) <= extents.X;
        var yn = Abs(ynx) <= extents.X;

        if ((xp || xn || yp || yn) == false)
            return default; // no intersections

        float ax, ay, bx, by;
        if (a > 0)
        {
            // positive slope means we group results in (x,y) and (-x,-y)
            ax = xp ? extents.X : ypx;
            ay = xp ? xpy : extents.Y;
            bx = xn ? -extents.X : ynx;
            by = xn ? xny : -extents.Y;
        }
        else
        {
            // negative slope means we group results in (x,-y) and (-x,y)
            ax = xp ? extents.X : ynx;
            ay = xp ? xpy : -extents.Y;
            bx = xn ? -extents.X : ypx;
            by = xn ? xny : extents.Y;
        }

        // if the points are very close, this means we hit a corner and we should return only one point
        if (Abs(ax - bx) + Abs(ay - by) < 0.000001f)
            return new(new(center.X + ax, center.Y + ay));

        // else, two points
        return new(
            new(center.X + ax, center.Y + ay),
            new(center.X + bx, center.Y + by));
    }

    /// <summary>Returns whether or not two discs overlap. Unlike circles, discs overlap even if one is smaller and is completely inside the other</summary>
    /// <param name="aPos">The position of the first disc</param>
    /// <param name="aRadius">The radius of the first disc</param>
    /// <param name="bPos">The position of the second disc</param>
    /// <param name="bRadius">The radius of the second disc</param>
    [MethodImpl(Inline)]
    public static bool DiscsOverlap(Vector2 aPos, float aRadius, Vector2 bPos, float bRadius) =>
        DistanceSquared(aPos, bPos) <= (aRadius + bRadius).Square();
}