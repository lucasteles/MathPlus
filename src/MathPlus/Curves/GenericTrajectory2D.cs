using System.Numerics;

public class GenericTrajectory2D
{
    public Vector2[] derivatives;

    public GenericTrajectory2D(params Vector2[] derivatives) => this.derivatives = derivatives;

    public Vector2 GetPosition(float time)
    {
        var pt = derivatives[0];
        for (var i = 1; i < derivatives.Length; i++)
        {
            var scale = MathF.Pow(time, i) / MathPlus.Mathf.Factorial((uint) i);
            pt += scale * derivatives[i];
        }

        return pt;
    }
}