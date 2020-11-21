using System;
using MathExpressions;
using MathExpressions.Terms;

namespace NumericalIntegration
{
    public static class Program
    {
        private static void Main()
        {
            const double eps = 1e-8;

            var f = T.Sin(T.Var(0));
            var a = 0;
            var b = Math.PI;

            new TrapezoidMethodSolver().Solve(f, a, b, eps);
            new ModifiedTrapezoidMethodSolver().Solve(f, a, b, eps);
            new SimpsonMethodSolver().Solve(f, a, b, eps);
            new Gauss3MethodSolver().Solve(f, a, b, eps);
        }
    }

    public class TrapezoidMethodSolver
    {
        public void Solve(ITerm f, double a, double b, double eps)
        {
            var s1 = f.Evaluate(a.AsMatrix()) + f.Evaluate(b.AsMatrix());
            for (var n = 1; n <= 10; n++)
            {
                var step = (b - a) / (n + 1);
                var sum = 0d;
                for (var i = 1; i <= n; i++)
                {
                    var value = a + i * step;
                    sum += f.Evaluate(value.AsMatrix());
                }

                var square = step * (s1 / 2 + sum);
                Console.WriteLine($"n = {n}, square = {square}");
            }
        }
    }

    public class ModifiedTrapezoidMethodSolver
    {
        public void Solve(ITerm f, double a, double b, double eps)
        {
            var firstGrad = f.Grad();
            var s1 = f.Evaluate(a.AsMatrix()) + f.Evaluate(b.AsMatrix());
            var s2 = firstGrad.Evaluate(a.AsMatrix()) - firstGrad.Evaluate(b.AsMatrix());

            for (var n = 2; n <= 10; n += 2)
            {
                var step = (b - a) / (n + 1);

                var sum = 0d;
                for (var i = 1; i <= n; i++)
                {
                    var arg = a + i * step;
                    sum += f.Evaluate(arg.AsMatrix());
                }

                var square = step * (1 * s1 / 2 + sum) + step * step * s2 / 12;

                Console.WriteLine($"n = {n}, square = {square}");
            }
        }
    }

    public class SimpsonMethodSolver
    {
        public void Solve(ITerm f, double a, double b, double eps)
        {
            var s1 = f.Evaluate(a.AsMatrix()) + f.Evaluate(b.AsMatrix());
            for (var n = 2; n <= 10; n += 2)
            {
                var step = (b - a) / (n + 1);

                var s2 = 0d;
                for (var i = 1; i <= n / 2; i++)
                {
                    var value = a + step * (2 * i - 1);
                    s2 += f.Evaluate(value.AsMatrix());
                }

                var s3 = 0d;
                for (var i = 1; i <= n / 2 - 1; i++)
                {
                    var value = a + step * (2 * i);
                    s3 += f.Evaluate(value.AsMatrix());
                }

                var square = step * (s1 + 8 * s2 + 4 * s3) / 6;

                Console.WriteLine($"n = {n}, square = {square}");
            }
        }
    }

    public class Gauss3MethodSolver
    {
        public void Solve(ITerm f, double a, double b, double eps)
        {
            for (var n = 1; n <= 10; n++)
            {
                var step = (b - a) / (n + 1);
                var square = 0d;
                for (var i = 0; i <= n; i++)
                {
                    var start = a + step * i;
                    var end = a + step * (i + 1);

                    var amplitude = (end - start) / 2;
                    var center = start + amplitude;

                    var x0 = center - Math.Sqrt(0.6) * amplitude;
                    var x1 = center;
                    var x2 = center + Math.Sqrt(0.6) * amplitude;
                    var a0 = 5d / 9;
                    var a1 = 8d / 9;
                    var a2 = a0;

                    square += amplitude * (a0 * f.Evaluate(x0.AsMatrix())
                                           + a1 * f.Evaluate(x1.AsMatrix())
                                           + a2 * f.Evaluate(x2.AsMatrix()));
                }

                Console.WriteLine($"n = {n}, square = {square}");
            }
        }
    }

    public static class MatrixExtensions
    {
        public static double[][] AsMatrix(this double x) => new[] {new[] {x}};
    }
}
