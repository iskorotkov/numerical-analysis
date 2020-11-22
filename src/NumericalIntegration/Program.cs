using System;
using System.Collections.Generic;
using MathExpressions;
using MathExpressions.Terms;

namespace NumericalIntegration
{
    public static class Program
    {
        private static void Main()
        {
            const double eps = 1e-8;

            var x = T.Var(0);
            var a = 1;
            var b = 2;
            
            // Variant 8
            var f = T.Exp(T.Value(5), x) * 6 * x + 3;
            // Variant 9
            //var f = T.Exp(T.Value(Math.E), x) * 6 * x + 3;

            new TrapezoidMethodSolver().Solve(f, a, b, eps);
            new ModifiedTrapezoidMethodSolver().Solve(f, a, b, eps);
            new SimpsonMethodSolver().Solve(f, a, b, eps);
            new Gauss3MethodSolver().Solve(f, a, b, eps);
        }
    }

    public class TrapezoidMethodSolver
    {
        private readonly RungeMethodEvaluator _rungeMethod = new();

        public void Solve(ITerm f, double a, double b, double eps)
        {
            Console.WriteLine("Trapezoid method");
            Console.WriteLine($"{"n",8}{"Square",16}{"Abs error",16}{"Error",16}");

            var s1 = f.Evaluate(a.AsMatrix()) + f.Evaluate(b.AsMatrix());
            double? previousResult = null;
            for (var n = 1;; n *= 2)
            {
                var step = (b - a) / (n + 1);
                var sum = 0d;
                for (var i = 1; i <= n; i++)
                {
                    var value = a + i * step;
                    sum += f.Evaluate(value.AsMatrix());
                }

                var square = step * (s1 / 2 + sum);
                _rungeMethod.Add(square);

                if (previousResult is { } p)
                {
                    var absError = Math.Abs(p - square) / 3;
                    var relativeError = absError / square;
                    Console.WriteLine($"{n,8}{square,16:g10}{absError,16:g10}{relativeError,16:g10}");

                    if (relativeError < eps)
                    {
                        var order = _rungeMethod.CalcApproximationOrder();
                        Console.WriteLine(order != null
                            ? $"Approximation order = {order}"
                            : "Can't calculate approximation order: need more iterations");
                        return;
                    }
                }
                else
                {
                    Console.WriteLine($"{n,8}{square,16:g10}{"-",16}{"-",16}");
                }

                previousResult = square;
            }
        }
    }

    public class ModifiedTrapezoidMethodSolver
    {
        private readonly RungeMethodEvaluator _rungeMethod = new();

        public void Solve(ITerm f, double a, double b, double eps)
        {
            Console.WriteLine("Trapezoid method (with splines)");
            Console.WriteLine($"{"n",8}{"Square",16}{"Abs error",16}{"Error",16}");

            var firstGrad = f.Grad();
            var s1 = f.Evaluate(a.AsMatrix()) + f.Evaluate(b.AsMatrix());
            var s2 = firstGrad.Evaluate(a.AsMatrix()) - firstGrad.Evaluate(b.AsMatrix());

            double? previousResult = null;
            for (var n = 2;; n *= 2)
            {
                var step = (b - a) / (n + 1);

                var sum = 0d;
                for (var i = 1; i <= n; i++)
                {
                    var arg = a + i * step;
                    sum += f.Evaluate(arg.AsMatrix());
                }

                var square = step * (1 * s1 / 2 + sum) + step * step * s2 / 12;
                _rungeMethod.Add(square);

                if (previousResult is { } p)
                {
                    var absError = Math.Abs(p - square) / 3;
                    var relativeError = absError / square;
                    Console.WriteLine($"{n,8}{square,16:g10}{absError,16:g10}{relativeError,16:g10}");

                    if (relativeError < eps)
                    {
                        var order = _rungeMethod.CalcApproximationOrder();
                        Console.WriteLine(order != null
                            ? $"Approximation order = {order}"
                            : "Can't calculate approximation order: need more iterations");
                        return;
                    }
                }
                else
                {
                    Console.WriteLine($"{n,8}{square,16:g10}{"-",16}{"-",16}");
                }

                previousResult = square;
            }
        }
    }

    public class SimpsonMethodSolver
    {
        private readonly RungeMethodEvaluator _rungeMethod = new();

        public void Solve(ITerm f, double a, double b, double eps)
        {
            Console.WriteLine("Simpson's method");
            Console.WriteLine($"{"n",8}{"Square",16}{"Abs error",16}{"Error",16}");

            var s1 = f.Evaluate(a.AsMatrix()) + f.Evaluate(b.AsMatrix());

            double? previousResult = null;
            for (var n = 2;; n *= 2)
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
                _rungeMethod.Add(square);

                if (previousResult is { } p)
                {
                    var absError = Math.Abs(p - square) / 15;
                    var relativeError = absError / square;
                    Console.WriteLine($"{n,8}{square,16:g10}{absError,16:g10}{relativeError,16:g10}");

                    if (relativeError < eps)
                    {
                        var order = _rungeMethod.CalcApproximationOrder();
                        Console.WriteLine(order != null
                            ? $"Approximation order = {order}"
                            : "Can't calculate approximation order: need more iterations");
                        return;
                    }
                }
                else
                {
                    Console.WriteLine($"{n,8}{square,16:g10}{"-",16}{"-",16}");
                }

                previousResult = square;
            }
        }
    }

    public class Gauss3MethodSolver
    {
        private readonly RungeMethodEvaluator _rungeMethod = new();

        public void Solve(ITerm f, double a, double b, double eps)
        {
            Console.WriteLine("Gauss-3 method");
            Console.WriteLine($"{"n",8}{"Square",16}{"Abs error",16}{"Error",16}");

            double? previousResult = null;
            for (var n = 1;; n *= 2)
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

                _rungeMethod.Add(square);

                if (previousResult is { } p)
                {
                    var absError = Math.Abs(p - square) / 15;
                    var relativeError = absError / square;
                    Console.WriteLine($"{n,8}{square,16:g10}{absError,16:g10}{relativeError,16:g10}");

                    if (relativeError < eps)
                    {
                        var order = _rungeMethod.CalcApproximationOrder();
                        Console.WriteLine(order != null
                            ? $"Approximation order = {order}"
                            : "Can't calculate approximation order: need more iterations");
                        return;
                    }
                }
                else
                {
                    Console.WriteLine($"{n,8}{square,16:g10}{"-",16}{"-",16}");
                }

                previousResult = square;
            }
        }
    }

    public class RungeMethodEvaluator
    {
        private readonly List<double> _history = new();

        public void Add(double result) => _history.Add(result);

        public double? CalcApproximationOrder()
        {
            if (_history.Count < 3)
            {
                return null;
            }

            return 1 * Math.Log((_history[^1] - _history[^3]) / (_history[^2] - _history[^3]) - 1) / Math.Log(0.5);
        }
    }

    public static class MatrixExtensions
    {
        public static double[][] AsMatrix(this double x) => new[] {new[] {x}};
    }
}
