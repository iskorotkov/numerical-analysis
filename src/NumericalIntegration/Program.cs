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
            var f = T.Exp(T.Value(5), x) - 6 * x + 3;
            // Variant 9
            // var f = T.Exp(T.Value(Math.E), x) - 6 * x + 3;

            var headerTemplate = $"{"n",8}{"h",8}{"Integral",16}{"Error",16}{"Order",16}";
            var valuesTemplate = "{0,8}{1,8:f4}{2,16:g8}{3,16:g8}{4,16:f4}";

            new TrapezoidMethodSolver().Solve(f, a, b, eps, headerTemplate, valuesTemplate);
            new ModifiedTrapezoidMethodSolver().Solve(f, a, b, eps, headerTemplate, valuesTemplate);
            new SimpsonMethodSolver().Solve(f, a, b, eps, headerTemplate, valuesTemplate);
            new Gauss3MethodSolver().Solve(f, a, b, eps, headerTemplate, valuesTemplate);
        }
    }

    public class TrapezoidMethodSolver
    {
        private readonly RungeMethodEvaluator _rungeMethod = new();

        public void Solve(ITerm f, double a, double b, double eps, string headerTemplate, string valuesTemplate)
        {
            Console.WriteLine("Trapezoid method");
            Console.WriteLine(headerTemplate);

            var s1 = f.Evaluate(a.AsMatrix()) + f.Evaluate(b.AsMatrix());
            double? previousResult = null;
            for (var n = 1;; n *= 2)
            {
                var step = (b - a) / n;
                var sum = 0d;
                for (var i = 1; i < n; i++)
                {
                    var arg = a + i * step;
                    sum += f.Evaluate(arg.AsMatrix());
                }

                var value = step * (s1 / 2 + sum);
                _rungeMethod.Add(value);

                if (previousResult is { } p)
                {
                    var absError = Math.Abs(p - value) / 3;
                    var error = absError / value;
                    object order = _rungeMethod.CalcApproximationOrder();
                    Console.WriteLine(valuesTemplate, n, step, value, error, order ?? "-");

                    if (Math.Abs(error) < eps)
                    {
                        Console.WriteLine($"Integral = {value}");
                        Console.WriteLine($"Function evaluated {n + 1} times");
                        break;
                    }
                }
                else
                {
                    Console.WriteLine(valuesTemplate, n, step, value, "-", "-");
                }

                previousResult = value;
            }
        }
    }

    public class ModifiedTrapezoidMethodSolver
    {
        private readonly RungeMethodEvaluator _rungeMethod = new();

        public void Solve(ITerm f, double a, double b, double eps, string headerTemplate, string valuesTemplate)
        {
            Console.WriteLine("Trapezoid method (with splines)");
            Console.WriteLine(headerTemplate);

            var firstGrad = f.Grad();
            var s1 = f.Evaluate(a.AsMatrix()) + f.Evaluate(b.AsMatrix());
            var s2 = firstGrad.Evaluate(a.AsMatrix()) - firstGrad.Evaluate(b.AsMatrix());

            double? previousResult = null;
            for (var n = 1;; n *= 2)
            {
                var step = (b - a) / n;

                var sum = 0d;
                for (var i = 1; i < n; i++)
                {
                    var arg = a + i * step;
                    sum += f.Evaluate(arg.AsMatrix());
                }

                var value = step * (s1 / 2 + sum) + step * step * s2 / 12;
                _rungeMethod.Add(value);

                if (previousResult is { } p)
                {
                    var absError = Math.Abs(p - value) / 3;
                    var error = absError / value;
                    object order = _rungeMethod.CalcApproximationOrder();
                    Console.WriteLine(valuesTemplate, n, step, value, error, order ?? "-");

                    if (Math.Abs(error) < eps)
                    {
                        Console.WriteLine($"Integral = {value}");
                        Console.WriteLine($"Function evaluated {n + 1} times");
                        break;
                    }
                }
                else
                {
                    Console.WriteLine(valuesTemplate, n, step, value, "-", "-");
                }

                previousResult = value;
            }
        }
    }

    public class SimpsonMethodSolver
    {
        private readonly RungeMethodEvaluator _rungeMethod = new();

        public void Solve(ITerm f, double a, double b, double eps, string headerTemplate, string valuesTemplate)
        {
            Console.WriteLine("Simpson's method");
            Console.WriteLine(headerTemplate);

            var s1 = f.Evaluate(a.AsMatrix()) + f.Evaluate(b.AsMatrix());

            double? previousResult = null;
            for (var n = 1;; n *= 2)
            {
                var step = (b - a) / n;

                var (subareas, substep) = (n * 2, step / 2);

                var s2 = 0d;
                for (var i = 1; i <= subareas / 2; i++)
                {
                    var arg = a + substep * (2 * i - 1);
                    s2 += f.Evaluate(arg.AsMatrix());
                }

                var s3 = 0d;
                for (var i = 1; i <= subareas / 2 - 1; i++)
                {
                    var arg = a + substep * (2 * i);
                    s3 += f.Evaluate(arg.AsMatrix());
                }

                var value = substep * (s1 + 4 * s2 + 2 * s3) / 3;
                _rungeMethod.Add(value);

                if (previousResult is { } p)
                {
                    var absError = Math.Abs(p - value) / 15;
                    var error = absError / value;
                    object order = _rungeMethod.CalcApproximationOrder();
                    Console.WriteLine(valuesTemplate, n, step, value, error, order ?? "-");

                    if (Math.Abs(error) < eps)
                    {
                        Console.WriteLine($"Integral = {value}");
                        Console.WriteLine($"Function evaluated {2 * n + 1} times");
                        break;
                    }
                }
                else
                {
                    Console.WriteLine(valuesTemplate, n, step, value, "-", "-");
                }

                previousResult = value;
            }
        }
    }

    public class Gauss3MethodSolver
    {
        private readonly RungeMethodEvaluator _rungeMethod = new();

        public void Solve(ITerm f, double a, double b, double eps, string headerTemplate, string valuesTemplate)
        {
            Console.WriteLine("Gauss-3 method");
            Console.WriteLine(headerTemplate);

            double? previousResult = null;
            var evaluatedTimes = 0;
            for (var n = 1;; n *= 2)
            {
                var step = (b - a) / n;
                var value = 0d;
                for (var i = 0; i < n; i++)
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

                    value += amplitude * (a0 * f.Evaluate(x0.AsMatrix())
                                          + a1 * f.Evaluate(x1.AsMatrix())
                                          + a2 * f.Evaluate(x2.AsMatrix()));
                    evaluatedTimes += 3;
                }

                _rungeMethod.Add(value);

                if (previousResult is { } p)
                {
                    var absError = Math.Abs(p - value) / 15;
                    var error = absError / value;
                    object order = _rungeMethod.CalcApproximationOrder();
                    Console.WriteLine(valuesTemplate, n, step, value, error, order ?? "-");

                    if (Math.Abs(error) < eps)
                    {
                        Console.WriteLine($"Integral = {value}");
                        Console.WriteLine($"Function evaluated {evaluatedTimes} times");
                        break;
                    }
                }
                else
                {
                    Console.WriteLine(valuesTemplate, n, step, value, "-", "-");
                }

                previousResult = value;
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

            return Math.Log((_history[^1] - _history[^3]) / (_history[^2] - _history[^3]) - 1) / Math.Log(0.5);
        }
    }

    public static class MatrixExtensions
    {
        public static double[][] AsMatrix(this double x) => new[] { new[] { x } };
    }
}
