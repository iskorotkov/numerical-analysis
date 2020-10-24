using System;
using System.Linq;
using MathExpressions;
using MathExpressions.Terms;
using SLAE.DirectMethods;

namespace NonLinearEquations
{
    internal static class Program
    {
        private static void Main()
        {
            const double epsilon = 1e-4;

            var values = new[]
            {
                new[] {-0.1},
                new[] {0.5}
            };

            var f1 = -T.Value(1) + T.Sin(T.Var(1) + 0.5) - T.Var(0);
            var f2 = -T.Cos(T.Var(0) - 2) - T.Var(1);

            var parameters = new GradientDescendParameters
            {
                Alpha = 1.0,
                Lambda = 0.5
            };

            var f = new[]
            {
                new[] {f1},
                new[] {f2}
            };
            var fGrad = new[]
            {
                new[] {f1.GradBy(new Var(0)), f1.GradBy(new Var(1))},
                new[] {f2.GradBy(new Var(0)), f2.GradBy(new Var(1))}
            };

            var fi1 = T.Var(0) + f1;
            var fi2 = T.Var(1) + f2;

            var fi = new[]
            {
                new[] {fi1},
                new[] {fi2}
            };
            var fiGrad = new[]
            {
                new[] {fi1.GradBy(new Var(0)), fi1.GradBy(new Var(1))},
                new[] {fi2.GradBy(new Var(0)), fi2.GradBy(new Var(1))}
            };

            var minF = f1 * f1 + f2 * f2;
            var minFGrad = new[]
            {
                new[] {minF.GradBy(new Var(0))},
                new[] {minF.GradBy(new Var(1))}
            };

            PrintHeader(values, parameters);

            var simpleIteration = new SimpleIterationSolver(epsilon);
            var newton = new NewtonSolver(epsilon);
            var gradientDescend = new GradientDescendSolver(epsilon, parameters);

            simpleIteration.Solve(f, fi, fiGrad, values.CreateCopy());
            newton.Solve(f, fGrad, values.CreateCopy());
            gradientDescend.Solve(minF, minFGrad, f, values.CreateCopy());
        }

        private static void PrintHeader(double[][] values, GradientDescendParameters parameters)
        {
            Console.WriteLine($"x:\n{values.ToPrettyString()}");
            Console.WriteLine($"\nParameters:\n\talpha = {parameters.Alpha}\n\tlambda = {parameters.Lambda}");
        }
    }

    public class SimpleIterationSolver
    {
        private readonly double _eps;

        public SimpleIterationSolver(double epsilon) => _eps = epsilon;

        public void Solve(ITerm[][] f, ITerm[][] fi, ITerm[][] fiGrad,
            double[][] values)
        {
            Console.WriteLine("\n=== Simple iteration method ===");

            var jakobi = fiGrad.Evaluate(values);
            Console.WriteLine("\nJakobi:");
            Console.WriteLine(jakobi.ToPrettyString());

            var q = jakobi.Norm3();
            Console.WriteLine($"\nJakobi norm = {q}\n");

            Console.WriteLine($"{"i",3}{"x",15}{"y",15}{"residual norm",15}{"f1",15}{"f2",15}{"jakobi norm",15}");

            var i = 0;
            double[][] xDiff;
            do
            {
                i++;

                var previousValues = values;
                values = fi.Evaluate(values);

                var fValues = f.Evaluate(values);
                var residualNorm = fValues.Norm3();
                var jakobiValues = fiGrad.Evaluate(values);
                var jakobiNorm = jakobiValues.Norm3();

                Console.WriteLine(
                    $"{i,3}{values[0][0],15:G6}{values[1][0],15:G6}{residualNorm,15:G6}{fValues[0][0],15:G6}{fValues[1][0],15:G6}{jakobiNorm,15:G6}");

                xDiff = previousValues;
                xDiff.Subtract(values);
            } while (!Converged(xDiff));
        }

        private bool Converged(double[][] xDiff) =>
            xDiff.All(row => row
                .Select(Math.Abs)
                .All(value => value < _eps));
    }

    public class GradientDescendParameters
    {
        public double Alpha { get; set; } = 1.0;
        public double Lambda { get; set; } = 0.5;
    }

    public class NewtonSolver
    {
        private readonly double _eps;

        public NewtonSolver(double epsilon) => _eps = epsilon;

        public void Solve(ITerm[][] f, ITerm[][] fGrad, double[][] values)
        {
            Console.WriteLine("\n=== Newton method ===\n");

            Console.WriteLine($"{"i",3}{"x",15}{"y",15}{"residual norm",15}{"f1",15}{"f2",15}");

            var i = 0;
            var fValues = f.Evaluate(values);
            double[][] previousFValues;
            do
            {
                i++;

                var gradValues = fGrad.Evaluate(values);
                var inverseGradValues = gradValues.CreateInverse();
                values.Subtract(inverseGradValues.Multiply(fValues));

                previousFValues = fValues;
                fValues = f.Evaluate(values);
                var residual = fValues.Norm3();

                Console.WriteLine(
                    $"{i,3}{values[0][0],15:G6}{values[1][0],15:G6}{residual,15:G6}{fValues[0][0],15:G6}{fValues[1][0],15:G6}");
            } while (!Converged(previousFValues));
        }

        private bool Converged(double[][] previousFValues) =>
            previousFValues.All(row => row
                .Select(Math.Abs)
                .All(value => value < _eps));
    }

    public class GradientDescendSolver
    {
        private readonly double _eps;
        private readonly GradientDescendParameters _parameters;

        public GradientDescendSolver(double epsilon, GradientDescendParameters parameters) =>
            (_eps, _parameters) = (epsilon, parameters);

        public void Solve(ITerm minF, ITerm[][] minFGrad, ITerm[][] f, double[][] values)
        {
            Console.WriteLine("\n=== Gradient descend method ===");

            var gradient = minFGrad.Evaluate(values);
            Console.WriteLine($"\nGradient vector:\n{gradient.ToPrettyString()}\n");

            Console.WriteLine($"{"i",3}{"x",15}{"y",15}{"alpha",15}{"residual",15}{"f1",15}{"f2",15}{"F",15}{"k",3}");

            var k = 0;
            var i = 0;
            double[][] xDiff;
            do
            {
                i++;

                var newValues = MakeIteration(minFGrad, values);

                while (minF.Evaluate(newValues) >= minF.Evaluate(values))
                {
                    k++;
                    _parameters.Alpha *= _parameters.Lambda;
                    newValues = MakeIteration(minFGrad, values);
                }

                var previousX = values;
                values = newValues;

                var fValues = f.Evaluate(values);
                var residual = fValues.Norm3();
                var minFValue = minF.Evaluate(values);

                Console.WriteLine(
                    $"{i,3}{values[0][0],15:G6}{values[1][0],15:G6}{_parameters.Alpha,15:G6}{residual,15:G6}{fValues[0][0],15:G6}{fValues[1][0],15:G6}{minFValue,15:G6}{k,3}");

                xDiff = previousX;
                xDiff.Subtract(values);
            } while (!Converged(xDiff));
        }

        private double[][] MakeIteration(ITerm[][] minFGrad, double[][] values)
        {
            var gradValues = minFGrad.Evaluate(values);
            gradValues.Multiply(_parameters.Alpha);

            var newValues = values.CreateCopy();
            newValues.Subtract(gradValues);
            return newValues;
        }

        private bool Converged(double[][] xDiff) => xDiff.Norm3() < _eps;
    }

    public static class Evaluation
    {
        public static double[][] Evaluate(this ITerm[][] matrix, double[][] values)
        {
            var copy = matrix.CreateCompatible<ITerm, double>();
            for (var row = 0; row < matrix.Rows(); row++)
            {
                for (var column = 0; column < matrix.Columns(); column++)
                {
                    copy[row][column] = matrix[row][column].Evaluate(values);
                }
            }

            return copy;
        }
    }

    public static class MatrixMath
    {
        public static double[][] CreateInverse(this double[][] a)
        {
            var l = a.CreateCompatible();
            var u = a.CreateCopy();
            var p = a.CreateSwapMatrix();

            // Make LU decomposition
            for (var row = 0; row < u.Rows(); row++)
            {
                var rowToSwap = u.RowToSwapWith(row);

                if (rowToSwap != row)
                {
                    u.SwapRows(row, rowToSwap);
                    l.SwapRows(row, rowToSwap);
                    p.SwapRows(row, rowToSwap);
                }

                l.FillColumn(u, row);

                u.SubtractRow(row);
                u.NormalizeRow(row);
            }

            var pa = a.CreateCopy();
            pa.ApplySwaps(p);

            // Find A^(-1)
            var pe = pa.CreateOneMatrix();
            pe.ApplySwaps(p);
            var y = l.FindY(pe);
            var invertedA = u.FindX(y);
            return invertedA;
        }

        public static void Multiply(this double[][] matrix, double value)
        {
            for (var i = 0; i < matrix.Rows(); i++)
            {
                for (var j = 0; j < matrix.Columns(); j++)
                {
                    matrix[i][j] *= value;
                }
            }
        }

        public static void Map(this double[][] matrix, Func<double, double> f)
        {
            for (var i = 0; i < matrix.Rows(); i++)
            {
                for (var j = 0; j < matrix.Columns(); j++)
                {
                    matrix[i][j] = f(matrix[i][j]);
                }
            }
        }
    }
}
