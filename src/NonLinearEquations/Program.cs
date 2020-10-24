using System;
using MathExpressions;
using MathExpressions.Terms;
using SLAE.DirectMethods;

namespace NonLinearEquations
{
    internal static class Program
    {
        private static void Main()
        {
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

            PrintHeader(values, parameters);

            f.SolveWithSimpleIterationMethod(fi, fiGrad, values.CreateCopy());

            f.SolveWithNewtonMethod(fGrad, values.CreateCopy());

            var minF = f1 * f1 + f2 * f2;
            var minFGrad = new[]
            {
                new[] {minF.GradBy(new Var(0))},
                new[] {minF.GradBy(new Var(1))}
            };

            minF.SolveWithGradientDescendMethod(minFGrad, f, parameters, values.CreateCopy());
        }

        private static void PrintHeader(double[][] values, GradientDescendParameters parameters)
        {
            Console.WriteLine($"x:\n{values.ToPrettyString()}");
            Console.WriteLine($"\nParameters:\n\talpha = {parameters.Alpha}\n\tlambda = {parameters.Lambda}");
        }
    }

    public class GradientDescendParameters
    {
        public double Alpha { get; set; } = 1.0;
        public double Lambda { get; set; } = 0.5;
    }

    public static class SimpleIterationMethod
    {
        public static void SolveWithSimpleIterationMethod(this ITerm[][] f, ITerm[][] fi, ITerm[][] fiGrad,
            double[][] values)
        {
            Console.WriteLine("\n=== Simple iteration method ===");

            var jakobi = fiGrad.Evaluate(values);
            Console.WriteLine("\nJakobi:");
            Console.WriteLine(jakobi.ToPrettyString());
            Console.WriteLine($"\nJakobi norm = {jakobi.Norm3()}\n");

            Console.WriteLine($"{"i",3}{"x",15}{"y",15}{"residual norm",15}{"f1",15}{"f2",15}{"jakobi norm",15}");

            for (var i = 1; i <= 20; i++)
            {
                values = fi.Evaluate(values);

                var fValues = f.Evaluate(values);
                var residualNorm = fValues.Norm3();
                var jakobiValues = fiGrad.Evaluate(values);
                var jakobiNorm = jakobiValues.Norm3();

                Console.WriteLine(
                    $"{i,3}{values[0][0],15:G6}{values[1][0],15:G6}{residualNorm,15:G6}{fValues[0][0],15:G6}{fValues[1][0],15:G6}{jakobiNorm,15:G6}");
            }
        }
    }

    public static class NewtonMethod
    {
        public static void SolveWithNewtonMethod(this ITerm[][] f, ITerm[][] fGrad, double[][] values)
        {
            Console.WriteLine("\n=== Newton method ===\n");

            Console.WriteLine($"{"i",3}{"x",15}{"y",15}{"residual norm",15}{"f1",15}{"f2",15}");

            var fValues = f.Evaluate(values);
            for (var i = 1; i <= 10; i++)
            {
                var gradValues = fGrad.Evaluate(values);
                var inverseGradValues = gradValues.CreateInverse();
                values.Subtract(inverseGradValues.Multiply(fValues));

                fValues = f.Evaluate(values);
                var residual = fValues.Norm3();

                Console.WriteLine(
                    $"{i,3}{values[0][0],15:G6}{values[1][0],15:G6}{residual,15:G6}{fValues[0][0],15:G6}{fValues[1][0],15:G6}");
            }
        }
    }

    public static class GradientDescendMethod
    {
        public static void SolveWithGradientDescendMethod(this ITerm minF, ITerm[][] minFGrad, ITerm[][] f,
            GradientDescendParameters parameters, double[][] values)
        {
            Console.WriteLine("\n=== Gradient descend method ===");

            var gradient = minFGrad.Evaluate(values);
            Console.WriteLine($"\nGradient vector:\n{gradient.ToPrettyString()}\n");

            Console.WriteLine($"{"i",3}{"x",15}{"y",15}{"alpha",15}{"residual",15}{"f1",15}{"f2",15}{"F",15}{"k",3}");

            var k = 0;
            for (var i = 1; i <= 20; i++)
            {
                var newValues = MakeIteration(minFGrad, parameters, values);

                while (minF.Evaluate(newValues) >= minF.Evaluate(values))
                {
                    k++;
                    parameters.Alpha *= parameters.Lambda;
                    newValues = minFGrad.MakeIteration(parameters, values);
                }

                values = newValues;
                var fValues = f.Evaluate(values);
                var residual = fValues.Norm3();
                var minFValue = minF.Evaluate(values);

                Console.WriteLine(
                    $"{i,3}{values[0][0],15:G6}{values[1][0],15:G6}{parameters.Alpha,15:G6}{residual,15:G6}{fValues[0][0],15:G6}{fValues[1][0],15:G6}{minFValue,15:G6}{k,3}");
            }
        }

        private static double[][] MakeIteration(this ITerm[][] minFGrad, GradientDescendParameters parameters,
            double[][] values)
        {
            var gradValues = minFGrad.Evaluate(values);
            gradValues.Multiply(parameters.Alpha);

            var newValues = values.CreateCopy();
            newValues.Subtract(gradValues);
            return newValues;
        }
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
