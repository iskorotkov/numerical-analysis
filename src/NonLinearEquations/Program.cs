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

            var x = T.Var(0);
            var y = T.Var(1);

            var f1 = -T.Value(1) + T.Sin(y + 0.5);
            var f2 = -T.Cos(x - 2);

            var f = new[]
            {
                new[] {f1},
                new[] {f2}
            };
            var fGrad = new[]
            {
                new[] {f1.GradBy((Var) x), f1.GradBy((Var) y)},
                new[] {f2.GradBy((Var) x), f2.GradBy((Var) y)}
            };

            f.SolveWithSimpleIterationMethod(fGrad, values.CreateCopy());

            f.SolveWithNewtonMethod(fGrad, values.CreateCopy());

            var minimizedF = f1 * f1 + f2 * f2;
            var minimizedGrad = new[]
            {
                new[] {minimizedF.GradBy((Var) x)},
                new[] {minimizedF.GradBy((Var) y)}
            };

            var parameters = new GradientDescendParameters
            {
                Alpha = 1.0,
                Lambda = 0.5
            };
            minimizedF.SolveWithGradientDescendMethod(minimizedGrad, parameters, values.CreateCopy());
        }
    }

    public class GradientDescendParameters
    {
        public double Alpha { get; set; } = 1.0;
        public double Lambda { get; set; } = 0.5;
    }

    public static class SimpleIterationMethod
    {
        public static void SolveWithSimpleIterationMethod(this ITerm[][] f, ITerm[][] fGrad,
            double[][] values)
        {
            var tau = fGrad.Evaluate(values);
            Console.WriteLine("tau:");
            Console.WriteLine(tau.ToPrettyString());
            Console.WriteLine($"Norm={tau.Norm3()}");

            Console.WriteLine("{0,5}|{1,15}|{2,15}|{3,15}|{4,15}|{5,15}",
                "i", "x", "y", "f1", "f2", "residual");

            for (var i = 1; i <= 10; i++)
            {
                var fValues = f.Evaluate(values);

                values.Add(tau.Multiply(fValues));

                fValues = f.Evaluate(values);
                var residual = fValues.Norm3();

                Console.WriteLine("{0,5}|{1,15:G6}|{2,15:G6}|{3,15:G6}|{4,15:G6}|{5,15:G6}",
                    i, values[0][0], values[1][0], fValues[0][0], fValues[1][0], residual);
            }
        }
    }

    public static class NewtonMethod
    {
        public static void SolveWithNewtonMethod(this ITerm[][] f, ITerm[][] fGrad, double[][] values)
        {
            Console.WriteLine("{0,5}|{1,15}|{2,15}|{3,15}|{4,15}|{5,15}",
                "i", "x", "y", "f1", "f2", "residual");

            for (var i = 1; i <= 10; i++)
            {
                var fValues = f.Evaluate(values);
                var gradValues = fGrad.Evaluate(values);
                var inverseGradValues = gradValues.CreateInverse();

                values.Subtract(inverseGradValues.Multiply(fValues));

                fValues = f.Evaluate(values);
                var residual = fValues.Norm3();

                Console.WriteLine("{0,5}|{1,15:G6}|{2,15:G6}|{3,15:G6}|{4,15:G6}|{5,15:G6}",
                    i, values[0][0], values[1][0], fValues[0][0], fValues[1][0], residual);
            }
        }
    }

    public static class GradientDescendMethod
    {
        public static void SolveWithGradientDescendMethod(this ITerm f, ITerm[][] fGrad,
            GradientDescendParameters parameters, double[][] values)
        {
            Console.WriteLine(fGrad.Evaluate(values).Norm3());

            Console.WriteLine("{0,5}|{1,15}|{2,15}|{3,15}|{4,15}",
                "i", "x", "y", "F", "alpha");

            for (var i = 1; i <= 100; i++)
            {
                var gradValues = fGrad.Evaluate(values);
                gradValues.Multiply(parameters.Alpha);
                var arg = values.CreateCopy();
                arg.Subtract(gradValues);

                while (f.Evaluate(arg) >= f.Evaluate(values))
                {
                    parameters.Alpha *= parameters.Lambda;

                    gradValues = fGrad.Evaluate(values);
                    gradValues.Multiply(parameters.Alpha);
                    arg = values.CreateCopy();
                    arg.Subtract(gradValues);
                }

                values = arg;
                var value = f.Evaluate(values);

                Console.WriteLine("{0,5}|{1,15:G6}|{2,15:G6}|{3,15:G6}|{4,15:G6}",
                    i, values[0][0], values[1][0], value, parameters.Alpha);
            }
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
