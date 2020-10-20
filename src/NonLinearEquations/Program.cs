﻿using System;
using SLAE.DirectMethods;
using static System.Math;

namespace NonLinearEquations
{
    internal static class Program
    {
        private static void Main()
        {
            var xk = new[]
            {
                new[] {0.4},
                new[] {1.0}
            };

            var parameters = new GradientDescendParameters
            {
                Alpha = 1.0,
                Lambda = 0.5
            };

            var f = new[]
            {
                new Fun2[] {(x, y) => -1 + Sin(y + 0.5)},
                new Fun2[] {(x, y) => -Cos(x - 2)}
            };

            var fGrad = new[]
            {
                new Fun2[] {(_, _) => 0d, (_, y) => Cos(y + 0.5)},
                new Fun2[] {(x, _) => Sin(x - 2), (_, _) => 0d}
            };

            Fun2 minimizedF = (x, y) => Pow(f[0][0](x, y), 2) + Pow(f[1][0](x, y), 2);
            var minimizedGrad = new[]
            {
                new Fun2[] {(x, y) => -2 * Cos(x - 2) * Sin(x - 2)},
                new Fun2[] {(x, y) => 2 * (Sin(y + 0.5) - 1) * Cos(y + 0.5)}
            };

            var tau = fGrad.Evaluate2(xk);

            f.SolveWithNewtonMethod(fGrad, xk.CreateCopy());
            f.SolveWithSimpleIterationMethod(tau, xk.CreateCopy());
            minimizedF.SolveWithGradientDescendMethod(minimizedGrad, parameters, xk.CreateCopy());
        }
    }

    public delegate double Fun2(double x, double y);

    public class GradientDescendParameters
    {
        public double Alpha { get; set; } = 1.0;
        public double Lambda { get; set; } = 0.5;
    }

    public static class SimpleIterationMethod
    {
        public static void SolveWithSimpleIterationMethod(this Fun2[][] f, double[][] tau, double[][] xk)
        {
            Console.WriteLine("{0,5}|{1,15}|{2,15}|{3,15}|{4,15}|{5,15}",
                "i", "x", "y", "f1", "f2", "residual");

            for (var i = 1; i <= 100; i++)
            {
                var fValues = f.Evaluate2(xk);

                xk.Add(tau.Multiply(fValues));

                fValues = f.Evaluate2(xk);
                var residual = fValues.Norm3();

                Console.WriteLine("{0,5}|{1,15:G6}|{2,15:G6}|{3,15:G6}|{4,15:G6}|{5,15:G6}",
                    i, xk[0][0], xk[1][0], fValues[0][0], fValues[1][0], residual);
            }
        }
    }

    public static class NewtonMethod
    {
        public static void SolveWithNewtonMethod(this Fun2[][] f, Fun2[][] fGrad, double[][] xk)
        {
            Console.WriteLine(fGrad.Evaluate2(xk).Norm3());

            Console.WriteLine("{0,5}|{1,15}|{2,15}|{3,15}|{4,15}|{5,15}",
                "i", "x", "y", "f1", "f2", "residual");

            for (var i = 1; i <= 10; i++)
            {
                var fValues = f.Evaluate2(xk);
                var gradValues = fGrad.Evaluate2(xk);
                var inverseGradValues = gradValues.CreateInverse();

                xk.Subtract(inverseGradValues.Multiply(fValues));

                fValues = f.Evaluate2(xk);
                var residual = fValues.Norm3();

                Console.WriteLine("{0,5}|{1,15:G6}|{2,15:G6}|{3,15:G6}|{4,15:G6}|{5,15:G6}",
                    i, xk[0][0], xk[1][0], fValues[0][0], fValues[1][0], residual);
            }
        }
    }

    public static class GradientDescendMethod
    {
        public static void SolveWithGradientDescendMethod(this Fun2 f, Fun2[][] fGrad, GradientDescendParameters parameters,
            double[][] xk)
        {
            Console.WriteLine(fGrad.Evaluate2(xk).Norm3());

            Console.WriteLine("{0,5}|{1,15}|{2,15}|{3,15}|{4,15}",
                "i", "x", "y", "F", "alpha");

            for (var i = 1; i <= 100; i++)
            {
                var gradValues = fGrad.Evaluate2(xk);
                gradValues.Multiply(parameters.Alpha);
                var arg = xk.CreateCopy();
                arg.Subtract(gradValues);

                while (f.Evaluate2(arg) >= f.Evaluate2(xk))
                {
                    parameters.Alpha *= parameters.Lambda;

                    gradValues = fGrad.Evaluate2(xk);
                    gradValues.Multiply(parameters.Alpha);
                    arg = xk.CreateCopy();
                    arg.Subtract(gradValues);
                }

                xk = arg;
                var value = f.Evaluate2(xk);

                Console.WriteLine("{0,5}|{1,15:G6}|{2,15:G6}|{3,15:G6}|{4,15:G6}",
                    i, xk[0][0], xk[1][0], value, parameters.Alpha);
            }
        }
    }

    public static class Evaluation
    {
        public static double[][] Evaluate2(this Fun2[][] matrix, double[][] values)
        {
            var copy = matrix.CreateCompatible<Fun2, double>();
            for (var row = 0; row < matrix.Rows(); row++)
            {
                for (var column = 0; column < matrix.Columns(); column++)
                {
                    copy[row][column] = matrix[row][column].Evaluate2(values);
                }
            }

            return copy;
        }

        public static double Evaluate2(this Fun2 f, double[][] values)
        {
            var (x, y) = (values[0][0], values[1][0]);
            return f(x, y);
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

            // Check that LU = PA
            var lu = l.Multiply(u);
            lu.Subtract(pa);

            // Find A^(-1)
            var pe = pa.CreateOneMatrix();
            pe.ApplySwaps(p);
            var y2 = l.FindY(pe);
            var invertedA = u.FindX(y2);
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
    }

    public static class MatrixComparison
    {
        public static bool GreaterOrEqual(this double[][] matrix, double[][] other)
        {
            for (var i = 0; i < matrix.Rows(); i++)
            {
                for (var j = 0; j < matrix.Columns(); j++)
                {
                    if (matrix[i][j] < other[i][j])
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        public static bool CloserThan(this double[][] matrix, double[][] other)
        {
            for (var i = 0; i < matrix.Rows(); i++)
            {
                for (var j = 0; j < matrix.Columns(); j++)
                {
                    if (Abs(matrix[i][j]) > Abs(other[i][j]))
                    {
                        return false;
                    }
                }
            }

            return true;
        }
    }
}
