using System;
using System.Linq;
using System.Text;

namespace SLAE.DirectMethods
{
    internal static class Program
    {
        private static void Main()
        {
            // var a = new[,]
            // {
            //     {-2.8, -9.0, 8.1, 3.0},
            //     {4.6, -9.4, -3.6, -9.5},
            //     {8.4, 2.6, 2.9, 3.6},
            //     {3.6, -1.6, 0.4, 2.4}
            // };

            var a = new double[][]
            {
                new[] {-8.8, 0.7, 4.5, 3.0},
                new[] {-2.0, -0.4, 9.7, 8.6},
                new[] {-9.7, -9.5, 3.8, -3.5},
                new[] {-7.1, 1.6, 1.4, 2.9}
            };

            var x = new double[][]
            {
                new[] {1.0},
                new[] {2.0},
                new[] {3.0},
                new[] {4.0}
            };

            var b = a.Multiply(x);

            Console.WriteLine(a.ToPrettyString());

            var u = a.CreateCopy();
            var l = a.CreateCompatible();
            var p = a.CreateSwapMatrix();
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

                Console.WriteLine($"U:\n{u.ToPrettyString()}");
                Console.WriteLine($"L:\n{l.ToPrettyString()}");
                Console.WriteLine($"P:\n{p.ToPrettyString()}");
            }

            var lu = l.Multiply(u);
            Console.WriteLine($"L*U:\n{lu.ToPrettyString()}");

            a.ApplySwaps(p);
            Console.WriteLine($"P*A:\n{a.ToPrettyString()}");

            lu.Subtract(a);
            Console.WriteLine($"L*U-P*A:\n{lu.ToPrettyString()}");

            var determinant = l.DiagonalMatrixDeterminant();
            Console.WriteLine($"Determinant:\n{determinant}");

            Console.ReadKey();
        }
    }

    public static class MatrixExtensions
    {
        public static double[][] Multiply(this double[][] x, double[][] y)
        {
            var rows = x.Rows();
            var columns = y.Columns();
            var result = new double[rows][];
            for (var row = 0; row < rows; row++)
            {
                result[row] = new double[columns];
                for (var column = 0; column < columns; column++)
                {
                    result[row][column] = x.MultiplyRowByColumn(row, y, column);
                }
            }

            return result;
        }

        public static void Subtract(this double[][] x, double[][] y)
        {
            for (var row = 0; row < x.Rows(); row++)
            {
                for (var column = 0; column < x.Columns(); column++)
                {
                    x[row][column] -= y[row][column];
                }
            }
        }

        public static double DiagonalMatrixDeterminant(this double[][] matrix)
        {
            var result = 1.0;
            for (var row = 0; row < matrix.Rows(); row++)
            {
                result *= matrix[row][row];
            }

            return result;
        }

        private static double MultiplyRowByColumn(this double[][] x, int row, double[][] y, int column)
        {
            var result = 0.0;
            var columns = x.Columns();
            for (var i = 0; i < columns; i++)
            {
                result += x[row][i] * y[i][column];
            }

            return result;
        }

        public static string ToPrettyString(this double[][] x)
        {
            var builder = new StringBuilder();
            for (var row = 0; row < x.Rows(); row++)
            {
                for (var column = 0; column < x.Columns(); column++)
                {
                    builder.Append(x[row][column]);
                    builder.Append(";\t");
                }

                builder.Append('\n');
            }

            return builder.ToString();
        }

        public static string ToPrettyString(this int[] matrix)
        {
            var builder = new StringBuilder();
            foreach (var row in matrix)
            {
                builder.Append(row);
                builder.Append(";\n");
            }

            return builder.ToString();
        }

        public static int RowToSwapWith(this double[][] matrix, int column)
        {
            var targetRow = column;
            for (var row = column + 1; row < matrix.Rows(); row++)
            {
                if (Math.Abs(matrix[targetRow][column]) < Math.Abs(matrix[row][column]))
                {
                    targetRow = row;
                }
            }

            return targetRow;
        }

        public static void SwapRows(this double[][] matrix, int rowA, int rowB)
        {
            for (var column = 0; column < Columns(matrix); column++)
            {
                var temp = matrix[rowA][column];
                matrix[rowA][column] = matrix[rowB][column];
                matrix[rowB][column] = temp;
            }
        }

        public static void SubtractRow(this double[][] matrix, int subtractedRow)
        {
            for (var row = subtractedRow + 1; row < matrix.Rows(); row++)
            {
                var multiplier = matrix[row][subtractedRow] / matrix[subtractedRow][subtractedRow];
                matrix[row][subtractedRow] = 0.0;
                for (var column = subtractedRow + 1; column < matrix.Columns(); column++)
                {
                    matrix[row][column] -= multiplier * matrix[subtractedRow][column];
                }
            }
        }

        public static void NormalizeRow(this double[][] matrix, int row)
        {
            for (var column = row + 1; column < matrix.Columns(); column++)
            {
                matrix[row][column] /= matrix[row][row];
            }

            matrix[row][row] = 1.0;
        }

        public static void SwapRows(this int[] matrix, int row1, int row2)
        {
            var temp = matrix[row1];
            matrix[row1] = matrix[row2];
            matrix[row2] = temp;
        }

        public static void FillColumn(this double[][] l, double[][] u, int column)
        {
            for (var row = column; row < l.Columns(); row++)
            {
                l[row][column] = u[row][column];
            }
        }

        public static void ApplySwaps(this double[][] matrix, int[] swaps)
        {
            var result = new double[matrix.Rows()][];
            for (var i = 0; i < matrix.Rows(); i++)
            {
                var index = swaps[i];
                result[i] = matrix[index];
            }

            for (var i = 0; i < matrix.Rows(); i++)
            {
                matrix[i] = result[i];
            }
        }

        public static int[] CreateSwapMatrix(this double[][] matrix) => Enumerable.Range(0, matrix.Rows()).ToArray();

        public static double[][] CreateCopy(this double[][] matrix)
        {
            var result = matrix.CreateCompatible();
            for (var row = 0; row < matrix.Rows(); row++)
            {
                for (var column = 0; column < matrix.Columns(); column++)
                {
                    result[row][column] = matrix[row][column];
                }
            }

            return result;
        }

        public static double[][] CreateCompatible(this double[][] matrix)
        {
            var result = new double[matrix.Rows()][];
            for (var row = 0; row < matrix.Rows(); row++)
            {
                result[row] = new double[matrix.Columns()];
            }

            return result;
        }

        public static int Rows(this double[][] x) => x.Length;

        public static int Columns(this double[][] x) => x[0].Length;
    }
}
