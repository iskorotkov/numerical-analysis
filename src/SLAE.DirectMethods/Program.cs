using System;
using System.IO;
using System.Linq;
using System.Text;

namespace SLAE.DirectMethods
{
    internal static class Program
    {
        private const string InputFile = "data/in.txt";
        private const string OutputFile = "data/out.txt";

        private static void Main()
        {
            // Read A
            double[][] a;
            using (var reader = new StreamReader(InputFile))
            {
                a = MatrixIO.ReadMatrix(reader);
            }

            // x
            var x = new[] {new[] {1.0}, new[] {2.0}, new[] {3.0}, new[] {4.0}};

            // Calculate b = A*x
            var b = a.Multiply(x);

            using TextWriter writer = new StreamWriter(OutputFile);

            // Print input values
            writer.WriteLine("Variant = 8b\n");
            writer.WriteLine($"x:\n{x.ToPrettyString()}\n");
            writer.WriteLine($"A:\n{a.ToPrettyString()}\n");
            writer.WriteLine($"b:\n{b.ToPrettyString()}\n");

            var l = a.CreateCompatible();
            var u = a.CreateCopy();
            var p = a.CreateSwapMatrix();
            var rank = l.Rows();

            // Make LU decomposition
            for (var row = 0; row < u.Rows(); row++)
            {
                var rowToSwap = u.RowToSwapWith(row);
                writer.WriteLine($"k = {row + 1}; m = {rowToSwap + 1};\n");

                if (Math.Abs(u[rowToSwap][rowToSwap]) < 1e-6)
                {
                    rank--;
                }

                if (rowToSwap != row)
                {
                    u.SwapRows(row, rowToSwap);
                    l.SwapRows(row, rowToSwap);
                    p.SwapRows(row, rowToSwap);
                }

                l.FillColumn(u, row);

                u.SubtractRow(row);
                u.NormalizeRow(row);

                writer.WriteLine($"U:\n{u.ToPrettyString()}\n");
                writer.WriteLine($"L:\n{l.ToPrettyString()}\n");
            }

            writer.WriteLine($"P:\n{p.ToPrettyString()}\n");
            writer.WriteLine($"Rank = {rank}\n");

            var pa = a.CreateCopy();
            pa.ApplySwaps(p);

            var pb = b.CreateCopy();
            pb.ApplySwaps(p);

            // Check that LU = PA
            var lu = l.Multiply(u);
            lu.Subtract(pa);
            writer.WriteLine($"L*U-P*A:\n{lu.ToPrettyString()}\n");

            // Print determinant
            var determinant = l.DiagonalMatrixDeterminant();
            writer.WriteLine($"Determinant = {determinant}\n");

            // Compute x
            var y = l.FindY(pb);
            var computedX = u.FindX(y);
            writer.WriteLine($"x:\n{computedX.ToPrettyString()}\n");

            // Find A^(-1)
            var e = pa.CreateOneMatrix();
            var y2 = l.FindY(e);
            var invertedA = u.FindX(y2);
            writer.WriteLine($"A^(-1):\n{invertedA.ToPrettyString()}\n");

            // Check that A*A^(-1) = E
            var aByInvertedA = pa.Multiply(invertedA);
            writer.WriteLine($"A*A^(-1):\n{aByInvertedA.ToPrettyString()}\n");

            // Find condition numbers of matrix A
            writer.WriteLine("Condition number of matrix A:");
            writer.WriteLine($"cond-1(A) = {pa.Cond1(invertedA)}");
            writer.WriteLine($"cond-2(A) = {pa.Cond2(invertedA)}");
            writer.WriteLine($"cond-3(A) = {pa.Cond3(invertedA)}\n");

            // Find A*x-b
            var dif = a.Multiply(computedX);
            dif.Subtract(b);
            writer.WriteLine($"A*x-b:\n{dif.ToPrettyString()}\n");

            // Find relative error
            var dx = computedX.CreateCopy();
            dx.Subtract(x);
            dx.Divide(x);
            writer.WriteLine($"Relative error:\n{dx.ToPrettyString()}");
        }
    }

    public static class MatrixIO
    {
        public static double[][] ReadMatrix(StreamReader reader)
        {
            var dimension = int.Parse(reader.ReadLine()!);
            var matrix = CreateEmptyMatrix(dimension, dimension);
            for (var row = 0; row < dimension; row++)
            {
                var line = reader.ReadLine();
                var values = line!.Split(';', ',', ' ', '\t')
                    .Where(s => !string.IsNullOrWhiteSpace(s))
                    .Select(double.Parse)
                    .ToArray();
                for (var column = 0; column < dimension; column++)
                {
                    matrix[row][column] = values[column];
                }
            }

            return matrix;
        }

        private static double[][] CreateEmptyMatrix(int rows, int columns)
        {
            var result = new double[rows][];
            for (var row = 0; row < rows; row++)
            {
                result[row] = new double[columns];
            }

            return result;
        }
    }

    public static class JakobiRotation
    {
        public static double[][] CreateRotatedDiagonalMatrix(this double[][] matrix, double eps)
        {
            while (matrix.Distance() > eps)
            {
                for (var i = 0; i < matrix.Rows() - 1; i++)
                {
                    for (var j = i + 1; j < matrix.Rows(); j++)
                    {
                        matrix = matrix.CreateRotatedMatrix(i, j, eps);
                    }
                }
            }

            return matrix;
        }

        private static double[][] CreateRotatedMatrix(this double[][] matrix, int i, int j, double eps)
        {
            if (Math.Abs(matrix[i][j]) >= eps)
            {
                double c;
                double s;
                if (Math.Abs(matrix[i][i] - matrix[j][j]) < 1e-6)
                {
                    const double theta = Math.PI / 4;
                    c = Math.Cos(theta);
                    s = Math.Sin(theta);
                }
                else
                {
                    var tau = (matrix[i][i] - matrix[j][j]) / 2 / matrix[i][j];
                    var t = Math.Sign(tau) / (Math.Abs(tau) + Math.Sqrt(1 + Math.Pow(tau, 2.0)));
                    c = 1 / Math.Sqrt(1 + Math.Pow(t, 2.0));
                    s = c * t;
                }

                var r = matrix.CreateRotationMatrix(i, j, c, s);
                var rTransposed = r.CreateTransposed();
                return rTransposed.Multiply(matrix).Multiply(r);
            }

            return matrix;
        }

        private static double[][] CreateRotationMatrix(this double[][] matrix, int i, int j, double c, double s)
        {
            var result = matrix.CreateOneMatrix();

            result[i][i] = result[j][j] = c;

            result[i][j] = -s;
            result[j][i] = s;

            return result;
        }

        private static double Distance(this double[][] matrix)
        {
            var result = 0.0;
            for (var row = 0; row < matrix.Rows(); row++)
            {
                for (var column = 0; column < matrix.Columns(); column++)
                {
                    if (row != column)
                    {
                        result += Math.Pow(matrix[row][column], 2.0);
                    }
                }
            }

            return result;
        }
    }

    public static class MatrixMath
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

        public static void Divide(this double[][] x, double[][] y)
        {
            for (var row = 0; row < x.Rows(); row++)
            {
                for (var column = 0; column < x.Columns(); column++)
                {
                    x[row][column] /= y[row][column];
                }
            }
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
    }

    public static class MatrixFormatting
    {
        public static string ToPrettyString(this double[][] x)
        {
            var builder = new StringBuilder();
            for (var row = 0; row < x.Rows(); row++)
            {
                for (var column = 0; column < x.Columns(); column++)
                {
                    builder.AppendFormat("{0,15:g6};", x[row][column]);
                }

                if (row != x.Rows() - 1)
                {
                    builder.Append('\n');
                }
            }

            return builder.ToString();
        }

        public static string ToPrettyString(this int[] matrix)
        {
            var builder = new StringBuilder();
            for (var row = 0; row < matrix.Length; row++)
            {
                builder.AppendFormat("{0,15:g6};", matrix[row] + 1);
                if (row != matrix.Length - 1)
                {
                    builder.Append('\n');
                }
            }

            return builder.ToString();
        }
    }

    public static class MatrixNorms
    {
        public static double Cond1(this double[][] matrix, double[][] invertedMatrix) =>
            matrix.Norm1() * invertedMatrix.Norm1();

        public static double Cond2(this double[][] matrix, double[][] invertedMatrix) =>
            matrix.Norm2() * invertedMatrix.Norm2();

        public static double Cond3(this double[][] matrix, double[][] invertedMatrix) =>
            matrix.Norm3() * invertedMatrix.Norm3();

        private static double Norm1(this double[][] matrix)
        {
            var result = 0.0;
            for (var row = 0; row < matrix.Rows(); row++)
            {
                var sum = 0.0;
                for (var column = 0; column < matrix.Columns(); column++)
                {
                    sum += Math.Abs(matrix[row][column]);
                }

                if (sum > result)
                {
                    result = sum;
                }
            }

            return result;
        }

        private static double Norm2(this double[][] matrix)
        {
            var result = 0.0;
            for (var column = 0; column < matrix.Rows(); column++)
            {
                var sum = 0.0;
                for (var row = 0; row < matrix.Columns(); row++)
                {
                    sum += Math.Abs(matrix[row][column]);
                }

                if (sum > result)
                {
                    result = sum;
                }
            }

            return result;
        }

        private static double Norm3(this double[][] matrix)
        {
            var transposedMatrixByMatrix = matrix.CreateTransposed().Multiply(matrix);
            var diagonalMatrix = transposedMatrixByMatrix.CreateRotatedDiagonalMatrix(1e-6);
            return FindLargestEigenvalue(diagonalMatrix);
        }

        private static double FindLargestEigenvalue(double[][] diagonalMatrix)
        {
            var result = 0.0;
            for (var row = 0; row < diagonalMatrix.Rows(); row++)
            {
                if (Math.Abs(diagonalMatrix[row][row]) > result)
                {
                    result = Math.Abs(diagonalMatrix[row][row]);
                }
            }

            return Math.Sqrt(result);
        }
    }

    public static class MatrixRows
    {
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
            for (var column = 0; column < matrix.Columns(); column++)
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
    }

    public static class MatrixCreation
    {
        public static double[][] CreateTransposed(this double[][] matrix)
        {
            var rows = matrix.Columns();
            var columns = matrix.Rows();
            var result = new double[rows][];
            for (var row = 0; row < rows; row++)
            {
                result[row] = new double[columns];
                for (var column = 0; column < columns; column++)
                {
                    result[row][column] = matrix[column][row];
                }
            }

            return result;
        }

        public static double[][] CreateOneMatrix(this double[][] matrix)
        {
            var result = matrix.CreateCompatible();
            for (var row = 0; row < matrix.Rows(); row++)
            {
                result[row][row] = 1.0;
            }

            return result;
        }

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

        public static int[] CreateSwapMatrix(this double[][] matrix) => Enumerable.Range(0, matrix.Rows()).ToArray();
    }

    public static class SLAESolving
    {
        public static double[][] FindY(this double[][] l, double[][] b)
        {
            var y = b.CreateCopy();
            for (var i = 0; i < b.Columns(); i++)
            {
                for (var row = 0; row < l.Rows(); row++)
                {
                    for (var column = 0; column < row; column++)
                    {
                        y[row][i] -= l[row][column] * y[column][i];
                    }

                    y[row][i] /= l[row][row];
                }
            }

            return y;
        }

        public static double[][] FindX(this double[][] u, double[][] y)
        {
            var x = y.CreateCopy();
            for (var i = 0; i < y.Columns(); i++)
            {
                for (var row = u.Rows() - 1; row >= 0; row--)
                {
                    for (var column = row + 1; column < u.Columns(); column++)
                    {
                        x[row][i] -= u[row][column] * x[column][i];
                    }
                }
            }

            return x;
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

        public static void FillColumn(this double[][] l, double[][] u, int column)
        {
            for (var row = column; row < l.Columns(); row++)
            {
                l[row][column] = u[row][column];
            }
        }
    }

    public static class MatrixExtensions
    {
        public static int Rows(this double[][] x) => x.Length;
        public static int Columns(this double[][] x) => x[0].Length;
    }
}