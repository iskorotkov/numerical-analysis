using System;
using System.IO;
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
            const string outputFileName = "../../../data/output.txt";

            // Variables that are used in equations (needed for correct auto differentiation)
            var x = T.Var(0);
            var y = T.Var(1);

            // sample
            // Functions f(x, y)
            var f1 = -T.Value(1) + T.Sin(T.Var(1) + 0.5) - T.Var(0);
            var f2 = -T.Cos(T.Var(0) - 2) - T.Var(1);
            // Functions fi(x, y) (for simple iteration method)
            var fi1 = -T.Value(1) + T.Sin(T.Var(1) + 0.5);
            var fi2 = -T.Cos(T.Var(0) - 2);
            // Approximated root of system of equations
            var values = new[]
            {
                new[] {-0.1},
                new[] {0.5}
            };

            // variant 8
            // Functions f(x, y)
            //var f1 = 2 * y - T.Cos(x + 1);
            //var f2 = x + T.Sin(y) + 0.4;
            // Functions fi(x, y) (for simple iteration method)
            //var fi1 = -T.Sin(y) - 0.4;
            //var fi2 = T.Cos(x + 1) / 2;
            // Approximated root of system of equations
            //var values = new[]
            //{
            //    new[] {-0.9},
            //    new[] {0.5}
            //};

            // variant 9
            // Functions f(x, y)
            //var f1 = T.Cos(x + 0.5) - y - 2;
            //var f2 = T.Sin(y) - 2 * x - 1;
            // Functions fi(x, y) (for simple iteration method)
            //var fi1 = (T.Sin(y) - 1) / 2;
            //var fi2 = T.Cos(x + 0.5) - 2;
            // Approximated root of system of equations
            //var values = new[]
            //{
            //    new[] {-0.9},
            //    new[] {-1.1}
            //};

            // Parameters for gradient descend method
            var gdParams = new GradientDescendSolver.Params
            {
                Alpha = 1.0,
                Lambda = 0.5
            };

            // Vector representation of f(x, y)
            var f = new[]
            {
                new[] {f1},
                new[] {f2}
            };
            // Vector representation of gradient of f(x, y)
            var fGrad = new[]
            {
                new[] {f1.GradBy(x), f1.GradBy(y)},
                new[] {f2.GradBy(x), f2.GradBy(y)}
            };

            // Vector representation of fi(x, y)
            var fi = new[]
            {
                new[] {fi1},
                new[] {fi2}
            };
            // Vector representation of gradient of fi(x, y)
            var fiGrad = new[]
            {
                new[] {fi1.GradBy(x), fi1.GradBy(y)},
                new[] {fi2.GradBy(x), fi2.GradBy(y)}
            };

            // Function F to minimize: F(x, y) = f1(x, y) * f1(x, y) + f2(x, y) * f2(x, y)
            var minF = f1 * f1 + f2 * f2;
            // Gradient of function F
            var minFGrad = new[]
            {
                new[] {minF.GradBy(x)},
                new[] {minF.GradBy(y)}
            };

            // Prepare output file
            using var outputFile = File.Create(outputFileName);
            using var writer = new StreamWriter(outputFile);
            PrintHeader(writer, values, gdParams);

            // Calculations
            new SimpleIterationSolver(writer, epsilon).Solve(f, fi, fiGrad, values.CreateCopy());
            new NewtonSolver(writer, epsilon).Solve(f, fGrad, values.CreateCopy());
            new GradientDescendSolver(writer, epsilon, gdParams).Solve(minF, minFGrad, f, values.CreateCopy());
        }

        // Print header with values of x and parameters for gradient descend method
        private static void PrintHeader(TextWriter writer, double[][] values, GradientDescendSolver.Params parameters)
        {
            writer.WriteLine($"x:\n{values.ToPrettyString()}");
            writer.WriteLine($"\nParameters:\n\talpha = {parameters.Alpha}\n\tlambda = {parameters.Lambda}");
        }
    }

    // Solver based on simple iterations method
    public class SimpleIterationSolver
    {
        private readonly StreamWriter _streamWriter;
        private readonly double _eps;

        // Create solver that sends output to streamWriter with given epsilon
        public SimpleIterationSolver(StreamWriter streamWriter, double epsilon)
        {
            _streamWriter = streamWriter;
            _eps = epsilon;
        }

        // Solve system f with provided approximate values
        public void Solve(ITerm[][] f, ITerm[][] fi, ITerm[][] fiGrad, double[][] values)
        {
            WriteMethodName();

            var jakobiValues = fiGrad.Evaluate(values);
            var jakobiNorm = jakobiValues.Norm3();

            WriteInfo(jakobiValues, jakobiNorm);

            if (jakobiNorm >= 1d)
            {
                _streamWriter.WriteLine("Simple iteration method can't be used with provided input values.");
                return;
            }

            WriteTableHeader();
            IterateUntilConverges(f, fi, fiGrad, values);
        }

        // Make iteration until method converges with given epsilon
        private void IterateUntilConverges(ITerm[][] f, ITerm[][] fi, ITerm[][] fiGrad, double[][] values)
        {
            var i = 0;
            double[][] xDiff;
            do
            {
                i++;

                var previousValues = values;
                values = fi.Evaluate(values);

                WriteIterationInfo(i, f, fiGrad, values);

                xDiff = previousValues;
                xDiff.Subtract(values);
            } while (!Converged(xDiff));
        }

        // Print info about iteration i
        private void WriteIterationInfo(int i, ITerm[][] f, ITerm[][] fiGrad, double[][] values)
        {
            var fValues = f.Evaluate(values);
            var residualNorm = fValues.Norm3();
            var jakobiValues = fiGrad.Evaluate(values);
            var jakobiNorm = jakobiValues.Norm3();

            _streamWriter.WriteLine(
                $"{i,3}{values[0][0],15:G6}{values[1][0],15:G6}{residualNorm,15:G6}{fValues[0][0],15:G6}{fValues[1][0],15:G6}{jakobiNorm,15:G6}");
        }

        private void WriteTableHeader()
        {
            _streamWriter.WriteLine($"{"i",3}{"x",15}{"y",15}{"residual norm",15}{"f1",15}{"f2",15}{"jakobi norm",15}");
        }

        // Print Jakobi norm value
        private void WriteInfo(double[][] jakobiValues, double jakobiNorm)
        {
            _streamWriter.WriteLine("\nJakobi:");
            _streamWriter.WriteLine(jakobiValues.ToPrettyString());
            _streamWriter.WriteLine($"\nJakobi norm = {jakobiNorm}\n");
        }

        private void WriteMethodName()
        {
            _streamWriter.WriteLine("\n=== Simple iteration method ===");
        }

        // Check whether the method converged
        private bool Converged(double[][] xDiff) => xDiff.Norm3() < _eps;
    }

    // Solver based on Newton method
    public class NewtonSolver
    {
        private readonly StreamWriter _streamWriter;
        private readonly double _eps;

        // Create solver that sends output to streamWriter with given epsilon
        public NewtonSolver(StreamWriter streamWriter, double epsilon)
        {
            _streamWriter = streamWriter;
            _eps = epsilon;
        }

        // Solve system f with provided approximate values
        public void Solve(ITerm[][] f, ITerm[][] fGrad, double[][] values)
        {
            WriteMethodName();
            WriteTableHeader();
            IterateUntilConverges(f, fGrad, values);
        }

        // Make iteration until method converges with given epsilon
        private void IterateUntilConverges(ITerm[][] f, ITerm[][] fGrad, double[][] values)
        {
            var i = 0;
            var fValues = f.Evaluate(values);
            double[][] previousFValues;
            do
            {
                i++;

                RecalculateX(fGrad, values, fValues);

                previousFValues = fValues;
                fValues = f.Evaluate(values);

                WriteIterationInfo(i, values, fValues);
            } while (!Converged(previousFValues));
        }

        // Recalculate x = x - f(x)/f'(x)
        private static void RecalculateX(ITerm[][] fGrad, double[][] x, double[][] fValues)
        {
            var gradValues = fGrad.Evaluate(x);
            var inverseGradValues = gradValues.CreateInverse();
            x.Subtract(inverseGradValues.Multiply(fValues));
        }

        // Write info about iteration i
        private void WriteIterationInfo(int i, double[][] values, double[][] fValues)
        {
            var residual = fValues.Norm3();

            _streamWriter.WriteLine(
                $"{i,3}{values[0][0],15:G6}{values[1][0],15:G6}{residual,15:G6}{fValues[0][0],15:G6}{fValues[1][0],15:G6}");
        }

        private void WriteTableHeader()
        {
            _streamWriter.WriteLine($"{"i",3}{"x",15}{"y",15}{"residual norm",15}{"f1",15}{"f2",15}");
        }

        private void WriteMethodName()
        {
            _streamWriter.WriteLine("\n=== Newton method ===\n");
        }

        // Check whether the method converged
        private bool Converged(double[][] previousFValues) => previousFValues.Norm3() < _eps;
    }

    // Solver based on gradient descend method
    public class GradientDescendSolver
    {
        // Set of parameters for gradient descend method
        public class Params
        {
            // Alpha controls speed of gradient descend
            public double Alpha { get; set; } = 1.0;

            // Lambda is a multiplier by which Alpha gets multiplier in order to reduce gradient descend speed
            public double Lambda { get; set; } = 0.5;
        }

        private readonly StreamWriter _streamWriter;
        private readonly double _eps;
        private readonly Params _parameters;

        // Create solver that sends output to streamWriter with given epsilon and provided parameters
        public GradientDescendSolver(StreamWriter streamWriter, double epsilon, Params parameters)
        {
            _streamWriter = streamWriter;
            _eps = epsilon;
            _parameters = parameters;
        }

        // Minimize function minF with provided approximate values
        public void Solve(ITerm minF, ITerm[][] minFGrad, ITerm[][] f, double[][] values)
        {
            WriteMethodName();
            WriteInfo(minFGrad, values);
            WriteTableHeader();
            IterateUntilConverges(minF, minFGrad, f, values);
        }

        // Make iteration until method converges with given epsilon
        private void IterateUntilConverges(ITerm minF, ITerm[][] minFGrad, ITerm[][] f, double[][] values)
        {
            var k = 0;
            var i = 0;
            double[][] xDiff;
            do
            {
                i++;

                var newValues = RecalculateValues(minFGrad, values);
                while (minF.Evaluate(newValues) >= minF.Evaluate(values))
                {
                    k++;
                    _parameters.Alpha *= _parameters.Lambda;
                    newValues = RecalculateValues(minFGrad, values);
                }

                var previousX = values;
                values = newValues;

                WriteIterationInfo(i, minF, f, values, k);

                xDiff = previousX;
                xDiff.Subtract(values);
            } while (!Converged(xDiff));
        }

        // Write info about iteration i
        private void WriteIterationInfo(int i, ITerm minF, ITerm[][] f, double[][] values, int k)
        {
            var fValues = f.Evaluate(values);
            var residual = fValues.Norm3();
            var minFValue = minF.Evaluate(values);

            _streamWriter.WriteLine(
                $"{i,3}{values[0][0],15:G6}{values[1][0],15:G6}{_parameters.Alpha,15:G6}{residual,15:G6}{fValues[0][0],15:G6}{fValues[1][0],15:G6}{minFValue,15:G6}{k,3}");
        }

        private void WriteTableHeader()
        {
            _streamWriter.WriteLine(
                $"{"i",3}{"x",15}{"y",15}{"alpha",15}{"residual",15}{"f1",15}{"f2",15}{"F",15}{"k",3}");
        }

        // Write value of gradient vector
        private void WriteInfo(ITerm[][] minFGrad, double[][] values)
        {
            var gradient = minFGrad.Evaluate(values);
            _streamWriter.WriteLine($"\nGradient vector:\n{gradient.ToPrettyString()}\n");
        }

        private void WriteMethodName()
        {
            _streamWriter.WriteLine("\n=== Gradient descend method ===");
        }

        // Recalculate value of matrix x
        private double[][] RecalculateValues(ITerm[][] minFGrad, double[][] values)
        {
            var gradValues = minFGrad.Evaluate(values);
            gradValues.Multiply(_parameters.Alpha);

            var newValues = values.CreateCopy();
            newValues.Subtract(gradValues);
            return newValues;
        }

        // Check whether the method converged
        private bool Converged(double[][] xDiff) => xDiff.Norm3() < _eps;
    }

    // Class that provides methods for function evaluation
    public static class Evaluation
    {
        // Evaluate function in every cell with provided input values and return matrix of the same dimensions with results
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

    // Class that provides operations on matrices
    public static class MatrixMath
    {
        // Create inverse matrix of matrix a
        public static double[][] CreateInverse(this double[][] a)
        {
            var l = a.CreateCompatible();
            var u = a.CreateCopy();
            var p = a.CreateSwapMatrix();

            MakeLupDecomposition(u, l, p);
            return FindInverseMatrix(p, a, l, u);
        }

        // Find inverse matrix using LU = PA
        private static double[][] FindInverseMatrix(int[] p, double[][] a, double[][] l, double[][] u)
        {
            var pa = a.CreateCopy();
            pa.ApplySwaps(p);

            var pe = pa.CreateOneMatrix();
            pe.ApplySwaps(p);

            return SolveEquation(l, u, pe);
        }

        // Solve equation AX = LUX = E and find X = A^(-1)
        private static double[][] SolveEquation(double[][] l, double[][] u, double[][] pe)
        {
            var y = l.FindY(pe);
            var invertedA = u.FindX(y);
            return invertedA;
        }

        // Make LUP decomposition for given matrix a (LU = PA)
        private static void MakeLupDecomposition(double[][] a, double[][] l, int[] p)
        {
            var u = a;
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
        }

        // Multiply matrix by value
        public static void Multiply(this double[][] matrix, double value) => matrix.Map(x => x * value);

        // Apply map function to each element of the matrix, returning new matrix as a result
        private static void Map(this double[][] matrix, Func<double, double> f)
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
