using System;
using System.IO;
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
            var parameters = new GradientDescendParameters
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
            PrintHeader(writer, values, parameters);

            // Calculations
            new SimpleIterationSolver(writer, epsilon).Solve(f, fi, fiGrad, values.CreateCopy());
            new NewtonSolver(writer, epsilon).Solve(f, fGrad, values.CreateCopy());
            new GradientDescendSolver(writer, epsilon, parameters).Solve(minF, minFGrad, f, values.CreateCopy());
        }

        // Print header with values of x and parameters for gradient descend method
        private static void PrintHeader(TextWriter writer, double[][] values,
            GradientDescendParameters parameters)
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
            _streamWriter.WriteLine("\n=== Simple iteration method ===");

            var jakobiValues = fiGrad.Evaluate(values);
            _streamWriter.WriteLine("\nJakobi:");
            _streamWriter.WriteLine(jakobiValues.ToPrettyString());

            var jakobiNorm = jakobiValues.Norm3();
            _streamWriter.WriteLine($"\nJakobi norm = {jakobiNorm}\n");

            if (jakobiNorm >= 1d)
            {
                _streamWriter.WriteLine("Simple iteration method can't be used with provided input values.");
                return;
            }

            _streamWriter.WriteLine($"{"i",3}{"x",15}{"y",15}{"residual norm",15}{"f1",15}{"f2",15}{"jakobi norm",15}");

            var i = 0;
            double[][] xDiff;
            do
            {
                i++;

                var previousValues = values;
                values = fi.Evaluate(values);

                var fValues = f.Evaluate(values);
                var residualNorm = fValues.Norm3();
                jakobiValues = fiGrad.Evaluate(values);
                jakobiNorm = jakobiValues.Norm3();

                _streamWriter.WriteLine(
                    $"{i,3}{values[0][0],15:G6}{values[1][0],15:G6}{residualNorm,15:G6}{fValues[0][0],15:G6}{fValues[1][0],15:G6}{jakobiNorm,15:G6}");

                xDiff = previousValues;
                xDiff.Subtract(values);
            } while (!Converged(xDiff));
        }

        // Check whether the method converged
        private bool Converged(double[][] xDiff) =>
            xDiff.All(row => row
                .Select(Math.Abs)
                .All(value => value < _eps));
    }

    // Set of parameters for gradient descend method
    public class GradientDescendParameters
    {
        // Alpha controls speed of gradient descend
        public double Alpha { get; set; } = 1.0;

        // Lambda is a multiplier by which Alpha gets multiplier in order to reduce gradient descend speed
        public double Lambda { get; set; } = 0.5;
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
            _streamWriter.WriteLine("\n=== Newton method ===\n");

            _streamWriter.WriteLine($"{"i",3}{"x",15}{"y",15}{"residual norm",15}{"f1",15}{"f2",15}");

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

                _streamWriter.WriteLine(
                    $"{i,3}{values[0][0],15:G6}{values[1][0],15:G6}{residual,15:G6}{fValues[0][0],15:G6}{fValues[1][0],15:G6}");
            } while (!Converged(previousFValues));
        }

        // Check whether the method converged
        private bool Converged(double[][] previousFValues) =>
            previousFValues.All(row => row
                .Select(Math.Abs)
                .All(value => value < _eps));
    }

    // Solver based on gradient descend method
    public class GradientDescendSolver
    {
        private readonly StreamWriter _streamWriter;
        private readonly double _eps;
        private readonly GradientDescendParameters _parameters;

        // Create solver that sends output to streamWriter with given epsilon and provided parameters
        public GradientDescendSolver(StreamWriter streamWriter, double epsilon, GradientDescendParameters parameters)
        {
            _streamWriter = streamWriter;
            _eps = epsilon;
            _parameters = parameters;
        }

        // Minimize function minF with provided approximate values
        public void Solve(ITerm minF, ITerm[][] minFGrad, ITerm[][] f, double[][] values)
        {
            _streamWriter.WriteLine("\n=== Gradient descend method ===");

            var gradient = minFGrad.Evaluate(values);
            _streamWriter.WriteLine($"\nGradient vector:\n{gradient.ToPrettyString()}\n");

            _streamWriter.WriteLine(
                $"{"i",3}{"x",15}{"y",15}{"alpha",15}{"residual",15}{"f1",15}{"f2",15}{"F",15}{"k",3}");

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

                var fValues = f.Evaluate(values);
                var residual = fValues.Norm3();
                var minFValue = minF.Evaluate(values);

                _streamWriter.WriteLine(
                    $"{i,3}{values[0][0],15:G6}{values[1][0],15:G6}{_parameters.Alpha,15:G6}{residual,15:G6}{fValues[0][0],15:G6}{fValues[1][0],15:G6}{minFValue,15:G6}{k,3}");

                xDiff = previousX;
                xDiff.Subtract(values);
            } while (!Converged(xDiff));
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
        private bool Converged(double[][] xDiff) => Math.Abs(xDiff.Norm3()) < _eps;
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
