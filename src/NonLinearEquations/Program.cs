using System;

using static System.Math;

namespace NonLinearEquations
{
    internal static class Program
    {
        private static void Main()
        {
            new SimpleIterationSolver();
            new NewtonSolver();
        }
    }

    public class SimpleIterationSolver
    {
        public SimpleIterationSolver()
        {
            var x = -10d;

            Func<double, double> f = x => -x + 1;

            var m1 = 1.0; // f.Derivative.Abs.Min;
            var M1 = 1.0; // f.Derivative.Abs.Max;
            var tau = 2 / (m1 + M1);

            for (var i = 0; i < 20; i++)
            {
                x = x + tau * Abs(f(x));
                Console.WriteLine($"i={i}, tau={tau}, x={x}");
            }
        }
    }

    public class NewtonSolver
    {
        public NewtonSolver()
        {
            var x = -10d;

            Func<double, double> f = x => -x + 1;
            Func<double, double> fDerivative = x => -1;

            for (var i = 0; i < 20; i++)
            {
                x = x - f(x) / fDerivative(x);
                Console.WriteLine($"i={i}, x={x}");
            }
        }
    }

    // TODO: Implement gradient descend method
    public class GradientDescend
    {
        public GradientDescend()
        {
        }
    }
}
