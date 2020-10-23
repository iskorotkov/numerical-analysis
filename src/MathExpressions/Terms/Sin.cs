using System;

namespace MathExpressions.Terms
{
    public readonly struct Sin : ITerm
    {
        private readonly ITerm _arg;

        public Sin(ITerm arg) => _arg = arg;
        public double Evaluate(double[][] x) => Math.Sin(_arg.Evaluate(x));
        public ITerm Grad() => new Product(new Cos(_arg), _arg.Grad());
        public ITerm GradBy(ITerm variable) => new Product(new Cos(_arg), _arg.GradBy(variable));

        public override string ToString() => $"sin({_arg})";
    }
}
