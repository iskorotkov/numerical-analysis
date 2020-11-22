using System;

namespace MathExpressions.Terms
{
    public readonly struct Logarithm : ITerm
    {
        private readonly ITerm _base;
        private readonly ITerm _arg;

        public Logarithm(ITerm arg, ITerm @base) => (_base, _arg) = (@base, arg);

        public double Evaluate(double[][] x) => Math.Log(_arg.Evaluate(x), _base.Evaluate(x));

        public ITerm Grad() => _arg.Grad() / (_arg * T.Log(_base, T.Value(Math.E)));

        public ITerm GradBy(ITerm v) => _arg.GradBy(v) / (_arg * T.Log(_base, T.Value(Math.E)));

        public override string ToString() => $"log_{_base}({_arg})";
    }
}
