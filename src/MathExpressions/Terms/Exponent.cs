using System;

namespace MathExpressions.Terms
{
    public readonly struct Exponent: ITerm
    {
        private readonly ITerm _base;
        private readonly ITerm _power;

        public Exponent(ITerm @base, ITerm power) => (_base, _power) = (@base, power);

        public double Evaluate(double[][] x) => Math.Pow(_base.Evaluate(x), _power.Evaluate(x));

        public ITerm Grad() => this * T.Log(_base, T.Value(Math.E)) * _power.Grad();

        public ITerm GradBy(ITerm v) => this * T.Log(_base, T.Value(Math.E)) * _power.GradBy(v);

        public override string ToString() => $"{_base}^{_power}";
    }
}
