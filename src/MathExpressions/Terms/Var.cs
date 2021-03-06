﻿namespace MathExpressions.Terms
{
    public readonly struct Var : ITerm
    {
        private readonly int _index;

        public Var(int index) => _index = index;

        public double Evaluate(double[][] x) => x[_index][0];
        public ITerm Grad() => new Value(1d);
        public ITerm GradBy(ITerm v) => new Value(Equals(v) ? 1d : 0d);

        public override string ToString() =>
            _index switch
            {
                0 => "x",
                1 => "y",
                2 => "z",
                3 => "w",
                _ => $"var{_index}"
            };
    }
}
