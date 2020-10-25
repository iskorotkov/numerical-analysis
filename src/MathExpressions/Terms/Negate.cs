namespace MathExpressions.Terms
{
    public readonly struct Negate : ITerm
    {
        private readonly ITerm _arg;

        public Negate(ITerm arg) => _arg = arg;

        public double Evaluate(double[][] x) => -_arg.Evaluate(x);
        public ITerm Grad() => -_arg.Grad();
        public ITerm GradBy(ITerm v) => -_arg.GradBy(v);
        public override string ToString() => $"-{_arg}";
    }
}
