namespace MathExpressions.Terms
{
    public readonly struct Fraction : ITerm
    {
        private readonly ITerm _leftArg;
        private readonly ITerm _rightArg;

        public Fraction(ITerm leftArg, ITerm rightArg) => (_leftArg, _rightArg) = (leftArg, rightArg);

        public double Evaluate(double[][] x) => _leftArg.Evaluate(x) / _rightArg.Evaluate(x);

        public ITerm Grad() =>
            new Fraction(
                new Diff(
                    new Product(_leftArg.Grad(), _rightArg),
                    new Product(_leftArg, _rightArg.Grad())),
                new Product(_rightArg, _rightArg));

        public ITerm GradBy(ITerm v) =>
            new Fraction(
                new Diff(
                    new Product(_leftArg.GradBy(v), _rightArg),
                    new Product(_leftArg, _rightArg.GradBy(v))),
                new Product(_rightArg, _rightArg));

        public override string ToString() => $"({_leftArg}) / ({_rightArg})";
    }
}