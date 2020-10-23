namespace MathExpressions.Terms
{
    public readonly struct Diff : ITerm
    {
        private readonly ITerm _leftArg;
        private readonly ITerm _rightArg;

        public Diff(ITerm leftArg, ITerm rightArg) => (_leftArg, _rightArg) = (leftArg, rightArg);

        public double Evaluate(double[][] x) => _leftArg.Evaluate(x) - _rightArg.Evaluate(x);
        public ITerm Grad() => new Diff(_leftArg.Grad(), _rightArg.Grad());
        public ITerm GradBy(Var v) => new Diff(_leftArg.GradBy(v), _rightArg.GradBy(v));

        public override string ToString() => $"{_leftArg} - {_rightArg}";
    }
}