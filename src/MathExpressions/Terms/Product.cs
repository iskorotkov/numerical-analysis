namespace MathExpressions.Terms
{
    public readonly struct Product : ITerm
    {
        private readonly ITerm _leftArg;
        private readonly ITerm _rightArg;

        public Product(ITerm leftArg, ITerm rightArg) => (_leftArg, _rightArg) = (leftArg, rightArg);

        public double Evaluate(double[][] x) => _leftArg.Evaluate(x) * _rightArg.Evaluate(x);

        public ITerm Grad() =>
            new Sum(
                new Product(_leftArg.Grad(), _rightArg),
                new Product(_leftArg, _rightArg.Grad()));

        public ITerm GradBy(ITerm v) =>
            new Sum(
                new Product(_leftArg.GradBy(v), _rightArg),
                new Product(_leftArg, _rightArg.GradBy(v)));

        public override string ToString() => $"({_leftArg}) * ({_rightArg})";
    }
}