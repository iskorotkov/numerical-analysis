using System.Globalization;

namespace MathExpressions.Terms
{
    public readonly struct Value : ITerm
    {
        private readonly double _value;

        public Value(double value) => _value = value;

        public double Evaluate(double[][] x) => _value;
        public ITerm Grad() => new Value(0d);
        public ITerm GradBy(ITerm v) => new Value(0d);

        public override string ToString() => _value.ToString(CultureInfo.InvariantCulture);
    }
}