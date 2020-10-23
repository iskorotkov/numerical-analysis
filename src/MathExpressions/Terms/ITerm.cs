namespace MathExpressions.Terms
{
    public interface ITerm
    {
        double Evaluate(double[][] x);
        ITerm Grad();
        ITerm GradBy(ITerm variable);

        public static ITerm operator +(ITerm left, ITerm right) => new Sum(left, right);
        public static ITerm operator -(ITerm left, ITerm right) => new Diff(left, right);
        public static ITerm operator *(ITerm left, ITerm right) => new Product(left, right);
        public static ITerm operator /(ITerm left, ITerm right) => new Fraction(left, right);

        public static ITerm operator +(double d, ITerm right) => new Sum(new Value(d), right);
        public static ITerm operator -(double d, ITerm right) => new Diff(new Value(d), right);
        public static ITerm operator *(double d, ITerm right) => new Product(new Value(d), right);
        public static ITerm operator /(double d, ITerm right) => new Fraction(new Value(d), right);

        public static ITerm operator +(ITerm left, double d) => new Sum(left, new Value(d));
        public static ITerm operator -(ITerm left, double d) => new Diff(left, new Value(d));
        public static ITerm operator *(ITerm left, double d) => new Product(left, new Value(d));
        public static ITerm operator /(ITerm left, double d) => new Fraction(left, new Value(d));

        public static ITerm operator +(ITerm term) => term;
        public static ITerm operator -(ITerm term) => new Diff(new Value(0), term);
    }
}