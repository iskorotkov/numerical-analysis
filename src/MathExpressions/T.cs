using MathExpressions.Terms;

namespace MathExpressions
{
    public static class T
    {
        public static ITerm Cos(ITerm t) => new Cos(t);
        public static ITerm Cos(double d) => new Cos(new Value(d));

        public static ITerm Sin(ITerm t) => new Sin(t);
        public static ITerm Sin(double d) => new Sin(new Value(d));

        public static ITerm Value(double d) => new Value(d);
        public static ITerm Var(int index) => new Var(index);
    }
}
