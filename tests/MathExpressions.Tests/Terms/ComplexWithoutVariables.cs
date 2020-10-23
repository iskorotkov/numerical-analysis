using System;
using MathExpressions.Terms;
using Shouldly;
using Xunit;

namespace MathExpressions.Tests.Terms
{
    public class ComplexWithoutVariables
    {
        private readonly ITerm _term;
        private const double Tolerance = 1e-10;

        public ComplexWithoutVariables()
        {
            _term = T.Sin(T.Cos(15)) + T.Sin(1) * T.Cos(2);
        }

        [Fact]
        public void EvaluatesCorrectly()
        {
            var expected = Math.Sin(Math.Cos(15)) + Math.Sin(1) * Math.Cos(2);

            _term.Evaluate(null).ShouldBe(expected, Tolerance);
            _term.Evaluate(new double[10][]).ShouldBe(expected, Tolerance);
        }

        [Fact]
        public void GradIsCorrect()
        {
            var expected = T.Cos(T.Cos(15)) * (-T.Sin(15) * T.Value(0)) +
                           (T.Cos(1) * T.Value(0) * T.Cos(2)
                            + T.Sin(1) * (-T.Sin(2) * T.Value(0)));

            var actual = _term.Grad();
            actual.ShouldBe(expected);
        }

        [Fact]
        public void GradEvaluatesToZero()
        {
            _term.Grad().Evaluate(null).ShouldBe(0d);
            _term.Grad().Evaluate(new double[10][]).ShouldBe(0d);
        }

        [Fact]
        public void GradByVarIsCorrect()
        {
            var expected = T.Cos(T.Cos(15)) * (-T.Sin(15) * T.Value(0)) +
                           (T.Cos(1) * T.Value(0) * T.Cos(2)
                            + T.Sin(1) * (-T.Sin(2) * T.Value(0)));
            var x = new Var(0);

            var actual = _term.GradBy(x);
            actual.ShouldBe(expected);
        }

        [Fact]
        public void GradByVarEvaluatesToZero()
        {
            var x = new Var(0);
            _term.GradBy(x).Evaluate(null).ShouldBe(0d);
            _term.GradBy(x).Evaluate(new double[10][]).ShouldBe(0d);
        }
    }
}
