using System;
using MathExpressions.Terms;
using Shouldly;
using Xunit;

namespace MathExpressions.Tests.Terms
{
    public class Exponents
    {
        private readonly ITerm _term = T.Exp(T.Sin(T.Var(0)), T.Cos(T.Var(0)));

        [Fact]
        public void EvaluatesCorrectly()
        {
            double[][] input = {new[] {1d}};
            var actual = _term.Evaluate(input);
            var expected = Math.Pow(Math.Sin(1), Math.Cos(1));
            actual.ShouldBe(expected);
        }

        [Fact]
        public void CalculatesGradCorrectly()
        {
            var actual = _term.Grad();
            var expected = _term * T.Log(T.Sin(T.Var(0)), T.Value(Math.E)) * (-T.Sin(T.Var(0)) * 1);
            actual.ShouldBe(expected);
        }

        [Fact]
        public void CalculatesGradByVarCorrectly()
        {
            var actual = _term.GradBy(T.Var(0));
            var expected = _term * T.Log(T.Sin(T.Var(0)), T.Value(Math.E)) * (-T.Sin(T.Var(0)) * 1);
            actual.ShouldBe(expected);
        }

        [Fact]
        public void CalculatesGradByVarAsZero()
        {
            var actual = _term.GradBy(T.Var(1));
            // Still evaluates to 0
            var expected = _term * T.Log(T.Sin(T.Var(0)), T.Value(Math.E)) * (-T.Sin(T.Var(0)) * 0);
            actual.ShouldBe(expected);
        }
    }
}
