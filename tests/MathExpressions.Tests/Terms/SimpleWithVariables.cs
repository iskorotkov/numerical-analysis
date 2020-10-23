using System;
using MathExpressions.Terms;
using Shouldly;
using Xunit;

namespace MathExpressions.Tests.Terms
{
    public class SimpleWithVariables
    {
        private readonly ITerm _term;
        private const double Tolerance = 1e-10;

        public SimpleWithVariables()
        {
            _term = T.Sin(new Var(0)) / T.Cos(new Var(1));
        }

        [Fact]
        public void EvaluatesCorrectly()
        {
            var expected = Math.Sin(2) / Math.Cos(3);
            var values = new[]
            {
                new[] {2d},
                new[] {3d}
            };

            Assert.Throws<NullReferenceException>(() => _term.Evaluate(null));
            Assert.Throws<NullReferenceException>(() => _term.Evaluate(new double[10][]));
            Assert.Throws<IndexOutOfRangeException>(() => _term.Evaluate(new[] {new[] {2d, 3d}}));
            _term.Evaluate(values).ShouldBe(expected, Tolerance);
        }

        [Fact]
        public void GradIsCorrect()
        {
            var expected = (T.Cos(new Var(0)) * T.Value(1) * T.Cos(new Var(1))
                            - (T.Sin(new Var(0)) * (-T.Sin(new Var(1)) * new Value(1))))
                           / (T.Cos(new Var(1)) * T.Cos(new Var(1)));

            var actual = _term.Grad();
            actual.ShouldBe(expected);
        }

        [Fact]
        public void GradEvaluatesCorrectly()
        {
            var expected = (Math.Cos(2) * Math.Cos(3) + Math.Sin(2) * Math.Sin(3)) / (Math.Cos(3) * Math.Cos(3));
            var values = new[]
            {
                new[] {2d},
                new[] {3d}
            };

            Assert.Throws<NullReferenceException>(() => _term.Grad().Evaluate(null));
            Assert.Throws<NullReferenceException>(() => _term.Grad().Evaluate(new double[10][]));
            Assert.Throws<IndexOutOfRangeException>(() => _term.Grad().Evaluate(new[] {new[] {2d, 3d}}));
            _term.Grad().Evaluate(values).ShouldBe(expected, Tolerance);
        }

        [Fact]
        public void GradByVarIsCorrect()
        {
            var xGrad = (T.Cos(new Var(0)) * T.Value(1) * T.Cos(new Var(1))
                         - (T.Sin(new Var(0)) * (-T.Sin(new Var(1)) * new Value(0))))
                        / (T.Cos(new Var(1)) * T.Cos(new Var(1)));

            var yGrad = (T.Cos(new Var(0)) * T.Value(0) * T.Cos(new Var(1))
                         - (T.Sin(new Var(0)) * (-T.Sin(new Var(1)) * new Value(1))))
                        / (T.Cos(new Var(1)) * T.Cos(new Var(1)));

            var byX = _term.GradBy(new Var(0));
            var byY = _term.GradBy(new Var(1));
            byX.ShouldBe(xGrad);
            byY.ShouldBe(yGrad);
        }

        [Fact]
        public void GradByVarEvaluatesCorrectly()
        {
            var xGrad = Math.Cos(2) / Math.Cos(3);
            var yGrad = Math.Sin(2) * Math.Sin(3) / (Math.Cos(3) * Math.Cos(3));
            var values = new[]
            {
                new[] {2d},
                new[] {3d}
            };

            Assert.Throws<NullReferenceException>(() => _term.Evaluate(null));
            Assert.Throws<NullReferenceException>(() => _term.Evaluate(new double[10][]));
            Assert.Throws<IndexOutOfRangeException>(() => _term.Evaluate(new[] {new[] {2d, 3d}}));
            _term.GradBy(new Var(0)).Evaluate(values).ShouldBe(xGrad, Tolerance);
            _term.GradBy(new Var(1)).Evaluate(values).ShouldBe(yGrad, Tolerance);
        }
    }
}
