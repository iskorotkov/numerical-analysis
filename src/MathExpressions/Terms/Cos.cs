﻿using System;

namespace MathExpressions.Terms
{
    public readonly struct Cos : ITerm
    {
        private readonly ITerm _arg;

        public Cos(ITerm arg) => _arg = arg;
        public double Evaluate(double[][] x) => Math.Cos(_arg.Evaluate(x));
        public ITerm Grad() => new Product(-T.Sin(_arg), _arg.Grad());
        public ITerm GradBy(ITerm variable) => new Product(-T.Sin(_arg), _arg.GradBy(variable));

        public override string ToString() => $"cos({_arg})";
    }
}