#pragma once
#include <vector>
#include <algorithm>
#include <utility>
#include "matrix.h"


class Polynomial
{
public:
	std::vector<double> a;
	int degree() const { return a.size() - 1; };
	Polynomial() { a = { 0 }; };
	double EvaluateAt(double x);
	void Display();
	Polynomial operator+(const Polynomial& p) const
	{
		Polynomial Res;
		Res.a.resize(std::max(a.size(), p.a.size()));
		for (int i = 0; i < Res.a.size(); i++)
		{
			Res.a[i] = 0.0;
			if (i < a.size())
				Res.a[i] += a[i];
			if (i < p.a.size())
				Res.a[i] += p.a[i];
		}
		while (IsZero(Res.a.back()))
			Res.a.pop_back();
		return Res;
	}
	Polynomial operator-(const Polynomial& p) const
	{
		Polynomial Res;
		Res.a.resize(std::max(a.size(), p.a.size()));
		for (int i = 0; i < Res.a.size(); i++)
		{
			Res.a[i] = 0.0;
			if (i < a.size())
				Res.a[i] += a[i];
			if (i < p.a.size())
				Res.a[i] -= p.a[i];
		}
		while (IsZero(Res.a.back()))
			Res.a.pop_back();
		return Res;
	}
	Polynomial operator*(const Polynomial& p) const
	{
		Polynomial Res;
		Res.a.resize(a.size()+p.a.size()-1);
		for (int i = 0; i < Res.a.size(); i++)
			Res.a[i] = 0.0;
		for (int i = 0; i < a.size(); i++)
		{
			for (int j = 0; j < p.a.size(); j++)
				Res.a[i + j] += a[i] * p.a[j];
		}
		return Res;
	}
};

Polynomial PolyInterpolation(std::vector<std::pair<double,double>> values); 
Polynomial PolyInterpolationLeadingOne(std::vector<std::pair<double, double>> values);