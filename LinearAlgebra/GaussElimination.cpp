#include "GaussElimination.h"

void SwapRows(Matrix& mx, int x, int y)
{
	if (x < 0 || x >= mx.RowNum() || y < 0 || y >= mx.RowNum())
		return;
	for (int i = 0; i < mx.ColNum(); i++)
		std::swap(mx.a[x][i], mx.a[y][i]);
}

void MultiplyRow(Matrix& mx, int r, double c)
{
	if (r < 0 || r >= mx.RowNum())
		return;
	for (int i = 0; i < mx.ColNum(); i++)
		mx.a[r][i] *= c;
}

void AddRowToOther(Matrix& mx, int x, int y, double c)
{
	if (x < 0 || x >= mx.RowNum() || y < 0 || y >= mx.RowNum())
		return;
	for (int i = 0; i < mx.ColNum(); i++)
		mx.a[y][i] += c * mx.a[x][i];
}

double GaussElimination(Matrix& mx)
{
	int x = 0;
	int y = 0;
	double res = 1.0;
	while (x < mx.RowNum() && y < mx.ColNum())
	{
		int x0 = x;
		while (x0 < mx.RowNum() && IsZero(mx.a[x0][y]))
			x0++;
		if (x0 == mx.RowNum())
		{
			y++;
			continue;
		}
		else
		{
			if (x != x0)
			{
				SwapRows(mx, x, x0);
				res *= -1.0;
			}
			for (int j = x + 1; j < mx.RowNum(); j++)
			{
				double c = -mx.a[j][y] / mx.a[x][y];
				AddRowToOther(mx, x, j, c);
			}
			x++;
			y++;
		}
	}
	return res;
}

double Determinant(Matrix& mx)
{
	if (mx.ColNum()!=mx.RowNum())
		return 0.0;
	Matrix mx_copy = mx;
	double res = GaussElimination(mx_copy);
	for (int i = 0; i < mx.RowNum(); i++)
		res *= mx_copy.a[i][i];
	return res;
}
