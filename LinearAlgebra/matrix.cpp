#include "matrix.h"
#include <algorithm>

std::string double_to_string(double x, int precision)
{
	std::string s=std::to_string(x);
	int ind = 0;
	while (ind < s.size() && s[ind] != '.')
		ind++;
	if (ind < s.size())
	{
		s = s.substr(0, std::min((int)s.size(), ind + 1 + precision));
		while (s.back() == '0')
			s.pop_back();
	}
	return s;
}

bool IsZero(double x)
{
	return (x<zero_error && x>-zero_error);
}

void Matrix::Display()
{
	for (int i = 0; i < RowNum(); i++)
	{
		std::string row;
		for (int j = 0; j < ColNum(); j++)
			row += double_to_string(a[i][j], double_precision)+'\t';
		row += '\n';
		std::cout << row;
	}
}
