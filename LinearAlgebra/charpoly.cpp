#include "charpoly.h"
#include "GaussElimination.h"

Polynomial CharacteristicPolynomial(Matrix& mx)
{
    if (mx.RowNum()!=mx.ColNum())
        throw std::invalid_argument("Characteristic polynomial is not defined for non-square matrices.");
    int n = mx.RowNum();
    std::vector<std::pair<double, double>> values;
    values.resize(n);
    Matrix mx_copy = mx * -1.0;
    for (int i = 0; i < n; i++)
    {
        values[i].first = i;
        values[i].second = Determinant(mx_copy);
        for (int j = 0; j < n; j++)
            mx_copy.a[j][j] += 1.0;
    }
    return PolyInterpolationLeadingOne(values);
}
