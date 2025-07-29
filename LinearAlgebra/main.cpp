#include <iostream>
#include "matrix.h"
#include "GaussElimination.h"
#include "polynomial.h"
#include "charpoly.h"
#include <random>

int main()
{
/*	std::cout << "Number of equations: ";
	int m = 0;
	std::cin >> m;
	std::vector<std::string> equs;
	for (int i = 0; i < m; i++)
	{
		std::string line;
		std::cout << "Equ. " << i + 1 << ": " ;
		std::cin >> line;
		equs.push_back(line);
	}
	Matrix mx= CreateAugmentedMatrix(equs);
	std::cout << "\nSolutions:\n";
	DisplayParametricSolution(mx);
	std::cout << "\nAugmented matrix:\n";
	mx.Display();
	std::cout << "\nGauss-Jordan elimination:\n";
	GaussJordanElimination(mx,true);
*/
	srand(time(NULL));
	int n = 10;
	Matrix mx(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			mx.a[i][j] = rand() % 2;
		}
	}
	mx.Display();
	Polynomial charpoly = CharacteristicPolynomial(mx);
	std::cout << "Characteristic polynomial: ";
	charpoly.Display();
	return 0;
}