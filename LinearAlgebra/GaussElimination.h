#pragma once
#include "matrix.h"

void SwapRows(Matrix& mx, int x, int y);				//swap rows x and y
void MultiplyRow(Matrix& mx, int r, double c);			//multiply every element of row r by c
void AddRowToOther(Matrix& mx, int x, int y, double c); //add c times row x to row y
double GaussElimination(Matrix& mx,bool bDisplaySteps);					//does Gauss elimination of mx and returns (-1)^s, where s is the number of swaps
double Determinant(Matrix& mx);							//returns determinant
void GaussJordanElimination(Matrix& mx, bool bDisplaySteps);


void DisplayParametricSolution(Matrix& mx);