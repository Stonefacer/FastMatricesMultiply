
#include "stdafx.h"
#include "StrassenAlgorithm.h"

void print(int ** matrix, size_t size) {
	if (size > 16) {
		return;
	}
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

int Check(int ** matrix1, int ** matrix2, size_t size) {
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			if (matrix1[i][j] != matrix2[i][j]) {
				return -1;
			}
		}
	}
	return 0;
}

/*
15	18	21
42	54	66
69	90	111
*/

int main(int argc, char ** args)
{
	size_t size = 512;
	if (argc > 1) {
		size = atoi(args[1]);
	}
	std::cout << "Size: " << size << std::endl;
	auto matrix1 = new int*[size];
	for (size_t i = 0; i < size; i++) {
		matrix1[i] = new int[size];
		for (size_t j = 0; j < size; j++) {
			matrix1[i][j] = (i*size + j) % 10 + 1;
		}
	}
	auto matrix2 = new int*[size];
	for (size_t i = 0; i < size; i++) {
		matrix2[i] = new int[size];
		for (size_t j = 0; j < size; j++) {
			matrix2[i][j] = (i*size + j) % 10 + 1;
		}
	}
	long long timeSpendFast = time(nullptr);
	auto resultMulFast = StrassenAlgorithm::FastMultiply(matrix1, matrix2, size);
	timeSpendFast = time(nullptr) - timeSpendFast;
	StrassenAlgorithm::Clear();
	print(matrix1, size);
	print(matrix2, size);
	std::cout << "MUL FAST Time:" << timeSpendFast << std::endl;
	std::cout << "Cache missed: " << StrassenAlgorithm::CacheMissed << " Cache found: " << StrassenAlgorithm::CacheFound << std::endl;
	print(resultMulFast, size);

	long long timeSpend = time(nullptr);
	auto result = StrassenAlgorithm::Multiply(matrix1, matrix2, size);
	timeSpend = time(nullptr) - timeSpend;
	print(matrix1, size);
	print(matrix2, size);
	std::cout << "MUL Time:" << timeSpend << std::endl;
	print(result, size);

	std::cout << (float)timeSpendFast / timeSpend << std::endl;
	std::cout << "Checking..." << std::endl;
	if(Check(result, resultMulFast, size)){
		std::cout << "FAILED!" << std::endl;
	}
	else {
		std::cout << "SUCCESS!" << std::endl;
	}
#if defined __DEBUG || defined _DEBUG || defined DEBUG
	system("pause");
#endif
	return 0;
}

