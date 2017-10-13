#include "stdafx.h"
#include "StrassenAlgorithm.h"

#define MAX_DEPTH 3
#define MAX_DEPTH_THREAD 10
#define THREADS_COUNT 8

//private

std::map<size_t, std::stack<int**>*> StrassenAlgorithm::_allocatedMatrices;
int StrassenAlgorithm::CacheMissed = 0;
int StrassenAlgorithm::CacheFound = 0;
int StrassenAlgorithm::RecursionDepth = 0;
int StrassenAlgorithm::ThreadsCount = 8;
std::mutex StrassenAlgorithm::_syncMutex;

int ** StrassenAlgorithm::AllocateMatrix(size_t size)
{
	int ** result = nullptr;
	_syncMutex.lock();
	if (_allocatedMatrices.count(size)) {
		auto stack = _allocatedMatrices[size];
		if (stack->size() != 0) {
			StrassenAlgorithm::CacheFound++;
			result = stack->top();
			stack->pop();
		}
	}
	if (result == nullptr) {
		StrassenAlgorithm::CacheMissed++;
		result = new int*[size];
		for (size_t i = 0; i < size; i++) {
			result[i] = new int[size];
			for (size_t j = 0; j < size; j++) {
				result[i][j] = 0;
			}
		}
	}
	_syncMutex.unlock();
	return result;
}

void StrassenAlgorithm::FreeMatrix(int ** matrix, size_t size)
{
	_syncMutex.lock();
	std::stack<int**> * stack;
	if (_allocatedMatrices.count(size)) {
		stack = _allocatedMatrices[size];
	}
	else {
		stack = new std::stack<int**>();
		_allocatedMatrices[size] = stack;
	}
	stack->push(matrix);
	_syncMutex.unlock();
}

int ** StrassenAlgorithm::Multiply(int ** matrix1, int ** matrix2, size_t size)
{
	auto result = AllocateMatrix(size);
	auto threads = new std::thread[ThreadsCount];
	for (size_t offset = 0; offset < ThreadsCount; offset++) {
		threads[offset] = std::thread([=, &matrix1, &matrix2, &result]() {
			auto localOffset = offset;
			for (size_t i = localOffset; i < size; i += ThreadsCount) {
				for (size_t j = 0; j < size; j++) {
					result[i][j] = 0;
					for (size_t k = 0; k < size; k++) {
						result[i][j] += matrix1[i][k] * matrix2[k][j];
					}
				}
			}
		});
	}
	for (size_t offset = 0; offset < ThreadsCount; offset++) {
		threads[offset].join();
	}
	return result;
}

int ** StrassenAlgorithm::MultiplySync(int ** matrix1, int ** matrix2,
	size_t minimalIndexXMatrix1, size_t minimalIndexYMatrix1,
	size_t minimalIndexXMatrix2, size_t minimalIndexYMatrix2, size_t size)
{
	auto result = AllocateMatrix(size);
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			result[i][j] = 0;
			for (size_t k = 0; k < size; k++) {
				result[i][j] += matrix1[i + minimalIndexXMatrix1][k + minimalIndexYMatrix1] *
					matrix2[k + minimalIndexXMatrix2][j + minimalIndexYMatrix2];
			}
		}
	}
	return result;
}

int ** StrassenAlgorithm::FastMultiply(int ** matrix1, int ** matrix2, size_t size)
{
	return MultiplySegment(matrix1, matrix1, 0, 0, 0, 0, size);
}

int ** StrassenAlgorithm::MultiplySegment(int ** matrix1, int ** matrix2,
	size_t minimalIndexXMatrix1, size_t minimalIndexYMatrix1,
	size_t minimalIndexXMatrix2, size_t minimalIndexYMatrix2, size_t size)
{
	if (size < 128) {
		return MultiplySync(matrix1, matrix2,
			minimalIndexXMatrix1, minimalIndexYMatrix1,
			minimalIndexXMatrix2, minimalIndexYMatrix2, size);
	}
	RecursionDepth++;
	auto half = size >> 1;
	int ** x1;
	int ** x2;
	int ** x3;
	int ** x4;
	int ** x5;
	int ** x6;
	int ** x7;

	if (RecursionDepth > MAX_DEPTH_THREAD) {
		//p1
		auto x1sum1 = Plus(matrix1, matrix1, minimalIndexXMatrix1, minimalIndexYMatrix1, half + minimalIndexXMatrix1, half + minimalIndexYMatrix1, half);
		auto x1sum2 = Plus(matrix2, matrix2, minimalIndexXMatrix2, minimalIndexYMatrix2, half + minimalIndexXMatrix2, half + minimalIndexYMatrix2, half);
		x1 = MultiplySegment(x1sum1, x1sum2, 0, 0, 0, 0, half);
		FreeMatrix(x1sum1, half);
		FreeMatrix(x1sum2, half);
		//p2
		auto x2sum1 = Plus(matrix1, matrix1, half + minimalIndexXMatrix1, minimalIndexYMatrix1, half + minimalIndexXMatrix1, half + minimalIndexYMatrix1, half);
		x2 = MultiplySegment(x2sum1, matrix2, 0, 0, minimalIndexXMatrix2, minimalIndexYMatrix2, half);
		FreeMatrix(x2sum1, half);
		//p3
		auto x3sum1 = Minus(matrix2, matrix2, minimalIndexXMatrix2, half + minimalIndexYMatrix2, half + minimalIndexXMatrix2, half + minimalIndexYMatrix2, half);
		x3 = MultiplySegment(matrix1, x3sum1, minimalIndexXMatrix1, minimalIndexYMatrix1, 0, 0, half);
		FreeMatrix(x3sum1, half);
		//p4
		auto x4diff1 = Minus(matrix2, matrix2, half + minimalIndexXMatrix2, minimalIndexYMatrix2, minimalIndexXMatrix2, minimalIndexYMatrix2, half);
		x4 = MultiplySegment(matrix1, x4diff1, half + minimalIndexXMatrix1, half + minimalIndexYMatrix1, 0, 0, half);
		FreeMatrix(x4diff1, half);
		//p5
		auto x5sum1 = Plus(matrix1, matrix1, minimalIndexXMatrix1, minimalIndexYMatrix1, minimalIndexXMatrix1, half + minimalIndexYMatrix1, half);
		x5 = MultiplySegment(x5sum1, matrix2, 0, 0, half + minimalIndexXMatrix2, half + minimalIndexYMatrix2, half);
		FreeMatrix(x5sum1, half);
		//p6
		auto x6diff1 = Minus(matrix1, matrix1, half + minimalIndexXMatrix1, minimalIndexYMatrix1, minimalIndexXMatrix1, minimalIndexYMatrix1, half);
		auto x6sum1 = Plus(matrix2, matrix2, minimalIndexXMatrix2, minimalIndexYMatrix2, minimalIndexXMatrix2, half + minimalIndexYMatrix2, half);
		x6 = MultiplySegment(x6diff1, x6sum1, 0, 0, 0, 0, half);
		FreeMatrix(x6diff1, half);
		FreeMatrix(x6sum1, half);
		//p7
		auto x7diff1 = Minus(matrix1, matrix1, minimalIndexXMatrix1, half + minimalIndexYMatrix1, half + minimalIndexXMatrix1, half + minimalIndexYMatrix1, half);
		auto x7sum1 = Plus(matrix2, matrix2, half + minimalIndexXMatrix2, minimalIndexYMatrix2, half + minimalIndexXMatrix2, half + minimalIndexYMatrix2, half);
		x7 = MultiplySegment(x7diff1, x7sum1, 0, 0, 0, 0, half);
		FreeMatrix(x7diff1, half);
		FreeMatrix(x7sum1, half);
	}
	else {
		//p1
		std::thread thread1([=, &matrix1, &matrix2, &x1]() {
			auto x1sum1 = Plus(matrix1, matrix1, minimalIndexXMatrix1, minimalIndexYMatrix1, half + minimalIndexXMatrix1, half + minimalIndexYMatrix1, half);
			auto x1sum2 = Plus(matrix2, matrix2, minimalIndexXMatrix2, minimalIndexYMatrix2, half + minimalIndexXMatrix2, half + minimalIndexYMatrix2, half);
			x1 = MultiplySegment(x1sum1, x1sum2, 0, 0, 0, 0, half);
			FreeMatrix(x1sum1, half);
			FreeMatrix(x1sum2, half);
		});
		//p2
		std::thread thread2([=, &matrix1, &matrix2, &x2]() {
			auto x2sum1 = Plus(matrix1, matrix1, half + minimalIndexXMatrix1, minimalIndexYMatrix1, half + minimalIndexXMatrix1, half + minimalIndexYMatrix1, half);
			x2 = MultiplySegment(x2sum1, matrix2, 0, 0, minimalIndexXMatrix2, minimalIndexYMatrix2, half);
			FreeMatrix(x2sum1, half);
		});
		//p3
		std::thread thread3([=, &matrix1, &matrix2, &x3]() {
			auto x3sum1 = Minus(matrix2, matrix2, minimalIndexXMatrix2, half + minimalIndexYMatrix2, half + minimalIndexXMatrix2, half + minimalIndexYMatrix2, half);
			x3 = MultiplySegment(matrix1, x3sum1, minimalIndexXMatrix1, minimalIndexYMatrix1, 0, 0, half);
			FreeMatrix(x3sum1, half);
		});
		//p4
		std::thread thread4([=, &matrix1, &matrix2, &x4]() {
			auto x4diff1 = Minus(matrix2, matrix2, half + minimalIndexXMatrix2, minimalIndexYMatrix2, minimalIndexXMatrix2, minimalIndexYMatrix2, half);
			x4 = MultiplySegment(matrix1, x4diff1, half + minimalIndexXMatrix1, half + minimalIndexYMatrix1, 0, 0, half);
			FreeMatrix(x4diff1, half);
		});
		//p5
		std::thread thread5([=, &matrix1, &matrix2, &x5]() {
			auto x5sum1 = Plus(matrix1, matrix1, minimalIndexXMatrix1, minimalIndexYMatrix1, minimalIndexXMatrix1, half + minimalIndexYMatrix1, half);
			x5 = MultiplySegment(x5sum1, matrix2, 0, 0, half + minimalIndexXMatrix2, half + minimalIndexYMatrix2, half);
			FreeMatrix(x5sum1, half);
		});
		//p6
		std::thread thread6([=, &matrix1, &matrix2, &x6]() {
			auto x6diff1 = Minus(matrix1, matrix1, half + minimalIndexXMatrix1, minimalIndexYMatrix1, minimalIndexXMatrix1, minimalIndexYMatrix1, half);
			auto x6sum1 = Plus(matrix2, matrix2, minimalIndexXMatrix2, minimalIndexYMatrix2, minimalIndexXMatrix2, half + minimalIndexYMatrix2, half);
			x6 = MultiplySegment(x6diff1, x6sum1, 0, 0, 0, 0, half);
			FreeMatrix(x6diff1, half);
			FreeMatrix(x6sum1, half);
		});
		//p7
		std::thread thread7([=, &matrix1, &matrix2, &x7]() {
			auto x7diff1 = Minus(matrix1, matrix1, minimalIndexXMatrix1, half + minimalIndexYMatrix1, half + minimalIndexXMatrix1, half + minimalIndexYMatrix1, half);
			auto x7sum1 = Plus(matrix2, matrix2, half + minimalIndexXMatrix2, minimalIndexYMatrix2, half + minimalIndexXMatrix2, half + minimalIndexYMatrix2, half);
			x7 = MultiplySegment(x7diff1, x7sum1, 0, 0, 0, 0, half);
			FreeMatrix(x7diff1, half);
			FreeMatrix(x7sum1, half);
		});

		thread1.join();
		thread2.join();
		thread3.join();
		thread4.join();
		thread5.join();
		thread6.join();
		thread7.join();
	}
	auto result = AllocateMatrix(size);
	// C1,1
	Plus(x1, x4, result, 0, 0, 0, 0, 0, 0, half);
	Minus(result, x5, result, 0, 0, 0, 0, 0, 0, half);
	Plus(result, x7, result, 0, 0, 0, 0, 0, 0, half);
	// C1,2
	Plus(x3, x5, result, 0, 0, 0, 0, 0, half, half);
	// C2,1
	Plus(x2, x4, result, 0, 0, 0, 0, half, 0, half);
	// C2,2
	Plus(x1, x3, result, 0, 0, 0, 0, half, half, half);
	Minus(result, x2, result, half, half, 0, 0, half, half, half);
	Plus(result, x6, result, half, half, 0, 0, half, half, half);

	FreeMatrix(x1, half);
	FreeMatrix(x2, half);
	FreeMatrix(x3, half);
	FreeMatrix(x4, half);
	FreeMatrix(x5, half);
	FreeMatrix(x6, half);
	FreeMatrix(x7, half);
	RecursionDepth--;
	return result;
}

int ** StrassenAlgorithm::Plus(int ** matrix1, int ** matrix2,
	size_t minimalIndexXMatrix1, size_t minimalIndexYMatrix1,
	size_t minimalIndexXMatrix2, size_t minimalIndexYMatrix2, size_t size)
{
	auto result = AllocateMatrix(size);
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			result[i][j] = matrix1[i + minimalIndexXMatrix1][j + minimalIndexYMatrix1] +
				matrix2[i + minimalIndexXMatrix2][j + minimalIndexYMatrix2];
		}
	}
	return result;
}

void StrassenAlgorithm::Plus(int ** matrix1, int ** matrix2, int ** result,
	size_t minimalIndexXMatrix1, size_t minimalIndexYMatrix1,
	size_t minimalIndexXMatrix2, size_t minimalIndexYMatrix2,
	size_t minimalIndexXResult, size_t minimalIndexYResult, size_t size)
{
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			result[i + minimalIndexXResult][j + minimalIndexYResult] =
				matrix1[i + minimalIndexXMatrix1][j + minimalIndexYMatrix1] +
				matrix2[i + minimalIndexXMatrix2][j + minimalIndexYMatrix2];
		}
	}
}

int ** StrassenAlgorithm::Minus(int ** matrix1, int ** matrix2,
	size_t minimalIndexXMatrix1, size_t minimalIndexYMatrix1,
	size_t minimalIndexXMatrix2, size_t minimalIndexYMatrix2, size_t size)
{
	auto result = AllocateMatrix(size);
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			result[i][j] = matrix1[i + minimalIndexXMatrix1][j + minimalIndexYMatrix1] -
				matrix2[i + minimalIndexXMatrix2][j + minimalIndexYMatrix2];
		}
	}
	return result;
}

void StrassenAlgorithm::Minus(int ** matrix1, int ** matrix2, int ** result,
	size_t minimalIndexXMatrix1, size_t minimalIndexYMatrix1,
	size_t minimalIndexXMatrix2, size_t minimalIndexYMatrix2,
	size_t minimalIndexXResult, size_t minimalIndexYResult, size_t size)
{
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			result[i + minimalIndexXResult][j + minimalIndexYResult] =
				matrix1[i + minimalIndexXMatrix1][j + minimalIndexYMatrix1] -
				matrix2[i + minimalIndexXMatrix2][j + minimalIndexYMatrix2];
		}
	}
}

void StrassenAlgorithm::Clear() {
	for (auto iter = _allocatedMatrices.begin(); iter != _allocatedMatrices.end(); iter++) {
		auto value = iter->second;
		auto size = iter->first;
		while (value->size() > 0)
		{
			int ** data = value->top();
			for (size_t i = 0; i < size; i++) {
				delete[] data[i];
			}
			delete[] data;
			value->pop();
		}
	}
}
