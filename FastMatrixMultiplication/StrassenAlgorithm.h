#pragma once
class StrassenAlgorithm
{
private:
	static std::map<size_t, std::stack<int**>*> _allocatedMatrices;
	static std::mutex _syncMutex;

	static int ** AllocateMatrix(size_t size);
	static void FreeMatrix(int ** matrix, size_t size);

public:
	static int CacheMissed;
	static int CacheFound;
	static int RecursionDepth;
	static int ThreadsCount;

	static int ** FastMultiply(int ** matrix1, int ** matrix2, size_t size);
	static int ** Multiply(int ** matrix1, int ** matrix2, size_t size);
	static int ** MultiplySync(int ** matrix1, int ** matrix2,
		size_t minimalIndexXMatrix1, size_t minimalIndexYMatrix1,
		size_t minimalIndexXMatrix2, size_t minimalIndexYMatrix2,
		size_t size);
	static int ** MultiplySegment(int ** matrix1, int ** matrix2,
		size_t minimalIndexXMatrix1, size_t minimalIndexYMatrix1,
		size_t minimalIndexXMatrix2, size_t minimalIndexYMatrix2,
		size_t size);
	static int ** Plus(int ** matrix1, int ** matrix2,
		size_t minimalIndexXMatrix1, size_t minimalIndexYMatrix1,
		size_t minimalIndexXMatrix2, size_t minimalIndexYMatrix2,
		size_t size);
	static void Plus(int ** matrix1, int ** matrix2, int ** result,
		size_t minimalIndexXMatrix1, size_t minimalIndexYMatrix1,
		size_t minimalIndexXMatrix2, size_t minimalIndexYMatrix2,
		size_t minimalIndexXResult, size_t minimalIndexYResult,
		size_t size);
	static int ** Minus(int ** matrix1, int ** matrix2,
		size_t minimalIndexXMatrix1, size_t minimalIndexYMatrix1,
		size_t minimalIndexXMatrix2, size_t minimalIndexYMatrix2,
		size_t size);
	static void Minus(int ** matrix1, int ** matrix2, int ** result,
		size_t minimalIndexXMatrix1, size_t minimalIndexYMatrix1,
		size_t minimalIndexXMatrix2, size_t minimalIndexYMatrix2,
		size_t minimalIndexXResult, size_t minimalIndexYResult,
		size_t size);
	static void Clear();
};

