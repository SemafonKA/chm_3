#include "pch.h"
#include "CppUnitTest.h"
#include "../chm_3)/IterSolvers.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTests
{
	TEST_CLASS(UnitTests)
	{
	public:
		
		TEST_METHOD(TestMethod1)
		{
			
		}
	};
}

namespace MatricesTests 
{
	TEST_CLASS(MatrixMults)
	{
	private:
		SparseMatrix M;
		std::vector<double> x;
		std::vector<double> f;
		std::vector<double> f_transp;

		void InitMatVecAns() {
			// SparseMatrix:
			// ----------------------
			// | 5 | 2 | 1 | 0 | 0  |
			// | 3 | 6 | 2 | 1 | 0  |
			// | 0 | 4 | 8 | 1 | 2  |
			// | 0 | 0 | 3 | 5 | 1  |
			// | 0 | 0 | 6 | 2 | 10 |
			// ----------------------

			// x:
			// { 1, 2, 3, 4, 5 }

			// f: 
			// { 12, 25, 46, 34, 76 }

			// f_transp:
			// { 11, 26, 71, 35, 60 }
			x = { 1, 2, 3, 4, 5 };
			f = { 12, 25, 46, 34, 76 };
			f_transp = { 11, 26, 71, 35, 60 };
			M.di = { 5, 6, 8, 5, 10 };
			M.ig = { 0, 0, 1, 3, 5, 7 };
			M.jg = { 0, 0, 1, 1, 2, 2, 3 };
			M.ggl = { 3, 0, 4, 0, 3, 6, 2 };
			M.ggu = { 2, 1, 2, 1, 1, 2, 1 };
		}

	public:
		TEST_METHOD(MultToVec) 
		{
			InitMatVecAns();
			auto ans = M.MultToVec(x);
			for (size_t i = 0; i < ans.size(); i++)
			{
				Assert::AreEqual(ans[i], f[i]);
			}
		}
		TEST_METHOD(TranspMultToVec) 
		{
			InitMatVecAns();
			auto ans = M.TranspMultToVec(x);
			for (size_t i = 0; i < ans.size(); i++)
			{
				Assert::AreEqual(ans[i], f_transp[i]);
			}
		}
	};
}