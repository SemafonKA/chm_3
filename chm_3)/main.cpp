#include <iostream>
#include <fstream>
#include <vector>
#include <format>

#include "Chrono_Timer.h"
#include "Matrix.h"

using namespace std;

namespace Vec {
   inline double Scalar(const vector<double>& l, const vector<double>& r) {
      if (l.size() != r.size()) throw runtime_error("Размеры векторов не совпадают");

      double res = 0.0;
      for (size_t i = 0; i < l.size(); i++)
      {
         res += l[i] * r[i];
      }

      return res;
   }
}

namespace IterSolvers {
   double minEps = 1e-8;
   size_t maxIter = 50000;

   int MSG_Assimetric(Matrix& A, vector<double>& f, vector<double>& x, double& eps, bool debugOutput = false) {
      f = A.TranspMultToVec(f);                             // f* = A^T * f
      vector<double> r = A.TranspMultToVec(A * x);          // r0 = A^T * A * x
      for (uint16_t i = 0; i < r.size(); i++) r[i] = f[i] - r[i];

      vector<double> z = r;               // z0
      vector<double> t;

      double rPrevScalar = Vec::Scalar(r, r);         // (r_k-1, r_k-1)
      double rScalar = 0;
      double a = 0;                       // alpha_k,
      double b = 0;                       // beta_k
      double normF = sqrt(Vec::Scalar(f, f));   // ||f||
      eps = DBL_MAX;

      size_t iter;
      for (iter = 1; iter < maxIter && eps > minEps; iter++)
      {
         t = A.TranspMultToVec(A * z);          // t = A^t * A * z_k-1
         a = rPrevScalar / Vec::Scalar(t, z);   // a_k = (r_k-1, r_k-1) / (t_k-1, z_k-1)

         for (uint16_t i = 0; i < x.size(); i++)      
         {
            x[i] += a * z[i];                         // x_k = x_k-1 + a * z_k-1
            r[i] -= a * t[i];                         // r_k = r_k-1 - a * t_k-1
         }
         rScalar = Vec::Scalar(r, r);
         b = rScalar / rPrevScalar;                   // b = (r_k, r_k) / (r_k-1, r_k-1)

         for (uint16_t i = 0; i < x.size(); i++)
         {
            z[i] = r[i] + b * z[i];                   // z_k = r_k + b * z_k-1
         }

         rPrevScalar = rScalar;
         eps = sqrt(rPrevScalar) / normF;

         // Выводим на то же место, что и раньше (со сдвигом каретки)
         if (debugOutput)
         {
            cout << format("\rИтерация: {0:<10} относительная невязка: {1:<15.3e}", iter, eps);
         }
         if (isinf(eps))
         {
            break;
         }
      }

      if (debugOutput)
      {
         cout << endl;
         if (isinf(eps))
         {
            cout << "Выход по переполнению метода" << endl << endl;
         }
         else if (iter > maxIter)
         {
            cout << "Выход по числу итераций" << endl << endl;
         }
         else
         {
            cout << "Выход по относительной невязке" << endl << endl;
         }
      }

      return iter;
   }
};

int main() {
   setlocale(LC_ALL, "ru-RU");

   uint16_t matrixSize;
   Matrix mat;
   vector<double> f;
   vector<double> x;

   auto kuslauF = ifstream("./iofiles/kuslau.txt");
   if (!kuslauF.is_open())
   {
      cout << "Файл ./iofiles/kuslau.txt отсутствует в директории" << endl;
      return 1;
   }
   kuslauF >> matrixSize >> IterSolvers::maxIter >> IterSolvers::minEps;
   kuslauF.close();

   try
   {
      mat = Matrix::ReadFromFiles(matrixSize, "./iofiles/ig.txt", "./iofiles/jg.txt", "./iofiles/ggl.txt", "./iofiles/ggu.txt", "./iofiles/di.txt");
      f = ReadVecFromFile(matrixSize, "./iofiles/pr.txt");
      x = ReadVecFromFile(matrixSize, "./iofiles/initX.txt");
   }
   catch (exception& e)
   {
      cout << e.what() << endl;
      return 1;
   }

   cout << "Все данные успешно считанны из файлов." << endl;
   cout << "Выберите метод для решения СЛАУ: " << endl;
   cout << "  1) МСГ для несимметричных матриц (без предобуславливания)" << endl;
   //cout << "  2) МСГ для нессиметричных матриц (диагональное предобуславливание)" << endl;

   int userCase;
   cin >> userCase;
   switch (userCase)
   {
      case 1:
      {
         cout << "Начало вычислений для метода МСГ для несимметричных матриц (без предобуславливания)" << endl << endl;
         Timer timer;
         double eps = 0;
         IterSolvers::MSG_Assimetric(mat, f, x, eps, true);
         timer.elapsed();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         break;
      }

      default:
         break;
   }


   auto outFile = ofstream("./iofiles/resultX.txt");
   outFile.precision(15);
   outFile.setf(std::ios_base::fixed);
   cout.precision(15);
   cout.setf(std::ios_base::fixed);
   for (auto& el : x)
   {
      outFile << el << endl;
      cout << el << endl;
   }
   outFile.close();

   return 0;
}