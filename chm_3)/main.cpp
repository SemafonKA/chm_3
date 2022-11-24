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
   size_t maxIter = 2000;

   namespace MSG_Assimetric {
      size_t Default(Matrix& A, vector<double>& f, vector<double>& x, double& eps, bool debugOutput = false) {
      vector<double> r = A * x;
      for (uint16_t i = 0; i < r.size(); i++) r[i] = f[i] - r[i]; // r0 = f - A * x
      r = A.TranspMultToVec(r);                             // r0 = A^t * (f - A * x)

      vector<double> z = r;               // z0
      vector<double> t;

      double rPrevScalar = Vec::Scalar(r, r);         // (r_k-1, r_k-1)
      double rScalar = 0;
      double a = 0;                       // alpha_k,
      double b = 0;                       // beta_k
      double normF = sqrt(Vec::Scalar(f, f));   // ||f||
      size_t size = x.size();
      eps = DBL_MAX;

      size_t iter;
      for (iter = 1; iter <= maxIter && eps > minEps; iter++)
      {
         t = A.TranspMultToVec(A * z);          // t = A^t * A * z_k-1
         a = rPrevScalar / Vec::Scalar(t, z);   // a_k = (r_k-1, r_k-1) / (t_k-1, z_k-1)

         for (uint16_t i = 0; i < size; i++)      
         {
            x[i] += a * z[i];                         // x_k = x_k-1 + a * z_k-1
            r[i] -= a * t[i];                         // r_k = r_k-1 - a * t_k-1
         }
         rScalar = Vec::Scalar(r, r);
         b = rScalar / rPrevScalar;                   // b = (r_k, r_k) / (r_k-1, r_k-1)

         for (uint16_t i = 0; i < size; i++)
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
   }

   namespace LOS {
      size_t resetIter = 10;
      size_t Default(Matrix& A, vector<double>& f, vector<double>& x, double& eps, bool debugOutput = false) {
         vector<double> r = A * x;
         for (uint16_t i = 0; i < r.size(); i++) r[i] = f[i] - r[i]; // r0 = f - A * x

         vector<double> z = r;      // z0
         vector<double> p = A * z;  // p0 = A * z0
         vector<double> Ar;         // A * r

         double ppScalar;
         double nev = Vec::Scalar(r, r);
         double ffScalar = Vec::Scalar(f, f);
         eps = DBL_MAX;
         double a;                  // alpha
         double b;                  // beta
         uint16_t size = x.size();
         size_t iter;

         for (iter = 1; iter <= maxIter && eps > minEps; iter++)
         {
            ppScalar = Vec::Scalar(p, p);     // (p_k-1, p_k-1)
            a = Vec::Scalar(p, r) / ppScalar;    // (p_k-1, r_k-1) / (p_k-1, p_k-1)

            for (uint16_t i = 0; i < size; i++)
            {
               x[i] += a * z[i];                   // [x_k] = [x_k-1] + a*z_k-1
               r[i] -= a * p[i];                   // [r_k] = [r_k-1] - a*p_k-1
            }

            Ar = A * r;                            // A * r_k
            b = -Vec::Scalar(p, Ar) / ppScalar; // b = - (p_k-1, A * r_k)

            for (uint16_t i = 0; i < size; i++)
            {
               z[i] = r[i] + b * z[i];             // [z_k] = r_k + b * [z_k-1]
               p[i] = Ar[i] + b * p[i];            // [p_k] = A * r_k + b * [p_k-1]
            }

            if (iter % resetIter == 0)    //iter% resetIter == 0
            {
               r = A * x;
               for (uint16_t i = 0; i < size; i++) r[i] = f[i] - r[i];
               z = r;
               p = A * z;
               nev = Vec::Scalar(r, r);
            }
            else
            {
               nev -= a * a * ppScalar;               // [(r_k, r_k)] = [(r_k-1, r_k-1)] - a * a * (p_k-1, p_k-1)
            }
            eps = sqrt(nev / ffScalar);

            // Выводим на то же место, что и раньше (со сдвигом каретки)
            if (debugOutput)
            {
               //cout << format("Итерация: {0:<10} относительная невязка: {1:<15.3e}\n", iter, eps);
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
   cout << "  3) ЛОС (без предобуславливания)" << endl;

   int userCase;
   cin >> userCase;
   switch (userCase)
   {
      case 1:
      {
         cout << "Начало вычислений для метода МСГ для несимметричных матриц (без предобуславливания)" << endl << endl;
         Timer timer;
         double eps = 0;
         IterSolvers::MSG_Assimetric::Default(mat, f, x, eps, true);
         timer.elapsed();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         break;
      }
      case 3:
      {
         cout << "Начало вычислений для метода ЛОС (без предобуславливания)" << endl << endl;
         Timer timer;
         double eps = 0;
         IterSolvers::LOS::Default(mat, f, x, eps, true);
         timer.elapsed();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         break;
      }
      default:
         break;
   }

   cout << "Полученное решение: " << endl;

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