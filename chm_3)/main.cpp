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

   // l or r may be similar vectors to ans
   inline void Mult(const vector<double>& l, const vector<double>& r, vector<double>& ans) {
      if (ans.size() != l.size() || ans.size() != r.size()) throw runtime_error("Ошибка: размеры векторов должны совпадать.");

      for (size_t i = 0; i < ans.size(); i++)
      {
         ans[i] = l[i] * r[i];
      }
   }
   inline vector<double> Mult(const vector<double>& l, const vector<double>& r) {
      if (r.size() != l.size()) throw runtime_error("Ошибка: размеры векторов должны совпадать.");
      vector<double> ans(l.size());

      for (size_t i = 0; i < ans.size(); i++)
      {
         ans[i] = l[i] * r[i];
      }
      return ans;
   }
}

namespace IterSolvers {
   double minEps = 1e-8;
   size_t maxIter = 2000;
   bool globalDebugOutput = false;

   vector<double>* _tmp1 = nullptr, * _tmp2 = nullptr,
      * _tmp3 = nullptr, * _tmp4 = nullptr, * _tmp5 = nullptr, * _tmp6 = nullptr;

   inline void VecInit(vector<double>*& vec, size_t size) {
      if (vec == nullptr)
      {
         vec = new vector<double>(size);
      }
      else if (vec->size() != size)
      {
         vec->resize(size);
      }
   }

   namespace MSG_Assimetric {
      inline void Init_Default(size_t size) {
         VecInit(_tmp1, size); // Массив для вектора r метода
         VecInit(_tmp2, size); // Массив для вектора z
         VecInit(_tmp3, size); // Массив для вектора t
         VecInit(_tmp4, size); // Массив для временного вектора
      }

      size_t Default(Matrix& A, vector<double>& f, vector<double>& x, double& eps, bool debugOutput = globalDebugOutput) {
         size_t size = x.size();
         Init_Default(size);

         vector<double>& r = *_tmp1;
         vector<double>& tmp = *_tmp4;
         A.MultToVec(x, tmp);
         for (uint16_t i = 0; i < size; i++) tmp[i] = f[i] - tmp[i]; // r0 = f - A * x
         A.TranspMultToVec(tmp, r);                                  // r0 = A^t * (f - A * x)

         vector<double>& z = *_tmp2;
         z = r;                        // z0
         vector<double>& t = *_tmp3;

         double rPrevScalar = Vec::Scalar(r, r);         // (r_k-1, r_k-1)
         double rScalar = 0;
         double a = 0;                       // alpha_k,
         double b = 0;                       // beta_k
         double normF = Vec::Scalar(f, f);   // ||f||
         eps = DBL_MAX;

         size_t iter;
         for (iter = 1; iter <= maxIter && eps > minEps; iter++)
         {
            A.MultToVec(z, tmp);
            A.TranspMultToVec(tmp, t);             // t = A^t * A * z_k-1
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
            eps = sqrt(rPrevScalar / normF);

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

         return iter - 1;
      }

      inline void Init_DiagPrecond(size_t size) {
         Init_Default(size);
         VecInit(_tmp5, size);      // Массив для вектора D
      }

      size_t DiagPrecond(Matrix& A, vector<double>& f, vector<double>& x, double& eps, bool debugOutput = globalDebugOutput) {
         size_t size = x.size();
         Init_DiagPrecond(size);

         vector<double>& D = *_tmp5;         // D = обратный корень от диагонали матрицы
         for (uint16_t i = 0; i < size; i++) D[i] = 1 / sqrt(A.di[i]);

         for (uint16_t i = 0; i < size; i++) x[i] /= D[i];     // local_x

         vector<double>& r = *_tmp1;             // r = U^-t * A^t * L^-t * L^-1 (f - A * x)
         vector<double>& tmp = *_tmp4;
         A.MultToVec(x, tmp);
         for (uint16_t i = 0; i < size; i++) tmp[i] = f[i] - tmp[i];
         Vec::Mult(D, tmp, tmp);
         Vec::Mult(D, tmp, tmp);
         A.TranspMultToVec(tmp, r);
         Vec::Mult(D, r, r);

         vector<double>& z = *_tmp2;
         z = r;

         vector<double>& t = *_tmp3;             // t = U^-1 * A^t * L^-t * L^-1 * A * U^-1 * z

         double rPrevScalar = Vec::Scalar(r, r);         // (r_k-1, r_k-1)
         double rScalar = 0;
         double a = 0;                       // alpha_k,
         double b = 0;                       // beta_k
         double normF = Vec::Scalar(f, f);   // ||f||
         eps = sqrt(rPrevScalar / normF);

         size_t iter;
         for (iter = 1; iter <= maxIter && eps > minEps; iter++)
         {
            Vec::Mult(D, z, t);
            A.MultToVec(t, tmp);
            Vec::Mult(D, tmp, tmp);
            Vec::Mult(D, tmp, tmp);
            A.TranspMultToVec(tmp, t);
            Vec::Mult(D, t, t);

            a = rPrevScalar / Vec::Scalar(t, z);         // a_k = (r_k-1, r_k-1) / (t_k-1, z_k-1)
            for (uint16_t i = 0; i < size; i++)
            {
               x[i] += a * z[i];                         // local_x_k = local_x_k-1 + a * z_k-1
               r[i] -= a * t[i];                         // r_k = r_k-1 - a * t_k-1
            }

            rScalar = Vec::Scalar(r, r);
            b = rScalar / rPrevScalar;                   // b = (r_k, r_k) / (r_k-1, r_k-1)

            for (uint16_t i = 0; i < size; i++)
            {
               z[i] = r[i] + b * z[i];                   // z_k = r_k + b * z_k-1
            }

            rPrevScalar = rScalar;
            eps = sqrt(rPrevScalar / normF);

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
         Vec::Mult(D, x, x);        // x = U^-1 * local_x

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

         return iter - 1;
      }
   }

   namespace LOS {
      size_t resetIter = 10;

      void Init_Default(size_t size) {
         VecInit(_tmp1, size); // Массив для вектора r метода
         VecInit(_tmp2, size); // Массив для вектора z
         VecInit(_tmp3, size); // Массив для вектора p
         VecInit(_tmp4, size); // Массив для вектора Ar
      }

      size_t Default(Matrix& A, vector<double>& f, vector<double>& x, double& eps, bool debugOutput = globalDebugOutput) {
         uint16_t size = x.size();
         Init_Default(size);

         vector<double>& r = *_tmp1;
         A.MultToVec(x, r);
         for (uint16_t i = 0; i < size; i++) r[i] = f[i] - r[i]; // r0 = f - A * x

         vector<double>& z = *_tmp2;         // z0
         z = r;
         vector<double>& p = *_tmp3;         // p0 = A * z0
         A.MultToVec(z, p);
         vector<double>& Ar = *_tmp4;        // A * r

         double ppScalar;
         double nev = Vec::Scalar(r, r);
         double ffScalar = Vec::Scalar(f, f);
         eps = nev / ffScalar;
         double a;                  // alpha
         double b;                  // beta
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

            A.MultToVec(r, Ar);      // A * r_k
            b = -Vec::Scalar(p, Ar) / ppScalar; // b = - (p_k-1, A * r_k) / (p_k-1, p_k-1)

            for (uint16_t i = 0; i < size; i++)
            {
               z[i] = r[i] + b * z[i];             // [z_k] = r_k + b * [z_k-1]
               p[i] = Ar[i] + b * p[i];            // [p_k] = A * r_k + b * [p_k-1]
            }

            if (iter % resetIter == 0)
            {
               A.MultToVec(x, r);
               for (uint16_t i = 0; i < size; i++) r[i] = f[i] - r[i];
               z = r;
               A.MultToVec(z, p);
            }
            nev = Vec::Scalar(r, r);
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

         return iter - 1;
      }

      void Init_DiagPrecond(size_t size) {
         VecInit(_tmp1, size); // Массив для вектора r метода
         VecInit(_tmp2, size); // Массив для вектора z
         VecInit(_tmp3, size); // Массив для вектора p
         VecInit(_tmp4, size); // Массив для вектора Ar
         VecInit(_tmp5, size); // Массив для вектора D
         VecInit(_tmp6, size); // Массив для вектора tmp
      }

      size_t DiagPrecond(Matrix& A, vector<double>& f, vector<double>& x, double& eps, bool debugOutput = globalDebugOutput) {
         uint16_t size = x.size();
         Init_DiagPrecond(size);

         vector<double>& D = *_tmp5;               // обратный корень от диагонали матрицы
         for (uint16_t i = 0; i < size; i++) D[i] = 1 / sqrt(A.di[i]);

         vector<double>& r = *_tmp1;               // r0 = L^-1 * (f - A * x)
         A.MultToVec(x, r);
         for (uint16_t i = 0; i < size; i++) r[i] = f[i] - r[i];
         Vec::Mult(D, r, r);

         vector<double>& z = *_tmp2;               // z0 = U^-1 * r
         Vec::Mult(D, r, z);

         vector<double>& p = *_tmp3;               // p0 = L^-1 * A * z0
         A.MultToVec(z, p);
         Vec::Mult(D, p, p);

         vector<double>& Ar = *_tmp4;              // Ar = L^-1 * A * U^-1 * r
         vector<double>& tmp = *_tmp6;

         double ppScalar;
         double nev = Vec::Scalar(r, r);
         double ffScalar = Vec::Scalar(f, f);
         eps = nev / ffScalar;
         double a;                  // alpha
         double b;                  // beta
         size_t iter;

         for (iter = 1; iter <= maxIter && eps > minEps; iter++)
         {
            ppScalar = Vec::Scalar(p, p);          // (p_k-1, p_k-1)
            a = Vec::Scalar(p, r) / ppScalar;      // (p_k-1, r_k-1) / (p_k-1, p_k-1)

            for (uint16_t i = 0; i < size; i++)
            {
               x[i] += a * z[i];                   // [x_k] = [x_k-1] + a*z_k-1
               r[i] -= a * p[i];                   // [r_k] = [r_k-1] - a*p_k-1
            }

            Vec::Mult(D, r, tmp);
            A.MultToVec(tmp, Ar);
            Vec::Mult(D, Ar, Ar);                  // Ar = L^-1 * A * U^-1 * r

            b = -Vec::Scalar(p, Ar) / ppScalar;    // b = - (p_k-1, L^-1 * A * U^-1 * r_k) / (p_k-1, p_k-1)
            Vec::Mult(D, r, tmp);                  // tmp = U^-1 * r_k

            for (uint16_t i = 0; i < size; i++)
            {
               z[i] = tmp[i] + b * z[i];            // [z_k] = U^-1 * r_k + b * [z_k-1]
               p[i] = Ar[i] + b * p[i];             // [p_k] = A * r_k + b * [p_k-1]
            }

            if (iter % resetIter == 0)
            {
               vector<double> r = A * x;           // r0 = L^-1 * (f - A * x)
               for (uint16_t i = 0; i < size; i++) r[i] = f[i] - r[i];
               Vec::Mult(D, r, r);

               vector<double> z = Vec::Mult(D, r); // z0 = U^-1 * r

               vector<double> p = A * z;           // p0 = L^-1 * A * z0
               Vec::Mult(D, p, p);
            }
            nev = Vec::Scalar(r, r);
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

         return iter - 1;
      }
   }

   void Destruct() {
      delete _tmp1, _tmp2, _tmp3, _tmp4, _tmp5, _tmp6;
      _tmp1 = _tmp2 = _tmp3 = _tmp4 = _tmp5 = _tmp6 = nullptr;
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

   IterSolvers::globalDebugOutput = false;

   cout << "Все данные успешно считанны из файлов." << endl;
   cout << "Выберите метод для решения СЛАУ: " << endl;
   cout << "  1) МСГ для несимметричных матриц (без предобуславливания)" << endl;
   cout << "  2) МСГ для нессиметричных матриц (диагональное предобуславливание)" << endl;
   cout << "  3) ЛОС (без предобуславливания)" << endl;
   cout << "  4) ЛОС (диагональное предобуславливание)" << endl;

   int userCase;
   cin >> userCase;
   switch (userCase)
   {
      case 1:
      {
         cout << "Начало вычислений для метода МСГ для несимметричных матриц (без предобуславливания)" << endl << endl;
         Timer timer;
         double eps = 0;
         size_t it = IterSolvers::MSG_Assimetric::Default(mat, f, x, eps);
         timer.elapsed();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         cout << "Количество итераций: " << it << endl;
         cout << "Относительная невязка: " << eps << endl;
         break;
      }
      case 2:
      {
         cout << "Начало вычислений для метода МСГ для несимметричных матриц (диагональное предобуславливание)" << endl << endl;
         Timer timer;
         double eps = 0;
         size_t it = IterSolvers::MSG_Assimetric::DiagPrecond(mat, f, x, eps);
         timer.elapsed();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         cout << "Количество итераций: " << it << endl;
         cout << "Относительная невязка: " << eps << endl;
         break;
      }
      case 3:
      {
         cout << "Начало вычислений для метода ЛОС (без предобуславливания)" << endl << endl;
         Timer timer;
         double eps = 0;
         size_t it = IterSolvers::LOS::Default(mat, f, x, eps);
         timer.elapsed();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         cout << "Количество итераций: " << it << endl;
         cout << "Относительная невязка: " << eps << endl;
         break;
      }
      case 4:
      {
         cout << "Начало вычислений для метода ЛОС (диагональное предобуславливание)" << endl << endl;
         Timer timer;
         double eps = 0;
         size_t it = IterSolvers::LOS::DiagPrecond(mat, f, x, eps);
         timer.elapsed();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         cout << "Количество итераций: " << it << endl;
         cout << "Относительная невязка: " << eps << endl;
         break;
      }
      default:
         break;
   }
   IterSolvers::Destruct();

   if (IterSolvers::globalDebugOutput)
   {
      cout << "Полученное решение: " << endl;
      cout.precision(15);
      cout.setf(std::ios_base::fixed);
   }

   auto outFile = ofstream("./iofiles/resultX.txt");
   outFile.precision(15);
   outFile.setf(std::ios_base::fixed);
   for (auto& el : x)
   {
      outFile << el << endl;
      if (IterSolvers::globalDebugOutput)
      {
         cout << el << endl;
      }
   }
   outFile.close();

   return 0;
}