#include <iostream>
#include <fstream>
#include <vector>
#include <format>

#include "Chrono_Timer.h"
#include "Matrix.h"

using namespace std;

namespace Vec {
   inline double Scalar(const vector<double>& l, const vector<double>& r) {
      if (l.size() != r.size()) throw runtime_error("������� �������� �� ���������");

      double res = 0.0;
      for (size_t i = 0; i < l.size(); i++)
      {
         res += l[i] * r[i];
      }

      return res;
   }

   // l or r may be similar vectors to ans
   inline void Mult(const vector<double>& l, const vector<double>& r, vector<double>& ans) {
      if (ans.size() != l.size() || ans.size() != r.size()) throw runtime_error("������: ������� �������� ������ ���������.");

      for (size_t i = 0; i < ans.size(); i++)
      {
         ans[i] = l[i] * r[i];
      }
   }
   inline vector<double> Mult(const vector<double>& l, const vector<double>& r) {
      if (r.size() != l.size()) throw runtime_error("������: ������� �������� ������ ���������.");
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

   namespace MSG_Assimetric {

      size_t Default(Matrix& A, vector<double>& f, vector<double>& x, double& eps, bool debugOutput = globalDebugOutput) {
         size_t size = x.size();

         vector<double> r = A * x;
         for (uint16_t i = 0; i < size; i++) r[i] = f[i] - r[i]; // r0 = f - A * x
         r = A.TranspMultToVec(r);                             // r0 = A^t * (f - A * x)

         vector<double> z = r;               // z0
         vector<double> t(size);

         double rPrevScalar = Vec::Scalar(r, r);         // (r_k-1, r_k-1)
         double rScalar = 0;
         double a = 0;                       // alpha_k,
         double b = 0;                       // beta_k
         double normF = Vec::Scalar(f, f);   // ||f||
         eps = DBL_MAX;

         size_t iter;
         for (iter = 1; iter <= maxIter && eps > minEps; iter++)
         {
            A.TranspMultToVec(A * z, t);          // t = A^t * A * z_k-1
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

            // ������� �� �� �� �����, ��� � ������ (�� ������� �������)
            if (debugOutput)
            {
               cout << format("\r��������: {0:<10} ������������� �������: {1:<15.3e}", iter, eps);
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
               cout << "����� �� ������������ ������" << endl << endl;
            }
            else if (iter > maxIter)
            {
               cout << "����� �� ����� ��������" << endl << endl;
            }
            else
            {
               cout << "����� �� ������������� �������" << endl << endl;
            }
         }

         return iter - 1;
      }

      size_t DiagPrecond(Matrix& A, vector<double>& f, vector<double>& x, double& eps, bool debugOutput = globalDebugOutput) {
         size_t size = x.size();

         vector<double> D(size);             // D = �������� ������ �� ��������� �������
         for (uint16_t i = 0; i < size; i++) D[i] = sqrt(A.di[i]);

         vector<double> local_x(size);
         Vec::Mult(D, x, local_x);
         for (uint16_t i = 0; i < size; i++) D[i] = 1 / D[i];

         vector<double> r(size);             // r = U^-t * A^t * L^-t * L^-1 (f - A * x)
         vector<double> tmp = A * x;
         for (uint16_t i = 0; i < size; i++) tmp[i] = f[i] - tmp[i];
         Vec::Mult(D, tmp, tmp);
         Vec::Mult(D, tmp, tmp);
         A.TranspMultToVec(tmp, r);
         Vec::Mult(D, r, r);

         vector<double> z = r;

         vector<double> t(size);             // t = U^-1 * A^t * L^-t * L^-1 * A * U^-1 * z

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

            a = rPrevScalar / Vec::Scalar(t, z);   // a_k = (r_k-1, r_k-1) / (t_k-1, z_k-1)
            for (uint16_t i = 0; i < size; i++)
            {
               local_x[i] += a * z[i];                   // local_x_k = local_x_k-1 + a * z_k-1
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

            // ������� �� �� �� �����, ��� � ������ (�� ������� �������)
            if (debugOutput)
            {
               cout << format("\r��������: {0:<10} ������������� �������: {1:<15.3e}", iter, eps);
            }
            if (isinf(eps))
            {
               break;
            }
         }
         Vec::Mult(D, local_x, x);        // x = U^-1 * local_x

         if (debugOutput)
         {
            cout << endl;
            if (isinf(eps))
            {
               cout << "����� �� ������������ ������" << endl << endl;
            }
            else if (iter > maxIter)
            {
               cout << "����� �� ����� ��������" << endl << endl;
            }
            else
            {
               cout << "����� �� ������������� �������" << endl << endl;
            }
         }

         return iter - 1;
      }
   }

   namespace LOS {
      size_t resetIter = 10;

      size_t Default(Matrix& A, vector<double>& f, vector<double>& x, double& eps, bool debugOutput = globalDebugOutput) {
         vector<double> r = A * x;
         for (uint16_t i = 0; i < r.size(); i++) r[i] = f[i] - r[i]; // r0 = f - A * x

         vector<double> z = r;      // z0
         vector<double> p = A * z;  // p0 = A * z0
         vector<double> Ar(x.size());         // A * r

         double ppScalar;
         double nev = Vec::Scalar(r, r);
         double ffScalar = Vec::Scalar(f, f);
         eps = nev / ffScalar;
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

            // ������� �� �� �� �����, ��� � ������ (�� ������� �������)
            if (debugOutput)
            {
               //cout << format("��������: {0:<10} ������������� �������: {1:<15.3e}\n", iter, eps);
               cout << format("\r��������: {0:<10} ������������� �������: {1:<15.3e}", iter, eps);
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
               cout << "����� �� ������������ ������" << endl << endl;
            }
            else if (iter > maxIter)
            {
               cout << "����� �� ����� ��������" << endl << endl;
            }
            else
            {
               cout << "����� �� ������������� �������" << endl << endl;
            }
         }

         return iter - 1;
      }

      size_t DiagPrecond(Matrix& A, vector<double>& f, vector<double>& x, double& eps, bool debugOutput = globalDebugOutput) {
         uint16_t size = x.size();

         vector<double> D(size);       // �������� ������ �� ��������� �������
         for (uint16_t i = 0; i < size; i++) D[i] = 1 / sqrt(A.di[i]);

         vector<double> r = A * x;     // r0 = L^-1 * (f - A * x)
         for (uint16_t i = 0; i < size; i++) r[i] = f[i] - r[i];
         Vec::Mult(D, r, r);

         vector<double> z = Vec::Mult(D, r);      // z0 = U^-1 * r

         vector<double> p = A * z;     // p0 = L^-1 * A * z0
         Vec::Mult(D, p, p);

         vector<double> Ar(size);      // Ar = L^-1 * A * U^-1 * r
         vector<double> tmp(size);

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

            // ������� �� �� �� �����, ��� � ������ (�� ������� �������)
            if (debugOutput)
            {
               //cout << format("��������: {0:<10} ������������� �������: {1:<15.3e}\n", iter, eps);
               cout << format("\r��������: {0:<10} ������������� �������: {1:<15.3e}", iter, eps);
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
               cout << "����� �� ������������ ������" << endl << endl;
            }
            else if (iter > maxIter)
            {
               cout << "����� �� ����� ��������" << endl << endl;
            }
            else
            {
               cout << "����� �� ������������� �������" << endl << endl;
            }
         }

         return iter - 1;
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
      cout << "���� ./iofiles/kuslau.txt ����������� � ����������" << endl;
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

   cout << "��� ������ ������� �������� �� ������." << endl;
   cout << "�������� ����� ��� ������� ����: " << endl;
   cout << "  1) ��� ��� �������������� ������ (��� ������������������)" << endl;
   cout << "  2) ��� ��� �������������� ������ (������������ ������������������)" << endl;
   cout << "  3) ��� (��� ������������������)" << endl;
   cout << "  4) ��� (������������ ������������������)" << endl;

   int userCase;
   cin >> userCase;
   switch (userCase)
   {
      case 1:
      {
         cout << "������ ���������� ��� ������ ��� ��� �������������� ������ (��� ������������������)" << endl << endl;
         Timer timer;
         double eps = 0;
         size_t it = IterSolvers::MSG_Assimetric::Default(mat, f, x, eps);
         timer.elapsed();
         cout << "����� �������� ������ �� " << timer.elapsedValue * 1000 << " ��" << endl << endl;
         cout << "���������� ��������: " << it << endl;
         cout << "������������� �������: " << eps << endl;
         break;
      }
      case 2:
      {
         cout << "������ ���������� ��� ������ ��� ��� �������������� ������ (������������ ������������������)" << endl << endl;
         Timer timer;
         double eps = 0;
         size_t it = IterSolvers::MSG_Assimetric::DiagPrecond(mat, f, x, eps);
         timer.elapsed();
         cout << "����� �������� ������ �� " << timer.elapsedValue * 1000 << " ��" << endl << endl;
         cout << "���������� ��������: " << it << endl;
         cout << "������������� �������: " << eps << endl;
         break;
      }
      case 3:
      {
         cout << "������ ���������� ��� ������ ��� (��� ������������������)" << endl << endl;
         Timer timer;
         double eps = 0;
         size_t it = IterSolvers::LOS::Default(mat, f, x, eps);
         timer.elapsed();
         cout << "����� �������� ������ �� " << timer.elapsedValue * 1000 << " ��" << endl << endl;
         cout << "���������� ��������: " << it << endl;
         cout << "������������� �������: " << eps << endl;
         break;
      }
      case 4:
      {
         cout << "������ ���������� ��� ������ ��� (������������ ������������������)" << endl << endl;
         Timer timer;
         double eps = 0;
         size_t it = IterSolvers::LOS::DiagPrecond(mat, f, x, eps);
         timer.elapsed();
         cout << "����� �������� ������ �� " << timer.elapsedValue * 1000 << " ��" << endl << endl;
         cout << "���������� ��������: " << it << endl;
         cout << "������������� �������: " << eps << endl;
         break;
      }
      default:
         break;
   }


   if (IterSolvers::globalDebugOutput)
   {
      cout << "���������� �������: " << endl;
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