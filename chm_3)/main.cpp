#include <iostream>
#include <fstream>

#include "Chrono_Timer.h"
#include "IterSolvers.h"

using namespace std;

int main() {
   setlocale(LC_ALL, "ru-RU");

   uint16_t matrixSize;
   Matrix mat;
   vector<double> f;
   vector<double> x;

   cout << "*** ��������� ��� ���������� ���� ����������� �������� ***" << endl;
   cout << "���������� ���������� ������ �� �����..." << endl;

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

   cout << "��� ������ ������� �������� �� ������." << endl << endl;
   cout << "�������� ����� ��� ������� ����: " << endl;
   cout << "  1) ��� ��� �������������� ������ (��� ������������������)" << endl;
   cout << "  2) ��� ��� ������������� ������ (������������ ������������������)" << endl;
   cout << "  3) ��� ��� ������������� ������ (�������� LU(sq)-������������������)" << endl;
   cout << "  4) ��� (��� ������������������)" << endl;
   cout << "  5) ��� (������������ ������������������)" << endl;
   cout << "  6) ��� (�������� LU(sq)-������������������)" << endl;

   int userCase;
   cin >> userCase;
   switch (userCase)
   {
      case 1:
      {
         cout << "������ ���������� ��� ������ ��� ��� �������������� ������ (��� ������������������)" << endl << endl;
         double eps = 0;
         IterSolvers::MSG_Assimetric::Init_Default(mat.Size());
         Timer timer;
         size_t it = IterSolvers::MSG_Assimetric::Default(mat, f, x, eps);
         timer.elapsed();
         IterSolvers::Destruct();
         cout << "����� �������� ������ �� " << timer.elapsedValue * 1000 << " ��" << endl << endl;
         cout << "���������� ��������: " << it << endl;
         cout << "������������� �������: " << eps << endl;
         break;
      }
      case 2:
      {
         cout << "������ ���������� ��� ������ ��� ��� �������������� ������ (������������ ������������������)" << endl << endl;
         double eps = 0;
         IterSolvers::MSG_Assimetric::Init_DiagPrecond(mat.Size());
         Timer timer;
         size_t it = IterSolvers::MSG_Assimetric::DiagPrecond(mat, f, x, eps);
         timer.elapsed();
         IterSolvers::Destruct();
         cout << "����� �������� ������ �� " << timer.elapsedValue * 1000 << " ��" << endl << endl;
         cout << "���������� ��������: " << it << endl;
         cout << "������������� �������: " << eps << endl;
         break;
      }
      case 3:
      {
         cout << "������ ���������� ��� ������ ��� ��� �������������� ������ (�������� LU(sq)-������������������)" << endl << endl;
         double eps = 0;
         IterSolvers::MSG_Assimetric::Init_LuPrecond(mat.Size(), mat.ggl.size());
         Timer timer;
         size_t it = IterSolvers::MSG_Assimetric::LuPrecond(mat, f, x, eps);
         timer.elapsed();
         IterSolvers::Destruct();
         cout << "����� �������� ������ �� " << timer.elapsedValue * 1000 << " ��" << endl << endl;
         cout << "���������� ��������: " << it << endl;
         cout << "������������� �������: " << eps << endl;
         break;
      }
      case 4:
      {
         cout << "������ ���������� ��� ������ ��� (��� ������������������)" << endl << endl;
         double eps = 0;
         IterSolvers::LOS::Init_Default(mat.Size());
         Timer timer;
         size_t it = IterSolvers::LOS::Default(mat, f, x, eps);
         timer.elapsed();
         IterSolvers::Destruct();
         cout << "����� �������� ������ �� " << timer.elapsedValue * 1000 << " ��" << endl << endl;
         cout << "���������� ��������: " << it << endl;
         cout << "������������� �������: " << eps << endl;
         break;
      }
      case 5:
      {
         cout << "������ ���������� ��� ������ ��� (������������ ������������������)" << endl << endl;
         double eps = 0;
         IterSolvers::LOS::Init_DiagPrecond(mat.Size());
         Timer timer;
         size_t it = IterSolvers::LOS::DiagPrecond(mat, f, x, eps);
         timer.elapsed();
         IterSolvers::Destruct();
         cout << "����� �������� ������ �� " << timer.elapsedValue * 1000 << " ��" << endl << endl;
         cout << "���������� ��������: " << it << endl;
         cout << "������������� �������: " << eps << endl;
         break;
      }
      case 6:
      {
         cout << "������ ���������� ��� ������ ��� (�������� LU(sq)-������������������)" << endl << endl;
         double eps = 0;
         IterSolvers::LOS::Init_LuPrecond(mat.Size(), mat.ggl.size());
         Timer timer;
         size_t it = IterSolvers::LOS::LuPrecond(mat, f, x, eps);
         timer.elapsed();
         IterSolvers::Destruct();
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