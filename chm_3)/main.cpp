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

   cout << "*** Программа для вычисления СЛАУ трёхшаговыми методами ***" << endl;
   cout << "Начинается считывание данных из файла..." << endl;

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

   cout << "Все данные успешно считанны из файлов." << endl << endl;
   cout << "Выберите метод для решения СЛАУ: " << endl;
   cout << "  1) МСГ для несимметричных матриц (без предобуславливания)" << endl;
   cout << "  2) МСГ для несиметричных матриц (диагональное предобуславливание)" << endl;
   cout << "  3) МСГ для несиметричных матриц (неполное LU(sq)-предобуславливание)" << endl;
   cout << "  4) ЛОС (без предобуславливания)" << endl;
   cout << "  5) ЛОС (диагональное предобуславливание)" << endl;
   cout << "  6) ЛОС (неполное LU(sq)-предобуславливание)" << endl;

   int userCase;
   cin >> userCase;
   switch (userCase)
   {
      case 1:
      {
         cout << "Начало вычислений для метода МСГ для несимметричных матриц (без предобуславливания)" << endl << endl;
         double eps = 0;
         IterSolvers::MSG_Assimetric::Init_Default(mat.Size());
         Timer timer;
         size_t it = IterSolvers::MSG_Assimetric::Default(mat, f, x, eps);
         timer.elapsed();
         IterSolvers::Destruct();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         cout << "Количество итераций: " << it << endl;
         cout << "Относительная невязка: " << eps << endl;
         break;
      }
      case 2:
      {
         cout << "Начало вычислений для метода МСГ для несимметричных матриц (диагональное предобуславливание)" << endl << endl;
         double eps = 0;
         IterSolvers::MSG_Assimetric::Init_DiagPrecond(mat.Size());
         Timer timer;
         size_t it = IterSolvers::MSG_Assimetric::DiagPrecond(mat, f, x, eps);
         timer.elapsed();
         IterSolvers::Destruct();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         cout << "Количество итераций: " << it << endl;
         cout << "Относительная невязка: " << eps << endl;
         break;
      }
      case 3:
      {
         cout << "Начало вычислений для метода МСГ для несимметричных матриц (неполное LU(sq)-предобуславливание)" << endl << endl;
         double eps = 0;
         IterSolvers::MSG_Assimetric::Init_LuPrecond(mat.Size(), mat);
         Timer timer;
         size_t it = IterSolvers::MSG_Assimetric::LuPrecond(mat, f, x, eps);
         timer.elapsed();
         IterSolvers::Destruct();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         cout << "Количество итераций: " << it << endl;
         cout << "Относительная невязка: " << eps << endl;
         break;
      }
      case 4:
      {
         cout << "Начало вычислений для метода ЛОС (без предобуславливания)" << endl << endl;
         double eps = 0;
         IterSolvers::LOS::Init_Default(mat.Size());
         Timer timer;
         size_t it = IterSolvers::LOS::Default(mat, f, x, eps);
         timer.elapsed();
         IterSolvers::Destruct();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         cout << "Количество итераций: " << it << endl;
         cout << "Относительная невязка: " << eps << endl;
         break;
      }
      case 5:
      {
         cout << "Начало вычислений для метода ЛОС (диагональное предобуславливание)" << endl << endl;
         double eps = 0;
         IterSolvers::LOS::Init_DiagPrecond(mat.Size());
         Timer timer;
         size_t it = IterSolvers::LOS::DiagPrecond(mat, f, x, eps);
         timer.elapsed();
         IterSolvers::Destruct();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         cout << "Количество итераций: " << it << endl;
         cout << "Относительная невязка: " << eps << endl;
         break;
      }
      case 6:
      {
         cout << "Начало вычислений для метода ЛОС (неполное LU(sq)-предобуславливание)" << endl << endl;
         double eps = 0;
         IterSolvers::LOS::Init_LuPrecond(mat.Size(), mat);
         Timer timer;
         size_t it = IterSolvers::LOS::LuPrecond(mat, f, x, eps);
         timer.elapsed();
         IterSolvers::Destruct();
         cout << "Метод закончил работу за " << timer.elapsedValue * 1000 << " мс" << endl << endl;
         cout << "Количество итераций: " << it << endl;
         cout << "Относительная невязка: " << eps << endl;
         break;
      }
      default:
         break;
   }

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