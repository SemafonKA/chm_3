#include "Matrix.h"

using namespace std;

vector<double> ReadVecFromFile(size_t size, const string& path) {
   vector<double> vec(size);
   auto file = ifstream(path);
   if (!file.is_open())
   {
      throw runtime_error("Файл " + path + " отсутствует в директории");
   }
   for (size_t i = 0; i < size; i++)
   {
      file >> vec[i];
   }
   file.close();

   return vec;
}

   /// <summary>
   /// Массив индексов строк/столбцов, вида 0, 0, 0 + k2, ..., 0+k2+...+kn, где ki - число элементов в i cтроке/столбце
   /// <para> Помимо этого первый элемент i строки можно найти как ggl[ig[i]] </para>
   /// <para> Пример массива для матрицы 3х3: </para>
   /// 
   /// <para> Матрица: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> ig: { 0, 0, 1, 2 } </para>
   /// </summary>

// Методы матрицы
   uint16_t Matrix::Size() const { return di.size(); }

   /// <summary>
   /// Умножение матрицы на вектор
   /// </summary>
   vector<double> Matrix::MultToVec(const vector<double>& right) const {
      if (right.size() != di.size()) throw runtime_error("Размеры матрицы и вектора не совпадают.");

      vector<double> result(right.size());

      for (uint16_t i = 0; i < result.size(); i++)
      {
         // Умножаем диагональ
         result[i] = di[i] * right[i];

         // Умножаем нижний и верхний треугольники
         for (uint32_t j = ig[i]; j < ig[i + 1]; j++)
         {
            result[i] += ggl[j] * right[jg[j]];
            result[jg[j]] += ggu[j] * right[i];
         }
      }

      return result;
   }

   vector<double> Matrix::operator*(const vector<double>& right) const {
      return MultToVec(right);
   }

   /// <summary>
   /// Умножение транспонированной матрицы на вектор
   /// </summary>
   vector<double> Matrix::TranspMultToVec(const vector<double>& right) {
      if (right.size() != di.size()) throw runtime_error("Размеры матрицы и вектора не совпадают.");

      vector<double> result(right.size());

      for (uint16_t i = 0; i < result.size(); i++)
      {
         // Умножаем диагональ
         result[i] = di[i] * right[i];

         // Умножаем нижний и верхний треугольники
         for (uint32_t j = ig[i]; j < ig[i + 1]; j++)
         {
            result[i] += ggu[j] * right[jg[j]];
            result[jg[j]] += ggl[j] * right[i];
         }
      }

      return result;
   }

   Matrix& Matrix::operator= (Matrix&& right) noexcept {
      ig = std::move(right.ig);
      jg = std::move(right.jg);
      ggl = std::move(right.ggl);
      ggu = std::move(right.ggu);
      di = std::move(right.di);

      return *this;
   }

// Конструкторы матрицы
   Matrix::Matrix() {}

   Matrix::Matrix(Matrix& right) :
      ig{ right.ig },
      jg{ right.jg },
      ggl{ right.ggl },
      ggu{ right.ggu },
      di{ right.di }
   {}

   // Конструктор перемещения (нужен для метода ReadFromFiles)
   Matrix::Matrix(Matrix&& right) noexcept
   {
      ig = std::move(right.ig);
      jg = std::move(right.jg);
      ggl = std::move(right.ggl);
      ggu = std::move(right.ggu);
      di = std::move(right.di);
   }

// Статические методы матрицы
   Matrix Matrix::ReadFromFiles(uint16_t matrixSize, const string& igP, const string& jgP, const string& gglP, const string& gguP, const string& diP) {
      Matrix mat;

      {
         mat.ig.resize(matrixSize + 1);
         auto igS = ifstream(igP);
         if (!igS.is_open()) throw runtime_error("Файл " + igP + " отсутствует в директории.");
         for (uint16_t i = 0; i <= matrixSize; i++)
         {
            igS >> mat.ig[i];
         }
         // Если массив ig в файле начинался с 1, то меняем его под наши параметры (под 0)
         if (mat.ig[0] == 1)
         {
            for (uint16_t i = 0; i <= matrixSize; i++)
            {
               mat.ig[i]--;
            }
         }
      }

      {
         auto jgS = ifstream(jgP);
         if (!jgS.is_open()) throw runtime_error("Файл " + jgP + " отсутствует в директории.");
         uint16_t var;
         while (!jgS.eof())
         {
            jgS >> var;
            mat.jg.push_back(var);
         }
      }
      try
      {
         mat.di = ReadVecFromFile(matrixSize, diP);
         mat.ggl = ReadVecFromFile(mat.jg.size(), gglP);
         mat.ggu = ReadVecFromFile(mat.jg.size(), gguP);
      }
      catch (exception& e)
      {
         throw e;
      }

      return mat;
   }
