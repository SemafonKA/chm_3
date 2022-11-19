#include "Matrix.h"

using namespace std;

vector<double> ReadVecFromFile(size_t size, const string& path) {
   vector<double> vec(size);
   auto file = ifstream(path);
   if (!file.is_open())
   {
      throw runtime_error("���� " + path + " ����������� � ����������");
   }
   for (size_t i = 0; i < size; i++)
   {
      file >> vec[i];
   }
   file.close();

   return vec;
}

   /// <summary>
   /// ������ �������� �����/��������, ���� 0, 0, 0 + k2, ..., 0+k2+...+kn, ��� ki - ����� ��������� � i c�����/�������
   /// <para> ������ ����� ������ ������� i ������ ����� ����� ��� ggl[ig[i]] </para>
   /// <para> ������ ������� ��� ������� 3�3: </para>
   /// 
   /// <para> �������: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> ig: { 0, 0, 1, 2 } </para>
   /// </summary>

// ������ �������
   uint16_t Matrix::Size() const { return di.size(); }

   /// <summary>
   /// ��������� ������� �� ������
   /// </summary>
   vector<double> Matrix::MultToVec(const vector<double>& right) const {
      if (right.size() != di.size()) throw runtime_error("������� ������� � ������� �� ���������.");

      vector<double> result(right.size());

      for (uint16_t i = 0; i < result.size(); i++)
      {
         // �������� ���������
         result[i] = di[i] * right[i];

         // �������� ������ � ������� ������������
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
   /// ��������� ����������������� ������� �� ������
   /// </summary>
   vector<double> Matrix::TranspMultToVec(const vector<double>& right) {
      if (right.size() != di.size()) throw runtime_error("������� ������� � ������� �� ���������.");

      vector<double> result(right.size());

      for (uint16_t i = 0; i < result.size(); i++)
      {
         // �������� ���������
         result[i] = di[i] * right[i];

         // �������� ������ � ������� ������������
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

// ������������ �������
   Matrix::Matrix() {}

   Matrix::Matrix(Matrix& right) :
      ig{ right.ig },
      jg{ right.jg },
      ggl{ right.ggl },
      ggu{ right.ggu },
      di{ right.di }
   {}

   // ����������� ����������� (����� ��� ������ ReadFromFiles)
   Matrix::Matrix(Matrix&& right) noexcept
   {
      ig = std::move(right.ig);
      jg = std::move(right.jg);
      ggl = std::move(right.ggl);
      ggu = std::move(right.ggu);
      di = std::move(right.di);
   }

// ����������� ������ �������
   Matrix Matrix::ReadFromFiles(uint16_t matrixSize, const string& igP, const string& jgP, const string& gglP, const string& gguP, const string& diP) {
      Matrix mat;

      {
         mat.ig.resize(matrixSize + 1);
         auto igS = ifstream(igP);
         if (!igS.is_open()) throw runtime_error("���� " + igP + " ����������� � ����������.");
         for (uint16_t i = 0; i <= matrixSize; i++)
         {
            igS >> mat.ig[i];
         }
         // ���� ������ ig � ����� ��������� � 1, �� ������ ��� ��� ���� ��������� (��� 0)
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
         if (!jgS.is_open()) throw runtime_error("���� " + jgP + " ����������� � ����������.");
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
