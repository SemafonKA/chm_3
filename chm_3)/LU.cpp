#include "LU.h"



/// <summary>
/// ����������� � ��������������� ������ ��� ����������
/// </summary>
/// <param name="diSize"> - ������ ���������,</param>
/// <param name="luSize"> - ������ �������� ������� � �������� ������������</param>
LU::LU(size_t diSize, size_t luSize) {
   Resize(diSize, luSize);
}


/// <summary>
/// ����������� � ����������� ��������� LU(sq)-���������� �� ������� mat
/// </summary>
/// <param name="mat"> - �������, �� ������� ���������� LU-����������, � ��������� ���� ������� � �������</param>
LU::LU(const Matrix& mat)
{
   MakeLuFor(mat);
}



/// <summary>
/// ��������� ������� mat � �������� LU(sq) - ����������
/// </summary>
/// <param name="mat"> - �������, ������� ��������� ���������. ��� �� ����� �������������� ��� ��������� �������� ������</param>
void LU::MakeLuFor(const Matrix& mat) {
   parent = &mat;
   if (di.size() != mat.di.size())
      di.resize(mat.di.size());
   if (ggl.size() != mat.ggl.size())
      ggl.resize(mat.ggl.size());
   if (ggu.size() != mat.ggu.size())
      ggu.resize(mat.ggu.size());

   const auto& ig = mat.ig;
   const auto& jg = mat.jg;

   di[0] = sqrt(mat.di[0]);
   for (size_t i = 1; i < mat.Size(); i++)
   {
      di[i] = 0;
      for (size_t j = ig[i]; j < ig[i + 1]; j++)
      {
         size_t k = ig[i];
         size_t v = ig[jg[j]];
         ggl[j] = ggu[j] = 0;
         while (k < j && v < ig[jg[j] + 1])
         {
            if (jg[k] > jg[v]) v++;
            else if (jg[k] < jg[v]) k++;
            else
            {
               ggl[j] += ggl[k] * ggu[v];
               ggu[j] += ggl[v] * ggu[k];
               k++;
               v++;
            }
         }
         ggl[j] = (mat.ggl[j] - ggl[j]) / di[jg[j]];
         ggu[j] = (mat.ggu[j] - ggu[j]) / di[jg[j]];

         di[i] += ggl[j] * ggu[j];
      }

      di[i] = sqrt(mat.di[i] - di[i]);
   }
}

void LU::Resize(size_t diSize, size_t luSize) {
   di.resize(diSize);
   ggl.resize(luSize);
   ggu.resize(luSize);
}



/// <summary>
/// ��������� ������ ������� L �� ������ vec. �������� ������ ��� ������ ������, �� ������ ������� LU
/// </summary>
/// <param name="vec"> - ������, �� ������� ����� ����������� ��������� ������� L</param>
/// <returns>������ � ����������� ������������ (���������� � ������)</returns>
std::vector<double> LU::LMultToVec(const std::vector<double>& vec) const
{
   std::vector<double> ans(vec.size());
   return LMultToVec(vec, ans);
}

/// <summary>
/// ��������� ������ ������� L �� ������ vec. ����� ������������ � ������ ans, �� ������ ������� LU
/// </summary>
/// <param name="vec"> - ������, �� ������� ����� ����������� ��������� ������� L;</param>
/// <param name="ans"> - ������, ���� ��������� ����� ��� ��������� ������ (������ ���������� �� vec!)</param>
/// <returns>������ �� ������ ans</returns>
std::vector<double>& LU::LMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const
{
   if (vec.size() != di.size()) throw std::runtime_error("������� ������� � ������� �� ���������.");
   if (vec.size() != ans.size()) throw std::runtime_error("������� ������� � ��������������� ������� �� ���������.");

   for (uint16_t i = 0; i < ans.size(); i++)
   {
      // �������� ���������
      ans[i] = di[i] * vec[i];

      // �������� ������ �����������
      for (uint32_t j = parent->ig[i]; j < parent->ig[i + 1]; j++)
      {
         ans[i] += ggl[j] * vec[parent->jg[j]];
      }
   }

   return ans;
}

std::vector<double> LU::LTranspMultToVec(const std::vector<double>& vec) const
{
   std::vector<double> ans(vec.size());
   return LTranspMultToVec(vec, ans);
}

std::vector<double>& LU::LTranspMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const
{
   if (vec.size() != di.size()) throw std::runtime_error("������� ������� � ������� �� ���������.");
   if (vec.size() != ans.size()) throw std::runtime_error("������� ������� � ��������������� ������� �� ���������.");

   for (uint16_t i = 0; i < ans.size(); i++)
   {
      // �������� ���������
      ans[i] = di[i] * vec[i];

      // �������� �� ������� ����������� � ������� �������
      for (uint32_t j = parent->ig[i]; j < parent->ig[i + 1]; j++)
      {
         ans[parent->jg[j]] += ggl[j] * vec[i];
      }
   }

   return ans;
}

/// <summary>
/// ��������� ������� ������� U �� ������ vec. �������� ������ ��� ������ ������, �� ������ ������� LU
/// </summary>
/// <param name="vec"> - ������, �� ������� ����� ����������� ��������� ������� U</param>
/// <returns>������ � ����������� ������������ (���������� � ������)</returns>
std::vector<double> LU::UMultToVec(const std::vector<double>& vec) const
{
   std::vector<double> ans(vec.size());
   return UMultToVec(vec, ans);
}

/// <summary>
/// ��������� ������� ������� U �� ������ vec. ����� ������������ � ������ ans, �� ������ ������� LU
/// </summary>
/// <param name="vec"> - ������, �� ������� ����� ����������� ��������� ������� U;</param>
/// <param name="ans"> - ������, ���� ��������� ����� ��� ��������� ������ (������ ���������� �� vec!)</param>
/// <returns>������ �� ������ ans</returns>
std::vector<double>& LU::UMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const
{
   if (vec.size() != di.size()) throw std::runtime_error("������� ������� � ������� �� ���������.");
   if (vec.size() != ans.size()) throw std::runtime_error("������� ������� � ��������������� ������� �� ���������.");

   for (uint16_t i = 0; i < ans.size(); i++)
   {
      // �������� ���������
      ans[i] = di[i] * vec[i];

      // �������� ������� �����������
      for (uint32_t j = parent->ig[i]; j < parent->ig[i + 1]; j++)
      {
         ans[parent->jg[j]] += ggu[j] * vec[i];
      }
   }

   return ans;
}

std::vector<double> LU::UTranspMultToVec(const std::vector<double>& vec) const
{
   std::vector<double> ans(vec.size());
   return UTranspMultToVec(vec, ans);
}

std::vector<double>& LU::UTranspMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const
{
   if (vec.size() != di.size()) throw std::runtime_error("������� ������� � ������� �� ���������.");
   if (vec.size() != ans.size()) throw std::runtime_error("������� ������� � ��������������� ������� �� ���������.");

   for (uint16_t i = 0; i < ans.size(); i++)
   {
      // �������� ���������
      ans[i] = di[i] * vec[i];

      // �������� ������ ����������� � ������� �������� ������������
      for (uint32_t j = parent->ig[i]; j < parent->ig[i + 1]; j++)
      {
         ans[i] += ggu[j] * vec[parent->jg[j]];
      }
   }

   return ans;
}



/// <summary>
/// ������� ���� ���� Lx = right. �� �������� ������ ��� ������ x, �� ������ ������� LU
/// </summary>
/// <param name="right"> - ������ ������ ����� ���������;</param>
/// <param name="x"> - ������, ���� ����� ������� �����. ������ ���� � ��� ���������� �������. ������ ���������� �� right!</param>
/// <returns>������ �� ������ x</returns>
std::vector<double>& LU::LSlauSolve(const std::vector<double>& right, std::vector<double>& x) const
{
   if (right.size() != di.size()) throw std::runtime_error("������� ������� � ������� �� ���������.");
   if (right.size() != x.size()) throw std::runtime_error("������� ������� � ��������������� ������� �� ���������.");

   size_t size = x.size();
   for (size_t i = 0; i < size; i++)
   {
      x[i] = 0;
      for (size_t j = parent->ig[i]; j < parent->ig[i + 1]; j++)
      {
         x[i] += x[parent->jg[j]] * ggl[j];
      }
      x[i] = (right[i] - x[i]) / di[i];
   }

   return x;
}

/// <summary>
/// ������� ���� ���� Lx = right. �������� ������ ��� ������ x, �� ������ ������� LU
/// </summary>
/// <param name="right"> - ������ ������ ����� ���������;</param>
/// <returns>���������� ������ x</returns>
std::vector<double> LU::LSlauSolve(const std::vector<double>& right) const
{
   std::vector<double> x(right.size());
   return LSlauSolve(right, x);
}

std::vector<double>& LU::LTranspSlauSolve(const std::vector<double>& right, std::vector<double>& x) const
{
   if (right.size() != di.size()) throw std::runtime_error("������� ������� � ������� �� ���������.");
   if (right.size() != x.size()) throw std::runtime_error("������� ������� � ��������������� ������� �� ���������.");

   size_t size = x.size();
   for (size_t i = 0; i < size; i++)
      x[i] = 0;

   for (size_t it = 0, i = size - 1; it < size; it++, i--)
   {
      x[i] = (right[i] - x[i]) / di[i];
      for (size_t j = parent->ig[i]; j < parent->ig[i + 1]; j++)
      {
         x[parent->jg[j]] += ggl[j] * x[i];
      }
   }

   return x;
}

std::vector<double> LU::LTranspSlauSolve(const std::vector<double>& right) const
{
   std::vector<double> x(right.size());
   return LTranspSlauSolve(right, x);
}


/// <summary>
/// ������� ���� ���� Ux = right. �� �������� ������ ��� ������ x, �� ������ ������� LU
/// </summary>
/// <param name="right"> - ������ ������ ����� ���������;</param>
/// <param name="x"> - ������, ���� ����� ������� �����. ������ ���� � ��� ���������� �������. ������ ���������� �� right!</param>
/// <returns>������ �� ������ x</returns>
std::vector<double>& LU::USlauSolve(const std::vector<double>& right, std::vector<double>& x) const
{
   if (right.size() != di.size()) throw std::runtime_error("������� ������� � ������� �� ���������.");
   if (right.size() != x.size()) throw std::runtime_error("������� ������� � ��������������� ������� �� ���������.");

   size_t size = x.size();
   for (size_t i = 0; i < size; i++) 
      x[i] = 0;

   for (size_t it = 0, i = size - 1; it < size; it++, i--)
   {
      x[i] = (right[i] - x[i]) / di[i];
      for (size_t j = parent->ig[i]; j < parent->ig[i + 1]; j++)
      {
         x[parent->jg[j]] += ggu[j] * x[i];
      }
   }

   return x;
}

/// <summary>
/// ������� ���� ���� Ux = right. �������� ������ ��� ������ x, �� ������ ������� LU
/// </summary>
/// <param name="right"> - ������ ������ ����� ���������;</param>
/// <returns>���������� ������ x</returns>
std::vector<double> LU::USlauSolve(const std::vector<double>& right) const
{
   std::vector<double> x(right.size());
   return USlauSolve(right, x);
}

std::vector<double>& LU::UTranspSlauSolve(const std::vector<double>& right, std::vector<double>& x) const
{
   if (right.size() != di.size()) throw std::runtime_error("������� ������� � ������� �� ���������.");
   if (right.size() != x.size()) throw std::runtime_error("������� ������� � ��������������� ������� �� ���������.");

   size_t size = x.size();
   for (size_t i = 0; i < size; i++)
   {
      x[i] = 0;
      for (size_t j = parent->ig[i]; j < parent->ig[i + 1]; j++)
      {
         x[i] += x[parent->jg[j]] * ggu[j];
      }
      x[i] = (right[i] - x[i]) / di[i];
   }

   return x;
}

std::vector<double> LU::UTranspSlauSolve(const std::vector<double>& right) const
{
   std::vector<double> x(right.size());
   return UTranspSlauSolve(right, x);
}
