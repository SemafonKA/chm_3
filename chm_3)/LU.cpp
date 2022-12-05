#include "LU.h"

/// <summary>
/// ����������� � ��������������� ������ ��� ����������
/// </summary>
/// <param name="diSize"> - ������ ���������,</param>
/// <param name="luSize"> - ������ �������� ������� � �������� ������������</param>
LU::LU(size_t diSize, size_t luSize) {
   di.resize(diSize);
   ggl.resize(luSize);
   ggu.resize(luSize);
}

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
