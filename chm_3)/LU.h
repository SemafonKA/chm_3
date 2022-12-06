#pragma once
#include <vector>
#include "Matrix.h"

// �������� ���������� LU(sq) ������� ������������ �������-����������� ������� Matrix
// �� ������ ������� �������, �� ���������� ������� �������� ������� (� ����� ��������� �� ��)
class LU {

// ���� ���������� ���������� ���������� LU
public:
   const Matrix* parent = nullptr;

   // ������ ������������ ��������� LU ����������. � ������ ������ ��������� L � U ���������
   std::vector<double> di;

   // ������ ��������� ������� ������������ L
   std::vector<double> ggl;

   // ������ ��������� �������� ������������ U
   std::vector<double> ggu;


// ���� �������� ������������� ������
public:

   /// <summary>
   /// ����������� � ��������������� ������ ��� ����������
   /// </summary>
   /// <param name="diSize"> - ������ ���������,</param>
   /// <param name="luSize"> - ������ �������� ������� � �������� ������������</param>
   LU(size_t diSize, size_t luSize);

   /// <summary>
   /// ����������� � ����������� ��������� LU(sq)-���������� �� ������� mat
   /// </summary>
   /// <param name="mat"> - �������, �� ������� ���������� LU-����������, � ��������� ���� ������� � �������</param>
   LU(const Matrix& mat);


// ���� �������� ������������� ������� ������
public:

   /// <summary>
   /// ��������� ������� mat � �������� LU(sq) - ����������
   /// </summary>
   /// <param name="mat"> - �������, ������� ��������� ���������. ��� �� ����� �������������� ��� ��������� �������� ������</param>
   void MakeLuFor(const Matrix& mat);

   /// <summary>
   /// ����� ��������� ������� ����������
   /// </summary>
   /// <param name="diSize"> - ������ ���������,</param>
   /// <param name="luSize"> - ������ �������� ������� � �������� ������������</param>
   void Resize(size_t diSize, size_t luSize);


// ��������� ������ �� ������

   /// <summary>
   /// ��������� ������ ������� L �� ������ vec. �������� ������ ��� ������ ������, �� ������ ������� LU
   /// </summary>
   /// <param name="vec"> - ������, �� ������� ����� ����������� ��������� ������� L</param>
   /// <returns>������ � ����������� ������������ (���������� � ������)</returns>
   std::vector<double> LMultToVec(const std::vector<double>& vec) const;

   /// <summary>
   /// ��������� ������ ������� L �� ������ vec. ����� ������������ � ������ ans, �� ������ ������� LU
   /// </summary>
   /// <param name="vec"> - ������, �� ������� ����� ����������� ��������� ������� L;</param>
   /// <param name="ans"> - ������, ���� ��������� ����� ��� ��������� ������ (������ ���������� �� vec!)</param>
   /// <returns>������ �� ������ ans</returns>
   std::vector<double>& LMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const;


   /// <summary>
   /// ��������� ������ ������� L^T �� ������ vec. �������� ������ ��� ������ ������, �� ������ ������� LU
   /// </summary>
   /// <param name="vec"> - ������, �� ������� ����� ����������� ��������� ������� L</param>
   /// <returns>������ � ����������� ������������ (���������� � ������)</returns>
   std::vector<double> LTranspMultToVec(const std::vector<double>& vec) const;

   /// <summary>
   /// ��������� ������ ������� L^T �� ������ vec. ����� ������������ � ������ ans, �� ������ ������� LU
   /// </summary>
   /// <param name="vec"> - ������, �� ������� ����� ����������� ��������� ������� L;</param>
   /// <param name="ans"> - ������, ���� ��������� ����� ��� ��������� ������ (������ ���������� �� vec!)</param>
   /// <returns>������ �� ������ ans</returns>
   std::vector<double>& LTranspMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const;


   /// <summary>
   /// ��������� ������� ������� U �� ������ vec. �������� ������ ��� ������ ������, �� ������ ������� LU
   /// </summary>
   /// <param name="vec"> - ������, �� ������� ����� ����������� ��������� ������� U</param>
   /// <returns>������ � ����������� ������������ (���������� � ������)</returns>
   std::vector<double> UMultToVec(const std::vector<double>& vec) const;

   /// <summary>
   /// ��������� ������� ������� U �� ������ vec. ����� ������������ � ������ ans, �� ������ ������� LU
   /// </summary>
   /// <param name="vec"> - ������, �� ������� ����� ����������� ��������� ������� U;</param>
   /// <param name="ans"> - ������, ���� ��������� ����� ��� ��������� ������ (������ ���������� �� vec!)</param>
   /// <returns>������ �� ������ ans</returns>
   std::vector<double>& UMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const;


   /// <summary>
   /// ��������� ������� ������� U^T �� ������ vec. �������� ������ ��� ������ ������, �� ������ ������� LU
   /// </summary>
   /// <param name="vec"> - ������, �� ������� ����� ����������� ��������� ������� U</param>
   /// <returns>������ � ����������� ������������ (���������� � ������)</returns>
   std::vector<double> UTranspMultToVec(const std::vector<double>& vec) const;

   /// <summary>
   /// ��������� ������� ������� U^T �� ������ vec. ����� ������������ � ������ ans, �� ������ ������� LU
   /// </summary>
   /// <param name="vec"> - ������, �� ������� ����� ����������� ��������� ������� U;</param>
   /// <param name="ans"> - ������, ���� ��������� ����� ��� ��������� ������ (������ ���������� �� vec!)</param>
   /// <returns>������ �� ������ ans</returns>
   std::vector<double>& UTranspMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const;


// ������� ���� � �������������� ������ � ������� ������ �����

   /// <summary>
   /// ������� ���� ���� Lx = right. �� �������� ������ ��� ������ x, �� ������ ������� LU
   /// </summary>
   /// <param name="right"> - ������ ������ ����� ���������;</param>
   /// <param name="x"> - ������, ���� ����� ������� �����. ������ ���� � ��� ���������� �������. ������ ���������� �� right!</param>
   /// <returns>������ �� ������ x</returns>
   std::vector<double>& LSlauSolve(const std::vector<double>& right, std::vector<double>& x) const;

   /// <summary>
   /// ������� ���� ���� Lx = right. �������� ������ ��� ������ x, �� ������ ������� LU
   /// </summary>
   /// <param name="right"> - ������ ������ ����� ���������;</param>
   /// <returns>���������� ������ x</returns>
   std::vector<double> LSlauSolve(const std::vector<double>& right) const;


   /// <summary>
   /// ������� ���� ���� L^T * x = right. �� �������� ������ ��� ������ x, �� ������ ������� LU
   /// </summary>
   /// <param name="right"> - ������ ������ ����� ���������;</param>
   /// <param name="x"> - ������, ���� ����� ������� �����. ������ ���� � ��� ���������� �������. ������ ���������� �� right!</param>
   /// <returns>������ �� ������ x</returns>
   std::vector<double>& LTranspSlauSolve(const std::vector<double>& right, std::vector<double>& x) const;

   /// <summary>
   /// ������� ���� ���� L^T * x = right. �������� ������ ��� ������ x, �� ������ ������� LU
   /// </summary>
   /// <param name="right"> - ������ ������ ����� ���������;</param>
   /// <returns>���������� ������ x</returns>
   std::vector<double> LTranspSlauSolve(const std::vector<double>& right) const;


   /// <summary>
   /// ������� ���� ���� Ux = right. �� �������� ������ ��� ������ x, �� ������ ������� LU
   /// </summary>
   /// <param name="right"> - ������ ������ ����� ���������;</param>
   /// <param name="x"> - ������, ���� ����� ������� �����. ������ ���� � ��� ���������� �������. ������ ���������� �� right!</param>
   /// <returns>������ �� ������ x</returns>
   std::vector<double>& USlauSolve(const std::vector<double>& right, std::vector<double>& x) const;

   /// <summary>
   /// ������� ���� ���� Ux = right. �������� ������ ��� ������ x, �� ������ ������� LU
   /// </summary>
   /// <param name="right"> - ������ ������ ����� ���������;</param>
   /// <returns>���������� ������ x</returns>
   std::vector<double> USlauSolve(const std::vector<double>& right) const;


   /// <summary>
   /// ������� ���� ���� U^T * x = right. �� �������� ������ ��� ������ x, �� ������ ������� LU
   /// </summary>
   /// <param name="right"> - ������ ������ ����� ���������;</param>
   /// <param name="x"> - ������, ���� ����� ������� �����. ������ ���� � ��� ���������� �������. ������ ���������� �� right!</param>
   /// <returns>������ �� ������ x</returns>
   std::vector<double>& UTranspSlauSolve(const std::vector<double>& right, std::vector<double>& x) const;

   /// <summary>
   /// ������� ���� ���� U^T * x = right. �������� ������ ��� ������ x, �� ������ ������� LU
   /// </summary>
   /// <param name="right"> - ������ ������ ����� ���������;</param>
   /// <returns>���������� ������ x</returns>
   std::vector<double> UTranspSlauSolve(const std::vector<double>& right) const;
};