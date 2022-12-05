#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <string>

std::vector<double> ReadVecFromFile(size_t size, const std::string& path);

/// <summary>
/// ����� �������� ������, ���������� � ����������� �������-���������� ����
/// <para> �������� �������� ��������� - double </para>
/// </summary>
class Matrix {
// ���������� �������
public:
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
   std::vector<uint32_t> ig;

   /// <summary>
   /// ������ �������� ��������/����� ��������� (������ ������ � ������������ ��������)
   /// <para> ������ ������� ��� ������� 3�3: </para>
   /// 
   /// <para> �������: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> jg: { 0, 1 } </para>
   /// </summary>
   std::vector<uint16_t> jg;

   /// <summary>
   /// ������ ��������� ������� ������������ �������
   /// <para> ������ ������� ��� ������� 3�3: </para>
   /// 
   /// <para> �������: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> ggl: { 3, 2 } </para>
   /// </summary>
   std::vector<double> ggl;

   /// <summary>
   /// ������ ��������� �������� ������������ �������
   /// <para> ������ ������� ��� ������� 3�3: </para>
   /// 
   /// <para> �������: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> ggu: { 2, 1 } </para>
   /// </summary>
   std::vector<double> ggu;

   /// <summary>
   /// ������ ��������� ��������� �������
   /// <para> ������ ������� ��� ������� 3�3: </para>
   /// 
   /// <para> �������: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> di: { 1, 8, 4 } </para>
   /// </summary>
   std::vector<double> di;

// ������ �������
public:
   uint16_t Size() const;

   /// <summary>
   /// ��������� ������� �� ������
   /// </summary>
   std::vector<double> MultToVec(const std::vector<double>& right) const;
   std::vector<double>& MultToVec(const std::vector<double>& right, std::vector<double>& result) const;

   std::vector<double> operator*(const std::vector<double>& right) const;

   /// <summary>
   /// ��������� ����������������� ������� �� ������
   /// </summary>
   std::vector<double> TranspMultToVec(const std::vector<double>& right) const;
   std::vector<double>& TranspMultToVec(const std::vector<double>& right, std::vector<double>& result) const;

   Matrix& operator= (Matrix&& right) noexcept;

// ������������ �������
public:
   Matrix();

   Matrix(Matrix& right);

   // ����������� ����������� (����� ��� ������ ReadFromFiles)
   Matrix(Matrix&& right) noexcept;

// ����������� ������ �������
public:
   static Matrix ReadFromFiles(uint16_t matrixSize, const std::string& igP, const std::string& jgP, const std::string& gglP, const std::string& gguP, const std::string& diP);
};