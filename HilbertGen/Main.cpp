#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;
using matrix = vector<vector<double>>;

#define FILE_OPEN(file,path) auto file = ofstream(path); file.precision(15); file.setf(std::ios::scientific);

const string diagFilePath = "./iofiles/di.txt";
const string gglFilePath = "./iofiles/ggl.txt";
const string gguFilePath = "./iofiles/ggu.txt";
const string igFilePath = "./iofiles/ig.txt";
const string jgFilePath = "./iofiles/jg.txt";
const string sizeFilePath = "./iofiles/matrix_size.txt";
const string vectorbFilePath = "./iofiles/pr.txt";
const string vectorxFilePath = "./iofiles/x.txt";

void getHilbertMatrix(matrix& mas) {
   for (int i = 0; i < mas.size(); i++)
   {
      for (int j = 0; j < mas.size(); j++)
      {
         mas[i][j] = 1.0 / (i + j + 1);
      }
   }
}

void getVectorX(vector<double>& x) {
   for (int i = 0; i < x.size(); i++)
   {
      x[i] = i + 1;
   }
}

void getVectorB(matrix& mas, vector<double>& x, vector<double>& b) {
   for (int i = 0; i < b.size(); i++)
   {
      for (int j = 0; j < x.size(); j++)
      {
         b[i] += mas[i][j] * x[j];
      }
   }
}

void fprintVector(vector<double>& b, const string& path) {
   FILE_OPEN(vecFile, path);

   for (auto elem : b)
   {
      vecFile << elem << " ";
   }
   vecFile << endl;
   vecFile.close();
}

void matrixPrint(matrix& mat) {
   FILE_OPEN(diagFile, diagFilePath);
   FILE_OPEN(gglFile, gglFilePath);
   FILE_OPEN(gguFile, gguFilePath);
   FILE_OPEN(igFile, igFilePath);
   FILE_OPEN(jgFile, jgFilePath);

   igFile << 0 << " " << 0;
   size_t allCount = 0;
   for (size_t i = 0; i < mat.size(); i++)
   {
      diagFile << mat[i][i] << " ";
      int count = 0;
      for (size_t j = 0; j < i; j++)
      {
         if (mat[i][j] != 0 || mat[j][i] != 0)
         {
            count++;
            gglFile << mat[i][j] << " ";
            gguFile << mat[j][i] << " ";
            jgFile << j << " ";
         }
      }
      gglFile << endl;
      gguFile << endl;
      jgFile << endl;
      allCount += count;
      igFile << allCount << " ";
   }
}

int main(int argc, char** argv) {
   setlocale(LC_ALL, "ru-RU");
   size_t matrixSize;

   cout << "Программа генераций матриц Гильберта произвольного размера в формате строчно-столбцовом разреженном" << endl << endl;
   if (argc < 2)
   {
      cout << "** Введите размер генерируемой матрицы: ";
      cin >> matrixSize;
      cout << endl;
   }
   else
   {
      matrixSize = atoi(argv[1]);
   }

   cout << "** Размер генерируемой матрицы: " << matrixSize << endl;
   cout << "Начало генерации..." << endl << endl;

   matrix mat(matrixSize);
   vector<double> x(matrixSize);
   vector<double> b(matrixSize);

   for (auto& elem : mat)
   {
      elem.resize(matrixSize);
   }
   getHilbertMatrix(mat);
   getVectorX(x);
   getVectorB(mat, x, b);

   cout << "Генерация завершена. Начало вывода в файлы..." << endl << endl;
   fprintVector(x, vectorxFilePath);
   fprintVector(b, vectorbFilePath);
   matrixPrint(mat);
   auto sizeFile = ofstream(sizeFilePath);
   sizeFile << matrixSize << endl;
   sizeFile.close();
   
   cout << "Вывод завершен. Программа завершается" << endl;
}