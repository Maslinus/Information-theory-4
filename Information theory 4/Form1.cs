using System;
using System.Collections.Generic;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO;
using MathNet.Numerics.LinearAlgebra;
using System.Globalization;
using System.Windows.Forms.DataVisualization.Charting;

namespace Information_theory_4
{
    public partial class Form1 : Form
    {
        List<Vector<double>> x = new List<Vector<double>>();
        List<double> y = new List<double>();
        List<double> max = new List<double>();
        int discharge = 7;

        public Form1()
        {
            InitializeComponent();
            CultureInfo.CurrentCulture = CultureInfo.GetCultureInfo("en-US");

            string textpath = Application.StartupPath + "\\БД Титаник.csv";
            string[] lines = File.ReadAllLines(textpath, Encoding.GetEncoding(1251));
            for (int i = 3; i < lines.Length; i++)
            {
                string[] numbers = lines[i].Split(new char[] { ',' });
                double[] doubles = new double[numbers.Length];
                doubles[0] = 1;
                for (int j = 1; j < discharge; j++)
                {
                    doubles[j] = Convert.ToDouble(numbers[j]);
                }
                Vector<double> vector = Vector<double>.Build.DenseOfArray(new double[7] { doubles[0], doubles[1], doubles[2], doubles[3], doubles[4], doubles[5], doubles[6] });
                x.Add(vector);
                y.Add(Convert.ToDouble(numbers[0]));
            }
            // из листа в масив
            double[][] data = new double[x.Count][];
            for (int i = 0; i < x.Count; i++)
            {
                data[i] = new double[7];
            }
            for (int i = 0; i < x.Count; i++)
            {
                for (int j = 0; j < discharge; j++)
                {
                    data[i][j] = x[i][j];
                }
            }
/*            PrintMatrix(data);
*/
            Vector<double> mean = CalculateMean(data);
            Vector<double> stdDev = CalculateStdDev(data, mean);

            // 2. Применение Z-нормализации
            double[][] normalizedData = ApplyZNormalization(data, mean, stdDev);

//            Вывод результатов
/*            Console.WriteLine("Исходные данные:");
            PrintMatrix(data);*/

/*            Console.WriteLine("\nСреднее значение:");
            PrintVector(mean);*/

/*            Console.WriteLine("\nСтандартное отклонение:");
            PrintVector(stdDev);*/

/*            Console.WriteLine("\nНормализованные данные:");
            PrintMatrix(normalizedData);*/

            for (int i = 0; i < x.Count; i++)
            {
                for (int j = 0; j < discharge; j++)
                {
                    if (double.IsNaN(normalizedData[i][j]))
                    {
                        x[i][j] = 1;
                    }
                    else
                    {
                        x[i][j] = normalizedData[i][j];
                    }
                }
            }
/*            foreach (var v in x)
            {
                PrintVector(v);
            }*/

            double lambda = 0.01;
            Vector<double> x_k = Vector<double>.Build.DenseOfArray(new double[7] { 1, -1, 0, 0, 0, -1, -1 });
            Vector<double> x_k_1 = Vector<double>.Build.DenseOfArray(new double[7] { 0, 0, 0, 0, 0, 0, 0 }); ;
            // начальные значения
            double e = 0.1; // точность
            int n = 0; // количетсво итераций
            Vector<double> E = Vector<double>.Build.DenseOfArray(new double[7] { e, e, e, e, e, e, e });
            List<double> res = new List<double>();
            int max_count = 30;
            Vector<double> a = Vector<double>.Build.Dense(discharge);

            Console.WriteLine("Метод градиента с постоянным шагом = " + lambda);
            //while (IsDifferenceBelowThreshold(x_k_1 - x_k, E))
            while (n < max_count)
            {
                x_k = x_k_1;
                x_k_1 = x_k + lambda * GradL(x_k); // находим новое
                n++; // наращиваем итерации
                res.Add(L(x_k_1));
                if (n == 3) a = x_k_1;
/*                Console.WriteLine(n + " " + x_k + " " + x_k_1);
                Console.WriteLine(L(x_k_1));*/
                max.Add(L(x_k_1));
            }

            Console.WriteLine("Минимум: (" + res.Min() + ")");
            Console.WriteLine("Точка минимума: (" + x_k_1 + ")");
            Console.WriteLine("Минимальное значение функции: " + L(x_k_1));
            Console.WriteLine("Количество итераций: " + n);

            Console.WriteLine(x_k_1[0] - x_k_1[1] * x[0][1] - x_k_1[2] * x[0][2] - x_k_1[3] * x[0][3] - x_k_1[4] * x[0][4] - x_k_1[5] * x[0][5] - x_k_1[6] * x[0][6]);
            Console.WriteLine(x_k_1[0] - x_k_1[1] * x[1][1] - x_k_1[2] * x[1][2] - x_k_1[3] * x[1][3] - x_k_1[4] * x[1][4] - x_k_1[5] * x[1][5] - x_k_1[6] * x[1][6]);
            Console.WriteLine(a.DotProduct(x[0]));
            Console.WriteLine(a.DotProduct(x[1]));

            chart1.Series.Add(new Series());
            chart1.Series[0].ChartType = SeriesChartType.Line;
            chart1.Series[0].Color = Color.Red;
            chart1.Series[0].IsVisibleInLegend = false;

            for (int i = 0; i < max_count; i++)
            {
                chart1.Series[0].Points.AddXY(i, max[i]);
            }            
        }

        public double L(Vector<double> theta)
        {
            double res = new double();
            for (int i = 0; i < x.Count; i++)
            {
                double s1 = Math.Pow(1 / (1 + Math.Pow(Math.E, -theta.DotProduct(x[i]))), y[i]);
                double s2 = Math.Pow(1 / (1 + Math.Pow(Math.E, theta.DotProduct(x[i]))), 1 - y[i]);
                res += Math.Log(s1 * s2);
            }
            return res;
        }

        public Vector<double> GradL(Vector<double> theta)
        {
            Vector<double> res = Vector<double>.Build.Dense(discharge);
            for (int i = 0; i < x.Count; i++)
            {
                res += (y[i] - (1 / (1 + Math.Pow(Math.E, -theta.DotProduct(x[i]))))) * x[i];
            }
            return res;
        }

        static Vector<double> CalculateMean(double[][] data)
        {
            var matrix = Matrix<double>.Build.DenseOfRowArrays(data);
            var columnSums = matrix.ColumnSums();
            return columnSums / matrix.RowCount;
        }

        static Vector<double> CalculateStdDev(double[][] data, Vector<double> mean)
        {
            var matrix = Matrix<double>.Build.DenseOfRowArrays(data);

            // Вычисление суммы квадратов разностей с учетом длины векторов
            var squaredDifferences = matrix.EnumerateRows()
                .Select(row => (row - mean).PointwisePower(2))
                .ToArray();

            // Построение матрицы суммы квадратов
            var sumSquaredDiff = Matrix<double>.Build.DenseOfRowArrays(squaredDifferences.Select(row => row.ToArray()).ToArray());

            // Вычисление стандартного отклонения
            var stdDev = sumSquaredDiff.ColumnSums().PointwisePower(0.5) / Math.Sqrt(matrix.RowCount);

            return stdDev;
        }

        static double[][] ApplyZNormalization(double[][] data, Vector<double> mean, Vector<double> stdDev)
        {
            var matrix = Matrix<double>.Build.DenseOfRowArrays(data);
            var normalizedMatrix = matrix.EnumerateRows().Select(row => (row - mean).PointwiseDivide(stdDev)).ToArray();

            return normalizedMatrix.Select(v => v.ToArray()).ToArray();
        }
    }
}
