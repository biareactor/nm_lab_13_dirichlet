#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "ui_mainwindow.h"

#include <vector>
#include <algorithm>

using vec = std::vector<double>;
using vecvec = std::vector<std::vector<double>>;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

bool MainWindow::idx_is_part_of_area(size_t i, size_t j)
{
    if (j <= m / 2)
        return true;
    else
        return i <= n / 2;
}

bool MainWindow::coord_is_part_of_area(double x, double y)
{
    if (y < (d-c) / 2)
        return true;
    else
        return x <= (b-a) / 2;
}

void MainWindow::init_v(vecvec& v, double h, double k)
{
    for (int j = 0; j <= m; j++)
    {
        double y = c + j * k;
        v[0][j] = mu1_test(y);
    }

    for (int i = 0; i <= n/2; i++)
    {
        double x = a + i * h;
        v[i][m] = mu2_test(x);
    }

    for (int j = m / 2; j <= m; j++)
    {
        double y = c + j * k;
        v[n/2][j] = mu3_test(y);
    }

    for (int i = n/2; i <= n; i++)
    {
        double x = a + i * h;
        v[i][m / 2] = mu4_test(x);
    }

    for (int j = 0; j <= m / 2; j++)
    {
        double y = c + j * k;
        v[n][j] = mu5_test(y);
    }

    for (int i = 0; i <= n; i++)
    {
        double x = a + i * h;
        v[i][0] = mu6_test(x);
    }
}

void MainWindow::init_f(vecvec& f, double h, double k)
{
    for (int j = 0; j <= m; j++)
    {
        for (int i = 0; i <= n; i++)
        {
            double x = a + i * h;
            double y = c + j * k;
            f[i][j] = f_test(x, y) * idx_is_part_of_area(i, j);
        }
    }
}

void MainWindow::fill_r(vecvec& r, const vecvec& v, const vecvec& f, double A, double h, double k, double h2, double k2)
{
    for (int j = 1; j < m; j++)
    {
        for (int i = 1; i < n; i++)
        {
            if (idx_is_part_of_area(i, j))
            {
                r[i][j] = A * v[i][j] + h2 * v[i + 1][j] + h2 * v[i - 1][j] + k2 * v[i][j + 1] + k2 * v[i][j - 1] - f[i][j];

                const double x = a + i * h;
                const double y = c + j * k;

                if (i == 1)       r[i][j] -= h2 * v[i-1][j] - mu1_test(y) * h2;
                if (j == m - 1)   r[i][j] -= k2 * v[i][j+1] - mu2_test(x) * k2;
                if (i == n/2 - 1) r[i][j] -= h2 * v[i+1][j] - mu3_test(y) * h2;
                if (j == m/2 - 1) r[i][j] -= k2 * v[i][j+1] - mu4_test(x) * k2;
                if (i == n - 1)   r[i][j] -= h2 * v[i+1][j] - mu5_test(y) * h2;
                if (j == 1)       r[i][j] -= k2 * v[i][j-1] - mu6_test(x) * k2;
            }
        }
    }
}

void MainWindow::fill_u_test(vecvec& u_test, size_t rows, size_t columns)
{
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < columns; j++)
        {
            const double h = (b - a) / n;
            const double k = (d - c) / m;
            double x = a + i * h;
            double y = c + j * k;
            u_test[i][j] = std::exp(1 - x * x - y * y);
        }
    }
}

double MainWindow::get_beta(double A, const vecvec& r, const vecvec& h, double h2, double k2)
{
    long double beta = 0, tmp_beta = 0;

    for (size_t j = 1; j < m; j++)
    {
        for (size_t i = 1; i < n; i++)
        {
            if (idx_is_part_of_area(i, j))
            {
                double Ah = A * h[i][j] + h2 * h[i+1][j] + h2 * h[i-1][j] + k2 * h[i][j+1] + k2 * h[i][j-1];

                if (i == 1)       Ah -= h2 * h[i-1][j];
                if (j == 1)       Ah -= k2 * h[i][j-1];
                if (i == n/2 - 1) Ah -= h2 * h[i+1][j];
                if (j == m/2 - 1) Ah -= k2 * h[i][j+1];
                if (i == n - 1)   Ah -= h2 * h[i+1][j];
                if (j == m - 1)   Ah -= k2 * h[i][j+1];

                beta += Ah * r[i][j];
                tmp_beta += Ah * h[i][j];
            }
        }
    }

    if (tmp_beta == 0) {
        beta = 0;
    }
    else beta /= tmp_beta;

    return beta;
}

double MainWindow::get_alpha(double A, const vecvec& r, const vecvec& h, double h2, double k2)
{
    long double alpha = 0,  tmp_alpha = 0;

    for (int j = 1; j < m; j++)
    {
        for (int i = 1; i < n; i++)
        {
            if (idx_is_part_of_area(i, j))
            {
                double Ah = A * h[i][j] + h2 * h[i + 1][j] + h2 * h[i - 1][j] + k2 * h[i][j + 1] + k2 * h[i][j - 1];

                if (i == 1)       Ah -= h2 * h[i-1][j];
                if (j == 1)       Ah -= k2 * h[i][j-1];
                if (i == n/2 - 1) Ah -= h2 * h[i+1][j];
                if (j == m/2 - 1) Ah -= k2 * h[i][j+1];
                if (i == n - 1)   Ah -= h2 * h[i+1][j];
                if (j == m - 1)   Ah -= k2 * h[i][j+1];

                alpha -= r[i][j] * h[i][j];
                tmp_alpha += Ah * h[i][j];
            }
        }
    }

    alpha /= tmp_alpha;

    return alpha;
}

double MainWindow::get_R(const vecvec& r)
{
    double R0 = 0;
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= m; j++)
        {
            if (idx_is_part_of_area(i, j))
                R0 += r[i][j] * r[i][j];
        }
    }

    return std::sqrt(R0);
}

void MainWindow::fill_table(QTableWidget* table, size_t columns, size_t rows, const vec& values)
{
    table->clear();
    table->setColumnCount(columns+1);
    table->setRowCount(rows+1);

    table->setItem(0, 0, new QTableWidgetItem("i/j"));

    for (size_t i = 0; i < rows; i++)
        table->setItem(i+1, 0, new QTableWidgetItem(QString::number(i)));

    for (size_t j = 0; j < columns; j++)
        table->setItem(0, j+1, new QTableWidgetItem(QString::number(j)));

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < columns; j++)
        {
            table->setItem(i+1, j+1, new QTableWidgetItem(QString::number(round(values[i][j] * 1e8) / 1e8)));
        }
    }
}

void MainWindow::on_pushButton_clicked()
{
    n = ui->n->text().toUInt();
    m = ui->m->text().toUInt();
    eps = ui->eps->text().toDouble();
    Nmax = ui->Nmax->text().toDouble();

    int s = 0;
    double eps_max = 0.0;
    double eps_cur = 0.0;
    double solve_err = 0.0;
    double solve_eps = 0.0;
    double A, k2, h2;
    double h = (b - a) / n;
    double k = (d - c) / m;

    double uv_max = 0.0;
    double x_max, y_max;

    std::vector<std::vector<double>> v(n + 1, std::vector<double>(m + 1, 0.0));
    std::vector<std::vector<double>> f(n + 1, std::vector<double>(m + 1, 0.0));
    std::vector<std::vector<double>> r(n + 1, std::vector<double>(m + 1, 0.0));
    std::vector<std::vector<double>> hs(n + 1, std::vector<double>(m + 1, 0.0));

    bool flag = false;

    h2 = (1 / (h * h));
    k2 = (1 / (k * k));
    A = -2 * (h2 + k2);

    init_v(v, h, k);
    init_f(f, h, k);
    fill_r(r, v, f, A, h, k, h2, k2);

    double R0 = get_R(r);

    const size_t columns = std::min<size_t>(100, n+1);
    const size_t rows = std::min<size_t>(100, m+1);

    fill_table(ui->v0_table, columns, rows, v);

    std::vector<std::vector<double>> v_old(n + 1, std::vector<double>(m + 1, 0.0));
    std::vector<std::vector<double>> h_old(n + 1, std::vector<double>(m + 1, 0.0));

    for(size_t steps = 0; !(eps_max <= eps || steps >= Nmax); steps++)
    {
        eps_max = 0.0;
        // (n, m)
        for (int j = 0; j <= m; j++)
        {
            int x_end = (j < m / 2) ? n : n / 2;

            for (int i = 0; i <= x_end; i++)
            {
                v_old[i][j] = v[i][j];
                h_old[i][j] = hs[i][j];
            }
        }

        fill_r(r, v_old, f, A, h, k, h2, k2);
        long double beta = get_beta(A, r, h_old, h2, k2);

        for (int j = 1; j < m; j++)
        {
            for (int i = 1; i < n; i++)
            {
                if (idx_is_part_of_area(i, j))
                    hs[i][j] = -r[i][j] + beta * h_old[i][j];
            }
        }

        long double alpha = get_alpha(A, r, hs, h2, k2);

        for (int j = 1; j < m; j++)
        {
            for (int i = 1; i < n; i++)
            {
                if (idx_is_part_of_area(i, j))
                {
                    v[i][j] = v_old[i][j] + alpha * hs[i][j];

                    eps_cur = std::abs(v_old[i][j] - v[i][j]);
                    if (eps_cur > eps_max) {
                        eps_max = eps_cur;
                    }
                }
            }
        }
    }

    double RN = get_R(r);
    std::vector<std::vector<double>> u_test(n + 1, std::vector<double>(m + 1, 0.0));
    fill_u_test(u_test, rows, columns);

    for (int j = 0; j <= m; j++)
    {
        for (int i = 0; i <= n; i++)
        {
            if (idx_is_part_of_area(i, j))
            {
                if (std::abs(u_test[i][j] - v[i][j]) > solve_err) {
                    solve_err = std::abs(u_test[i][j] - v[i][j]);
                }
            }
        }
    }

    for (int j = 0; j <= m; j++)
    {
        for (int i = 0; i <= n; i++)
        {
            if (idx_is_part_of_area(i, j))
            {
                if (std::abs(u_test[i][j] - v[i][j]) > uv_max) {
                    uv_max = std::abs(u_test[i][j] - v[i][j]);
                    x_max = a + i * h;
                    y_max = c + j * k;
                }
            }
        }
    }

    fill_table(ui->u_table, columns, rows, u_test);
    fill_table(ui->vn_table, columns, rows, v);

    std::vector<std::vector<double>> uv(n + 1, std::vector<double>(m + 1, 0.0));

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < columns; j++)
        {
            uv[i][j] = u_test[i][j]-v[i][j];
        }
    }

}

//System::Void buttonCalculateTest_Click(System::Object^ sender, System::EventArgs^ e) {



//   this->label13->Text = ("Макс. разность двух решений |U-V|\n достигается в точке (x, y):" + "   (" + Convert::ToString(round(x_max * 1e8) / 1e8) + ", " + Convert::ToString(round(y_max * 1e8) / 1e8) + ")");


//   dataGridView5->Rows->Clear();
//   dataGridView5->Columns->Clear();
//   this->dataGridView5->Columns->Add("", "j / i");
//   for (int j = 0; j <= std::min(newm, 100); j++) {
//       this->dataGridView5->Rows->Add(Convert::ToString(j));
//   }

//   for (int i = 0; i <= std::min(newn, 100); i++) {
//       this->dataGridView5->Columns->Add("", Convert::ToString(i));
//   }

//   for (int i = 0; i <= std::min(100, newn); i++)
//   {
//       int y_end = (i < std::min(100, newn) / 2) ? std::min(100, newm) : std::min(100, newm) / 2;

//       for (int j = 0; j <= y_end; j++)
//       {
//           this->dataGridView5->Rows[j]->Cells[i + 1]->Value = round(u_test(i, j, n, m) * 1e8) / 1e8;
//       }
//   }

//   dataGridView6->Rows->Clear();
//   dataGridView6->Columns->Clear();
//   this->dataGridView6->Columns->Add("", "j / i");
//   for (int j = 0; j <= std::min(newm, 100); j++) {
//       this->dataGridView6->Rows->Add(Convert::ToString(j));
//   }

//   for (int i = 0; i <= std::min(newn, 100); i++) {
//       this->dataGridView6->Columns->Add("", Convert::ToString(i));
//   }

//   for (int i = 0; i <= std::min(100, newn); i++)
//   {
//       int y_end = (i < std::min(100, newn) / 2) ? std::min(100, newm) : std::min(100, newm) / 2;

//       for (int j = 0; j <= y_end; j++)
//       {
//           this->dataGridView6->Rows[j]->Cells[i + 1]->Value = round(v[i][j] * 1e8) / 1e8;
//       }
//   }

//   dataGridView7->Rows->Clear();
//   dataGridView7->Columns->Clear();
//   this->dataGridView7->Columns->Add("", "j / i");
//   for (int j = 0; j <= std::min(newm, 100); j++) {
//       this->dataGridView7->Rows->Add(Convert::ToString(j));
//   }

//   for (int i = 0; i <= std::min(newn, 100); i++) {
//       this->dataGridView7->Columns->Add("", Convert::ToString(i));
//   }

//   for (int i = 0; i <= std::min(100, newn); i++)
//   {
//       int y_end = (i < std::min(100, newn) / 2) ? std::min(100, newm) : std::min(100, newm) / 2;

//       for (int j = 0; j <= y_end; j++)
//       {
//           this->dataGridView7->Rows[j]->Cells[i + 1]->Value = round(std::abs(u_test(i, j, n, m) - v[i][j]) * 1e8) / 1e8;
//       }
//   }

//   //if (radioButtonDrawU->Checked) {
//   chart1->ChartAreas[0]->RecalculateAxesScale();
//   double x1, y1, z1;
//   double xstep1 = (b - a) / n, ystep1 = (d - c) / m;
//   double max1 = u_test(a, c, n, m), min1 = max1;
//   Color Color2;

//   for (x1 = a; x1 <= b; x1 += xstep1)
//   {
//       int y_end = (x1 < b / 2) ? d : d / 2;

//       for (y1 = c; y1 <= y_end; y1 += ystep1)
//       {
//           int i = int((x1 - a) / xstep1);
//           int j = int((y1 - c) / ystep1);
//           z1 = u_test(i, j, n, m);
//           if (z1 > max1) max1 = z1;
//           else if (z1 < min1) min1 = z1;
//       }
//   }

//   for (int i = 0; i <= n; i++)
//   {
//       int y_end = (i < n / 2) ? m : m / 2;

//       Series^ Series1 = gcnew Series();
//       Series1->ChartType = SeriesChartType::Spline;
//       Series1->Color = Color::Blue;
//       chart1->Series->Add(Series1);

//       for (int j = 0; j <= y_end; j++)
//       {
//           x1 = a + i * xstep1;
//           y1 = c + j * ystep1;
//           z1 = u_test(i, j, n, m);
//           chart1->Series[i]->Points->AddXY((double)y1, z1);
//           Color2 = calculateColorPurple(min1, z1, max1);
//           chart1->Series[i]->Points[j]->Color = Color2;
//           chart1->Series[i]->Points[j]->BackSecondaryColor = Color2;
//       }
//   }

//   //}

//   //if (radioButtonDrawV->Checked) {
//   chart8->ChartAreas[0]->RecalculateAxesScale();
//   double x8, y8, z8;
//   double xstep8 = (b - a) / n, ystep8 = (d - c) / m;
//   double max8 = v[0][0], min8 = max8;
//   Color Color8;

//   for (x8 = a; x8 <= b; x8 += xstep8)
//   {
//       int y_end = (x8 < b / 2) ? d : d / 2;

//       for (y8 = c; y8 <= y_end; y8 += ystep8)
//       {
//           int i = int((x8 - a) / xstep8);
//           int j = int((y8 - c) / ystep8);
//           z8 = v[i][j];
//           if (z8 > max8) max8 = z8;
//           else if (z8 < min8) min8 = z8;
//       }
//   }

//   for (int i = 0; i <= n; i++)
//   {
//       int y_end = (i < n / 2) ? m : m / 2;

//       Series^ Series1 = gcnew Series();
//       Series1->ChartType = SeriesChartType::Spline;
//       Series1->Color = Color::Blue;
//       chart8->Series->Add(Series1);

//       for (int j = 0; j <= y_end; j++)
//       {
//           x8 = a + i * xstep8;
//           y8 = c + j * ystep8;
//           z8 = v[i][j];
//           chart8->Series[i]->Points->AddXY((double)y8, z8);
//           Color8 = calculateColorBlue(min8, z8, max8);
//           chart8->Series[i]->Points[j]->Color = Color8;
//           chart8->Series[i]->Points[j]->BackSecondaryColor = Color8;
//       }
//   }
//   //}

//   //if (radioButtonDrawUV->Checked) {
//   chart9->ChartAreas[0]->RecalculateAxesScale();
//   double x9, y9, z9;
//   double xstep9 = (b - a) / n, ystep9 = (d - c) / m;
//   double max9 = u_test(a, c, n, m) - v[0][0], min9 = max9;
//   Color Color9;

//   for (x9 = a; x9 <= b; x9 += xstep9)
//   {
//       int y_end = (x9 < b / 2) ? d : d / 2;

//       for (y9 = c; y9 <= y_end; y9 += ystep9)
//       {
//           int i = int((x9 - a) / xstep9);
//           int j = int((y9 - c) / ystep9);
//           z9 = u_test(i, j, n, m) - v[i][j];
//           if (z9 > max9) max9 = z9;
//           else if (z9 < min9) min9 = z9;
//       }
//   }

//   for (int i = 0; i <= n; i++)
//   {
//       int y_end = (i < n / 2) ? m : m / 2;

//       Series^ Series1 = gcnew Series();
//       Series1->ChartType = SeriesChartType::Spline;
//       Series1->Color = Color::Blue;
//       chart9->Series->Add(Series1);

//       for (int j = 0; j <= y_end; j++)
//       {
//           x9 = a + i * xstep9;
//           y9 = c + j * ystep9;
//           z9 = u_test(i, j, n, m) - v[i][j];
//           chart9->Series[i]->Points->AddXY((double)y9, z9);
//           Color9 = calculateColorRed(min9, z9, max9);
//           chart9->Series[i]->Points[j]->Color = Color9;
//           chart9->Series[i]->Points[j]->BackSecondaryColor = Color9;
//       }
//   }
//   //}
//}

//System::Void buttonCalculateMain_Click(System::Object^ sender, System::EventArgs^ e) {
//   n = Convert::ToInt32(textBox8->Text);
//   m = Convert::ToInt32(textBox7->Text);
//   eps = Convert::ToDouble(textBox6->Text);
//   Nmax = Convert::ToInt32(textBox5->Text);

//   int n_2 = 2 * n;
//   int m_2 = 2 * m;

//   int s = 0, s2 = 0;
//   double eps_max = 0.0, eps_max2 = 0.0;
//   double eps_cur = 0.0, eps_curr2 = 0.0;
//   double solve_eps = 0.0;
//   double A, k2, h2, A_2, k2_2, h2_2;
//   double h = (b - a) / n;
//   double k = (d - c) / m;
//   double h_2 = (b - a) / n_2;
//   double k_2 = (d - c) / m_2;

//   double v2v_max = 0.0;
//   double x_max, y_max;

//   std::vector<std::vector<double>> v(n + 1, std::vector<double>(m + 1, 0.0));
//   std::vector<std::vector<double>> f(n + 1, std::vector<double>(m + 1, 0.0));
//   std::vector<std::vector<double>> r(n + 1, std::vector<double>(m + 1, 0.0));

//   std::vector<std::vector<double>> v2(n_2 + 1, std::vector<double>(m_2 + 1, 0.0));
//   std::vector<std::vector<double>> f2(n_2 + 1, std::vector<double>(m_2 + 1, 0.0));
//   std::vector<std::vector<double>> r2(n_2 + 1, std::vector<double>(m_2 + 1, 0.0));


//   double R_N_2 = 0.0, R2_N_2 = 0.0;
//   double R_0_2 = 0.0, R2_0_2 = 0.0, R_0_inf = 0.0, R2_0_inf = 0.0;
//   bool flag = false, flag2 = false;

//   h2 = (1 / (h * h));
//   k2 = (1 / (k * k));
//   A = -2 * (h2 + k2);

//   h2_2 = (1 / (h_2 * h_2));
//   k2_2 = (1 / (k_2 * k_2));
//   A_2 = -2 * (h2_2 + k2_2);


//   if (radioButtonXMain->Checked) {
//       for (int j = 1; j < m; j++)
//           for (int i = 1; i < n; i++) {
//               v[i][j] = interpolationXTest(j, m);
//           }

//       for (int j = 1; j < m_2; j++)
//           for (int i = 1; i < n_2; i++) {
//               v2[i][j] = interpolationXTest(j, m_2);
//           }
//   }

//   if (radioButtonYMain->Checked) {
//       for (int j = 1; j < m; j++)
//           for (int i = 1; i < n; i++) {
//               v[i][j] = interpolationYTest(i, n);
//           }

//       for (int j = 1; j < m_2; j++)
//           for (int i = 1; i < n_2; i++) {
//               v2[i][j] = interpolationYTest(i, n_2);
//           }
//   }

//   // initialize u1(y), u2(y) for (n, m)
//   for (int j = 0; j <= m; j++) {
//       double y = c + j * k;
//       v[0][j] = mu1(y);
//       v[n][j] = mu2(y);
//   }

//   // initialize u3(x), u4(x) for (n, m)
//   for (int i = 0; i <= n; i++) {
//       double x = a + i * h;
//       v[i][0] = mu3(x);
//       v[i][m] = mu4(x);
//   }

//   // initialize f(x,y) for (n, m)
//   for (int j = 0; j <= m; j++)
//       for (int i = 0; i <= n; i++) {
//           double x = a + i * h;
//           double y = c + j * k;
//           f[i][j] = f_main(x, y);
//       }
//   // initialize u1(y), u2(y) for (2n, 2m)
//   for (int j = 0; j <= m_2; j++) {
//       double y = c + j * k_2;
//       v2[0][j] = mu1(y);
//       v2[n_2][j] = mu2(y);
//   }

//   // initialize u3(x), u4(x) for (2n, 2m)
//   for (int i = 0; i <= n_2; i++) {
//       double x = a + i * h_2;
//       v2[i][0] = mu3(x);
//       v2[i][m_2] = mu4(x);
//   }

//   // initialize f(x,y) for (2n, 2m)
//   for (int j = 0; j <= m_2; j++)
//       for (int i = 0; i <= n_2; i++) {
//           double x = a + i * h_2;
//           double y = c + j * k_2;
//           f2[i][j] = f_main(x, y);
//       }

//   for (int j = 1; j < m; j++)
//       for (int i = 1; i < n; i++) {
//           r[i][j] = A * v[i][j] + h2 * v[i + 1][j] + h2 * v[i - 1][j] + k2 * v[i][j + 1] + k2 * v[i][j - 1] - f[i][j];
//           if (i == 1) {
//               r[i][j] -= h2 * v[i - 1][j] - mu1(c + j * k) * h2;
//           }
//           if (j == 1) {
//               r[i][j] -= k2 * v[i][j - 1] - mu3(a + i * h) * k2;
//           }
//           if (i == n - 1) {
//               r[i][j] -= h2 * v[i + 1][j] - mu2(c + j * k) * h2;
//           }
//           if (j == m - 1) {
//               r[i][j] -= k2 * v[i][j + 1] - mu4(a + i * h) * k2;
//           }
//       }

//   for (int i = 0; i < n + 1; i++)
//       for (int j = 0; j < m + 1; j++) {
//           R_0_inf = std::max(R_0_inf, std::abs(r[i][j]));
//           R_0_2 += r[i][j] * r[i][j];
//       }
//   R_0_2 = std::sqrt(R_0_2);
//   // end of calculating R0 (n,m)

//   // calculate R0 (2n,2m)
//   for (int j = 1; j < m_2; j++)
//       for (int i = 1; i < n_2; i++) {
//           r2[i][j] = A_2 * v2[i][j] + h2_2 * v2[i + 1][j] + h2_2 * v2[i - 1][j] + k2_2 * v2[i][j + 1] + k2_2 * v2[i][j - 1] - f2[i][j];
//           if (i == 1) {
//               r2[i][j] -= h2_2 * v2[i - 1][j] - mu1(c + j * k_2) * h2_2;
//           }
//           if (j == 1) {
//               r2[i][j] -= k2_2 * v2[i][j - 1] - mu3(a + i * h_2) * k2_2;
//           }
//           if (i == n - 1) {
//               r2[i][j] -= h2_2 * v2[i + 1][j] - mu2(c + j * k_2) * h2_2;
//           }
//           if (j == m - 1) {
//               r2[i][j] -= k2_2 * v2[i][j + 1] - mu4(a + i * h_2) * k2_2;
//           }
//       }

//   for (int i = 0; i < n_2 + 1; i++)
//       for (int j = 0; j < m_2 + 1; j++) {
//           R2_0_inf = std::max(R2_0_inf, std::abs(r2[i][j]));
//           R2_0_2 += r2[i][j] * r2[i][j];
//       }
//   R2_0_2 = std::sqrt(R2_0_2);
//   // end of calculating R0 (2n,2m)

//   //this->chart2->Series->Clear();
//   //this->chart3->Series->Clear();
//   this->chart4->Series->Clear();
//   this->chart5->Series->Clear();
//   this->chart6->Series->Clear();

//   //if (radioButtonDrawV0Main->Checked) {
//   //chart2->ChartAreas[0]->RecalculateAxesScale();
//   //double x1, y1, z1;
//   //double xstep1 = (b - a) / n, ystep1 = (d - c) / m;
//   //double max1 = v[0][0], min1 = max1;
//   //Color Color1;
//   //for (x1 = a; x1 <= b; x1 += xstep1)
//      // for (y1 = c; y1 <= d; y1 += ystep1) {
//         //  int i = int((x1 - a) / xstep1);
//         //  int j = int((y1 - c) / ystep1);
//         //  z1 = v[i][j];
//         //  if (z1 > max1) max1 = z1;
//         //  else if (z1 < min1) min1 = z1;
//      // }
//   //for (int i = 0; i <= n; i++) {
//      // Series^ Series1 = gcnew Series();
//      // Series1->ChartType = SeriesChartType::Spline;
//      // Series1->Color = Color::Blue;
//      // //chart2->Series->Add(Series1);
//      // for (int j = 0; j <= m; j++) {
//         //  x1 = a + i * xstep1;
//         //  y1 = c + j * ystep1;
//         //  z1 = v[i][j];
//         //  //chart2->Series[i]->Points->AddXY((double)y1, z1);
//         //  Color1 = calculateColorBlue(min1, z1, max1);
//         //  chart2->Series[i]->Points[j]->Color = Color1;
//         //  chart2->Series[i]->Points[j]->BackSecondaryColor = Color1;
//      // }
//   //}
//   //}

//   //if (radioButtonDrawV20->Checked) {
//   //chart3->ChartAreas[0]->RecalculateAxesScale();
//   //double x2, y2, z2;
//  // double xstep2 = (b - a) / n_2, ystep2 = (d - c) / m_2;
//   //double max2 = v2[0][0], min2 = max2;
//   //Color Color2;
//   //for (x2 = a; x2 <= b; x2 += xstep2)
//      // for (y2 = c; y2 <= d; y2 += ystep2) {
//         //  int i = int((x2 - a) / xstep2);
//         //  int j = int((y2 - c) / ystep2);
//         //  z2 = v2[i][j];
//         //  if (z2 > max2) max2 = z2;
//         //  else if (z2 < min2) min2 = z2;
//      // }
//   //for (int i = 0; i <= n_2; i++) {
//      // Series^ Series1 = gcnew Series();
//      // Series1->ChartType = SeriesChartType::Spline;
//      // Series1->Color = Color::Blue;
//      // //chart3->Series->Add(Series1);
//      // for (int j = 0; j <= m_2; j++) {
//         //  x2 = a + i * xstep2;
//         //  y2 = c + j * ystep2;
//         //  z2 = v2[i][j];
//         //  //chart3->Series[i]->Points->AddXY((double)y2, z2);
//         //  Color2 = calculateColorPurple(min2, z2, max2);
//         //  //chart3->Series[i]->Points[j]->Color = Color2;
//         //  //chart3->Series[i]->Points[j]->BackSecondaryColor = Color2;
//      // }
//   //}
//   //}

//   std::vector<std::vector<double>> v_old(n + 1, std::vector<double>(m + 1, 0.0)),
//       v_old2(n_2 + 1, std::vector<double>(m_2 + 1, 0.0));
//   std::vector<std::vector<double>> h_old(n + 1, std::vector<double>(m + 1, 0.0));
//   std::vector<std::vector<double>> hs(n + 1, std::vector<double>(m + 1, 0.0));
//   std::vector<std::vector<double>> h_old2(n_2 + 1, std::vector<double>(m_2 + 1, 0.0));
//   std::vector<std::vector<double>> hs2(n_2 + 1, std::vector<double>(m_2 + 1, 0.0));

//   // Метод сопряженных градиентов
//   while (!flag) {
//       eps_max = 0.0;
//       // (n, m)
//       for (int j = 0; j < m + 1; j++)
//           for (int i = 0; i < n + 1; i++) {
//               v_old[i][j] = v[i][j];
//               h_old[i][j] = hs[i][j];
//           }
//       long double alpha = 0, beta = 0;
//       long double tmp_beta = 0, tmp_alpha = 0;

//       for (int j = 1; j < m; j++) {
//           for (int i = 1; i < n; i++) {
//               r[i][j] = A * v_old[i][j] + h2 * v_old[i + 1][j] + h2 * v_old[i - 1][j] +
//                   k2 * v_old[i][j + 1] + k2 * v_old[i][j - 1] - f[i][j];

//               double Ah_old = (A * h_old[i][j] + h2 * h_old[i + 1][j] + h2 * h_old[i - 1][j] +
//                   k2 * h_old[i][j + 1] + k2 * h_old[i][j - 1]);
//               if (i == 1) {
//                   r[i][j] -= h2 * v_old[i - 1][j] - mu1(c + j * k) * h2;
//                   Ah_old -= h2 * h_old[i - 1][j];
//               }
//               if (j == 1) {
//                   r[i][j] -= k2 * v_old[i][j - 1] - mu3(a + i * h) * k2;
//                   Ah_old -= k2 * h_old[i][j - 1];
//               }
//               if (i == n - 1) {
//                   r[i][j] -= h2 * v_old[i + 1][j] - mu2(c + j * k) * h2;
//                   Ah_old -= h2 * h_old[i + 1][j];
//               }
//               if (j == m - 1) {
//                   r[i][j] -= k2 * v_old[i][j + 1] - mu4(a + i * h) * k2;
//                   Ah_old -= k2 * h_old[i][j + 1];
//               }

//               beta += Ah_old * r[i][j];
//               tmp_beta += Ah_old * h_old[i][j];
//           }
//       }

//       if (tmp_beta == 0)
//           beta = 0;
//       else
//           beta /= tmp_beta;

//       for (int j = 1; j < m; ++j)
//           for (int i = 1; i < n; ++i)
//               hs[i][j] = -r[i][j] + beta * h_old[i][j];

//       for (int j = 1; j < m; ++j) {
//           for (int i = 1; i < n; ++i) {
//               double Ah = (A * hs[i][j] + h2 * hs[i + 1][j] + h2 * hs[i - 1][j] +
//                   k2 * hs[i][j + 1] + k2 * hs[i][j - 1]);
//               if (i == 1) {
//                   Ah -= h2 * hs[i - 1][j];
//               }
//               if (j == 1) {
//                   Ah -= k2 * hs[i][j - 1];
//               }
//               if (i == n - 1) {
//                   Ah -= h2 * hs[i + 1][j];
//               }
//               if (j == m - 1) {
//                   Ah -= k2 * hs[i][j + 1];
//               }
//               alpha -= r[i][j] * hs[i][j];
//               tmp_alpha += Ah * hs[i][j];
//           }
//       }
//       alpha /= tmp_alpha;

//       for (int j = 1; j < m; ++j) {
//           for (int i = 1; i < n; ++i) {
//               v[i][j] = v_old[i][j] + alpha * hs[i][j];

//               eps_cur = std::abs(v_old[i][j] - v[i][j]);
//               if (eps_cur > eps_max) {
//                   eps_max = eps_cur;
//               }
//           }
//       }

//       s++;
//       if (eps_max <= eps || s >= Nmax) {
//           flag = true;
//       }
//   }
//   //==========================================================================
//   while (!flag2) {
//       eps_max2 = 0.0;
//       // (n, m)
//       for (int j = 0; j < m_2 + 1; j++)
//           for (int i = 0; i < n_2 + 1; i++) {
//               v_old2[i][j] = v2[i][j];
//               h_old2[i][j] = hs2[i][j];
//           }
//       long double alpha = 0, beta = 0;
//       long double tmp_beta = 0, tmp_alpha = 0;

//       for (int j = 1; j < m_2; j++) {
//           for (int i = 1; i < n_2; i++) {
//               r2[i][j] = A_2 * v_old2[i][j] + h2_2 * v_old2[i + 1][j] + h2_2 * v_old2[i - 1][j] +
//                   k2_2 * v_old2[i][j + 1] + k2_2 * v_old2[i][j - 1] - f2[i][j];
//               double Ah_old = (A_2 * h_old2[i][j] + h2_2 * h_old2[i + 1][j] + h2_2 * h_old2[i - 1][j] +
//                   k2_2 * h_old2[i][j + 1] + k2_2 * h_old2[i][j - 1]);
//               if (i == 1) {
//                   r2[i][j] -= h2_2 * v_old2[i - 1][j] - mu1(c + j * k_2) * h2_2;
//                   Ah_old -= h2_2 * h_old2[i - 1][j];
//               }
//               if (j == 1) {
//                   r2[i][j] -= k2_2 * v_old2[i][j - 1] - mu3(a + i * h_2) * k2_2;
//                   Ah_old -= k2_2 * h_old2[i][j - 1];
//               }
//               if (i == n_2 - 1) {
//                   r2[i][j] -= h2_2 * v_old2[i + 1][j] - mu2(c + j * k_2) * h2_2;
//                   Ah_old -= h2_2 * h_old2[i + 1][j];
//               }
//               if (j == m_2 - 1) {
//                   r2[i][j] -= k2_2 * v_old2[i][j + 1] - mu4(a + i * h_2) * k2_2;
//                   Ah_old -= k2_2 * h_old2[i][j + 1];
//               }

//               beta += Ah_old * r2[i][j];
//               tmp_beta += Ah_old * h_old2[i][j];
//           }
//       }

//       if (tmp_beta == 0) {
//           beta = 0;
//       }
//       else beta /= tmp_beta;

//       for (int j = 1; j < m_2; ++j)
//           for (int i = 1; i < n_2; ++i)
//               hs2[i][j] = -r2[i][j] + beta * h_old2[i][j];


//       for (int j = 1; j < m_2; ++j) {
//           for (int i = 1; i < n_2; ++i) {
//               double Ah = (A_2 * hs2[i][j] + h2_2 * hs2[i + 1][j] + h2_2 * hs2[i - 1][j] +
//                   k2_2 * hs2[i][j + 1] + k2_2 * hs2[i][j - 1]);
//               if (i == 1) {
//                   Ah -= h2_2 * hs2[i - 1][j];
//               }
//               if (j == 1) {
//                   Ah -= k2_2 * hs2[i][j - 1];
//               }
//               if (i == n_2 - 1) {
//                   Ah -= h2_2 * hs2[i + 1][j];
//               }
//               if (j == m_2 - 1) {
//                   Ah -= k2_2 * hs2[i][j + 1];
//               }
//               alpha -= r2[i][j] * hs2[i][j];
//               tmp_alpha += Ah * hs2[i][j];
//           }
//       }
//       alpha /= tmp_alpha;

//       for (int j = 1; j < m_2; ++j) {
//           for (int i = 1; i < n_2; ++i) {
//               v2[i][j] = v_old2[i][j] + alpha * hs2[i][j];

//               eps_curr2 = std::abs(v_old2[i][j] - v2[i][j]);
//               if (eps_curr2 > eps_max2) {
//                   eps_max2 = eps_curr2;
//               }
//           }
//       }

//       s2++;
//       if (eps_max2 <= eps || s2 >= Nmax) {
//           flag2 = true;
//       }
//   }
//   //==========================================

//   R_N_2 = 0.0;
//   for (int i = 0; i < n + 1; i++)
//       for (int j = 0; j < m + 1; j++) {
//           R_N_2 += r[i][j] * r[i][j];
//       }
//   R_N_2 = std::sqrt(R_N_2);

//   R2_N_2 = 0.0;
//   for (int i = 0; i < n_2 + 1; i++)
//       for (int j = 0; j < m_2 + 1; j++) {
//           R2_N_2 += r2[i][j] * r2[i][j];
//       }
//   R2_N_2 = std::sqrt(R2_N_2);

//   for (int j = 0; j < m + 1; j++)
//       for (int i = 0; i < n + 1; i++) {
//           if (std::abs(v2[2 * i][2 * j] - v[i][j]) > solve_eps) {
//               solve_eps = std::abs(v2[2 * i][2 * j] - v[i][j]);
//           }

//       }

//   this->label27->Text = "Метод:    метод сопряженных градиентов";
//   this->label25->Text = ("Количество итераций:" + "   " + Convert::ToString(s));
//   this->label24->Text = ("Достиг. точность метода:" + "   " + Convert::ToString(eps_max));
//   this->label23->Text = ("Невязка СЛАУ на нач. приближении R0_2:" + "   " + Convert::ToString(round(R_0_2 * 1e10) / 1e10));
//   this->label14->Text = ("Невязка СЛАУ на нач. приближении R0_inf:" + "   " + Convert::ToString(round(R_0_inf * 1e10) / 1e10));

//   this->label22->Text = ("Схема на сетке решена с невязкой Rs_2:" + "   " + Convert::ToString(R_N_2));
//   this->label30->Text = ("Достиг. точность решения:" + "   " + Convert::ToString(solve_eps));

//   this->label32->Text = "Метод:    метод сопряженных градиентов";
//   this->label40->Text = ("Количество итераций:" + "   " + Convert::ToString(s2));
//   this->label38->Text = ("Достиг. точность метода:" + "   " + Convert::ToString(eps_max2));
//   this->label39->Text = ("Невязка СЛАУ на нач. приближении R2_0_2:" + "   " + Convert::ToString(round(R2_0_2 * 1e10) / 1e10));
//   this->label15->Text = ("Невязка СЛАУ на нач. приближении R2_0_inf:" + "   " + Convert::ToString(round(R2_0_inf * 1e10) / 1e10));
//   this->label37->Text = ("Схема на сетке решена с невязкой R2_s_2:" + "   " + Convert::ToString(R2_N_2));

//   for (int j = 0; j <= m; j++)
//       for (int i = 0; i <= n; i++) {
//           if (std::abs(v2[2 * i][2 * j] - v[i][j]) > v2v_max) {
//               v2v_max = std::abs(v2[2 * i][2 * j] - v[i][j]);
//               x_max = a + i * h;
//               y_max = c + j * k;
//           }
//       }

//   //this->label21->Text = ("Макс. разность двух решений |V2-V|:" + "   " + Convert::ToString(round(v2v_max * 1e10) / 1e10));
//   this->label20->Text = ("Макс. разность двух решений |V2-V|\n достигается в точке (x, y):" + "   (" +
//       Convert::ToString(round(x_max * 1e8) / 1e8) + ", " + Convert::ToString(round(y_max * 1e10) / 1e10) + ")");

//   int newm = m;
//   int newn = n;

//   //if (this->radioButtonV_main->Checked) {
//   dataGridView2->Rows->Clear();
//   dataGridView2->Columns->Clear();
//   this->dataGridView2->Columns->Add("", "j / i");
//   for (int j = 0; j <= std::min(newm, 100); j++) {
//       this->dataGridView2->Rows->Add(Convert::ToString(j));
//   }

//   for (int i = 0; i <= std::min(newn, 100); i++) {
//       this->dataGridView2->Columns->Add("", Convert::ToString(i));
//   }

//   for (int j = 0; j <= std::min(newm, 100); j++)
//       for (int i = 0; i <= std::min(newn, 100); i++)
//           this->dataGridView2->Rows[j]->Cells[i + 1]->Value = round(v[i][j] * 1e8) / 1e8;
//   //}

//   //if (this->radioButtonV2_main->Checked) {
//   dataGridView3->Rows->Clear();
//   dataGridView3->Columns->Clear();
//   this->dataGridView3->Columns->Add("", "j / i");
//   for (int j = 0; j <= std::min(m_2, 100); j++) {
//       this->dataGridView3->Rows->Add(Convert::ToString(j));
//   }

//   for (int i = 0; i <= std::min(n_2, 100); i++) {
//       this->dataGridView3->Columns->Add("", Convert::ToString(i));
//   }

//   for (int i = 0; i <= std::min(n_2, 100); i++)
//       for (int j = 0; j <= std::min(m_2, 100); j++)
//           this->dataGridView3->Rows[j]->Cells[i + 1]->Value = round(v2[i][j] * 1e8) / 1e8;
//   //}

//   //if (this->radioButtonV2V_main->Checked) {
//   dataGridView4->Rows->Clear();
//   dataGridView4->Columns->Clear();
//   this->dataGridView4->Columns->Add("", "j / i");
//   for (int j = 0; j <= std::min(newm, 100); j++) {
//       this->dataGridView4->Rows->Add(Convert::ToString(j));
//   }

//   for (int i = 0; i <= std::min(newn, 100); i++) {
//       this->dataGridView4->Columns->Add("", Convert::ToString(i));
//   }

//   for (int i = 0; i <= std::min(newn, 100); i++)
//       for (int j = 0; j <= std::min(newm, 100); j++)
//           this->dataGridView4->Rows[j]->Cells[i + 1]->Value = round(std::abs(v2[2 * i][2 * j] - v[i][j]) * 1e8) / 1e8;
//   //}


//   //if (radioButtonDrawVMain->Checked) {
//   chart4->ChartAreas[0]->RecalculateAxesScale();
//   double x3, y3, z3;
//   double xstep3 = (b - a) / n, ystep3 = (d - c) / m;
//   double max3 = v[0][0], min3 = max3;
//   Color Color3;
//   for (x3 = a; x3 <= b; x3 += xstep3)
//       for (y3 = c; y3 <= d; y3 += ystep3) {
//           int i = int((x3 - a) / xstep3);
//           int j = int((y3 - c) / ystep3);
//           z3 = v[i][j];
//           if (z3 > max3) max3 = z3;
//           else if (z3 < min3) min3 = z3;
//       }
//   for (int i = 0; i <= n; i++) {
//       Series^ Series1 = gcnew Series();
//       Series1->ChartType = SeriesChartType::Spline;
//       Series1->Color = Color::Blue;
//       chart4->Series->Add(Series1);
//       for (int j = 0; j <= m; j++) {
//           x3 = a + i * xstep3;
//           y3 = c + j * ystep3;
//           z3 = v[i][j];
//           chart4->Series[i]->Points->AddXY((double)y3, z3);
//           Color3 = calculateColorBlue(min3, z3, max3);
//           chart4->Series[i]->Points[j]->Color = Color3;
//           chart4->Series[i]->Points[j]->BackSecondaryColor = Color3;
//       }
//   }
//   //}

//   //if (radioButtonDrawV2->Checked) {
//   chart5->ChartAreas[0]->RecalculateAxesScale();
//   double x4, y4, z4;
//   double xstep4 = (b - a) / n_2, ystep4 = (d - c) / m_2; //øàãè
//   double max4 = v2[0][0], min4 = max4;
//   Color Color4;
//   for (x4 = a; x4 <= b; x4 += xstep4)
//       for (y4 = c; y4 <= d; y4 += ystep4) {
//           int i = int((x4 - a) / xstep4);
//           int j = int((y4 - c) / ystep4);
//           z4 = v2[i][j];
//           if (z4 > max4) max4 = z4;
//           else if (z4 < min4) min4 = z4;
//       }
//   for (int i = 0; i <= n_2; i++) {
//       Series^ Series1 = gcnew Series();
//       Series1->ChartType = SeriesChartType::Spline;
//       Series1->Color = Color::Blue;
//       chart5->Series->Add(Series1);
//       for (int j = 0; j <= m_2; j++) {
//           x4 = a + i * xstep4;
//           y4 = c + j * ystep4;
//           z4 = v2[i][j];
//           chart5->Series[i]->Points->AddXY((double)y4, z4);
//           Color4 = calculateColorPurple(min4, z4, max4);
//           chart5->Series[i]->Points[j]->Color = Color4;
//           chart5->Series[i]->Points[j]->BackSecondaryColor = Color4;
//       }
//   }
//   //}

//   //if (radioButtonDrawVV2->Checked) {
//   chart6->ChartAreas[0]->RecalculateAxesScale();
//   double x5, y5, z5;
//   double xstep5 = (b - a) / n, ystep5 = (d - c) / m;
//   double max5 = v[0][0] - v2[0][0], min5 = max5;
//   Color Color5;
//   for (x5 = a; x5 <= b; x5 += xstep5)
//       for (y5 = c; y5 <= d; y5 += ystep5) {
//           int i = int((x5 - a) / xstep5);
//           int j = int((y5 - c) / ystep5);
//           z5 = v[i][j] - v2[2 * i][2 * j];
//           if (z5 > max5) max5 = z5;
//           else if (z5 < min5) min5 = z5;
//       }
//   for (int i = 0; i <= n; i++) {
//       Series^ Series1 = gcnew Series();
//       Series1->ChartType = SeriesChartType::Spline;
//       Series1->Color = Color::Blue;
//       chart6->Series->Add(Series1);
//       for (int j = 0; j <= m; j++) {
//           x5 = a + i * xstep5;
//           y5 = c + j * ystep5;
//           z5 = v[i][j] - v2[2 * i][2 * j];
//           chart6->Series[i]->Points->AddXY((double)y5, z5);
//           Color5 = calculateColorRed(min5, z5, max5);
//           chart6->Series[i]->Points[j]->Color = Color5;
//           chart6->Series[i]->Points[j]->BackSecondaryColor = Color5;
//       }
//   }
//   // }
//}*/

MainWindow::~MainWindow()
{
    delete ui;
}
