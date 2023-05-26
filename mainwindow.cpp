#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "ui_mainwindow.h"

#include <vector>
#include <algorithm>
#include <QFile>

using vec = std::vector<double>;
using vecvec = std::vector<std::vector<double>>;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

bool MainWindow::idx_is_part_of_area_include(size_t i, size_t j)
{
    if (j <= m / 4 || j >= m * 3 / 4)
        return true;
    else
        return i <= n / 4 || i >= n * 3 / 4;
}

bool MainWindow::idx_is_part_of_area(size_t i, size_t j)
{
    if (j < m / 4 || j > m * 3 / 4)
        return true;
    else
        return i < n / 4 || i > n * 3 / 4;
}

void MainWindow::fill_u_test(vecvec& u_test_vec)
{
    for (int j = 0; j < m + 1; j++)
       for (int i = 0; i < n + 1; i++) {
           if (idx_is_part_of_area_include(i,j))
           {
               u_test_vec[i][j] = u_test(i, j, n, m);
           }

       }
}

void MainWindow::fill_table(QTableWidget* table, size_t columns, size_t rows, const vecvec& values)
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

bool MainWindow::save_to_file(const vecvec& values, const QString& filename)
{
    bool success = false;

    QFile data(filename);
    if (data.open(QIODevice::ReadWrite))
    {
        QTextStream output(&data);

        for (size_t i = 0; i <= n; i++)
            output << a + i*h << ',';
        output << '\n';

        for (size_t j = 0; j <= m; j++)
            output << c + j*k << ',';
        output << '\n';

        for (size_t i = 0; i <= n; i++)
        {
            for (size_t j = 0; j <= m; j++)
            {
                output << values[i][j] << ',';
            }
        }
        output << '\n';

        success = true;
    }

    data.close();
    return success;
}

void MainWindow::on_pushButton_clicked()
{
    n = ui->n->text().toUInt();
    m = ui->m->text().toUInt();
    eps = ui->eps->text().toDouble();
    Nmax = ui->Nmax->text().toUInt();

    int s = 0;
    double eps_max = 0.0;
    double eps_cur = 0.0;
    double solve_err = 0.0;
    double solve_eps = 0.0;
    h = (b - a) / n;
    k = (d - c) / m;

    double uv_max = 0.0;
    double x_max, y_max;

    std::vector<std::vector<double>> v(n + 1, std::vector<double>(m + 1, 0.0));
    std::vector<std::vector<double>> f(n + 1, std::vector<double>(m + 1, 0.0));
    std::vector<std::vector<double>> r(n + 1, std::vector<double>(m + 1, 0.0));
    std::vector<std::vector<double>> hs(n + 1, std::vector<double>(m + 1, 0.0));

    double R_N_2 = 0.0;
    double R_0_2 = 0.0;
    double R_0_inf = 0.0;
    bool flag = false;

    h2 = (1 / (h * h));
    k2 = (1 / (k * k));
    A = -2 * (h2 + k2);

    for (size_t j = 0; j <= m; j++)
    {
        double y = c + j * k;
        v[0][j] = mu1_test(y);
        v[n][j] = mu2_test(y);
    }
    for (size_t i = 0; i <= n; i++)
    {
        double x = a + i * h;
        v[i][0] = mu3_test(x);
        v[i][m] = mu4_test(x);
    }
    for (size_t j = m / 4; j <= m * 3 / 4; j++)
    {
        double y = c + j * k;
        v[n/4][j] = mu5_test(y);
        v[n*3/4][j] = mu6_test(y);
    }

    for (size_t i = n/4; i <= n*3/4; i++)
    {
        double x = a + i * h;
        v[i][m/4] = mu7_test(x);
        v[i][m*3/4] = mu8_test(x);
    }


    const size_t columns = std::min<size_t>(101, n+1);
    const size_t rows = std::min<size_t>(101, m+1);

    fill_table(ui->v0_table, columns, rows, v);
    save_to_file(v, "v0");

    std::vector<std::vector<double>> u_test_vec(n + 1, std::vector<double>(m + 1, 0.0));
    fill_u_test(u_test_vec);
    fill_table(ui->u_table, columns, rows, u_test_vec);
    save_to_file(u_test_vec, "u");

    for (int j = 0; j <= m; j++)
       for (int i = 0; i <= n; i++) {
           if (idx_is_part_of_area_include(i,j))
           {
               double x = a + i * h;
               double y = c + j * k;
               f[i][j] = f_test(x, y);
           }

       }

    //невязка
    for (int j = 1; j < m; j++)
       for (int i = 1; i < n; i++) {
           if (idx_is_part_of_area(i,j))
           {
               r[i][j] = A * v[i][j] + h2 * v[i + 1][j] + h2 * v[i - 1][j] + k2 * v[i][j + 1] + k2 * v[i][j - 1] - f[i][j];

               if (i == 1)         r[i][j] -= h2 * v[i - 1][j] - mu1_test(c + j * k) * h2;
               if (i == n - 1)     r[i][j] -= h2 * v[i + 1][j] - mu2_test(c + j * k) * h2;

               if (j == 1)         r[i][j] -= k2 * v[i][j - 1] - mu3_test(a + i * h) * k2;
               if (j == m - 1)     r[i][j] -= k2 * v[i][j + 1] - mu4_test(a + i * h) * k2;

               if (i == n*3/4 + 1) r[i][j] -= h2 * v[i - 1][j] - mu6_test(c + j * k) * h2;
               if (i == n/4 - 1)   r[i][j] -= h2 * v[i + 1][j] - mu5_test(c + j * k) * h2;

               if (j == m*3/4 + 1) r[i][j] -= k2 * v[i][j - 1] - mu8_test(a + i * h) * k2;
               if (j == m/4 - 1)   r[i][j] -= k2 * v[i][j + 1] - mu7_test(a + i * h) * k2;
           }

       }

    //евклидова норма невязки для начального приближения
    for (int i = 0; i < n + 1; i++)
       for (int j = 0; j < m + 1; j++) {
           if (idx_is_part_of_area_include(i,j))
           {
               R_0_inf = std::max(R_0_inf, std::abs(r[i][j]));
               R_0_2 += r[i][j] * r[i][j];
           }

       }
    R_0_2 = std::sqrt(R_0_2);

    std::vector<std::vector<double>> v_old(n + 1, std::vector<double>(m + 1, 0.0));
    std::vector<std::vector<double>> h_old(n + 1, std::vector<double>(m + 1, 0.0));
    while (!flag) {
       eps_max = 0.0;
       // (n, m)
       for (int j = 0; j < m + 1; j++)
           for (int i = 0; i < n + 1; i++) {
               if (idx_is_part_of_area_include(i,j))
               {
                   v_old[i][j] = v[i][j];
                   h_old[i][j] = hs[i][j];
               }

           }
       long double alpha = 0, beta = 0;
       long double tmp_beta = 0, tmp_alpha = 0;

       for (int j = 1; j < m; j++) {
           for (int i = 1; i < n; i++) {
               if (idx_is_part_of_area(i,j))
               {
                   r[i][j] = A * v_old[i][j] + h2 * v_old[i + 1][j] + h2 * v_old[i - 1][j] +
                       k2 * v_old[i][j + 1] + k2 * v_old[i][j - 1] - f[i][j];

                   if (i == 1)         r[i][j] -= h2 * v[i - 1][j] - mu1_test(c + j * k) * h2;
                   if (i == n - 1)     r[i][j] -= h2 * v[i + 1][j] - mu2_test(c + j * k) * h2;

                   if (j == 1)         r[i][j] -= k2 * v[i][j - 1] - mu3_test(a + i * h) * k2;
                   if (j == m - 1)     r[i][j] -= k2 * v[i][j + 1] - mu4_test(a + i * h) * k2;

                   if (i == n*3/4 + 1) r[i][j] -= h2 * v[i - 1][j] - mu6_test(c + j * k) * h2;
                   if (i == n/4 - 1)   r[i][j] -= h2 * v[i + 1][j] - mu5_test(c + j * k) * h2;

                   if (j == m*3/4 + 1) r[i][j] -= k2 * v[i][j - 1] - mu8_test(a + i * h) * k2;
                   if (j == m/4 - 1)   r[i][j] -= k2 * v[i][j + 1] - mu7_test(a + i * h) * k2;


                   double Ah_old = (A * h_old[i][j] + h2 * h_old[i + 1][j] + h2 * h_old[i - 1][j] +
                       k2 * h_old[i][j + 1] + k2 * h_old[i][j - 1]);

                   if (i == 1)         Ah_old -= h2 * h_old[i - 1][j];
                   if (j == 1)         Ah_old -= k2 * h_old[i][j - 1];

                   if (i == n - 1)     Ah_old -= h2 * h_old[i + 1][j];
                   if (j == m - 1)     Ah_old -= k2 * h_old[i][j + 1];

                   if (i == n*3/4 + 1) Ah_old -= h2 * h_old[i - 1][j];
                   if (i == n/4 - 1)   Ah_old -= h2 * h_old[i + 1][j];

                   if (j == m*3/4 + 1) Ah_old -= k2 * h_old[i][j - 1];
                   if (j == m/4 - 1)   Ah_old -= k2 * h_old[i][j + 1];

                   beta += Ah_old * r[i][j];
                   tmp_beta += Ah_old * h_old[i][j];
               }

           }
       }

       if (tmp_beta == 0) {
           beta = 0;
       }
       else beta /= tmp_beta;

       for (int j = 1; j < m; ++j)
           for (int i = 1; i < n; ++i)
               if (idx_is_part_of_area(i,j))
               {
                   hs[i][j] = -r[i][j] + beta * h_old[i][j];
               }



       for (int j = 1; j < m; ++j) {
           for (int i = 1; i < n; ++i) {
               if (idx_is_part_of_area(i,j))
               {
                   double Ah = (A * hs[i][j] + h2 * hs[i + 1][j] + h2 * hs[i - 1][j] +
                       k2 * hs[i][j + 1] + k2 * hs[i][j - 1]);

                   if (i == 1)         Ah -= h2 * hs[i - 1][j];
                   if (j == 1)         Ah -= k2 * hs[i][j - 1];

                   if (i == n - 1)     Ah -= h2 * hs[i + 1][j];
                   if (j == m - 1)     Ah -= k2 * hs[i][j + 1];

                   if (i == n*3/4 + 1) Ah -= h2 * hs[i - 1][j];
                   if (i == n/4 - 1)   Ah -= h2 * hs[i + 1][j];

                   if (j == m*3/4 + 1) Ah -= k2 * hs[i][j - 1];
                   if (j == m/4 - 1)   Ah -= k2 * hs[i][j + 1];

                   alpha -= r[i][j] * hs[i][j];
                   tmp_alpha += Ah * hs[i][j];
               }

           }
       }
       alpha /= tmp_alpha;

       for (int j = 1; j < m; ++j) {
           for (int i = 1; i < n; ++i) {
               if (idx_is_part_of_area(i,j))
               {
                   v[i][j] = v_old[i][j] + alpha * hs[i][j];

                   eps_cur = std::abs(v_old[i][j] - v[i][j]);
                   if (eps_cur > eps_max) {
                       eps_max = eps_cur;
                   }
               }

           }
       }

       s++;
       if (eps_max <= eps || s >= Nmax) {
           flag = true;
       }
    }
    // End

    fill_table(ui->vn_table, columns, rows, v);
    save_to_file(v, "vN");
    std::vector<std::vector<double>> uv(n + 1, std::vector<double>(m + 1, 0.0));

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < columns; j++)
        {
            uv[i][j] = u_test_vec[i][j]-v[i][j];
        }
    }

    fill_table(ui->uv_table, columns, rows, uv);
    save_to_file(uv, "uv");

    R_N_2 = 0.0;
    for (int i = 0; i < n + 1; i++)
       for (int j = 0; j < m + 1; j++) {
           if (idx_is_part_of_area_include(i,j))
           {
               R_N_2 += r[i][j] * r[i][j];
           }

       }
    R_N_2 = std::sqrt(R_N_2);

    for (int j = 0; j < m + 1; j++)
       for (int i = 0; i < n + 1; i++) {
           if (idx_is_part_of_area_include(i,j))
           {
               if (std::abs(u_test(i, j, n, m) - v[i][j]) > uv_max) {
                   uv_max = std::abs(u_test(i, j, n, m) - v[i][j]);
                   x_max = a + i * h;
                   y_max = c + j * k;
               }
           }

       }

    ui->eps_max->setText(QString::number(eps_max));
    ui->N->setText(QString::number(s));
    ui->UV->setText(QString::number(uv_max));
    ui->UV_coord->setText("(" + QString::number(x_max) + ", " + QString::number(y_max) + ")");
    ui->R0->setText(QString::number(R_0_2));
    ui->RN->setText(QString::number(R_N_2));
}

MainWindow::~MainWindow()
{
    delete ui;
}

//void MainWindow::on_pushButton_clicked()
//{
//    n = ui->n->text().toUInt();
//    m = ui->m->text().toUInt();
//    eps = ui->eps->text().toDouble();
//    Nmax = ui->Nmax->text().toUInt();

//    h = (b - a) / n;
//    k = (d - c) / m;
//    h2 = (1 / (h * h));
//    k2 = (1 / (k * k));
//    A = -2 * (h2 + k2);

//    double eps_max = 0.0;
//    double eps_cur = 0.0;

//    double uv_max = 0.0;
//    double x_max, y_max;

//    std::vector<std::vector<double>> v(n + 1, std::vector<double>(m + 1, 0.0));
//    std::vector<std::vector<double>> f(n + 1, std::vector<double>(m + 1, 0.0));
//    std::vector<std::vector<double>> r(n + 1, std::vector<double>(m + 1, 0.0));
//    std::vector<std::vector<double>> hs(n + 1, std::vector<double>(m + 1, 0.0));

//    init_v(v);
//    init_f(f);
//    fill_r(r, v, f);

//    double R0 = get_R(r);

//    const size_t columns = std::min<size_t>(101, n+1);
//    const size_t rows = std::min<size_t>(101, m+1);

//    fill_table(ui->v0_table, columns, rows, v);

//    std::vector<std::vector<double>> v_old(n + 1, std::vector<double>(m + 1, 0.0));
//    std::vector<std::vector<double>> h_old(n + 1, std::vector<double>(m + 1, 0.0));

//    size_t steps = 0;
//    for(; steps < Nmax; steps++)
//    {
//        eps_max = 0.0;

//        for (size_t j = 0; j <= m; j++)
//        {
//            for (size_t i = 0; i <= n; i++)
//            {
//                if (idx_is_part_of_area_include(i, j))
//                {
//                    v_old[i][j] = v[i][j];
//                    h_old[i][j] = hs[i][j];
//                }
//            }
//        }

//        fill_r(r, v_old, f);
//        long double beta = get_beta(r, h_old);

//        for (size_t j = 1; j < m; j++)
//        {
//            for (size_t i = 1; i < n; i++)
//            {
//                if (idx_is_part_of_area(i, j))
//                    hs[i][j] = -r[i][j] + beta * h_old[i][j];
//            }
//        }

//        long double alpha = get_alpha(r, hs);

//        for (size_t j = 1; j < m; j++)
//        {
//            for (size_t i = 1; i < n; i++)
//            {
//                if (idx_is_part_of_area(i, j))
//                {
//                    v[i][j] = v_old[i][j] + alpha * hs[i][j];

//                    eps_cur = std::abs(v_old[i][j] - v[i][j]);
//                    if (eps_cur > eps_max) {
//                        eps_max = eps_cur;
//                    }
//                }
//            }
//        }

//        if (eps_max < eps)
//            break;
//    }

//    double RN = get_R(r);
//    std::vector<std::vector<double>> u_test(n + 1, std::vector<double>(m + 1, 0.0));
//    fill_u_test(u_test);

//    for (size_t j = 0; j <= m; j++)
//    {
//        for (size_t i = 0; i <= n; i++)
//        {
//            if (idx_is_part_of_area_include(i, j))
//            {
//                if (std::abs(u_test[i][j] - v[i][j]) > uv_max) {
//                    uv_max = std::abs(u_test[i][j] - v[i][j]);
//                    x_max = a + i * h;
//                    y_max = c + j * k;
//                }
//            }
//        }
//    }

//    fill_table(ui->u_table, columns, rows, u_test);
//    fill_table(ui->vn_table, columns, rows, v);

//    std::vector<std::vector<double>> uv(n + 1, std::vector<double>(m + 1, 0.0));

//    for (size_t i = 0; i < rows; i++)
//    {
//        for (size_t j = 0; j < columns; j++)
//        {
//            uv[i][j] = u_test[i][j]-v[i][j];
//        }
//    }

//    fill_table(ui->uv_table, columns, rows, uv);
//    ui->eps_max->setText(QString::number(eps_max));
//    ui->N->setText(QString::number(steps));
//    ui->UV->setText(QString::number(uv_max));
//    ui->UV_coord->setText("(" + QString::number(x_max) + ", " + QString::number(y_max) + ")");
//    ui->R0->setText(QString::number(R0));
//    ui->RN->setText(QString::number(RN));
//}

void MainWindow::init_v(vecvec& v)
{
//    for (size_t j = 0; j <= m; j++)
//    {
//        double y = c + j * k;
//        v[0][j] = mu1_test(y);
//    }

//    //for (size_t i = 0; i <= n/2; i++)
//    for (size_t i = 0; i <= n; i++)
//    {
//        double x = a + i * h;
//        v[i][m] = mu2_test(x);
//    }

////    for (size_t j = m / 2; j <= m; j++)
////    {
////        double y = c + j * k;
////        v[n/2][j] = mu3_test(y);
////    }

////    for (size_t i = n/2; i <= n; i++)
////    {
////        double x = a + i * h;
////        v[i][m / 2] = mu4_test(x);
////    }

//    //for (size_t j = 0; j <= m / 2; j++)
//    for (size_t j = 0; j <= m; j++)
//    {
//        double y = c + j * k;
//        v[n][j] = mu5_test(y);
//    }

//    for (size_t i = 0; i <= n; i++)
//    {
//        double x = a + i * h;
//        v[i][0] = mu6_test(x);
//    }
}

void MainWindow::init_f(vecvec& f)
{
    for (size_t j = 0; j <= m; j++)
    {
        for (size_t i = 0; i <= n; i++)
        {
            if (idx_is_part_of_area_include(i, j))
            {
                double x = a + i * h;
                double y = c + j * k;
                f[i][j] = f_test(x, y);
            }
        }
    }
}

void MainWindow::fill_r(vecvec& r, const vecvec& v, const vecvec& f)
{
    for (size_t j = 1; j < m; j++)
    {
        for (size_t i = 1; i < n; i++)
        {
            if (idx_is_part_of_area(i, j))
            {
                r[i][j] = A * v[i][j] + h2 * v[i + 1][j] + h2 * v[i - 1][j] + k2 * v[i][j + 1] + k2 * v[i][j - 1] - f[i][j];

                const double x = a + i * h;
                const double y = c + j * k;

                if (i == 1)       r[i][j] -= h2 * v[i-1][j] - mu1_test(y) * h2;
                if (j == m - 1)   r[i][j] -= k2 * v[i][j+1] - mu2_test(x) * k2;
//                if (i == n/2 - 1) r[i][j] -= h2 * v[i+1][j] - mu3_test(y) * h2;
//                if (j == m/2 - 1) r[i][j] -= k2 * v[i][j+1] - mu4_test(x) * k2;
//                if (i == n - 1)   r[i][j] -= h2 * v[i+1][j] - mu5_test(y) * h2;
//                if (j == 1)       r[i][j] -= k2 * v[i][j-1] - mu6_test(x) * k2;
            }
        }
    }
}

long double MainWindow::get_beta(const vecvec& r, const vecvec& h)
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
                if (j == m - 1)   Ah -= k2 * h[i][j+1];
//                if (i == n/2 - 1) Ah -= h2 * h[i+1][j];
//                if (j == m/2 - 1) Ah -= k2 * h[i][j+1];
                if (i == n - 1)   Ah -= h2 * h[i+1][j];
                if (j == 1)       Ah -= k2 * h[i][j-1];

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

long double MainWindow::get_alpha(const vecvec& r, const vecvec& h)
{
    long double alpha = 0,  tmp_alpha = 0;

    for (size_t j = 1; j < m; j++)
    {
        for (size_t i = 1; i < n; i++)
        {
            if (idx_is_part_of_area(i, j))
            {
                double Ah = A * h[i][j] + h2 * h[i + 1][j] + h2 * h[i - 1][j] + k2 * h[i][j + 1] + k2 * h[i][j - 1];

                if (i == 1)       Ah -= h2 * h[i-1][j];
                if (j == m - 1)   Ah -= k2 * h[i][j+1];
//                if (i == n/2 - 1) Ah -= h2 * h[i+1][j];
//                if (j == m/2 - 1) Ah -= k2 * h[i][j+1];
                if (i == n - 1)   Ah -= h2 * h[i+1][j];
                if (j == 1)       Ah -= k2 * h[i][j-1];

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
    for (size_t i = 0; i <= n; i++)
    {
        for (size_t j = 0; j <= m; j++)
        {
            if (idx_is_part_of_area_include(i, j))
                R0 += r[i][j] * r[i][j];
        }
    }

    return std::sqrt(R0);
}
