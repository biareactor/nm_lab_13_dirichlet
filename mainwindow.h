#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <float.h>
#include <utility>
#include <QTableWidget>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_clicked();

private:
    Ui::MainWindow *ui;

private:
    using vec = std::vector<double>;
    using vecvec = std::vector<std::vector<double>>;

    bool idx_is_part_of_area(size_t i, size_t j);

    void init_v(vecvec& v, double h, double k);
    void init_f(vecvec& f, double h, double k);

    void fill_r(vecvec& r, const vecvec& v, const vecvec& f, double A, double h, double k, double h2, double k2);
    void fill_u_test(vecvec& u_test);
    void fill_table(QTableWidget* table, size_t columns, size_t rows, const vecvec& values);

    long double get_beta(double A, const vecvec& r, const vecvec& h, double h2, double k2);
    long double get_alpha(double A, const vecvec& r, const vecvec& h, double h2, double k2);
    double get_R(const vecvec& r);

    double mu1_test(double y) {
        double res = std::exp(-y * y); //x = a = -1
        return res;
    }

    double mu2_test(double x) {
        double res = std::exp(-x * x); //y = d = 1
        return res;
    }

    double mu3_test(double y) {
        double res = std::exp(1 - y * y); //x = a + (b-a)/2 = 0
        return res;
    }

    double mu4_test(double x) {
        double res = std::exp(1 - x * x); //y = c + (d-c)/2 = 0
        return res;
    }

    double mu5_test(double y) {
        double res = std::exp(-y * y); //x = b = 1
        return res;
    }

    double mu6_test(double x) {
        double res = std::exp(-x * x); //y = c = -1
        return res;
    }

    double f_test(double x, double y) {
        return 4 * (1 - x * x - y * y) * std::exp(1 - x * x - y * y);
    }

private:
    size_t Nmax;
    double eps;
    size_t n, m;
    const double a = -1.0;
    const double b = 1.0;
    const double c = -1.0;
    const double d = 1.0;
};

#endif // MAINWINDOW_H
