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

    bool idx_is_part_of_area_include(size_t i, size_t j);
    bool idx_is_part_of_area(size_t i, size_t j);

    void init_v(vecvec& v);
    void init_f(vecvec& f);

    void fill_r(vecvec& r, const vecvec& v, const vecvec& f);
    void fill_u_test(vecvec& u_test);
    void fill_table(QTableWidget* table, size_t columns, size_t rows, const vecvec& values);

    long double get_beta(const vecvec& r, const vecvec& h);
    long double get_alpha(const vecvec& r, const vecvec& h);
    double get_R(const vecvec& r);

//    double mu1_test(double y) {
//        double res = std::exp(-y * y); //x = a = -1
//        return res;
//    }

//    double mu2_test(double x) {
//        double res = std::exp(-x * x); //y = d = 1
//        return res;
//    }

//    double mu3_test(double y) {
//        double res = std::exp(1 - y * y); //x = a + (b-a)/2 = 0
//        return res;
//    }

//    double mu4_test(double x) {
//        double res = std::exp(1 - x * x); //y = c + (d-c)/2 = 0
//        return res;
//    }

//    double mu5_test(double y) {
//        double res = std::exp(-y * y); //x = b = 1
//        return res;
//    }

//    double mu6_test(double x) {
//        double res = std::exp(-x * x); //y = c = -1
//        return res;
//    }

//    double f_test(double x, double y) {
//        return -4 * (1 - x * x - y * y) * std::exp(1 - x * x - y * y);
//    }
    double u_test(double i, double j, int _n, int _m) {
        const double h = (b - a) / _n;
        const double k = (d - c) / _m;
        double x = a + i * h;
        double y = c + j * k;
        return std::exp(1 - x * x - y * y);
    }

    double mu1_test(double y) {
        return  std::exp(-y * y);
    }

    double mu2_test(double y) {
        return std::exp(-y * y);
    }

    double mu3_test(double x) {
        return std::exp(-x * x);
    }

    double mu4_test(double x) {
        return std::exp(-x * x);
    }

    double mu5_test(double y) {
        double res = std::exp(1 - 0.25 - y * y); //x = a + (b-a)/4 = -0.5
        return res;
    }

    double mu6_test(double y) {
        double res = std::exp(1 - 0.25 - y * y); //x = a + 3(b-a)/4 = 0.5
        return res;
    }

    double mu7_test(double x) {
        double res = std::exp(1 - 0.25 - x * x); //y = c + (d-c)/4 = -0.5
        return res;
    }
    double mu8_test(double x) {
        double res = std::exp(1 - 0.25 - x * x); //y = c + 3(d-c)/4 = 0.5
        return res;
    }

    double f_test(double x, double y) {
        return -4 * (1 - x * x - y * y) * std::exp(1 - x * x - y * y);
    }

    bool save_to_file(const vecvec& values, const QString& filename);

private:
    size_t Nmax;
    double eps;
    size_t n, m;
    const double a = -1.0;
    const double b = 1.0;
    const double c = -1.0;
    const double d = 1.0;
    double A;
    double h;
    double k;
    double h2;
    double k2;
};

#endif // MAINWINDOW_H
